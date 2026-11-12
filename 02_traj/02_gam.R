# Lifespan GAM per APA site with individual-normalised profile weights.
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
})

# --- paths ---

root <- Sys.getenv("APA_ROOT")
DATA <- file.path(root, "data/02_traj")
APA_MATRIX <- file.path(root, "data/01_quant/apa_combat_donor_mean.tsv")
COV_FILE <- file.path(DATA, "00_cov.tsv")
SUPPORT_FILE <- file.path(DATA, "01_support.tsv")
PHENO_GROUP <- file.path(root, "raw/phenotype.group")

# --- settings ---

K_LADDER <- c(20L, 30L, 40L, 60L)
K_INDEX_P_THR <- 0.05
EDF_MARGIN <- 1.0
Q_THR <- 0.05
MAX_SITES <- as.integer(Sys.getenv("APA_GAM_MAX_SITES", "0"))

# --- inputs ---

cov <- fread(COV_FILE)
cov[, age_model := log10(pmax(as.numeric(age) + 1, 1e-6))]
support <- fread(SUPPORT_FILE)
grid_model <- support$age_model
supported <- as.logical(support$is_supported)
GRID_SIZE <- nrow(support)

apa_dt <- if (MAX_SITES > 0L) {
  fread(APA_MATRIX, nrows = MAX_SITES, showProgress = FALSE)
} else {
  fread(APA_MATRIX, showProgress = FALSE)
}
site_ids <- as.character(apa_dt[[1]])
mat <- as.matrix(apa_dt[, -1, with = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- site_ids
samples <- intersect(cov$sample, colnames(apa_dt)[-1])
stopifnot(length(samples) >= 3)
mat <- mat[, samples, drop = FALSE]
cov <- cov[match(samples, cov$sample)]

# --- site annotation ---

parts <- strsplit(site_ids, "|", fixed = TRUE)
ann <- data.table(
  site_id = site_ids,
  transcript_id = sapply(parts, function(x) if (length(x) >= 1) x[1] else ""),
  gene_id = sapply(parts, function(x) if (length(x) >= 2) x[2] else ""),
  locus = sapply(parts, function(x) if (length(x) >= 3) x[3] else ""),
  strand = sapply(parts, function(x) if (length(x) >= 4) x[4] else "")
)
if (file.exists(PHENO_GROUP) && any(ann$gene_id == "")) {
  pg <- fread(PHENO_GROUP, header = FALSE, fill = TRUE,
              col.names = c("site_id", "gene_id_pg"))
  pg <- unique(pg[, .(site_id = as.character(site_id), gene_id_pg = as.character(gene_id_pg))])
  ann <- merge(ann, pg, by = "site_id", all.x = TRUE, sort = FALSE)
  ann[gene_id == "", gene_id := gene_id_pg]
  ann[, gene_id_pg := NULL]
}

# --- model design ---

model_df <- as.data.frame(cov[, .(age_model, sex_binary, RIN_z, profile_weight)])
grid_df <- data.frame(age_model = grid_model)
grid_df$sex_binary <- weighted.mean(model_df$sex_binary, model_df$profile_weight)
grid_df$RIN_z <- 0
covar_rhs <- "sex_binary + RIN_z"

# --- fit every site ---

n_sites <- nrow(mat)
site_stats <- vector("list", n_sites)
traj_z <- matrix(NA_real_, n_sites, GRID_SIZE)
traj_logit <- matrix(NA_real_, n_sites, GRID_SIZE)

for (i in seq_len(n_sites)) {
  dat <- model_df
  dat$y <- as.numeric(mat[i, ])
  if (anyNA(dat$y)) dat$y[is.na(dat$y)] <- mean(dat$y, na.rm = TRUE)

  null_fit <- tryCatch(
    gam(as.formula(sprintf("y ~ %s", covar_rhs)), data = dat, weights = profile_weight, method = "REML"),
    error = function(e) NULL
  )
  best <- NULL
  if (!is.null(null_fit)) {
    for (k in K_LADDER) {
      full_form <- as.formula(sprintf("y ~ %s + s(age_model, k = %d, bs = 'tp')", covar_rhs, k))
      fit <- tryCatch(gam(full_form, data = dat, weights = profile_weight, method = "REML"),
                      error = function(e) NULL)
      if (is.null(fit)) next

      s_tab <- summary(fit)$s.table
      edf <- as.numeric(s_tab[1, "edf"])
      kc <- tryCatch(
        k.check(fit, subsample = 5000, n.rep = 400),
        error = function(e) matrix(NA_real_, 1, 4,
          dimnames = list(NULL, c("k'", "edf", "k-index", "p-value")))
      )
      k_p <- suppressWarnings(as.numeric(kc[1, "p-value"]))
      k_idx <- suppressWarnings(as.numeric(kc[1, "k-index"]))

      dev <- tryCatch(anova(null_fit, fit, test = "F"), error = function(e) NULL)
      age_p <- if (!is.null(dev)) suppressWarnings(as.numeric(dev$`Pr(>F)`[2])) else
        as.numeric(s_tab[1, "p-value"])
      if (!is.finite(age_p)) age_p <- as.numeric(s_tab[1, "p-value"])

      best <- list(k = k, edf = edf, k_p = k_p, k_idx = k_idx, age_p = age_p, fit = fit)
      k_ok <- is.finite(k_p) && k_p > K_INDEX_P_THR && !(is.finite(edf) && edf >= k - EDF_MARGIN)
      if (k_ok) break
    }
  }

  if (is.null(best)) {
    site_stats[[i]] <- data.table(
      site_id = site_ids[i], age_pvalue = NA_real_, age_qvalue = NA_real_,
      smooth_edf = NA_real_, k_index_pvalue = NA_real_, k_index = NA_real_,
      final_k = NA_integer_, delta_logit = NA_real_, amplitude_logit = NA_real_
    )
    next
  }

  pred_logit <- as.numeric(predict(best$fit, newdata = grid_df))

  dat_z <- dat
  weighted_y_mean <- weighted.mean(dat$y, dat$profile_weight)
  weighted_y_sd <- sqrt(weighted.mean((dat$y - weighted_y_mean)^2, dat$profile_weight))
  dat_z$y <- (dat$y - weighted_y_mean) / weighted_y_sd
  full_form_k <- as.formula(sprintf("y ~ %s + s(age_model, k = %d, bs = 'tp')", covar_rhs, best$k))
  fit_z <- tryCatch(gam(full_form_k, data = dat_z, weights = profile_weight, method = "REML"),
                    error = function(e) NULL)
  pred_z <- if (!is.null(fit_z)) as.numeric(predict(fit_z, newdata = grid_df)) else
    as.numeric(scale(pred_logit))

  sup_idx <- which(supported)
  if (length(sup_idx) >= 2) {
    y_sup <- pred_logit[sup_idx]
    delta <- y_sup[length(y_sup)] - y_sup[1]
    ampl <- max(y_sup) - min(y_sup)
  } else {
    delta <- ampl <- NA_real_
  }

  traj_z[i, ] <- pred_z
  traj_logit[i, ] <- pred_logit
  site_stats[[i]] <- data.table(
    site_id = site_ids[i], age_pvalue = best$age_p, age_qvalue = NA_real_,
    smooth_edf = best$edf, k_index_pvalue = best$k_p, k_index = best$k_idx,
    final_k = best$k, delta_logit = delta, amplitude_logit = ampl
  )
}

# --- multiple testing and outputs ---

stats_dt <- rbindlist(site_stats, fill = TRUE)
stats_dt[, age_qvalue := p.adjust(age_pvalue, method = "BH")]
stats_dt[, association_pass := age_qvalue < Q_THR]
stats_dt <- merge(ann, stats_dt, by = "site_id", all.y = TRUE, sort = FALSE)

all_names <- sprintf("grid_%03d", 0:(GRID_SIZE - 1L))
colnames(traj_z) <- all_names
colnames(traj_logit) <- all_names

fwrite(stats_dt, file.path(DATA, "02_sites.tsv"), sep = "\t")
fwrite(data.table(site_id = site_ids, as.data.table(traj_z)), file.path(DATA, "traj_z.csv"), sep = ",")
fwrite(data.table(site_id = site_ids, as.data.table(traj_logit)), file.path(DATA, "traj_logit.csv"), sep = ",")
