#!/usr/bin/env Rscript
# Cell-type-adjusted GAM age models for each directional APA event (t4).
# Reads 02_traj tables and apa_after_combat.tsv from data/04_qtl/; writes T4_celltype_adj_events.tsv to data/02_traj/.
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(readxl)
})

# --- config ---
out <- file.path(Sys.getenv("APA_ROOT"), "data", "02_traj")
k_ladder <- c(20L, 30L, 40L, 60L)
ilr_names <- paste0("cellprop_ILR", 1:6)
ilr_files <- c(prenatal = "covariates_with_cell_proportions.pre.txt",
               postnatal = "covariates_with_cell_proportions.postnatal.txt")

# --- loading data ---
sites <- fread(file.path(out, "02_sites_clustered.tsv"))[
  directional_pass == TRUE & cluster_id %in% c(1L, 2L)]
tx_map <- fread(file.path(Sys.getenv("APA_ROOT"), "ref", "gencode.v40.tx2gene.tsv"),
                select = c("tx_name", "gene_name"))
tx_map[, tx_base := sub("\\.[0-9]+$", "", tx_name)]
sites[, tx_base := sub("\\.[0-9]+$", "", transcript_id)]
sites <- merge(sites, unique(tx_map[, .(tx_base, gene_symbol = gene_name)],
                                by = "tx_base"), by = "tx_base", all.x = TRUE)
sites[is.na(gene_symbol) | gene_symbol == "", gene_symbol := gene_id]
cov <- fread(file.path(out, "00_cov.tsv"))
ilr <- rbindlist(lapply(names(ilr_files), function(stage) {
  tab <- fread(file.path(out, ilr_files[[stage]]))
  setnames(tab, 1L, "covariate")
  values <- as.matrix(tab[match(ilr_names, covariate), -1L, with = FALSE])
  storage.mode(values) <- "numeric"
  result <- as.data.table(t(values))
  setnames(result, ilr_names)
  result[, sample := names(tab)[-1L]]
  result
}))
cov <- merge(cov, ilr, by = "sample")
stopifnot(nrow(cov) == 697L, uniqueN(cov$sample) == 697L,
          all(is.finite(as.matrix(cov[, ..ilr_names]))))
cov[, profile_weight := 1 / .N, by = individual_id]
cov[, age_model := log10(pmax(as.numeric(age) + 1, 1e-6))]

apa <- fread(file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl", "apa_after_combat.tsv"), showProgress = FALSE)
setnames(apa, 1L, "site_id")
samples <- intersect(cov$sample, names(apa))
stopifnot(length(samples) == 697L)
cov <- cov[match(samples, sample)]
apa <- apa[match(sites$site_id, site_id)]
stopifnot(nrow(apa) == nrow(sites), !anyNA(apa$site_id),
          identical(apa$site_id, sites$site_id))
mat <- as.matrix(apa[, ..samples])
storage.mode(mat) <- "numeric"
rm(apa)
model_df <- as.data.frame(cov[, c("age_model", "sex_binary", "RIN_z",
                                  "profile_weight", ilr_names), with = FALSE])
support <- fread(file.path(out, "01_support.tsv"))
grid <- data.frame(age_model = support$age_model,
                   sex_binary = weighted.mean(cov$sex_binary, cov$profile_weight),
                   RIN_z = 0)
for (name in ilr_names) {
  grid[[name]] <- weighted.mean(cov[[name]], cov$profile_weight)
}
supported <- which(support$is_supported == TRUE)
stopifnot(length(supported) >= 2L)

# --- cell-adjusted age models ---
covariates <- paste(c("sex_binary", "RIN_z", ilr_names), collapse = " + ")
null_formula <- as.formula(paste("y ~", covariates))
null_frame <- model_df
results <- vector("list", nrow(sites))
for (i in seq_len(nrow(sites))) {
  dat <- null_frame
  dat$y <- mat[i, ]
  if (anyNA(dat$y)) dat$y[is.na(dat$y)] <- mean(dat$y, na.rm = TRUE)
  null_fit <- tryCatch(gam(null_formula, data = dat,
                           weights = profile_weight, method = "REML"),
                       error = function(e) NULL)
  best <- NULL
  if (!is.null(null_fit)) {
    for (k in k_ladder) {
      full_formula <- as.formula(sprintf("y ~ %s + s(age_model, k = %d, bs = 'tp')",
                                         covariates, k))
      fit <- tryCatch(gam(full_formula, data = dat,
                          weights = profile_weight, method = "REML"),
                      error = function(e) NULL)
      if (is.null(fit)) next
      edf <- as.numeric(summary(fit)$s.table[1L, "edf"])
      kc <- tryCatch(k.check(fit, subsample = 5000, n.rep = 400),
                     error = function(e) NULL)
      k_p <- if (is.null(kc)) NA_real_ else as.numeric(kc[1L, "p-value"])
      dev <- tryCatch(anova(null_fit, fit, test = "F"), error = function(e) NULL)
      age_p <- if (is.null(dev)) NA_real_ else as.numeric(dev$`Pr(>F)`[2L])
      if (!is.finite(age_p)) age_p <- as.numeric(summary(fit)$s.table[1L, "p-value"])
      best <- list(fit = fit, k = k, edf = edf, k_p = k_p, age_p = age_p)
      if (is.finite(k_p) && k_p > 0.05 && edf < k - 1) break
    }
  }
  if (is.null(best)) {
    results[[i]] <- data.table(age_p_adjusted = NA_real_, delta_adjusted = NA_real_,
                               final_k_adjusted = NA_integer_)
  } else {
    pred <- as.numeric(predict(best$fit, newdata = grid))
    results[[i]] <- data.table(age_p_adjusted = best$age_p,
                               delta_adjusted = pred[tail(supported, 1L)] - pred[supported[1L]],
                               final_k_adjusted = best$k)
  }
}

# --- t4 event summary ---
t4 <- cbind(sites[, .(site_id, gene_id, gene_symbol, cluster_id, age_p_original = age_pvalue,
                      age_q_original = age_qvalue, delta_original = delta_logit)],
            rbindlist(results))
t4[, age_q_adjusted := p.adjust(age_p_adjusted, "BH"), by = cluster_id]
t4[, age_q_all_directional := p.adjust(age_p_adjusted, "BH")]
t4[, `:=`(retained = age_q_adjusted < 0.05 &
           sign(delta_adjusted) == sign(delta_original),
           n_matched_samples = nrow(cov))]
mito <- suppressWarnings(read_excel(
  paste0(Sys.getenv("APA_HPC_HOME"), "/project/.ref/enrichment/Human.MitoCarta3.0.xls"),
  sheet = "A Human MitoCarta3.0"))
oxphos <- unique(mito$Symbol[grepl("OXPHOS > OXPHOS subunits",
                                    mito$MitoCarta3.0_MitoPathways, fixed = TRUE)])
proteasome <- fread(paste0(Sys.getenv("APA_HPC_HOME"), "/project/.ref/enrichment/HGNC_Proteasome.tsv"))
t4[, `:=`(c2_oxphos = cluster_id == 2L & gene_symbol %in% oxphos,
           c2_proteasome = cluster_id == 2L &
             gene_symbol %in% proteasome$`Approved symbol`)]
pairs <- fread(file.path(out, "T1_shared_isoform_pairs.tsv"),
               select = c("site_id_c1", "site_id_c2"))
t4[, shared_pair_event := site_id %in% c(pairs$site_id_c1, pairs$site_id_c2)]
fwrite(t4, file.path(out, "T4_celltype_adj_events.tsv"), sep = "\t")