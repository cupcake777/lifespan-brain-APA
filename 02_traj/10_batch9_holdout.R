#!/usr/bin/env Rscript
# Batch9 holdout validation of the discovery APA trajectory classes.
# Trajectories and clusters are fit on training samples only; batch9 holdout samples are scored afterwards.

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(parallel)
  library(sva)
  library(softImpute)
  library(flashier)
  library(ebnm)
  library(irlba)
  library(cluster)
  library(bluster)
})

# --- config ---

APA_ROOT <- Sys.getenv("APA_ROOT")
OUT_DIR <- file.path(APA_ROOT, "data", "02_traj")
RAW_ALL_FILE <- Sys.getenv("RAW_ALL_FILE", file.path(APA_ROOT, "raw", "all.apa.count"))
BATCH9_FILE <- Sys.getenv("BATCH9_FILE", file.path(APA_ROOT, "raw", "batch9.apa.count"))
META_FILE <- file.path(APA_ROOT, "data", "04_qtl", "CHB.all.meta")
QC_MATRIX_FILE <- file.path(APA_ROOT, "data", "04_qtl", "apa_after_combat.tsv")
ORIGINAL_CLUSTER_FILE <- file.path(OUT_DIR, "02_sites_clustered.tsv")

OUT_DESIGN <- file.path(OUT_DIR, "10_batch9_holdout_design.tsv")
OUT_EVENTS <- file.path(OUT_DIR, "10_batch9_holdout_events.tsv")
OUT_CLUSTERS <- file.path(OUT_DIR, "10_batch9_holdout_clusters.tsv")
OUT_NULLS <- file.path(OUT_DIR, "10_batch9_holdout_nulls.tsv")
OUT_STABILITY <- file.path(OUT_DIR, "10_batch9_holdout_discovery_stability.tsv")
OUT_REPLICATES <- file.path(OUT_DIR, "10_batch9_replicate_qc.tsv")

PDUI_EPSILON <- 0.001
GRID_SIZE <- 128L
INNER_WINDOW <- 0.04
OUTER_WINDOW <- 0.08
INNER_MIN <- 3L
OUTER_MIN <- 5L
K_LADDER <- c(20L, 30L, 40L, 60L)
K_INDEX_P_THR <- 0.05
EDF_MARGIN <- 1.0
AGE_Q_THR <- 0.05
NULL_QUANTILE <- 0.95
NULL_PERM_ITER <- as.integer(Sys.getenv("APA_HOLDOUT_NULL_PERM", "200"))
NULL_SITES_PER <- as.integer(Sys.getenv("APA_HOLDOUT_NULL_SITES", "500"))
GENE_BOOT_ITER <- as.integer(Sys.getenv("APA_HOLDOUT_GENE_BOOT", "2000"))
AGE_PERM_ITER <- as.integer(Sys.getenv("APA_HOLDOUT_AGE_PERM", "200"))
CLUSTER_BOOT_ITER <- as.integer(Sys.getenv("APA_HOLDOUT_CLUSTER_BOOT", "100"))
MAX_SITES <- as.integer(Sys.getenv("APA_HOLDOUT_MAX_SITES", "0"))
WORKERS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
MIN_CLUSTER_N <- 30L
MIN_TEMPLATE_RANGE <- 0.05
MIN_TEMPLATE_SD <- 0.02
MIN_TEMPLATE_RESIDUAL_SS <- 0.02
MAX_TEMPLATE_VIF <- 10
MIN_VALIDATION_OBS <- 20L
SEED <- 20260717L

set.seed(SEED)

# --- helpers ---

choose2 <- function(x) x * (x - 1) / 2

fit_template_rows <- function(y, template, nuisance) {
  nuisance_qr <- qr(nuisance)
  nuisance <- nuisance[, nuisance_qr$pivot[seq_len(nuisance_qr$rank)], drop = FALSE]
  nuisance_inv <- solve(crossprod(nuisance))
  residualize <- function(x) {
    x - (x %*% nuisance %*% nuisance_inv) %*% t(nuisance)
  }
  y_res <- residualize(y)
  t_res <- residualize(template)
  template_ss <- rowSums(t_res^2)
  beta <- rowSums(t_res * y_res) / template_ss
  residual <- y_res - beta * t_res
  h_nuisance <- rowSums((nuisance %*% nuisance_inv) * nuisance)
  leverage <- sweep(t_res^2, 1L, template_ss, "/") +
    matrix(h_nuisance, nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
  hc3_term <- t_res * residual / pmax(1 - leverage, 1e-8)
  se <- sqrt(rowSums(hc3_term^2) / template_ss^2)
  df <- ncol(y) - ncol(nuisance) - 1L
  y_ss <- rowSums(y_res^2)
  partial_r <- rowSums(t_res * y_res) / sqrt(template_ss * y_ss)
  centered_ss <- rowSums((template - rowMeans(template))^2)
  data.table(
    beta = beta,
    se_hc3 = se,
    ci_low = beta - qt(0.975, df) * se,
    ci_high = beta + qt(0.975, df) * se,
    p_hc3 = 2 * pt(abs(beta / se), df = df, lower.tail = FALSE),
    partial_r = partial_r,
    template_residual_ss = template_ss,
    template_vif = centered_ss / template_ss,
    residual_df = df
  )
}

fit_kmeans <- function(x, k) {
  kmeans(x, centers = k, nstart = 8L, iter.max = 500L, algorithm = "Lloyd")
}

cluster_gene_summary <- function(events, value_columns) {
  events[, lapply(.SD, median, na.rm = TRUE),
         by = .(cluster_id, gene_id), .SDcols = value_columns]
}

bootstrap_cluster_median <- function(values) {
  boot <- replicate(
    GENE_BOOT_ITER,
    median(sample(values, length(values), replace = TRUE), na.rm = TRUE)
  )
  c(median = median(values, na.rm = TRUE),
    low = unname(quantile(boot, 0.025, na.rm = TRUE)),
    high = unname(quantile(boot, 0.975, na.rm = TRUE)))
}

signflip_p <- function(values) {
  observed <- median(values, na.rm = TRUE)
  null <- replicate(
    GENE_BOOT_ITER,
    median(values * sample(c(-1, 1), length(values), replace = TRUE),
           na.rm = TRUE)
  )
  (1 + sum(null >= observed)) / (GENE_BOOT_ITER + 1)
}

# --- fixed sample split ---

meta <- fread(META_FILE)
setnames(meta, 1L, "sample")
meta[, sample := as.character(sample)]

qc_samples <- names(fread(QC_MATRIX_FILE, nrows = 0L))[-1L]
batch9_samples <- names(fread(BATCH9_FILE, nrows = 0L))[-seq_len(4L)]
holdout_samples <- batch9_samples[!grepl("_rep$", batch9_samples)]
bridge_samples <- batch9_samples[grepl("_rep$", batch9_samples)]
holdout_samples <- intersect(holdout_samples, qc_samples)
bridge_samples <- intersect(bridge_samples, qc_samples)
correction_samples <- setdiff(qc_samples, holdout_samples)
discovery_samples <- correction_samples[!grepl("_rep$", correction_samples)]

meta_correction <- meta[match(correction_samples, sample)]
meta_discovery <- meta[match(discovery_samples, sample)]
meta_holdout <- meta[match(holdout_samples, sample)]

meta_discovery[, `:=`(
  age_model = log10(pmax(as.numeric(age) + 1, 1e-6)),
  sex_binary = as.integer(tolower(sex) == "male"),
  RIN_z = as.numeric(scale(RIN))
)]
meta_holdout[, `:=`(
  age_model = log10(pmax(as.numeric(age) + 1, 1e-6)),
  sex_binary = as.integer(tolower(sex) == "male"),
  RIN_z = as.numeric(scale(RIN))
)]

analysis_scope <- if (MAX_SITES > 0L) "smoke_subset" else "full"

# --- training-only preprocessing ---

raw_header <- names(fread(RAW_ALL_FILE, nrows = 0L))
raw <- fread(
  RAW_ALL_FILE,
  select = c(raw_header[seq_len(4L)], correction_samples),
  nrows = if (MAX_SITES > 0L) MAX_SITES else -1L,
  showProgress = FALSE
)
parts <- tstrsplit(raw$Gene, "|", fixed = TRUE)
raw[, site_id := paste(parts[[1]], parts[[2]], Loci, parts[[4]], sep = "|")]
autosomal <- !parts[[3]] %chin% c("chrX", "chrY", "chrM")
n_raw_loci_loaded <- nrow(raw)

pdui <- as.matrix(raw[, ..correction_samples])
storage.mode(pdui) <- "numeric"
rownames(pdui) <- raw$site_id
row_na <- rowMeans(is.na(pdui))
row_sd <- apply(pdui, 1L, sd, na.rm = TRUE)
keep <- row_na <= 0.5 & is.finite(row_sd) & row_sd >= 0.001 & autosomal
pdui <- pdui[keep, , drop = FALSE]
site_ids <- rownames(pdui)
n_training_qc_loci <- nrow(pdui)

message(sprintf(
  "training preprocessing: %d loci x %d profiles; holdout=%d, bridge=%d",
  nrow(pdui), ncol(pdui), length(holdout_samples), length(bridge_samples)
))

lambda <- if (nrow(pdui) < 1500L) 5 else 30
soft_fit <- softImpute(
  as(pdui, "Incomplete"), rank.max = min(50L, min(dim(pdui)) - 1L),
  lambda = lambda, type = "svd"
)
soft_complete <- complete(pdui, soft_fit)
min_sd <- min(apply(soft_complete, 1L, sd))
pca_init <- irlba(soft_complete, nv = min(10L, min(dim(soft_complete)) - 1L))
flash_fit <- flash_init(pdui, S = min_sd, var_type = 1) |>
  flash_factors_init(
    list(d = pca_init$d, u = pca_init$u, v = pca_init$v),
    ebnm_fn = ebnm_point_laplace
  ) |>
  flash_backfit(maxiter = 300)

pdui_imputed <- pdui
missing <- is.na(pdui_imputed)
pdui_imputed[missing] <- fitted(flash_fit)[missing]
remaining <- is.na(pdui_imputed)
if (any(remaining)) {
  means <- rowMeans(pdui_imputed, na.rm = TRUE)
  pdui_imputed[remaining] <- means[row(pdui_imputed)[remaining]]
}
logit <- log2(
  pmin(pmax(pdui_imputed, PDUI_EPSILON), 1 - PDUI_EPSILON) /
    (1 - pmin(pmax(pdui_imputed, PDUI_EPSILON), 1 - PDUI_EPSILON))
)

combat <- suppressWarnings(ComBat(
  dat = logit,
  batch = factor(meta_correction$batch),
  mod = model.matrix(~ LifeStage, data = meta_correction),
  par.prior = TRUE,
  ref.batch = "batch9"
))
discovery <- combat[, match(discovery_samples, correction_samples), drop = FALSE]
rm(raw, pdui, soft_complete, flash_fit, pdui_imputed, logit, combat)
gc()

# --- training-only support profile ---

grid_model <- seq(
  min(meta_discovery$age_model), max(meta_discovery$age_model),
  length.out = GRID_SIZE
)
grid_age <- 10^grid_model - 1
age_span <- diff(range(meta_discovery$age_model))
inner_count <- rowSums(
  abs(outer(grid_model, meta_discovery$age_model, "-")) <=
    INNER_WINDOW * age_span
)
outer_count <- rowSums(
  abs(outer(grid_model, meta_discovery$age_model, "-")) <=
    OUTER_WINDOW * age_span
)
supported <- inner_count >= INNER_MIN & outer_count >= OUTER_MIN
support <- data.table(
  grid_index = seq_len(GRID_SIZE) - 1L,
  age_model = grid_model,
  grid_age = grid_age,
  local_n_inner = inner_count,
  local_n = outer_count,
  is_supported = supported
)
support[, segment_id := rleid(is_supported)]
support[is_supported == FALSE, segment_id := NA_integer_]
support[is_supported == TRUE, segment_id := match(segment_id, unique(segment_id))]
supported_index <- which(supported)

# --- training-only continuous-age trajectories ---

model_df <- as.data.frame(meta_discovery[, .(
  age_model, sex_binary, RIN_z
)])
grid_df <- data.frame(
  age_model = grid_model,
  sex_binary = mean(model_df$sex_binary),
  RIN_z = 0
)

n_sites <- nrow(discovery)
stats_rows <- vector("list", n_sites)
trajectory <- matrix(NA_real_, n_sites, GRID_SIZE)
trajectory_z <- matrix(NA_real_, n_sites, GRID_SIZE)
chunk_size <- max(100L, WORKERS * 250L)
for (start in seq.int(1L, n_sites, by = chunk_size)) {
  index <- start:min(n_sites, start + chunk_size - 1L)
  result <- mclapply(index, function(i) {
    dat <- model_df
    dat$y <- as.numeric(discovery[i, ])
    null_fit <- tryCatch(
      gam(y ~ sex_binary + RIN_z, data = dat, method = "REML"),
      error = function(e) NULL
    )
    if (is.null(null_fit)) return(NULL)
    best <- NULL
    for (k in K_LADDER) {
      fit <- tryCatch(
        gam(
          y ~ sex_binary + RIN_z + s(age_model, k = k, bs = "tp"),
          data = dat, method = "REML"
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      s_tab <- summary(fit)$s.table
      kc <- tryCatch(
        k.check(fit, subsample = 5000, n.rep = 400),
        error = function(e) matrix(
          NA_real_, 1L, 4L,
          dimnames = list(NULL, c("k'", "edf", "k-index", "p-value"))
        )
      )
      dev <- tryCatch(anova(null_fit, fit, test = "F"), error = function(e) NULL)
      age_p <- if (!is.null(dev)) as.numeric(dev$`Pr(>F)`[2]) else
        as.numeric(s_tab[1, "p-value"])
      best <- list(
        fit = fit,
        k = k,
        edf = as.numeric(s_tab[1, "edf"]),
        k_p = as.numeric(kc[1, "p-value"]),
        k_index = as.numeric(kc[1, "k-index"]),
        age_p = age_p
      )
      k_ok <- is.finite(best$k_p) && best$k_p > K_INDEX_P_THR &&
        !(is.finite(best$edf) && best$edf >= k - EDF_MARGIN)
      if (k_ok) break
    }
    if (is.null(best)) return(NULL)
    pred <- as.numeric(predict(best$fit, newdata = grid_df))
    y_mean <- mean(dat$y)
    y_sd <- sd(dat$y)
    dat$y <- (dat$y - y_mean) / y_sd
    fit_z <- tryCatch(
      gam(
        y ~ sex_binary + RIN_z + s(age_model, k = best$k, bs = "tp"),
        data = dat, method = "REML"
      ),
      error = function(e) NULL
    )
    pred_z <- if (is.null(fit_z)) rep(NA_real_, GRID_SIZE) else
      as.numeric(predict(fit_z, newdata = grid_df))
    list(
      age_p = best$age_p,
      edf = best$edf,
      k_p = best$k_p,
      k_index = best$k_index,
      final_k = best$k,
      pred = pred,
      pred_z = pred_z
    )
  }, mc.cores = WORKERS, mc.preschedule = TRUE)
  for (j in seq_along(index)) {
    i <- index[j]
    x <- result[[j]]
    if (is.null(x)) {
      stats_rows[[i]] <- data.table(
        site_id = site_ids[i], age_pvalue = NA_real_, smooth_edf = NA_real_,
        k_index_pvalue = NA_real_, k_index = NA_real_, final_k = NA_integer_
      )
      next
    }
    trajectory[i, ] <- x$pred
    trajectory_z[i, ] <- x$pred_z
    stats_rows[[i]] <- data.table(
      site_id = site_ids[i],
      age_pvalue = x$age_p,
      smooth_edf = x$edf,
      k_index_pvalue = x$k_p,
      k_index = x$k_index,
      final_k = x$final_k
    )
  }
  message(sprintf("trajectory fits: %d / %d", max(index), n_sites))
}

site_stats <- rbindlist(stats_rows)
site_stats[, age_qvalue := p.adjust(age_pvalue, method = "BH")]
site_stats[, association_pass := age_qvalue < AGE_Q_THR]
rownames(trajectory) <- rownames(trajectory_z) <- site_ids

# --- training-only dynamic null ---

eligible <- which(site_stats$association_pass)
null_sites_per <- min(NULL_SITES_PER, length(eligible))
null_formula <- y ~ sex_binary + RIN_z + s(age_model, k = 20, bs = "tp")

null_data <- rbindlist(mclapply(
  seq_len(NULL_PERM_ITER),
  function(iteration) {
    set.seed(SEED + iteration)
    dat <- model_df
    dat$age_model <- sample(dat$age_model)
    selected <- sample(eligible, null_sites_per)
    metric <- vapply(selected, function(i) {
      dat$y <- as.numeric(discovery[i, ])
      fit <- gam(null_formula, data = dat, method = "REML")
      pred <- as.numeric(predict(fit, newdata = grid_df))[supported_index]
      c(delta = abs(pred[length(pred)] - pred[1L]),
        amplitude = max(pred) - min(pred))
    }, numeric(2L))
    data.table(
      delta = metric["delta", ], amplitude = metric["amplitude", ]
    )
  },
  mc.cores = WORKERS, mc.preschedule = FALSE
))
delta_threshold <- as.numeric(quantile(null_data$delta, NULL_QUANTILE))
amplitude_threshold <- as.numeric(quantile(null_data$amplitude, NULL_QUANTILE))

supported_traj <- trajectory[, supported_index, drop = FALSE]
site_stats[, `:=`(
  delta_logit = supported_traj[, ncol(supported_traj)] - supported_traj[, 1L],
  amplitude_logit = apply(supported_traj, 1L, function(x) max(x) - min(x))
)]
site_stats[, directional_pass :=
  association_pass & abs(delta_logit) >= delta_threshold]

# --- training-only cluster discovery ---

directional_index <- which(site_stats$directional_pass)
cluster_input <- trajectory_z[directional_index, supported_index, drop = FALSE]
pca <- prcomp(cluster_input, center = TRUE, scale. = FALSE, rank. = 12L)
pc <- pca$x[, seq_len(min(12L, ncol(pca$x))), drop = FALSE]
k_range <- 2L:min(8L, nrow(pc) - 1L)
fits <- setNames(lapply(k_range, function(k) fit_kmeans(pc, k)), k_range)
silhouette_index <- sample(nrow(pc), min(5000L, nrow(pc)))
k_diagnostics <- rbindlist(lapply(k_range, function(k) {
  labels <- fits[[as.character(k)]]$cluster
  data.table(
    k = k,
    sample_silhouette = mean(silhouette(
      labels[silhouette_index], dist(pc[silhouette_index, , drop = FALSE])
    )[, "sil_width"]),
    min_cluster_size = min(tabulate(labels))
  )
}))
k_diagnostics[, size_gate_pass := min_cluster_size >= MIN_CLUSTER_N]
selected_k <- k_diagnostics[size_gate_pass == TRUE][
  order(-sample_silhouette, k), k
][1]
if (!is.finite(selected_k)) {
  selected_k <- k_diagnostics[order(-sample_silhouette, k), k][1]
}

base_labels <- fits[[as.character(selected_k)]]$cluster
median_delta <- tapply(
  site_stats$delta_logit[directional_index], base_labels,
  median, na.rm = TRUE
)
label_rank <- rank(-median_delta, ties.method = "first")
cluster_labels <- as.integer(label_rank[as.character(base_labels)])

cluster_consensus <- numeric(length(cluster_labels))
for (b in seq_len(CLUSTER_BOOT_ITER)) {
  set.seed(SEED + 10000L + b)
  boot_index <- sample(
    seq_len(nrow(pc)), max(round(0.8 * nrow(pc)), selected_k + 1L)
  )
  fit <- fit_kmeans(pc[boot_index, , drop = FALSE], selected_k)
  distance <- sapply(seq_len(selected_k), function(j) {
    rowSums((pc - matrix(
      fit$centers[j, ], nrow = nrow(pc), ncol = ncol(pc), byrow = TRUE
    ))^2)
  })
  boot_labels <- max.col(-distance)
  mapping <- apply(table(boot_labels[boot_index], base_labels[boot_index]),
                   1L, which.max)
  cluster_consensus <- cluster_consensus +
    (mapping[as.character(boot_labels)] == base_labels)
}
cluster_consensus <- cluster_consensus / CLUSTER_BOOT_ITER

parts <- tstrsplit(site_stats$site_id[directional_index], "|", fixed = TRUE)
assignments <- data.table(
  site_id = site_stats$site_id[directional_index],
  transcript_id = parts[[1]],
  gene_id = parts[[2]],
  delta_logit = site_stats$delta_logit[directional_index],
  amplitude_logit = site_stats$amplitude_logit[directional_index],
  age_pvalue = site_stats$age_pvalue[directional_index],
  age_qvalue = site_stats$age_qvalue[directional_index],
  cluster_id = cluster_labels,
  cluster_label = paste0("C", cluster_labels),
  cluster_consensus = cluster_consensus
)

# --- discovery stability against the full-data result ---

original <- fread(
  ORIGINAL_CLUSTER_FILE,
  select = c("site_id", "directional_pass", "cluster_id")
)[directional_pass == TRUE & !is.na(cluster_id)]
common <- merge(
  assignments[, .(site_id, holdout_cluster = cluster_id)],
  original[, .(site_id, original_cluster = cluster_id)],
  by = "site_id"
)

tab <- table(common$holdout_cluster, common$original_cluster)
sum_cells <- sum(choose2(tab))
sum_rows <- sum(choose2(rowSums(tab)))
sum_cols <- sum(choose2(colSums(tab)))
total <- choose2(sum(tab))
expected <- sum_rows * sum_cols / total
denominator <- (sum_rows + sum_cols) / 2 - expected
cluster_ari <- if (denominator == 0) 1 else (sum_cells - expected) / denominator

stability <- data.table(
  analysis_scope = analysis_scope,
  training_selected_k = selected_k,
  n_training_directional = nrow(assignments),
  n_original_directional = nrow(original),
  n_common_directional = nrow(common),
  event_jaccard = nrow(common) /
    length(union(assignments$site_id, original$site_id)),
  cluster_ari = cluster_ari,
  median_cluster_consensus = median(assignments$cluster_consensus)
)

# --- raw batch9 holdout outcomes ---

batch9 <- fread(BATCH9_FILE, showProgress = FALSE)
batch9_parts <- tstrsplit(batch9$Gene, "|", fixed = TRUE)
batch9[, site_id := paste(
  batch9_parts[[1]], batch9_parts[[2]], Loci, batch9_parts[[4]], sep = "|"
)]
batch9 <- batch9[match(assignments$site_id, site_id)]
holdout_pdui <- as.matrix(batch9[, ..holdout_samples])
storage.mode(holdout_pdui) <- "numeric"
n_observed <- rowSums(is.finite(holdout_pdui))
holdout_mean <- rowMeans(holdout_pdui, na.rm = TRUE)
holdout_pdui[!is.finite(holdout_pdui)] <-
  holdout_mean[row(holdout_pdui)[!is.finite(holdout_pdui)]]
holdout_y <- log2(
  pmin(pmax(holdout_pdui, PDUI_EPSILON), 1 - PDUI_EPSILON) /
    (1 - pmin(pmax(holdout_pdui, PDUI_EPSILON), 1 - PDUI_EPSILON))
)

template <- matrix(
  NA_real_, nrow(assignments), nrow(meta_holdout),
  dimnames = list(assignments$site_id, meta_holdout$sample)
)
holdout_segment <- rep(NA_integer_, nrow(meta_holdout))
for (segment in sort(unique(na.omit(support$segment_id)))) {
  grid_i <- which(support$segment_id == segment)
  sample_i <- which(
    meta_holdout$age >= min(support$grid_age[grid_i]) &
      meta_holdout$age <= max(support$grid_age[grid_i])
  )
  holdout_segment[sample_i] <- segment
  # linear interpolation of the training trajectory onto holdout ages
  seg_grid_age <- support$grid_age[grid_i]
  target_age <- meta_holdout$age[sample_i]
  interp_weights <- matrix(0, nrow = length(seg_grid_age), ncol = length(target_age))
  for (j in seq_along(target_age)) {
    upper <- which(seg_grid_age >= target_age[j])[1]
    lower <- tail(which(seg_grid_age <= target_age[j]), 1L)
    if (!length(lower) || is.na(lower)) lower <- 1L
    if (!length(upper) || is.na(upper)) upper <- length(seg_grid_age)
    if (lower == upper) {
      interp_weights[lower, j] <- 1
    } else {
      w <- (target_age[j] - seg_grid_age[lower]) /
        (seg_grid_age[upper] - seg_grid_age[lower])
      interp_weights[lower, j] <- 1 - w
      interp_weights[upper, j] <- w
    }
  }
  template[, sample_i] <- trajectory[
    directional_index, grid_i, drop = FALSE
  ] %*% interp_weights
}

supported_samples <- which(!is.na(holdout_segment))
holdout_y <- holdout_y[, supported_samples, drop = FALSE]
template <- template[, supported_samples, drop = FALSE]
holdout_meta_supported <- meta_holdout[supported_samples]
holdout_segment <- holdout_segment[supported_samples]

global_nuisance <- model.matrix(
  ~ sex_binary + RIN_z, data = holdout_meta_supported
)
global_fit <- fit_template_rows(holdout_y, template, global_nuisance)
setnames(global_fit, names(global_fit), paste0(names(global_fit), "_global"))

segment_fits <- list()
for (segment in sort(unique(holdout_segment))) {
  sample_i <- which(holdout_segment == segment)
  nuisance <- model.matrix(
    ~ sex_binary + RIN_z, data = holdout_meta_supported[sample_i]
  )
  fit <- fit_template_rows(
    holdout_y[, sample_i, drop = FALSE],
    template[, sample_i, drop = FALSE],
    nuisance
  )
  segment_fits[[as.character(segment)]] <- fit
}

# mismatched template: pair each site with a nearby site from the same cluster
wrong_index <- integer(nrow(assignments))
for (cluster_id in sort(unique(assignments$cluster_id))) {
  index <- which(assignments$cluster_id == cluster_id)
  x <- scale(as.matrix(assignments[index, .(delta_logit, amplitude_logit)]))
  x[!is.finite(x)] <- 0
  score <- prcomp(x, center = FALSE, scale. = FALSE, rank. = 1L)$x[, 1L]
  n_groups <- min(100L, max(1L, floor(length(index) / 30L)))
  groups <- ceiling(rank(score, ties.method = "first") /
    length(index) * n_groups)
  for (group in unique(groups)) {
    local <- which(groups == group)
    if (length(local) == 1L) {
      candidate <- setdiff(seq_along(index), local)
      distance <- rowSums((x[candidate, , drop = FALSE] -
        matrix(x[local, ], nrow = length(candidate), ncol = ncol(x),
               byrow = TRUE))^2)
      wrong_index[index[local]] <- index[candidate[which.min(distance)]]
    } else {
      distance <- as.matrix(dist(x[local, , drop = FALSE]))
      diag(distance) <- Inf
      wrong_index[index[local]] <- index[local[max.col(-distance)]]
    }
  }
}
wrong_fit <- fit_template_rows(
  holdout_y, template[wrong_index, ], global_nuisance
)

events <- cbind(assignments, global_fit)
events[, `:=`(
  n_observed = n_observed,
  n_supported_samples = length(supported_samples),
  template_range = apply(template, 1L, function(x) max(x) - min(x)),
  template_sd = apply(template, 1L, sd),
  mismatched_site_id = assignments$site_id[wrong_index],
  mismatched_partial_r = wrong_fit$partial_r,
  same_minus_mismatched = partial_r_global - wrong_fit$partial_r
)]

segment_names <- c("prenatal", "postnatal")
for (segment in sort(unique(holdout_segment))) {
  label <- segment_names[segment]
  fit <- segment_fits[[as.character(segment)]]
  events[, (paste0(c("beta_", "ci_low_", "ci_high_", "p_hc3_"), label)) :=
    list(fit$beta, fit$ci_low, fit$ci_high, fit$p_hc3)]
}

events[, estimable :=
  n_observed >= MIN_VALIDATION_OBS &
    template_range >= MIN_TEMPLATE_RANGE &
    template_sd >= MIN_TEMPLATE_SD &
    template_residual_ss_global >= MIN_TEMPLATE_RESIDUAL_SS &
    template_vif_global <= MAX_TEMPLATE_VIF]
events[, event_class := fcase(
  estimable == FALSE, "not_estimable",
  ci_high_global < 0, "reversal",
  ci_low_global > 0 &
    ((is.finite(ci_high_prenatal) & ci_high_prenatal < 0) |
       (is.finite(ci_high_postnatal) & ci_high_postnatal < 0)),
    "segment_inconsistent",
  ci_low_global > 0 & same_minus_mismatched > 0,
    "phenotype_preserved",
  ci_low_global > 0, "cluster_shape_only",
  default = "non_replication"
)]

# --- gene-block cluster inference ---

estimable_events <- events[estimable == TRUE]
gene <- cluster_gene_summary(
  estimable_events,
  c("beta_global", "beta_prenatal", "beta_postnatal",
    "same_minus_mismatched")
)

cluster_rows <- lapply(sort(unique(events$cluster_id)), function(cluster_id) {
  current_cluster <- cluster_id
  g <- gene[cluster_id == current_cluster]
  global <- bootstrap_cluster_median(g$beta_global)
  prenatal <- bootstrap_cluster_median(g$beta_prenatal)
  postnatal <- bootstrap_cluster_median(g$beta_postnatal)
  mismatch <- bootstrap_cluster_median(g$same_minus_mismatched)
  data.table(
    cluster_id = cluster_id,
    cluster_label = paste0("C", cluster_id),
    n_training_events = events[cluster_id == current_cluster, .N],
    n_estimable_events = estimable_events[cluster_id == current_cluster, .N],
    n_estimable_genes = nrow(g),
    median_beta = global["median"],
    median_beta_ci_low = global["low"],
    median_beta_ci_high = global["high"],
    fraction_gene_beta_positive = mean(g$beta_global > 0),
    prenatal_median_beta = prenatal["median"],
    prenatal_ci_low = prenatal["low"],
    prenatal_ci_high = prenatal["high"],
    postnatal_median_beta = postnatal["median"],
    postnatal_ci_low = postnatal["low"],
    postnatal_ci_high = postnatal["high"],
    same_minus_mismatched = mismatch["median"],
    mismatch_ci_low = mismatch["low"],
    mismatch_ci_high = mismatch["high"],
    beta_signflip_p = signflip_p(g$beta_global),
    mismatched_template_p = signflip_p(g$same_minus_mismatched)
  )
})
clusters <- rbindlist(cluster_rows)
clusters[, beta_signflip_holm_p := p.adjust(beta_signflip_p, method = "holm")]
clusters[, mismatched_template_holm_p :=
  p.adjust(mismatched_template_p, method = "holm")]

# --- age-permutation null ---

strata <- interaction(
  holdout_meta_supported$LifeStage,
  holdout_meta_supported$sex,
  drop = TRUE
)
age_null <- matrix(
  NA_real_, nrow = nrow(clusters), ncol = AGE_PERM_ITER,
  dimnames = list(clusters$cluster_id, NULL)
)
for (b in seq_len(AGE_PERM_ITER)) {
  set.seed(SEED + 20000L + b)
  perm <- seq_along(supported_samples)
  for (stratum in levels(strata)) {
    index <- which(strata == stratum)
    perm[index] <- if (length(index) > 1L) sample(index) else index
  }
  fit <- fit_template_rows(
    holdout_y, template[, perm, drop = FALSE], global_nuisance
  )
  perm_event <- copy(events[, .(cluster_id, gene_id, estimable)])
  perm_event[, beta_global := fit$beta]
  perm_gene <- cluster_gene_summary(perm_event[estimable == TRUE], "beta_global")
  for (i in seq_len(nrow(clusters))) {
    age_null[i, b] <- perm_gene[
      cluster_id == clusters$cluster_id[i], median(beta_global, na.rm = TRUE)
    ]
  }
}

clusters[, age_permutation_p := vapply(seq_len(.N), function(i) {
  (1 + sum(age_null[i, ] >= median_beta[i], na.rm = TRUE)) /
    (AGE_PERM_ITER + 1)
}, numeric(1))]
clusters[, age_permutation_holm_p :=
  p.adjust(age_permutation_p, method = "holm")]

clusters[, result_class := fcase(
  median_beta_ci_high < 0, "reversal",
  median_beta_ci_low > 0 &
    (prenatal_ci_high < 0 | postnatal_ci_high < 0),
    "segment_inconsistent",
  median_beta_ci_low > 0 & beta_signflip_holm_p < 0.05 &
    age_permutation_holm_p < 0.05 &
    mismatched_template_holm_p < 0.05,
    "phenotype_specific_preservation",
  median_beta_ci_low > 0 & beta_signflip_holm_p < 0.05 &
    age_permutation_holm_p < 0.05,
    "cluster_level_preservation_only",
  median_beta_ci_low > 0 & age_permutation_holm_p >= 0.05,
    "stage_displacement_only",
  (prenatal_ci_low > 0 & postnatal_ci_high >= 0) |
    (postnatal_ci_low > 0 & prenatal_ci_high >= 0),
    "segment_specific_preservation",
  default = "non_replication"
)]

nulls <- rbindlist(list(
  clusters[, .(
    cluster_id,
    null_type = "age_permutation",
    observed = median_beta,
    null_median = apply(age_null, 1L, median, na.rm = TRUE),
    pvalue = age_permutation_p,
    holm_p = age_permutation_holm_p
  )],
  clusters[, .(
    cluster_id,
    null_type = "within_cluster_mismatched_template",
    observed = same_minus_mismatched,
    null_median = 0,
    pvalue = mismatched_template_p,
    holm_p = mismatched_template_holm_p
  )]
))

# --- technical replicate reproducibility ---

rep_originals <- sub("_rep$", "", bridge_samples)
rep_columns <- unique(c(rep_originals, bridge_samples))
rep_raw <- fread(
  RAW_ALL_FILE,
  select = c(raw_header[seq_len(4L)], rep_columns),
  showProgress = FALSE
)
rep_parts <- tstrsplit(rep_raw$Gene, "|", fixed = TRUE)
rep_raw[, site_id := paste(
  rep_parts[[1]], rep_parts[[2]], Loci, rep_parts[[4]], sep = "|"
)]
rep_raw <- rep_raw[match(assignments$site_id, site_id)]

replicate_qc <- rbindlist(lapply(seq_along(bridge_samples), function(i) {
  original <- rep_originals[i]
  repeated <- bridge_samples[i]
  x <- rep_raw[[original]]
  y <- rep_raw[[repeated]]
  finite <- is.finite(x) & is.finite(y)
  data.table(
    original_sample = original,
    batch9_replicate = repeated,
    original_batch = meta[sample == original, batch][1],
    life_stage = meta[sample == original, LifeStage][1],
    n_common_events = sum(finite),
    pearson = cor(x[finite], y[finite], method = "pearson"),
    spearman = cor(x[finite], y[finite], method = "spearman"),
    median_abs_pdui_difference = median(abs(x[finite] - y[finite]))
  )
}))

# --- provenance and outputs ---

k_diagnostics[, `:=`(
  analysis_scope = analysis_scope,
  selected = k == selected_k
)]

design <- data.table(
  metric = c(
    "analysis_scope", "raw_loci_loaded", "training_qc_loci",
    "correction_profiles", "discovery_individuals", "bridge_profiles",
    "holdout_individuals", "holdout_prenatal", "holdout_postnatal",
    "holdout_supported_samples", "support_grid_points",
    "support_segments", "association_pass_sites", "directional_sites",
    "selected_k", "delta_threshold", "amplitude_threshold",
    "null_permutations", "null_sites_per_permutation",
    "gene_boot_iterations", "age_permutation_iterations",
    "median_replicate_pearson", "median_replicate_spearman",
    "holdout_outcome_used_for_selection"
  ),
  value = as.character(c(
    analysis_scope, n_raw_loci_loaded, n_training_qc_loci,
    length(correction_samples), length(discovery_samples), length(bridge_samples),
    length(holdout_samples), sum(meta_holdout$LifeStage == "Prenatal"),
    sum(meta_holdout$LifeStage == "Postnatal"), length(supported_samples),
    sum(supported), uniqueN(na.omit(support$segment_id)),
    sum(site_stats$association_pass), nrow(assignments), selected_k,
    delta_threshold, amplitude_threshold, NULL_PERM_ITER, null_sites_per,
    GENE_BOOT_ITER, AGE_PERM_ITER, median(replicate_qc$pearson),
    median(replicate_qc$spearman), FALSE
  ))
)

events[, `:=`(
  analysis_scope = analysis_scope,
  gene_boot_iter = GENE_BOOT_ITER,
  age_permutation_iter = AGE_PERM_ITER
)]
clusters[, `:=`(
  analysis_scope = analysis_scope,
  gene_boot_iter = GENE_BOOT_ITER,
  age_permutation_iter = AGE_PERM_ITER,
  null_training_permutations = NULL_PERM_ITER,
  null_training_sites = null_sites_per
)]
nulls[, `:=`(
  analysis_scope = analysis_scope,
  gene_boot_iter = GENE_BOOT_ITER,
  age_permutation_iter = AGE_PERM_ITER
)]

fwrite(design, OUT_DESIGN, sep = "\t")
fwrite(events, OUT_EVENTS, sep = "\t")
fwrite(clusters, OUT_CLUSTERS, sep = "\t")
fwrite(nulls, OUT_NULLS, sep = "\t")
fwrite(rbindlist(list(
  data.table(record_type = "summary", stability),
  data.table(record_type = "k_diagnostic", k_diagnostics)
), fill = TRUE), OUT_STABILITY, sep = "\t")
fwrite(replicate_qc, OUT_REPLICATES, sep = "\t")

print(clusters[, .(
  cluster_label, n_estimable_events, n_estimable_genes,
  median_beta, median_beta_ci_low, median_beta_ci_high,
  prenatal_median_beta, postnatal_median_beta,
  age_permutation_holm_p, mismatched_template_holm_p, result_class
)])
message("batch9 holdout validation complete")
