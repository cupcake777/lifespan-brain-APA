# Dynamic APA sites from individual-block permutation nulls and Support Mask sensitivity.
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(parallel)
})

# --- paths ---

root <- Sys.getenv("APA_ROOT")
DATA <- file.path(root, "data/02_traj")
APA_MATRIX <- file.path(DATA, "apa_after_combat_donor_mean.tsv")

# --- settings ---

N_PERM_PRIMARY <- as.integer(Sys.getenv("APA_DYNAMIC_N_PERM", "1000"))
N_PERM_SENSITIVITY <- as.integer(Sys.getenv("APA_DYNAMIC_SENS_N_PERM", "250"))
SITES_PER <- as.integer(Sys.getenv("APA_DYNAMIC_SITES_PER", "500"))
MAX_MATRIX_ROWS <- as.integer(Sys.getenv("APA_DYNAMIC_MAX_ROWS", "0"))
MAX_MASKS <- 5L
BASIS_K <- 20L
NULL_QUANTILES <- c(0.90, 0.95, 0.99)
PRIMARY_QUANTILE <- 0.95
PRIMARY_SEED <- 2026L
SENSITIVITY_SEEDS <- c(2027L, 2028L)
WORKERS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# --- Support Masks ---

calibration <- fread(file.path(DATA, "01_support_calibration.tsv"))
prediction <- fread(file.path(DATA, "01_support_prediction.tsv"))
selection <- fread(file.path(DATA, "01_support_selection_probability.tsv"))

mask_candidates <- merge(
  calibration,
  prediction[, .(candidate_id, prediction_support_better, support_to_nonsupport_rmse_ratio)],
  by = "candidate_id", all.x = TRUE
)
mask_candidates <- mask_candidates[is_baseline == TRUE | (design_plateau == TRUE & prediction_support_better == TRUE)]
setorder(mask_candidates, -is_baseline, parameter_distance, candidate_id)
mask_candidates <- mask_candidates[!duplicated(mask_signature)][seq_len(min(.N, MAX_MASKS))]
if (!any(mask_candidates$is_baseline)) {
  mask_candidates <- rbind(calibration[is_baseline == TRUE], mask_candidates, fill = TRUE)[
    seq_len(min(.N, MAX_MASKS))]
}

selected_ids <- mask_candidates$candidate_id
mask_list <- setNames(lapply(selected_ids, function(id) {
  selection[candidate_id == id][order(grid_index)]$full_supported
}), selected_ids)
baseline_id <- mask_candidates[is_baseline == TRUE, candidate_id][1]
fwrite(mask_candidates, file.path(DATA, "01_support_selected_masks.tsv"), sep = "\t")

# --- model inputs ---

cov <- fread(file.path(DATA, "00_cov.tsv"))
cov[, age_model := log10(pmax(as.numeric(age) + 1, 1e-6))]
sites <- fread(file.path(DATA, "02_sites.tsv"))
traj <- fread(file.path(DATA, "traj_logit.csv"), showProgress = FALSE)

grid_cols <- grep("^grid_", names(traj), value = TRUE)
grid_model <- selection[candidate_id == baseline_id][order(grid_index)]$age_model
traj_mat <- as.matrix(traj[, ..grid_cols])
storage.mode(traj_mat) <- "numeric"
rownames(traj_mat) <- traj$site_id
rm(traj)

apa <- if (MAX_MATRIX_ROWS > 0L) {
  fread(APA_MATRIX, nrows = MAX_MATRIX_ROWS, showProgress = FALSE)
} else {
  fread(APA_MATRIX, showProgress = FALSE)
}
setnames(apa, 1L, "site_id")
samples <- intersect(cov$sample, names(apa))
cov <- cov[match(samples, sample)]
mat <- as.matrix(apa[, ..samples])
storage.mode(mat) <- "numeric"
rownames(mat) <- apa$site_id
rm(apa)

# --- observed metrics across masks ---

observed_metrics <- rbindlist(lapply(selected_ids, function(id) {
  supported_index <- which(mask_list[[id]])
  supported_traj <- traj_mat[, supported_index, drop = FALSE]
  data.table(
    candidate_id = id,
    site_id = rownames(traj_mat),
    delta_logit = supported_traj[, ncol(supported_traj)] - supported_traj[, 1L],
    amplitude_logit = apply(supported_traj, 1L, function(x) max(x) - min(x))
  )
}))

# --- individual-block permutation null ---

form <- y ~ sex_binary + RIN_z + s(age_model, k = BASIS_K, bs = "tp")
grid_df <- data.frame(
  sex_binary = weighted.mean(cov$sex_binary, cov$profile_weight),
  RIN_z = 0,
  age_model = grid_model
)
individual_age <- unique(cov[, .(individual_id, age_model)])
eligible_ids <- intersect(sites[association_pass == TRUE, site_id], rownames(mat))

run_seed_group <- function(seed_base, n_perm) {
  run_permutation <- function(random_seed) {
    set.seed(random_seed)
    permuted_age <- sample(individual_age$age_model)
    names(permuted_age) <- individual_age$individual_id
    model_df <- as.data.frame(cov[, .(sex_binary, RIN_z, profile_weight)])
    model_df$age_model <- permuted_age[cov$individual_id]
    selected_sites <- sample(eligible_ids, min(SITES_PER, length(eligible_ids)))
    metrics <- vapply(selected_sites, function(id) {
      dat <- model_df
      dat$y <- as.numeric(mat[id, ])
      if (anyNA(dat$y)) {
        dat$y[is.na(dat$y)] <- weighted.mean(dat$y, dat$profile_weight, na.rm = TRUE)
      }
      fit <- tryCatch(gam(form, data = dat, weights = profile_weight, method = "REML"),
                      error = function(e) NULL)
      if (is.null(fit)) return(rep(NA_real_, length(mask_list) * 2L))
      pred <- as.numeric(predict(fit, newdata = grid_df))
      unlist(lapply(mask_list, function(mask) {
        supported_pred <- pred[mask]
        c(
          delta = abs(supported_pred[length(supported_pred)] - supported_pred[1L]),
          amplitude = max(supported_pred) - min(supported_pred)
        )
      }), use.names = FALSE)
    }, numeric(length(mask_list) * 2L))

    rbindlist(lapply(seq_along(selected_ids), function(j) {
      data.table(
        candidate_id = selected_ids[j],
        delta = metrics[(j - 1L) * 2L + 1L, ],
        amplitude = metrics[(j - 1L) * 2L + 2L, ]
      )
    }))
  }

  random_seeds <- seed_base * 100000L + seq_len(n_perm)
  if (WORKERS > 1L && .Platform$OS.type != "windows") {
    null_list <- mclapply(random_seeds, run_permutation, mc.cores = WORKERS, mc.preschedule = FALSE)
  } else {
    null_list <- lapply(random_seeds, run_permutation)
  }
  null_data <- rbindlist(null_list)[is.finite(delta) & is.finite(amplitude)]
  null_data[, rbindlist(lapply(NULL_QUANTILES, function(q) {
    data.table(
      seed = seed_base,
      n_perm = n_perm,
      sites_per = SITES_PER,
      n_null = .N,
      quantile = q,
      delta_threshold = as.numeric(quantile(delta, q)),
      amplitude_threshold = as.numeric(quantile(amplitude, q))
    )
  })), by = candidate_id]
}

null_summary <- rbindlist(c(
  list(run_seed_group(PRIMARY_SEED, N_PERM_PRIMARY)),
  lapply(SENSITIVITY_SEEDS, run_seed_group, n_perm = N_PERM_SENSITIVITY)
))
fwrite(null_summary, file.path(DATA, "03_null.tsv"), sep = "\t")

# --- primary and neighbouring-mask site flags ---

primary_thresholds <- null_summary[seed == PRIMARY_SEED & quantile == PRIMARY_QUANTILE]
mask_flags <- merge(observed_metrics, primary_thresholds[, .(
  candidate_id, delta_threshold, amplitude_threshold
)], by = "candidate_id")
association_ids <- sites[association_pass == TRUE, site_id]
mask_flags[, `:=`(
  association_pass = site_id %chin% association_ids,
  directional_pass = site_id %chin% association_ids &
    is.finite(delta_logit) & abs(delta_logit) >= delta_threshold,
  amplitude_pass = site_id %chin% association_ids &
    is.finite(amplitude_logit) & amplitude_logit >= amplitude_threshold
)]
mask_flags[, `:=`(
  transient_pass = amplitude_pass & !directional_pass,
  dynamic_any_pass = directional_pass | amplitude_pass
)]
fwrite(mask_flags, file.path(DATA, "02_mask_site_flags.tsv"), sep = "\t")

primary_flags <- mask_flags[candidate_id == baseline_id]
primary_flags <- primary_flags[match(sites$site_id, site_id)]
sites[, `:=`(
  delta_logit = primary_flags$delta_logit,
  amplitude_logit = primary_flags$amplitude_logit,
  dynamic_pass = primary_flags$directional_pass,
  delta_logit_threshold = primary_flags$delta_threshold,
  directional_pass = primary_flags$directional_pass,
  amplitude_pass = primary_flags$amplitude_pass,
  transient_pass = primary_flags$transient_pass,
  dynamic_any_pass = primary_flags$dynamic_any_pass,
  amplitude_logit_threshold = primary_flags$amplitude_threshold
)]
fwrite(sites, file.path(DATA, "02_sites.tsv"), sep = "\t")

# --- seed and Support Mask robustness ---

baseline_directional <- mask_flags[candidate_id == baseline_id, directional_pass]
mask_summary <- rbindlist(lapply(seq_len(nrow(null_summary)), function(i) {
  threshold <- null_summary[i]
  metrics <- observed_metrics[candidate_id == threshold$candidate_id]
  directional <- metrics$site_id %chin% association_ids &
    abs(metrics$delta_logit) >= threshold$delta_threshold
  amplitude <- metrics$site_id %chin% association_ids &
    metrics$amplitude_logit >= threshold$amplitude_threshold
  data.table(
    candidate_id = threshold$candidate_id,
    seed = threshold$seed,
    quantile = threshold$quantile,
    delta_threshold = threshold$delta_threshold,
    amplitude_threshold = threshold$amplitude_threshold,
    n_directional = sum(directional),
    n_amplitude = sum(amplitude),
    n_transient = sum(amplitude & !directional),
    n_union = sum(directional | amplitude),
    directional_jaccard_vs_primary = {
      union_n <- sum(directional | baseline_directional)
      if (union_n == 0L) 1 else sum(directional & baseline_directional) / union_n
    }
  )
}))
fwrite(mask_summary, file.path(DATA, "02_mask_dynamic_sensitivity.tsv"), sep = "\t")

threshold_sensitivity <- mask_summary[candidate_id == baseline_id & seed == PRIMARY_SEED]
fwrite(threshold_sensitivity, file.path(DATA, "02_dynamic_threshold_sensitivity.tsv"), sep = "\t")
