# Held-out donor-block prediction check for the Support Profile.
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(parallel)
})

# --- paths ---

root <- Sys.getenv("APA_ROOT")
DATA <- file.path(root, "data/02_traj")
APA_MATRIX <- file.path(root, "data/01_quant/apa_combat_donor_mean.tsv")
COV_FILE <- file.path(DATA, "00_cov.tsv")
SITES_FILE <- file.path(DATA, "02_sites.tsv")
CALIBRATION_FILE <- file.path(DATA, "01_support_calib.tsv")
SELECTION_FILE <- file.path(DATA, "01_support_selection_probability.tsv")

# --- settings ---

PREDICTION_SITES <- as.integer(Sys.getenv("APA_PREDICTION_SITES", "200"))
MAX_MATRIX_ROWS <- as.integer(Sys.getenv("APA_PREDICTION_MAX_ROWS", "0"))
N_FOLDS <- 5L
BOUNDARY_GRID_DISTANCE <- 2L
WORKERS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
SEED <- 20260715L

# --- helpers ---

weighted_rmse <- function(error, weight) sqrt(sum(weight * error^2) / sum(weight))

# --- inputs ---

cov <- fread(COV_FILE)
cov[, age_model := log10(pmax(as.numeric(age) + 1, 1e-6))]
sites <- fread(SITES_FILE)
calibration <- fread(CALIBRATION_FILE)
selection <- fread(SELECTION_FILE)

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

# --- seeded site and donor folds ---

set.seed(SEED)
sites[, prediction_stratum := fifelse(association_pass == TRUE, "age_associated", "not_age_associated")]
selected_sites <- sites[, .SD[sample(.N, min(.N, ceiling(PREDICTION_SITES / 2)))],
                        by = prediction_stratum]
if (nrow(selected_sites) > PREDICTION_SITES) {
  selected_sites <- selected_sites[sample(.N, PREDICTION_SITES)]
}
selected_sites <- selected_sites[site_id %chin% rownames(mat)]

individual_ids <- unique(cov$individual_id)
fold_id <- sample(rep(seq_len(N_FOLDS), length.out = length(individual_ids)))
names(fold_id) <- individual_ids
cov[, fold := fold_id[individual_id]]

# --- donor-block held-out prediction ---

prediction_rows <- mclapply(seq_len(nrow(selected_sites)), function(i) {
  site <- selected_sites[i]
  y <- as.numeric(mat[site$site_id, ])
  if (anyNA(y)) y[is.na(y)] <- weighted.mean(y, cov$profile_weight, na.rm = TRUE)
  y_mean <- weighted.mean(y, cov$profile_weight)
  y_sd <- sqrt(weighted.mean((y - y_mean)^2, cov$profile_weight))
  k <- max(5L, min(as.integer(site$final_k), 20L, uniqueN(cov$age_model) - 1L))

  rbindlist(lapply(seq_len(N_FOLDS), function(fold) {
    train <- cov$fold != fold
    test <- !train
    dat <- as.data.frame(cov[train, .(age_model, sex_binary, RIN_z, profile_weight)])
    dat$y <- y[train]
    fit <- tryCatch(
      gam(y ~ sex_binary + RIN_z + s(age_model, k = k, bs = "tp"),
          data = dat, weights = profile_weight, method = "REML"),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    newdata <- as.data.frame(cov[test, .(age_model, sex_binary, RIN_z)])
    pred <- as.numeric(predict(fit, newdata = newdata))
    data.table(
      site_id = site$site_id,
      prediction_stratum = site$prediction_stratum,
      sample = cov$sample[test],
      individual_id = cov$individual_id[test],
      profile_weight = cov$profile_weight[test],
      age_model = cov$age_model[test],
      fold = fold,
      standardized_error = (y[test] - pred) / y_sd
    )
  }), fill = TRUE)
}, mc.cores = WORKERS, mc.preschedule = TRUE)

errors <- rbindlist(prediction_rows, fill = TRUE)

# --- candidate-specific prediction reliability ---

plateau_ids <- calibration[design_plateau == TRUE, candidate_id]
selection <- selection[candidate_id %chin% plateau_ids]

summary_rows <- lapply(plateau_ids, function(id) {
  candidate_grid <- selection[candidate_id == id][order(grid_index)]
  mask <- candidate_grid$full_supported
  sample_grid <- vapply(cov$age_model, function(age_value) {
    which.min(abs(candidate_grid$age_model - age_value))
  }, integer(1))
  supported_distance <- vapply(sample_grid, function(grid_i) {
    min(abs(which(mask) - grid_i))
  }, integer(1))
  sample_class <- fifelse(
    mask[sample_grid], "supported",
    fifelse(supported_distance <= BOUNDARY_GRID_DISTANCE, "boundary", "unsupported")
  )
  support_class <- sample_class[match(errors$sample, cov$sample)]

  candidate_errors <- copy(errors)
  candidate_errors[, `:=`(candidate_id = id, support_class = support_class)]
  class_summary <- candidate_errors[, .(
    n_profiles = .N,
    n_individuals = uniqueN(individual_id),
    weighted_rmse = weighted_rmse(standardized_error, profile_weight),
    weighted_mae = sum(profile_weight * abs(standardized_error)) / sum(profile_weight)
  ), by = support_class]

  non_supported <- candidate_errors[support_class != "supported"]
  supported <- candidate_errors[support_class == "supported"]
  supported_rmse <- weighted_rmse(supported$standardized_error, supported$profile_weight)
  nonsupported_rmse <- weighted_rmse(non_supported$standardized_error, non_supported$profile_weight)
  boundary_rmse <- class_summary[support_class == "boundary", weighted_rmse]
  unsupported_rmse <- class_summary[support_class == "unsupported", weighted_rmse]
  if (length(boundary_rmse) == 0L) boundary_rmse <- NA_real_
  if (length(unsupported_rmse) == 0L) unsupported_rmse <- NA_real_
  data.table(
    candidate_id = id,
    n_sites = uniqueN(candidate_errors$site_id),
    supported_rmse = supported_rmse,
    boundary_rmse = boundary_rmse,
    unsupported_rmse = unsupported_rmse,
    nonsupported_rmse = nonsupported_rmse,
    support_to_nonsupport_rmse_ratio = supported_rmse / nonsupported_rmse,
    prediction_support_better = supported_rmse < nonsupported_rmse,
    class_metrics = paste(
      sprintf("%s:n=%d,rmse=%.4f", class_summary$support_class,
              class_summary$n_individuals, class_summary$weighted_rmse),
      collapse = ";"
    )
  )
})

prediction_summary <- rbindlist(summary_rows, fill = TRUE)
prediction_summary <- merge(
  calibration[, .(candidate_id, is_baseline, parameter_distance,
                  n_supported, segment_boundaries, mask_group_size, design_plateau)],
  prediction_summary,
  by = "candidate_id"
)
setorder(prediction_summary, parameter_distance, candidate_id)

fwrite(errors, file.path(DATA, "01_support_prediction_errors.tsv"), sep = "\t")
fwrite(prediction_summary, file.path(DATA, "01_support_prediction.tsv"), sep = "\t")
