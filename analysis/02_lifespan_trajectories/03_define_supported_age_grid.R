# Support Profile calibration from unique-individual counts, LOO and subsampling.
suppressPackageStartupMessages(library(data.table))

# --- paths ---

root <- Sys.getenv("APA_ROOT")
DATA <- file.path(root, "data/02_traj")

# --- settings ---

GRID_SIZE <- 128L
INNER_WINDOWS <- c(0.03, 0.04, 0.05)
OUTER_WINDOWS <- c(0.06, 0.08, 0.10)
INNER_MINIMUMS <- 2:4
OUTER_MINIMUMS <- 4:6
BASELINE <- list(inner_window = 0.04, outer_window = 0.08, inner_min = 3L, outer_min = 5L)
SUBSAMPLE_ITER <- 500L
SUBSAMPLE_FRACTION <- 0.90
SEED <- 20260715L
LOO_MIN_JACCARD <- 0.95
SUBSAMPLE_Q05_MIN_JACCARD <- 0.95
MIN_POINT_RETENTION <- 0.90

# --- helper functions ---

jaccard <- function(a, b) {
  union_n <- sum(a | b)
  if (union_n == 0L) return(1)
  sum(a & b) / union_n
}

# --- individual ages ---

cov <- fread(file.path(DATA, "00_cov.tsv"))
cov[, age_model := log10(pmax(as.numeric(age) + 1, 1e-6))]
individuals <- cov[, .(age = unique(age), age_model = unique(age_model)), by = individual_id]
stopifnot(lengths(individuals$age) == 1L, lengths(individuals$age_model) == 1L)
individuals[, `:=`(age = unlist(age), age_model = unlist(age_model))]

grid_model <- seq(min(individuals$age_model), max(individuals$age_model), length.out = GRID_SIZE)
grid_age <- 10^grid_model - 1
age_range <- diff(range(individuals$age_model))

# --- local individual-count matrices ---

all_windows <- sort(unique(c(INNER_WINDOWS, OUTER_WINDOWS)))
count_membership <- setNames(lapply(all_windows, function(window_fraction) {
  abs(outer(grid_model, individuals$age_model, "-")) <= window_fraction * age_range
}), sprintf("%.2f", all_windows))

set.seed(SEED)
n_subsample <- round(nrow(individuals) * SUBSAMPLE_FRACTION)
subsample_membership <- replicate(SUBSAMPLE_ITER, {
  selected <- sample.int(nrow(individuals), n_subsample, replace = FALSE)
  seq_len(nrow(individuals)) %in% selected
})
storage.mode(subsample_membership) <- "integer"

# --- predeclared support candidates ---

candidates <- CJ(
  inner_window = INNER_WINDOWS,
  outer_window = OUTER_WINDOWS,
  inner_min = INNER_MINIMUMS,
  outer_min = OUTER_MINIMUMS
)
candidates[, candidate_id := sprintf("I%02dN%d_O%02dN%d",
                                     round(inner_window * 100), inner_min,
                                     round(outer_window * 100), outer_min)]
candidates[, is_baseline :=
  inner_window == BASELINE$inner_window & outer_window == BASELINE$outer_window &
  inner_min == BASELINE$inner_min & outer_min == BASELINE$outer_min]
candidates[, parameter_distance :=
  abs(match(inner_window, INNER_WINDOWS) - match(BASELINE$inner_window, INNER_WINDOWS)) +
  abs(match(outer_window, OUTER_WINDOWS) - match(BASELINE$outer_window, OUTER_WINDOWS)) +
  abs(inner_min - BASELINE$inner_min) + abs(outer_min - BASELINE$outer_min)]

calibration_rows <- vector("list", nrow(candidates))
selection_rows <- vector("list", nrow(candidates))
mask_signatures <- character(nrow(candidates))

for (i in seq_len(nrow(candidates))) {
  x <- candidates[i]
  inner_members <- count_membership[[sprintf("%.2f", x$inner_window)]]
  outer_members <- count_membership[[sprintf("%.2f", x$outer_window)]]
  inner_count <- rowSums(inner_members)
  outer_count <- rowSums(outer_members)
  full_mask <- inner_count >= x$inner_min & outer_count >= x$outer_min

  inner_sub <- inner_members %*% subsample_membership
  outer_sub <- outer_members %*% subsample_membership
  subsample_masks <- inner_sub >= x$inner_min & outer_sub >= x$outer_min
  selection_probability <- rowMeans(subsample_masks)
  subsample_jaccard <- apply(subsample_masks, 2L, jaccard, b = full_mask)

  loo_jaccard <- vapply(seq_len(nrow(individuals)), function(j) {
    loo_mask <- (inner_count - inner_members[, j]) >= x$inner_min &
      (outer_count - outer_members[, j]) >= x$outer_min
    jaccard(loo_mask, full_mask)
  }, numeric(1))

  block <- rleid(full_mask)
  segments <- data.table(grid_age, mask = full_mask, block)[mask == TRUE, .(
    start_age = min(grid_age), end_age = max(grid_age), n_grid = .N
  ), by = block]
  full_retention <- selection_probability[full_mask]
  mask_signatures[i] <- paste(as.integer(full_mask), collapse = "")

  calibration_rows[[i]] <- data.table(
    candidate_id = x$candidate_id,
    inner_window = x$inner_window,
    outer_window = x$outer_window,
    inner_min = x$inner_min,
    outer_min = x$outer_min,
    is_baseline = x$is_baseline,
    parameter_distance = x$parameter_distance,
    n_supported = sum(full_mask),
    n_components = nrow(segments),
    n_singleton_components = sum(segments$n_grid == 1L),
    n_short_components_le2 = sum(segments$n_grid <= 2L),
    segment_boundaries = paste(sprintf("%.3f..%.3f", segments$start_age, segments$end_age),
                               collapse = ";"),
    loo_min_jaccard = min(loo_jaccard),
    loo_median_jaccard = median(loo_jaccard),
    loo_n_changed = sum(loo_jaccard < 1),
    subsample_min_jaccard = min(subsample_jaccard),
    subsample_q05_jaccard = as.numeric(quantile(subsample_jaccard, 0.05)),
    subsample_median_jaccard = median(subsample_jaccard),
    min_supported_point_retention = if (length(full_retention)) min(full_retention) else NA_real_
  )

  selection_rows[[i]] <- data.table(
    candidate_id = x$candidate_id,
    grid_index = seq_along(grid_model) - 1L,
    age_model = grid_model,
    grid_age = grid_age,
    inner_count = inner_count,
    outer_count = outer_count,
    full_supported = full_mask,
    subsample_selection_probability = selection_probability
  )
}

calibration <- rbindlist(calibration_rows)
selection <- rbindlist(selection_rows)
calibration[, mask_signature := mask_signatures]
calibration[, mask_group_size := .N, by = mask_signature]
calibration[, design_stable :=
  loo_min_jaccard >= LOO_MIN_JACCARD &
  subsample_q05_jaccard >= SUBSAMPLE_Q05_MIN_JACCARD &
  min_supported_point_retention >= MIN_POINT_RETENTION &
  n_components == 2L & n_singleton_components == 0L]
calibration[, design_plateau := design_stable & mask_group_size >= 2L]
setorder(calibration, parameter_distance, candidate_id)

# --- baseline Support Profile ---

baseline_id <- calibration[is_baseline == TRUE, candidate_id]
baseline <- selection[candidate_id == baseline_id]
setnames(baseline,
         c("inner_count", "outer_count", "full_supported"),
         c("local_n_inner", "local_n", "is_supported"))
baseline[, `:=`(
  inner_window = BASELINE$inner_window,
  outer_window = BASELINE$outer_window,
  inner_min = BASELINE$inner_min,
  outer_min = BASELINE$outer_min
)]
baseline <- baseline[, .(
  grid_index, age_model, grid_age, local_n, local_n_inner, is_supported,
  inner_window, outer_window, inner_min, outer_min, subsample_selection_probability
)]

fwrite(baseline, file.path(DATA, "01_support.tsv"), sep = "\t")
fwrite(calibration, file.path(DATA, "01_support_calibration.tsv"), sep = "\t")
fwrite(selection, file.path(DATA, "01_support_selection_probability.tsv"), sep = "\t")
