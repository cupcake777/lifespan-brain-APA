#!/usr/bin/env Rscript
# External C1/C2 pattern replication in BrainSpan: age association and matched-window direction.
# Discovery C1/C2 membership is frozen; no k is selected in the external data.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(splines)
})

# --- config ---

APA_ROOT <- Sys.getenv("APA_ROOT")
OUT_DIR <- file.path(APA_ROOT, "data", "02_traj")
BRAINSPAN_EVENT_FILE <- file.path(OUT_DIR, "08_brainspan_event_replication.tsv")
BRAINSPAN_MATRIX_FILE <- file.path(APA_ROOT, "raw", "all.apa.brainspan.count")
BRAINSPAN_META_FILE <- file.path(APA_ROOT, "raw", "all.apa.brainspan.meta")
SUPPORT_FILE <- file.path(OUT_DIR, "01_support.tsv")
DISCOVERY_TRAJECTORY_FILE <- file.path(OUT_DIR, "traj_logit.csv")

OUT_AGE <- file.path(OUT_DIR, "12_external_age_association.tsv")
OUT_EVENT <- file.path(OUT_DIR, "12_external_c1c2_event_replication.tsv")
OUT_SUMMARY <- file.path(OUT_DIR, "12_external_c1c2_trend_sum.tsv")
OUT_PLOT <- file.path(OUT_DIR, "12_external_c1c2_traj_plot.tsv")

AGE_FDR <- 0.05
BOOT_ITER <- 5000L
BOOT_SEED <- 20260720L
BRAINSPAN_SMOOTH_DF <- 2L
PDUI_EPSILON <- 0.001

# --- helpers ---

parse_field <- function(ids, pos) {
  vapply(strsplit(as.character(ids), "|", fixed = TRUE), function(x) {
    if (length(x) >= pos) x[pos] else ""
  }, character(1))
}

standardize_rows <- function(x) {
  centered <- x - rowMeans(x)
  row_sd <- sqrt(rowSums(centered^2) / (ncol(x) - 1L))
  keep <- is.finite(row_sd) & row_sd > 1e-8
  list(matrix = centered[keep, , drop = FALSE] / row_sd[keep], keep = keep)
}

# --- supported external samples ---

support_all <- fread(SUPPORT_FILE)[order(grid_index)]
support <- support_all[is_supported == TRUE]
support[, segment_id := cumsum(c(TRUE, diff(grid_index) != 1L))]

discovery_trajectory <- fread(DISCOVERY_TRAJECTORY_FILE)
discovery_site_id <- discovery_trajectory$site_id
discovery_trajectory[, site_id := NULL]
discovery_trajectory <- as.matrix(discovery_trajectory)

brainspan_meta <- fread(BRAINSPAN_META_FILE)
brainspan_meta[, `:=`(age = as.numeric(age), segment_id = NA_integer_)]
for (segment in unique(support$segment_id)) {
  age_range <- range(support[segment_id == segment, grid_age])
  brainspan_meta[age >= age_range[1] & age <= age_range[2], segment_id := segment]
}
brainspan_meta <- brainspan_meta[!is.na(segment_id)]
brainspan_meta[, sex_binary := as.integer(sex == "Male")]

# --- BrainSpan age-associated C1/C2 replication ---

brainspan_events <- fread(BRAINSPAN_EVENT_FILE)[mapped_qc == TRUE]
brainspan_matrix <- fread(BRAINSPAN_MATRIX_FILE, showProgress = FALSE)
raw_sample_columns <- setdiff(
  names(brainspan_matrix),
  c("Gene", "fit_value", "Predicted_Proximal_APA", "Loci")
)
brainspan_matrix[, brainspan_site_id := paste(
  parse_field(Gene, 1L), parse_field(Gene, 2L), Loci,
  parse_field(Gene, 4L), sep = "|"
)]
for (sample_name in raw_sample_columns) {
  pdui <- as.numeric(brainspan_matrix[[sample_name]])
  pdui <- pmin(pmax(pdui, PDUI_EPSILON), 1 - PDUI_EPSILON)
  set(brainspan_matrix, j = sample_name,
      value = log2(pdui / (1 - pdui)))
}
brainspan_matrix[, c("Gene", "fit_value", "Predicted_Proximal_APA", "Loci") := NULL]
brainspan_matrix <- brainspan_matrix[
  match(brainspan_events$brainspan_site_id, brainspan_site_id)
]
brainspan_samples <- brainspan_meta$sample
brainspan_y <- as.matrix(brainspan_matrix[, ..brainspan_samples])
storage.mode(brainspan_y) <- "numeric"
brainspan_y[!is.finite(brainspan_y)] <- rowMeans(
  brainspan_y, na.rm = TRUE
)[row(brainspan_y)[!is.finite(brainspan_y)]]

brainspan_reduced <- model.matrix(~ sex_binary, data = brainspan_meta)
brainspan_age_terms <- model.matrix(
  ~ factor(segment_id) - 1, data = brainspan_meta
)
for (segment in sort(unique(brainspan_meta$segment_id))) {
  sample_index <- which(brainspan_meta$segment_id == segment)
  basis <- ns(brainspan_meta$age[sample_index], df = BRAINSPAN_SMOOTH_DF)
  block <- matrix(0, nrow(brainspan_meta), ncol(basis))
  block[sample_index, ] <- basis
  colnames(block) <- paste0("segment", segment, "_age", seq_len(ncol(basis)))
  brainspan_age_terms <- cbind(brainspan_age_terms, block)
}

# age-term test: add independent age columns to the reduced model, then compare RSS
residualized <- qr.resid(qr(brainspan_reduced), brainspan_age_terms)
qr_terms <- qr(residualized)
keep <- sort(qr_terms$pivot[seq_len(qr_terms$rank)])
terms <- residualized[, keep, drop = FALSE]
colnames(terms) <- colnames(brainspan_age_terms)[keep]
brainspan_full <- cbind(brainspan_reduced, terms)

age_columns <- setdiff(colnames(brainspan_full), colnames(brainspan_reduced))
fit <- eBayes(lmFit(brainspan_y, brainspan_full))
test <- topTable(
  fit, coef = age_columns, number = Inf, sort.by = "none"
)

reduced_cross <- brainspan_y %*% brainspan_reduced
full_cross <- brainspan_y %*% brainspan_full
rss_reduced <- rowSums(brainspan_y^2) - rowSums(
  (reduced_cross %*% solve(crossprod(brainspan_reduced))) * reduced_cross
)
rss_full <- rowSums(brainspan_y^2) - rowSums(
  (full_cross %*% solve(crossprod(brainspan_full))) * full_cross
)

brainspan_age <- data.table(
  age_df = length(age_columns),
  residual_df = ncol(brainspan_y) - ncol(brainspan_full),
  age_f = test$F,
  age_p = test$P.Value,
  age_q = p.adjust(test$P.Value, method = "BH"),
  age_partial_r2 = pmax(rss_reduced - rss_full, 0) / rss_reduced
)
brainspan_age[, `:=`(
  cohort = "BrainSpan",
  site_id = brainspan_events$site_id,
  gene_id = brainspan_events$gene_id,
  discovery_cluster_id = brainspan_events$cluster_id,
  n_samples = ncol(brainspan_y),
  age_selected = age_q <= AGE_FDR
)]
brainspan_selected <- brainspan_age$age_selected
brainspan_events <- brainspan_events[brainspan_selected]
brainspan_y <- brainspan_y[brainspan_selected, , drop = FALSE]
brainspan_events[, `:=`(
  external_age_f = brainspan_age[brainspan_selected]$age_f,
  external_age_p = brainspan_age[brainspan_selected]$age_p,
  external_age_q = brainspan_age[brainspan_selected]$age_q,
  external_age_partial_r2 = brainspan_age[brainspan_selected]$age_partial_r2
)]

brainspan_grid <- unlist(lapply(
  sort(unique(brainspan_meta$segment_id)),
  function(segment) {
    age <- brainspan_meta[segment_id == segment, age]
    seq(min(age), max(age), length.out = 30L)
  }
))
brainspan_trajectory <- matrix(
  NA_real_, nrow(brainspan_y), length(brainspan_grid)
)
grid_start <- 1L
for (segment in sort(unique(brainspan_meta$segment_id))) {
  sample_index <- which(brainspan_meta$segment_id == segment)
  prediction_age <- brainspan_grid[grid_start:(grid_start + 29L)]
  seg_y <- brainspan_y[, sample_index, drop = FALSE]
  seg_age <- brainspan_meta$age[sample_index]
  seg_sex <- brainspan_meta$sex_binary[sample_index]
  basis <- ns(seg_age, df = BRAINSPAN_SMOOTH_DF)
  design <- cbind(intercept = 1, sex = seg_sex, basis)
  qr_design <- qr(design)
  design <- design[, sort(qr_design$pivot[seq_len(qr_design$rank)]), drop = FALSE]
  prediction_basis <- predict(basis, prediction_age)
  prediction <- cbind(intercept = 1, sex = mean(seg_sex), prediction_basis)
  prediction <- prediction[, colnames(design), drop = FALSE]
  coefficients <- seg_y %*% design %*% solve(crossprod(design))
  brainspan_trajectory[, grid_start:(grid_start + 29L)] <- coefficients %*% t(prediction)
  grid_start <- grid_start + 30L
}

# interpolate the discovery trajectory onto the BrainSpan grid
interp_weights <- matrix(0, nrow = length(support_all$grid_age), ncol = length(brainspan_grid))
for (j in seq_along(brainspan_grid)) {
  upper <- which(support_all$grid_age >= brainspan_grid[j])[1]
  lower <- tail(which(support_all$grid_age <= brainspan_grid[j]), 1)
  if (lower == upper) {
    interp_weights[lower, j] <- 1
  } else {
    upper_weight <- (brainspan_grid[j] - support_all$grid_age[lower]) /
      (support_all$grid_age[upper] - support_all$grid_age[lower])
    interp_weights[lower, j] <- 1 - upper_weight
    interp_weights[upper, j] <- upper_weight
  }
}
brainspan_discovery <- discovery_trajectory[
  match(brainspan_events$site_id, discovery_site_id), , drop = FALSE
] %*% interp_weights

# --- C1/C2 event and trend summary ---

cohort <- "BrainSpan"
events <- brainspan_events
external <- brainspan_trajectory
discovery <- brainspan_discovery
age <- brainspan_grid
segment_id <- rep(sort(unique(brainspan_meta$segment_id)), each = 30L)

external_scaled <- standardize_rows(external)
discovery_scaled <- standardize_rows(discovery)
keep <- external_scaled$keep & discovery_scaled$keep
events <- events[keep]
external_z <- standardize_rows(external[keep, , drop = FALSE])$matrix
discovery_z <- standardize_rows(discovery[keep, , drop = FALSE])$matrix

event_result <- events[, .(
  cohort,
  site_id,
  gene_id,
  discovery_cluster_id = cluster_id,
  external_age_f,
  external_age_p,
  external_age_q,
  external_age_partial_r2
)]
event_result[, `:=`(
  discovery_matched_delta_z =
    discovery_z[, ncol(discovery_z)] - discovery_z[, 1L],
  external_delta_z = external_z[, ncol(external_z)] - external_z[, 1L],
  shape_pearson = rowSums(discovery_z * external_z) /
    (ncol(external_z) - 1L),
  shape_spearman = vapply(seq_len(nrow(external_z)), function(i) {
    suppressWarnings(cor(
      discovery_z[i, ], external_z[i, ], method = "spearman"
    ))
  }, numeric(1))
)]
event_result[, `:=`(
  discovery_cluster_direction_concordant =
    sign(external_delta_z) == ifelse(discovery_cluster_id == 1L, 1, -1),
  matched_window_direction_concordant =
    sign(external_delta_z) == sign(discovery_matched_delta_z)
)]

summary_rows <- list()
plot_rows <- list()
set.seed(BOOT_SEED)
for (cluster in 1:2) {
  index <- which(events$cluster_id == cluster)
  gene_index <- split(index, events$gene_id[index])
  bootstrap_result <- replicate(BOOT_ITER, {
    sampled <- sample(names(gene_index), length(gene_index), replace = TRUE)
    sampled_index <- unlist(gene_index[sampled], use.names = FALSE)
    c(
      median_delta = median(
        external_z[sampled_index, ncol(external_z)] -
          external_z[sampled_index, 1L]
      ),
      colMeans(external_z[sampled_index, , drop = FALSE])
    )
  })
  bootstrap_delta <- bootstrap_result[1L, ]
  bootstrap_template <- bootstrap_result[-1L, , drop = FALSE]
  external_template <- colMeans(external_z[index, , drop = FALSE])
  discovery_template <- colMeans(discovery_z[index, , drop = FALSE])
  event_subset <- event_result[discovery_cluster_id == cluster]
  expected_sign <- ifelse(cluster == 1L, 1, -1)

  summary_rows[[cluster]] <- data.table(
    cohort = cohort,
    discovery_cluster_id = cluster,
    n_events = length(index),
    n_genes = length(gene_index),
    discovery_template_delta_z =
      discovery_template[length(discovery_template)] - discovery_template[1L],
    external_template_delta_z =
      external_template[length(external_template)] - external_template[1L],
    median_external_delta_z = median(event_subset$external_delta_z),
    median_external_delta_ci_low = quantile(bootstrap_delta, 0.025),
    median_external_delta_ci_high = quantile(bootstrap_delta, 0.975),
    expected_direction_fraction = mean(
      event_subset$discovery_cluster_direction_concordant
    ),
    matched_window_direction_fraction = mean(
      event_subset$matched_window_direction_concordant
    ),
    event_shape_pearson_median = median(event_subset$shape_pearson),
    event_shape_spearman_median = median(event_subset$shape_spearman),
    aggregate_template_pearson = cor(
      discovery_template, external_template, method = "pearson"
    ),
    aggregate_template_spearman = cor(
      discovery_template, external_template, method = "spearman"
    ),
    expected_direction_bootstrap_p =
      (sum(sign(bootstrap_delta) != expected_sign) + 1) /
        (BOOT_ITER + 1)
  )

  plot_rows[[cluster]] <- data.table(
    cohort = cohort,
    discovery_cluster_id = cluster,
    age = age,
    segment_id = segment_id,
    discovery_z = discovery_template,
    external_z = external_template,
    external_ci_low = apply(
      bootstrap_template, 1L, quantile, probs = 0.025
    ),
    external_ci_high = apply(
      bootstrap_template, 1L, quantile, probs = 0.975
    )
  )
}

event_replication <- event_result
trend_summary <- rbindlist(summary_rows)
trajectory_plot_data <- rbindlist(plot_rows)

# --- save results ---

age_association <- brainspan_age
age_association[, `:=`(
  brainspan_matrix_file = BRAINSPAN_MATRIX_FILE,
  brainspan_matrix_scale = "raw PDUI transformed to log2(PDUI/(1-PDUI))",
  brainspan_source_build = "hg19; events lifted to GRCh38 before matching"
)]
setcolorder(age_association, c(
  "cohort", "site_id", "gene_id", "discovery_cluster_id", "n_samples",
  "age_df", "residual_df", "age_f", "age_p", "age_q",
  "age_partial_r2", "age_selected", "brainspan_matrix_file",
  "brainspan_matrix_scale", "brainspan_source_build"
))

fwrite(age_association, OUT_AGE, sep = "\t")
fwrite(event_replication, OUT_EVENT, sep = "\t")
fwrite(trend_summary, OUT_SUMMARY, sep = "\t")
fwrite(trajectory_plot_data, OUT_PLOT, sep = "\t")

print(age_association[, .(
  n_mapped = .N,
  n_age_selected = sum(age_selected),
  age_selected_fraction = mean(age_selected)
), by = cohort])
print(trend_summary)
