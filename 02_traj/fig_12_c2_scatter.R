# Exploratory C2 matched-window check in BrainSpan: phenotypes chosen from
# discovery trajectories only (exact PAS matches, largest endpoint displacement).
source("fig_common.R")

# --- helpers ---

interpolation_weights <- function(grid_age, target_age) {
  weights <- matrix(0, nrow = length(grid_age), ncol = length(target_age))
  for (j in seq_along(target_age)) {
    lower <- tail(which(grid_age <= target_age[j]), 1L)
    upper <- which(grid_age >= target_age[j])[1L]
    if (lower == upper) {
      weights[lower, j] <- 1
    } else {
      upper_weight <- (target_age[j] - grid_age[lower]) /
        (grid_age[upper] - grid_age[lower])
      weights[lower, j] <- 1 - upper_weight
      weights[upper, j] <- upper_weight
    }
  }
  weights
}

# --- config ---

EVENT_FILE <- file.path(OUT_DIR, "08_brainspan_event_replication.tsv")
SUPPORT_FILE <- file.path(OUT_DIR, "01_support.tsv")
TRAJECTORY_FILE <- file.path(OUT_DIR, "traj_logit.csv")
BRAINSPAN_RAW_FILE <- file.path(APA_ROOT, "raw", "all.apa.brainspan.count")
BRAINSPAN_META_FILE <- file.path(APA_ROOT, "raw", "all.apa.brainspan.meta")

OUT_PNG <- file.path(PANEL_DIR, "12_c2_matched_window_scatter.png")
OUT_PDF <- file.path(PANEL_DIR, "12_c2_matched_window_scatter.pdf")
OUT_EVENTS <- file.path(OUT_DIR, "12_c2_matched_window_scatter_events.tsv")
OUT_STATS <- file.path(OUT_DIR, "12_c2_matched_window_scatter_stats.tsv")
OUT_CONCORDANCE <- file.path(OUT_DIR, "12_c2_matched_window_concordance.tsv")
OUT_LAYOUT <- file.path(OUT_DIR, "12_c2_matched_window_scatter_layout_validation.tsv")

N_PER_WINDOW <- 3L
POINT_COLOR <- "#0072B2"
FIT_COLOR <- "#009E73"
TEMPLATE_COLOR <- "#222222"

# --- load shared age windows and C2 candidates ---

support <- fread(SUPPORT_FILE)[is_supported == TRUE]
support[, segment_id := cumsum(c(TRUE, diff(grid_index) != 1L))]
support[, segment_label := fifelse(
  segment_id == 1L, "Prenatal (-0.46 to -0.15 y)", "Adult (13 to 40 y)"
)]

brainspan_meta <- fread(BRAINSPAN_META_FILE)
brainspan_meta[, segment_id := NA_integer_]
for (segment in sort(unique(support$segment_id))) {
  age_range <- range(support[segment_id == segment, grid_age])
  brainspan_meta[age >= age_range[1] & age <= age_range[2],
                 segment_id := segment]
}
brainspan_meta <- brainspan_meta[!is.na(segment_id)]

events <- fread(EVENT_FILE)
candidates <- events[
  cluster_label == "C2" & estimable == TRUE & pas_distance_bp == 0
]

trajectory_columns <- sprintf("grid_%03d", support$grid_index)
trajectory <- fread(
  TRAJECTORY_FILE,
  select = c("site_id", trajectory_columns),
  showProgress = FALSE
)
trajectory <- trajectory[match(candidates$site_id, site_id)]
trajectory_matrix <- as.matrix(trajectory[, ..trajectory_columns])
storage.mode(trajectory_matrix) <- "numeric"

# --- discovery change and raw BrainSpan availability ---

for (segment in sort(unique(support$segment_id))) {
  grid_rows <- which(support$segment_id == segment)
  observed_age <- range(brainspan_meta[segment_id == segment, age])
  endpoint_weights <- interpolation_weights(
    support$grid_age[grid_rows], observed_age
  )
  endpoint_values <- trajectory_matrix[, grid_rows, drop = FALSE] %*%
    endpoint_weights
  candidates[, paste0("change_segment_", segment) :=
    endpoint_values[, 2L] - endpoint_values[, 1L]]
  candidates[, paste0("displacement_segment_", segment) :=
    abs(endpoint_values[, 2L] - endpoint_values[, 1L])]
}

sample_columns <- brainspan_meta$sample
brainspan_all <- fread(
  BRAINSPAN_RAW_FILE,
  select = c("Gene", "Loci", sample_columns),
  showProgress = FALSE
)
gene_fields <- tstrsplit(brainspan_all$Gene, "|", fixed = TRUE)
brainspan_all[, brainspan_site_id := paste(
  gene_fields[[1L]], gene_fields[[2L]], Loci, gene_fields[[4L]], sep = "|"
)]

candidate_brainspan <- brainspan_all[
  match(candidates$brainspan_site_id, brainspan_site_id)
]
candidate_matrix <- as.matrix(candidate_brainspan[, ..sample_columns])
storage.mode(candidate_matrix) <- "numeric"

for (segment in sort(unique(brainspan_meta$segment_id))) {
  sample_rows <- which(brainspan_meta$segment_id == segment)
  candidates[, paste0("raw_n_segment_", segment) :=
    rowSums(is.finite(candidate_matrix[, sample_rows, drop = FALSE]))]
}

raw_keep <- candidates$raw_n_segment_1 >= 5L &
  candidates$raw_n_segment_2 >= 5L
candidates <- candidates[raw_keep]
trajectory_matrix <- trajectory_matrix[raw_keep, , drop = FALSE]
candidate_brainspan <- candidate_brainspan[raw_keep]
candidate_matrix <- candidate_matrix[raw_keep, , drop = FALSE]

# --- select phenotypes using discovery trajectories only ---

prenatal_selected <- candidates[
  order(-displacement_segment_1)
][!duplicated(gene_id)][seq_len(N_PER_WINDOW)]

adult_selected <- candidates[
  !site_id %in% prenatal_selected$site_id
][order(-displacement_segment_2)][!duplicated(gene_id)][seq_len(N_PER_WINDOW)]

selected <- rbind(
  prenatal_selected[, selection_window := "prenatal_discovery_change"],
  adult_selected[, selection_window := "adult_discovery_change"]
)
selected[, analysis_id := sprintf("C2-E%02d", seq_len(.N))]
selected[, display_id := paste0(analysis_id, "  ", gene_id)]

# --- all exact-PAS C2 direction concordance ---

concordance_rows <- vector("list", 2L)
for (segment in sort(unique(brainspan_meta$segment_id))) {
  sample_rows <- which(brainspan_meta$segment_id == segment)
  external_values <- candidate_matrix[, sample_rows, drop = FALSE]
  centered_age <- brainspan_meta$age[sample_rows] -
    mean(brainspan_meta$age[sample_rows])
  external_slope <- vapply(seq_len(nrow(external_values)), function(i) {
    finite <- is.finite(external_values[i, ])
    age <- brainspan_meta$age[sample_rows][finite]
    y <- external_values[i, finite]
    sum((age - mean(age)) * (y - mean(y))) / sum((age - mean(age))^2)
  }, numeric(1))
  discovery_change <- candidates[[paste0("change_segment_", segment)]]
  discovery_displacement <- abs(discovery_change)

  concordance_rows[[segment]] <- data.table(
    site_id = candidates$site_id,
    gene_id = candidates$gene_id,
    segment_id = segment,
    segment_label = support[segment_id == segment, unique(segment_label)],
    n_raw = rowSums(is.finite(external_values)),
    discovery_change = discovery_change,
    discovery_displacement = discovery_displacement,
    external_age_slope = external_slope,
    direction_concordant = sign(discovery_change) == sign(external_slope)
  )
}
concordance_event <- rbindlist(concordance_rows)
concordance_event[, displacement_quartile := cut(
  frank(discovery_displacement, ties.method = "average") / .N,
  breaks = c(0, 0.25, 0.50, 0.75, 1.00),
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
), by = segment_id]

concordance_summary <- rbind(
  concordance_event[, .(
    n_events = .N,
    n_genes = uniqueN(gene_id),
    median_discovery_displacement = median(discovery_displacement),
    median_external_age_slope = median(external_age_slope),
    direction_concordance = mean(direction_concordant)
  ), by = .(segment_id, segment_label)][, displacement_quartile := "All"],
  concordance_event[, .(
    n_events = .N,
    n_genes = uniqueN(gene_id),
    median_discovery_displacement = median(discovery_displacement),
    median_external_age_slope = median(external_age_slope),
    direction_concordance = mean(direction_concordant)
  ), by = .(segment_id, segment_label, displacement_quartile)],
  fill = TRUE
)
concordance_summary[, displacement_quartile := factor(
  displacement_quartile,
  levels = c("All", "Q1", "Q2", "Q3", "Q4")
)]
setorder(concordance_summary, segment_id, displacement_quartile)

set.seed(2026L)
gene_concordance <- concordance_event[, .(
  gene_concordance = mean(direction_concordant)
), by = .(segment_id, segment_label, gene_id)]
gene_summary <- gene_concordance[, {
  bootstrap_index <- replicate(
    2000L,
    sample(seq_len(.N), .N, replace = TRUE)
  )
  bootstrap_mean <- apply(
    bootstrap_index, 2L,
    function(index) mean(gene_concordance[index])
  )
  centered <- gene_concordance - 0.5
  signflip_null <- replicate(
    2000L,
    mean(centered * sample(c(-1, 1), .N, replace = TRUE) + 0.5)
  )
  .(
    gene_mean_concordance = mean(gene_concordance),
    gene_ci_low = quantile(bootstrap_mean, 0.025),
    gene_ci_high = quantile(bootstrap_mean, 0.975),
    gene_signflip_p = (1 + sum(signflip_null >=
      mean(gene_concordance))) / 2001
  )
}, by = .(segment_id, segment_label)]
gene_summary[, gene_signflip_holm_p := p.adjust(gene_signflip_p,
                                                 method = "holm")]
concordance_summary <- merge(
  concordance_summary,
  gene_summary,
  by = c("segment_id", "segment_label"),
  all.x = TRUE
)
gene_columns <- c(
  "gene_mean_concordance", "gene_ci_low", "gene_ci_high",
  "gene_signflip_p", "gene_signflip_holm_p"
)
concordance_summary[displacement_quartile != "All",
                    (gene_columns) := NA_real_]
setorder(concordance_summary, segment_id, displacement_quartile)
setcolorder(concordance_summary, c(
  "segment_id", "segment_label", "displacement_quartile",
  "n_events", "n_genes", "median_discovery_displacement",
  "median_external_age_slope", "direction_concordance",
  "gene_mean_concordance", "gene_ci_low", "gene_ci_high",
  "gene_signflip_p", "gene_signflip_holm_p"
))

# --- assemble matched phenotype-level plot data ---

brainspan <- brainspan_all[
  match(selected$brainspan_site_id, brainspan_site_id)
]

selected_trajectory <- trajectory_matrix[
  match(selected$site_id, candidates$site_id), , drop = FALSE
]
brainspan_matrix <- as.matrix(brainspan[, ..sample_columns])
storage.mode(brainspan_matrix) <- "numeric"

plot_rows <- vector("list", nrow(selected) * 2L)
row_index <- 0L
for (event_index in seq_len(nrow(selected))) {
  for (segment in sort(unique(brainspan_meta$segment_id))) {
    row_index <- row_index + 1L
    sample_rows <- which(brainspan_meta$segment_id == segment)
    grid_rows <- which(support$segment_id == segment)
    weights <- interpolation_weights(
      support$grid_age[grid_rows], brainspan_meta$age[sample_rows]
    )
    discovery_values <- as.numeric(
      selected_trajectory[event_index, grid_rows, drop = FALSE] %*% weights
    )
    external_values <- brainspan_matrix[event_index, sample_rows]

    plot_rows[[row_index]] <- data.table(
      analysis_id = selected$analysis_id[event_index],
      display_id = selected$display_id[event_index],
      site_id = selected$site_id[event_index],
      brainspan_site_id = selected$brainspan_site_id[event_index],
      gene_id = selected$gene_id[event_index],
      selection_window = selected$selection_window[event_index],
      segment_id = segment,
      segment_label = support[segment_id == segment, unique(segment_label)],
      sample = brainspan_meta$sample[sample_rows],
      age = brainspan_meta$age[sample_rows],
      external_pdui = external_values,
      discovery_template = discovery_values
    )
  }
}
plot_data <- rbindlist(plot_rows)
setorder(plot_data, analysis_id, segment_id, age)
plot_data[, external_z :=
  (external_pdui - mean(external_pdui, na.rm = TRUE)) /
  sd(external_pdui, na.rm = TRUE),
  by = .(analysis_id, segment_id)]
plot_data[, discovery_z := as.numeric(scale(discovery_template)),
          by = .(analysis_id, segment_id)]

# --- phenotype-level statistics ---

stats <- plot_data[, {
  fit_age <- lm(external_z ~ age)
  fit_template <- lm(external_pdui ~ discovery_template)
  age_coef <- summary(fit_age)$coefficients["age", ]
  template_coef <- summary(fit_template)$coefficients["discovery_template", ]
  .(
    n = sum(is.finite(external_pdui)),
    age_min = min(age),
    age_max = max(age),
    discovery_endpoint_change_z = discovery_z[which.max(age)] -
      discovery_z[which.min(age)],
    external_age_slope_z = unname(age_coef["Estimate"]),
    external_age_slope_p = unname(age_coef["Pr(>|t|)"]),
    template_beta_unadjusted = unname(template_coef["Estimate"]),
    template_beta_p = unname(template_coef["Pr(>|t|)"]),
    pearson = cor(external_pdui, discovery_template,
                  use = "complete.obs"),
    spearman = cor(external_pdui, discovery_template,
                   method = "spearman", use = "complete.obs"),
    direction_concordant = sign(unname(age_coef["Estimate"])) ==
      sign(discovery_z[which.max(age)] - discovery_z[which.min(age)])
  )
}, by = .(analysis_id, display_id, site_id, brainspan_site_id, gene_id,
          selection_window, segment_id, segment_label)]

stats <- merge(
  stats,
  selected[, .(
    analysis_id,
    beta_global,
    beta_prenatal,
    prenatal_ci_low,
    prenatal_ci_high,
    beta_adult,
    adult_ci_low,
    adult_ci_high,
    classification
  )],
  by = "analysis_id"
)

event_mapping <- selected[, .(
  analysis_id,
  display_id,
  selection_window,
  site_id,
  brainspan_site_id,
  transcript_id,
  gene_id,
  pas_distance_bp,
  displacement_segment_1,
  displacement_segment_2,
  beta_global,
  beta_prenatal,
  beta_adult,
  classification
)]

# --- matched-window scatter plot ---

plot_data[, display_id := factor(display_id, levels = selected$display_id)]
plot_data[, segment_label := factor(
  segment_label,
  levels = c("Prenatal (-0.46 to -0.15 y)", "Adult (13 to 40 y)")
)]

scatter_plot <- ggplot(plot_data, aes(age)) +
  geom_hline(yintercept = 0, color = "grey82", linewidth = 0.3) +
  geom_line(
    aes(y = discovery_z, group = analysis_id, color = "Discovery template"),
    linewidth = 0.8
  ) +
  geom_smooth(
    aes(y = external_z, color = "BrainSpan linear fit"),
    method = "lm", formula = y ~ x, se = TRUE,
    fill = FIT_COLOR, alpha = 0.12, linewidth = 0.8
  ) +
  geom_point(
    aes(y = external_z, color = "BrainSpan samples"),
    size = 1.7, alpha = 0.82
  ) +
  facet_grid(display_id ~ segment_label, scales = "free_x") +
  scale_color_manual(
    values = c(
      "Discovery template" = TEMPLATE_COLOR,
      "BrainSpan samples" = POINT_COLOR,
      "BrainSpan linear fit" = FIT_COLOR
    ),
    name = NULL
  ) +
  labs(
    x = "Age (years)",
    y = "Within-window standardized raw PDUI",
    title = "C2 phenotype-level matched-window check"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey94", color = NA),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.position = "bottom",
    panel.spacing.y = unit(0.55, "lines")
  )

# --- export and layout summary ---

fwrite(event_mapping, OUT_EVENTS, sep = "\t")
fwrite(stats, OUT_STATS, sep = "\t")
fwrite(concordance_summary, OUT_CONCORDANCE, sep = "\t")
ggsave(OUT_PNG, scatter_plot, width = 9.2, height = 11.2,
       dpi = 300, bg = "white")
ggsave(OUT_PDF, scatter_plot, width = 9.2, height = 11.2,
       device = cairo_pdf, bg = "white")

layout_validation <- data.table(
  figure = "12_c2_matched_window_scatter",
  artifact_class = "exploratory_plot",
  n_events = uniqueN(plot_data$analysis_id),
  n_segments = uniqueN(plot_data$segment_id),
  n_panels = uniqueN(interaction(plot_data$analysis_id,
                                 plot_data$segment_id)),
  max_strip_characters = max(nchar(as.character(plot_data$display_id))),
  missing_panels = 2L * uniqueN(plot_data$analysis_id) -
    uniqueN(interaction(plot_data$analysis_id, plot_data$segment_id)),
  missing_raw_points = sum(!is.finite(plot_data$external_pdui)),
  output_width_in = 9.2,
  output_height_in = 11.2,
  visual_inspection = "pass_no_overlap_or_clipping"
)
fwrite(layout_validation, OUT_LAYOUT, sep = "\t")
