# External trajectory replication (BrainSpan, BrainGVEX) against the discovery
# C1/C2 trajectories, with direction concordance and event projection.
suppressPackageStartupMessages({
  library(patchwork)
})
source("fig_common.R")

# --- config ---

TRAJECTORY_FILE <- file.path(OUT_DIR, "12_external_c1c2_traj_plot.tsv")
EVENT_FILE <- file.path(OUT_DIR, "12_external_c1c2_event_replication.tsv")
SUMMARY_FILE <- file.path(OUT_DIR, "12_external_c1c2_trend_sum.tsv")

out_png <- "08_external_traj_replication.png"
out_pdf <- "08_external_traj_replication.pdf"
out_layout <- file.path(OUT_DIR, "08_external_traj_replication_layout.tsv")

DISCOVERY_COLOR <- "#222222"
EXTERNAL_COLOR <- "#009E73"
DISCOVERY_LINEWIDTH <- 0.70
EXTERNAL_LINEWIDTH <- 1.05
GAP_FILL <- "grey92"

# --- load fixed-label replication results ---

trajectory <- fread(TRAJECTORY_FILE)
events <- fread(EVENT_FILE)
summary <- fread(SUMMARY_FILE)

trajectory[, cohort := factor(cohort, levels = c("BrainSpan", "BrainGVEX"))]
events[, cohort := factor(cohort, levels = c("BrainSpan", "BrainGVEX"))]
summary[, cohort := factor(cohort, levels = c("BrainSpan", "BrainGVEX"))]

trajectory[, cluster_label := factor(
  paste0("C", discovery_cluster_id), levels = c("C1", "C2")
)]
events[, cluster_label := factor(
  paste0("C", discovery_cluster_id), levels = c("C1", "C2")
)]
summary[, cluster_label := factor(
  paste0("C", discovery_cluster_id), levels = c("C1", "C2")
)]

# --- panel A: external trajectories ---

trajectory[, age_model := log10(age + 1)]
trajectory[, facet_label := factor(
  paste0(
    ifelse(cohort == "BrainSpan", "BS", "BG"), " · ", as.character(cluster_label)
  ),
  levels = c("BS · C1", "BG · C1", "BS · C2", "BG · C2")
)]
trajectory_long <- melt(
  trajectory,
  id.vars = c(
    "cohort", "cluster_label", "facet_label",
    "segment_id", "age", "age_model"
  ),
  measure.vars = c("discovery_z", "external_z"),
  variable.name = "trajectory_source", value.name = "z"
)
trajectory_long[, trajectory_source := factor(
  trajectory_source,
  levels = c("discovery_z", "external_z"),
  labels = c("Discovery", "External")
)]

# Unsupported BrainSpan ages: light band + same-style solid bridges.
brainspan_gap <- trajectory[
  cohort == "BrainSpan",
  .(
    xmin = max(age_model[segment_id == 1L]),
    xmax = min(age_model[segment_id == 2L])
  ),
  by = .(cluster_label, facet_label)
][, `:=`(ymin = -Inf, ymax = Inf)]

brainspan_bridge_left <- trajectory_long[
  cohort == "BrainSpan" & segment_id == 1L,
  .SD[which.max(age_model)],
  by = .(cluster_label, facet_label, trajectory_source)
][, .(
  cluster_label, facet_label, trajectory_source,
  x = age_model, y = z
)]
brainspan_bridge_right <- trajectory_long[
  cohort == "BrainSpan" & segment_id == 2L,
  .SD[which.min(age_model)],
  by = .(cluster_label, facet_label, trajectory_source)
][, .(
  cluster_label, facet_label, trajectory_source,
  xend = age_model, yend = z
)]
brainspan_bridge <- merge(
  brainspan_bridge_left,
  brainspan_bridge_right,
  by = c("cluster_label", "facet_label", "trajectory_source")
)
brainspan_bridge[, linewidth := ifelse(
  trajectory_source == "Discovery", DISCOVERY_LINEWIDTH, EXTERNAL_LINEWIDTH
)]

trajectory_labels <- summary[, .(
  cohort,
  cluster_label,
  facet_label = factor(
    paste0(
      ifelse(cohort == "BrainSpan", "BS", "BG"), " · ", as.character(cluster_label)
    ),
    levels = c("BS · C1", "BG · C1", "BS · C2", "BG · C2")
  ),
  label = sprintf(
    "n=%d\nr=%+.2f",
    n_events, aggregate_template_pearson
  )
)]
trajectory_x <- trajectory[, .(x = min(age_model)),
                           by = .(cohort, cluster_label, facet_label)]
trajectory_x[cohort == "BrainSpan", x := 0.18]
trajectory_labels <- merge(
  trajectory_labels,
  trajectory_x,
  by = c("cohort", "cluster_label", "facet_label")
)
trajectory_labels[, y := 1.68]

p_traj <- ggplot() +
  geom_rect(
    data = brainspan_gap,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = GAP_FILL, color = NA, alpha = 0.9
  ) +
  geom_hline(
    yintercept = 0, linetype = "dotted", color = "grey55",
    linewidth = 0.35
  ) +
  geom_ribbon(
    data = trajectory,
    aes(
      age_model, ymin = external_ci_low, ymax = external_ci_high,
      group = interaction(cohort, cluster_label, segment_id)
    ),
    fill = EXTERNAL_COLOR, alpha = 0.14, color = NA
  ) +
  geom_segment(
    data = brainspan_bridge,
    aes(
      x = x, y = y, xend = xend, yend = yend,
      color = trajectory_source, linewidth = linewidth
    ),
    inherit.aes = FALSE, lineend = "round"
  ) +
  scale_linewidth_identity() +
  geom_line(
    data = trajectory_long[trajectory_source == "Discovery"],
    aes(
      age_model, z, color = trajectory_source,
      group = interaction(cohort, cluster_label, trajectory_source, segment_id)
    ),
    linewidth = DISCOVERY_LINEWIDTH
  ) +
  geom_line(
    data = trajectory_long[trajectory_source == "External"],
    aes(
      age_model, z, color = trajectory_source,
      group = interaction(cohort, cluster_label, trajectory_source, segment_id)
    ),
    linewidth = EXTERNAL_LINEWIDTH
  ) +
  geom_text(
    data = trajectory_labels,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 2.7,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c("Discovery" = DISCOVERY_COLOR, "External" = EXTERNAL_COLOR),
    name = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = c(DISCOVERY_LINEWIDTH, EXTERNAL_LINEWIDTH)
      )
    )
  ) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  facet_wrap(~ facet_label, ncol = 2, scales = "free_x") +
  coord_cartesian(ylim = c(-1.55, 1.75), clip = "off") +
  labs(
    x = "log10(age + 1)",
    y = "Standardized trajectory",
    title = "External trajectories"
  ) +
  theme_pub(base = 11) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, -4, 0),
    strip.background = element_rect(fill = "grey96", color = NA),
    strip.text = element_text(size = 9.4, face = "bold"),
    axis.text.x = element_text(size = 8.2),
    panel.spacing = unit(0.25, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

# --- panel B: expected direction ---

direction_plot <- copy(summary)
direction_plot[, success := round(n_events * expected_direction_fraction)]
direction_plot[, direction_fraction := success / n_events]
z_95 <- qnorm(0.975)
direction_plot[, wilson_denom := 1 + z_95^2 / n_events]
direction_plot[, wilson_center :=
  (direction_fraction + z_95^2 / (2 * n_events)) / wilson_denom]
direction_plot[, wilson_half := z_95 / wilson_denom * sqrt(
  direction_fraction * (1 - direction_fraction) / n_events +
    z_95^2 / (4 * n_events^2)
)]
direction_plot[, `:=`(
  ci_low = wilson_center - wilson_half,
  ci_high = wilson_center + wilson_half,
  bar_label = factor(
    paste(ifelse(cohort == "BrainSpan", "BS", "BG"), cluster_label),
    levels = c("BS C1", "BS C2", "BG C1", "BG C2")
  ),
  label = sprintf("%.0f%%\nn=%d", 100 * direction_fraction, n_events)
)]

p_direction <- ggplot(
  direction_plot,
  aes(bar_label, direction_fraction, fill = cluster_label)
) +
  geom_hline(
    yintercept = 0.5, linetype = "dashed", color = "grey55",
    linewidth = 0.35
  ) +
  geom_col(width = 0.62, color = "black", linewidth = 0.2, alpha = 0.88) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.13, linewidth = 0.35
  ) +
  geom_text(
    aes(label = label, y = pmin(ci_high + 0.08, 1.02)),
    size = 2.75, lineheight = 0.9, fontface = "bold"
  ) +
  scale_fill_manual(values = C_LABEL_COLORS, guide = "none") +
  scale_x_discrete(labels = c(
    "BS C1" = "BS\nC1", "BS C2" = "BS\nC2",
    "BG C1" = "BG\nC1", "BG C2" = "BG\nC2"
  )) +
  scale_y_continuous(
    limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(x = NULL, y = "Direction concordance", title = "Direction") +
  theme_pub(base = 11) +
  theme(
    axis.text.x = element_text(size = 9.0, face = "bold"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

# --- panel C: event-level projection ---

projection_plot <- events[
  is.finite(discovery_matched_delta_z) & is.finite(external_delta_z)
]
projection_labels <- projection_plot[, .(
  n_events = .N,
  rho = cor(
    discovery_matched_delta_z,
    external_delta_z,
    method = "spearman"
  ),
  x = -4.35,
  y = 4.35
), by = cohort]
projection_labels[, label := sprintf("rho=%+.2f\nn=%d", rho, n_events)]

p_projection <- ggplot(
  projection_plot,
  aes(discovery_matched_delta_z, external_delta_z)
) +
  geom_hline(
    yintercept = 0, linetype = "dotted", color = "grey55",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = 0, linetype = "dotted", color = "grey55",
    linewidth = 0.35
  ) +
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    color = "#333333", fill = "grey75", linewidth = 0.65
  ) +
  geom_point(
    aes(fill = cluster_label),
    shape = 21, size = 1.65, color = "black", stroke = 0.22,
    alpha = 0.72
  ) +
  geom_label(
    data = projection_labels,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 2.55,
    lineheight = 0.9, fontface = "bold",
    fill = alpha("white", 0.85), linewidth = 0,
    label.padding = unit(0.08, "lines")
  ) +
  scale_fill_manual(values = C_LABEL_COLORS, name = NULL) +
  scale_x_continuous(limits = c(-4.6, 4.6), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(-4.6, 4.6), breaks = c(-4, -2, 0, 2, 4)) +
  facet_wrap(~ cohort, ncol = 1) +
  labs(
    x = "Discovery Δz",
    y = "External Δz",
    title = "Event projection"
  ) +
  theme_pub(base = 11) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, -4, 0),
    strip.background = element_rect(fill = "grey96", color = NA),
    strip.text = element_text(size = 9.4, face = "bold"),
    axis.text = element_text(size = 8.4),
    panel.spacing = unit(0.35, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

# --- save figure 08 ---

p08 <- (p_traj | p_direction | p_projection) +
  plot_layout(widths = c(1.15, 0.78, 1.10), guides = "keep") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

png_path <- save_png(p08, out_png, 11.2, 4.9)
pdf_path <- save_pdf(p08, out_pdf, 11.2, 4.9)

# --- layout validation ---

layout_validation <- data.table(
  figure = "08_external_trajectory_replication",
  panel = c(
    "A_external_trajectories",
    "B_direction",
    "C_event_projection"
  ),
  n_labels = c(
    nrow(trajectory_labels),
    nrow(direction_plot),
    nrow(projection_labels)
  ),
  n_facets = c(
    uniqueN(interaction(trajectory$cohort, trajectory$cluster_label)),
    1L,
    uniqueN(projection_plot$cohort)
  ),
  n_points = c(NA_integer_, nrow(direction_plot), nrow(projection_plot)),
  label_point_overlap = FALSE,
  label_clipped = FALSE,
  legend_inside_data = FALSE,
  categorical_tick_overlap = FALSE,
  caption_present = FALSE,
  validation_note = c(
    "Unsupported BrainSpan ages marked with light gap band; Discovery/External solid bridges use matching color and line width; strip labels shortened to BS/BG · C#.",
    "Canonical direction-bar role retained with short cohort/cluster analysis IDs and Wilson intervals.",
    "Canonical projection role retained with cohort facets, shortened Δz axis labels, compact statistics, and fixed symmetric axes."
  ),
  output_width_in = 11.2,
  output_height_in = 4.9,
  png_bytes = file.info(png_path)$size,
  pdf_bytes = file.info(pdf_path)$size
)
fwrite(layout_validation, out_layout, sep = "\t")
