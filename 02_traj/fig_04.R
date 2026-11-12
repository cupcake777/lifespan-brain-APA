#!/usr/bin/env Rscript
# Supported-grid trajectory inference from converged Lloyd clusters.
# Class medians with IQR ribbons, sampled raw lines, unsupported ages shaded.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("fig_common.R")

TRAJ_SAMPLE_PER_CLUSTER <- 55L

support <- fread(file.path(OUT_DIR, "01_support.tsv"))
sites <- fread(file.path(OUT_DIR, "02_sites_clustered.tsv"))
clustered <- sites[directional_pass == TRUE & !is.na(cluster_id)]
traj_z <- fread(file.path(INTERMEDIATE_DIR, "traj_z.csv"), showProgress = FALSE)

traj_cols <- grep("^grid_", names(traj_z), value = TRUE)
grid_idx <- as.integer(sub("grid_", "", traj_cols))
x_vals <- support$age_model[grid_idx + 1L]

traj_dyn <- merge(
  clustered[, .(site_id, cluster_id)],
  traj_z[, c("site_id", traj_cols), with = FALSE],
  by = "site_id"
)

# --- trajectory data: class medians and IQR ribbons ---

ribbon_dt <- rbindlist(lapply(sort(unique(traj_dyn$cluster_id)), function(k) {
  mat <- as.matrix(traj_dyn[cluster_id == k, ..traj_cols])
  data.table(
    cluster_id = k,
    age_model = x_vals,
    med = apply(mat, 2, median, na.rm = TRUE),
    q25 = apply(mat, 2, quantile, 0.25, na.rm = TRUE),
    q75 = apply(mat, 2, quantile, 0.75, na.rm = TRUE)
  )
}))
ribbon_dt[, cluster_label := cluster_label(cluster_id)]
ribbon_dt[, cluster_label := factor(cluster_label, levels = c("C1 lengthening", "C2 shortening"))]

line_sites <- clustered[, .SD[sample(.N, min(.N, TRAJ_SAMPLE_PER_CLUSTER))], by = cluster_id]
line_dt <- merge(
  line_sites[, .(site_id, cluster_id)],
  traj_z[, c("site_id", traj_cols), with = FALSE],
  by = "site_id"
)
line_long <- melt(
  line_dt,
  id.vars = c("site_id", "cluster_id"),
  measure.vars = traj_cols,
  variable.name = "grid",
  value.name = "z"
)
line_long[, age_model := x_vals[match(grid, traj_cols)]]
line_long[, cluster_label := cluster_label(cluster_id)]
line_long[, cluster_label := factor(cluster_label, levels = c("C1 lengthening", "C2 shortening"))]

label_dt <- ribbon_dt[age_model == max(age_model), .(
  cluster_id,
  cluster_label,
  age_model = age_model + 0.035,
  med,
  label = paste0("C", cluster_id)
)]

y_min <- min(-1.0, floor(min(ribbon_dt$q25, na.rm = TRUE) * 4) / 4 - 0.05)
y_max <- max(1.0, ceiling(max(ribbon_dt$q75, na.rm = TRUE) * 4) / 4 + 0.05)
unsupported_region <- support[is_supported == FALSE, .(
  xmin = min(age_model),
  xmax = max(age_model),
  label_x = mean(range(age_model))
)]

# --- plot ---

p04 <- ggplot() +
  geom_rect(data = unsupported_region,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey55", alpha = 0.16) +
  geom_hline(yintercept = 0, color = "#9CA3AF", linewidth = 0.35, linetype = "22") +
  geom_line(
    data = line_long,
    aes(age_model, z, group = site_id, color = cluster_label),
    alpha = 0.018,
    linewidth = 0.20
  ) +
  geom_ribbon(
    data = ribbon_dt,
    aes(age_model, ymin = q25, ymax = q75, fill = cluster_label),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    data = ribbon_dt,
    aes(age_model, med, color = cluster_label),
    linewidth = 1.15,
    lineend = "round"
  ) +
  geom_text(
    data = label_dt,
    aes(age_model, med, label = label, color = cluster_label),
    hjust = 0,
    vjust = 0.5,
    size = 3.35,
    fontface = "bold"
  ) +
  scale_color_manual(values = TRAJ_COLORS, guide = "none") +
  scale_fill_manual(values = TRAJ_COLORS, guide = "none") +
  scale_x_continuous(breaks = age_breaks, labels = age_labels, expand = c(0, 0)) +
  coord_cartesian(
    xlim = c(min(x_vals), max(x_vals) + 0.10),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  labs(x = "Age (years)", y = "Trajectory Z-score", title = "APA trajectory classes") +
  theme_pub(base = 10.8) +
  theme(
    axis.text = element_text(size = 9.8, face = "bold", color = "black"),
    axis.title = element_text(size = 11.5, face = "bold", color = "black"),
    panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.28),
    legend.position = "none",
    plot.title = element_text(size = 13.2, face = "bold"),
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  )

save_png(p04, "04_pattern.png", 6.2, 4.1)
save_pdf(p04, "04_pattern.pdf", 6.2, 4.1)

# --- layout validation ---

layout_validation <- rbindlist(list(
  clustered[, .(
    figure = "04_pattern",
    panel = paste0("C", cluster_id),
    cluster_full_label = cluster_label(cluster_id),
    n_dynamic_sites = .N,
    n_raw_lines_drawn = min(.N, TRAJ_SAMPLE_PER_CLUSTER),
    n_age_grid = length(traj_cols),
    y_min = y_min,
    y_max = y_max,
    max_label_chars = max(nchar(paste0("C", cluster_id))),
    validation_note = "All 128 fitted points are shown; the unsupported interval is shaded and excluded from clustering."
  ), by = cluster_id],
  data.table(
    figure = "04_pattern",
    panel = "axis",
    cluster_full_label = NA_character_,
    n_dynamic_sites = nrow(clustered),
    n_raw_lines_drawn = nrow(line_sites),
    n_age_grid = length(traj_cols),
    y_min = y_min,
    y_max = y_max,
    max_label_chars = max(nchar(age_labels)),
    validation_note = "Zero baseline is drawn; right plot margin reserves space for short endpoint labels."
  )
), fill = TRUE)
fwrite(layout_validation, file.path(OUT_DIR, "04_pattern_layout_validation.tsv"), sep = "\t")
