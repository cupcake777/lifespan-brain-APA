#!/usr/bin/env Rscript
# Supported-grid trajectory heatmap from converged Lloyd clusters.
# Full-grid context is shaded; unsupported ages are grey.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
source("fig_common.R")

HEATMAP_BINS_PER_CLUSTER <- 90L

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

# --- heatmap data: ordered trajectory quantile bins ---

site_order <- rbindlist(lapply(sort(unique(traj_dyn$cluster_id)), function(k) {
  sub <- traj_dyn[cluster_id == k]
  mat <- as.matrix(sub[, ..traj_cols])
  med <- apply(mat, 2, median, na.rm = TRUE)
  data.table(
    site_id = sub$site_id,
    cluster_id = k,
    rho_to_median = apply(mat, 1, safe_spearman, b = med)
  )
}))
site_order <- site_order[order(cluster_id, -rho_to_median, site_id)]
site_order[, site_rank := seq_len(.N), by = cluster_id]
site_order[, heat_bin := pmin(
  HEATMAP_BINS_PER_CLUSTER,
  pmax(1L, ceiling(site_rank / .N * HEATMAP_BINS_PER_CLUSTER))
), by = cluster_id]

heat_dt <- merge(
  site_order[, .(site_id, cluster_id, heat_bin)],
  traj_dyn,
  by = c("site_id", "cluster_id")
)

heat_long <- melt(
  heat_dt,
  id.vars = c("site_id", "cluster_id", "heat_bin"),
  measure.vars = traj_cols,
  variable.name = "grid",
  value.name = "z"
)
heat_long[, age_model := x_vals[match(grid, traj_cols)]]

heat_plot <- heat_long[, .(
  z = median(z, na.rm = TRUE),
  n_sites = .N
), by = .(cluster_id, heat_bin, age_model)]
heat_plot[, cluster_label := cluster_label(cluster_id)]
heat_plot[, cluster_label := factor(cluster_label, levels = c("C1 lengthening", "C2 shortening"))]
heat_plot[, cluster_short := factor(paste0("C", cluster_id), levels = c("C1", "C2"))]

# --- plot ---

tile_width <- median(diff(sort(unique(heat_plot$age_model))), na.rm = TRUE)
unsupported_region <- support[is_supported == FALSE, .(
  xmin = min(age_model) - tile_width / 2,
  xmax = max(age_model) + tile_width / 2
)]

p03 <- ggplot(heat_plot, aes(age_model, heat_bin, fill = z)) +
  geom_tile(width = tile_width, height = 1) +
  geom_rect(data = unsupported_region,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey35", alpha = 0.18) +
  facet_grid(cluster_short ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_vline(xintercept = 0, color = "#4B5563", linewidth = 0.26, alpha = 0.55) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"),
    values = rescale(c(-2.5, -1, 0, 1, 2.5)),
    limits = c(-2.5, 2.5),
    oob = squish,
    name = "Trajectory Z-score",
    na.value = "grey92",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3.0, "cm"),
      barheight = unit(0.22, "cm")
    )
  ) +
  scale_x_continuous(breaks = age_breaks, labels = age_labels, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Age (years)", y = NULL, title = "APA trajectory classes") +
  theme_pub(base = 10.5) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 9.6, face = "bold"),
    axis.title.x = element_text(size = 11.2, face = "bold"),
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 10.0, angle = 0, face = "bold", color = "#111827"),
    strip.placement = "outside",
    panel.spacing = unit(0.18, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 9.4, face = "bold"),
    legend.text = element_text(size = 8.8, face = "bold"),
    plot.title = element_text(size = 13.0, face = "bold"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

save_png(p03, "03_pattern_heatmap.png", 5.8, 5.2)
save_pdf(p03, "03_pattern_heatmap.pdf", 5.8, 5.2)

# --- layout validation ---

layout_validation <- heat_plot[, .(
  figure = "03_pattern_heatmap",
  panel = unique(as.character(cluster_short)),
  cluster_full_label = unique(as.character(cluster_label)),
  n_dynamic_sites = uniqueN(site_order[cluster_id == .BY$cluster_id]$site_id),
  n_heatmap_bins = uniqueN(heat_bin),
  n_age_grid = uniqueN(age_model),
  max_label_chars = max(nchar(as.character(cluster_short))),
  colorbar = "bottom horizontal",
  validation_note = "Raw site rows are compressed to ordered quantile-bin median rows; visible class labels are short IDs and full labels are recorded in this table."
), by = cluster_id]
fwrite(layout_validation, file.path(OUT_DIR, "03_pattern_heatmap_layout_validation.tsv"), sep = "\t")
