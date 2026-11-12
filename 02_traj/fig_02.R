#!/usr/bin/env Rscript
# Dynamic APA site filter map: age q-value against delta logit(PDUI).
# Sites are coloured by filter gate: Dynamic, FDR only, Effect only, or Neither.
source("fig_common.R")

SCATTER_SAMPLE <- 45000L
sites <- fread(file.path(OUT_DIR, "02_sites.tsv"))
thr_values <- sort(unique(na.omit(sites$delta_logit_threshold)))
thr <- median(thr_values)

# --- filter gates ---

sites[, q_pass := age_qvalue <= FDR_THR]
sites[, delta_pass := abs(delta_logit) >= delta_logit_threshold]
sites[, filter_gate := fifelse(q_pass == TRUE & delta_pass == TRUE, "Dynamic",
                        fifelse(q_pass == TRUE & delta_pass == FALSE, "FDR only",
                        fifelse(q_pass == FALSE & delta_pass == TRUE, "Effect only",
                                "Neither")))]
sites[, filter_gate := factor(filter_gate,
                              levels = c("Neither", "Effect only", "FDR only", "Dynamic"))]

plot_dt <- sites[is.finite(delta_logit) & is.finite(age_qvalue)]
full_counts <- plot_dt[, .(n = .N), by = filter_gate]
full_counts[, pct := round(n / sum(n) * 100, 1)]

if (nrow(plot_dt) > SCATTER_SAMPLE) {
  plot_dt <- plot_dt[sample(.N, SCATTER_SAMPLE)]
}
plot_dt[, filter_gate := factor(filter_gate, levels = levels(sites$filter_gate))]

x_lim <- quantile(sites$delta_logit, c(0.005, 0.995), na.rm = TRUE)
y_values <- -log10(pmax(sites$age_qvalue, .Machine$double.xmin))
y_cap <- min(quantile(y_values, 0.995, na.rm = TRUE), 60)

# --- plot ---

# quadrant boundaries
x_cut <- thr          # ~0.4
y_cut <- -log10(FDR_THR)  # ~1.3

p02 <- ggplot(plot_dt, aes(delta_logit, -log10(pmax(age_qvalue, .Machine$double.xmin)))) +
  # subtle quadrant shading, only the Dynamic quadrants
  annotate("rect", xmin = -Inf, xmax = -thr, ymin = y_cut, ymax = Inf,
           fill = "#D55E00", alpha = 0.03) +
  annotate("rect", xmin = thr, xmax = Inf, ymin = y_cut, ymax = Inf,
           fill = "#D55E00", alpha = 0.03) +
  # threshold lines
  geom_vline(xintercept = c(-thr, thr), color = "#1a1a1a", linetype = "dashed", linewidth = 0.55) +
  geom_hline(yintercept = y_cut, color = "#1a1a1a", linetype = "dotted", linewidth = 0.55) +
  # points
  geom_point(aes(color = filter_gate), alpha = 0.15, size = 0.6, stroke = 0) +
  # axis lines on top
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  # labels, each inside its quadrant with an arrow to the cluster

  # "Effect only": lower-left quadrant (|delta| < thr, qval > FDR)
  annotate("segment",
           x = -x_cut * 3.5, y = y_cut + 3.5,
           xend = -x_cut * 1.8, yend = y_cut + 1.5,
           color = "#0072B2", linewidth = 0.4,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  annotate("text",
           x = -x_cut * 3.7, y = y_cut + 4.5,
           label = paste0("Effect only\nn = ", comma(full_counts[filter_gate == "Effect only", n]),
                          " (", full_counts[filter_gate == "Effect only", pct], "%)"),
           color = "#0072B2", size = 3.1, hjust = 1, fontface = "bold") +

  # "Neither": sparse mid-left area, arrow down to the gray sub-FDR band
  annotate("segment",
           x = x_lim[1] + 1.2, y = y_cap * 0.17,
           xend = x_lim[1] + 1.2, yend = y_cut + 0.6,
           color = "#999999", linewidth = 0.4,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  annotate("text",
           x = x_lim[1] + 0.5, y = y_cap * 0.22,
           label = paste0("Neither\nn = ", comma(full_counts[filter_gate == "Neither", n]),
                          " (", full_counts[filter_gate == "Neither", pct], "%)"),
           color = "#4D4D4D", size = 3.1, hjust = 0, fontface = "bold") +

  # "FDR only": upper-center (|delta| < thr, qval < FDR)
  annotate("segment",
           x = -x_cut * 1.5, y = y_cap * 0.55,
           xend = -x_cut * 0.5, yend = y_cap * 0.38,
           color = "#009E73", linewidth = 0.4,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  annotate("text",
           x = -x_cut * 1.7, y = y_cap * 0.58,
           label = paste0("FDR only\nn = ", comma(full_counts[filter_gate == "FDR only", n]),
                          " (", full_counts[filter_gate == "FDR only", pct], "%)"),
           color = "#009E73", size = 3.1, hjust = 1, fontface = "bold") +

  # "Dynamic": upper-right quadrant (|delta| >= thr, qval < FDR)
  annotate("segment",
           x = x_cut * 3.5, y = y_cap * 0.92,
           xend = x_cut * 2.0, yend = y_cap * 0.72,
           color = "#D55E00", linewidth = 0.4,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  annotate("text",
           x = x_cut * 3.7, y = y_cap * 0.95,
           label = paste0("Dynamic\nn = ", comma(full_counts[filter_gate == "Dynamic", n]),
                          " (", full_counts[filter_gate == "Dynamic", pct], "%)"),
           color = "#D55E00", size = 3.1, hjust = 0, fontface = "bold") +

  scale_color_manual(values = c("Neither" = "#BDBDBD",
                                "Effect only" = "#0072B2",
                                "FDR only" = "#009E73",
                                "Dynamic" = "#D55E00"),
                     name = NULL,
                     guide = guide_legend(
                       override.aes = list(size = 3, alpha = 0.85))) +
  coord_cartesian(xlim = x_lim, ylim = c(0, y_cap)) +
  labs(x = expression(Delta * "logit(PDUI)"),
       y = expression(-log[10]("age q-value")),
       title = "Dynamic APA site filter") +
  theme_pub(base = 11) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 10, color = "#333333", face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    plot.title = element_text(size = 12.5, face = "bold", margin = margin(b = 6)),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, face = "bold")
  )

save_png(p02, "02_dynamic_site_filter_map.png", 7.0, 5.0)
save_pdf(p02, "02_dynamic_site_filter_map.pdf", 7.0, 5.0)

# --- layout validation ---

label_positions <- data.table(
  label = c("Effect only", "Neither", "FDR only", "Dynamic"),
  x = c(-x_cut * 3.7, x_lim[1] + 0.5, -x_cut * 1.7, x_cut * 3.7),
  y = c(y_cut + 4.5, y_cap * 0.22, y_cap * 0.58, y_cap * 0.95)
)
layout_validation <- data.table(
  panel = "02_dynamic_site_filter_map",
  n_total = nrow(sites),
  n_plotted = nrow(plot_dt),
  label_count = nrow(label_positions),
  legend_items = nlevels(plot_dt$filter_gate),
  label_clipped = any(
    label_positions$x < x_lim[1] | label_positions$x > x_lim[2] |
      label_positions$y < 0 | label_positions$y > y_cap
  ),
  label_axis_overlap = any(label_positions$y < y_cap * 0.05),
  annotation_background = "none",
  legend_inside_data = FALSE,
  categorical_tick_overlap = FALSE
)
fwrite(layout_validation, file.path(OUT_DIR, "02_dynamic_site_filter_map_layout_validation.tsv"), sep = "\t")
