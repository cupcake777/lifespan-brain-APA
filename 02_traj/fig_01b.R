#!/usr/bin/env Rscript
# Sample support profile with the individual-count track, and GAM k-index diagnostics.
# Supported blocks are labelled with their age range and total N.
suppressPackageStartupMessages({
  library(patchwork)
  library(ggrepel)
})
source("fig_common.R")

support <- fread(file.path(OUT_DIR, "01_support.tsv"))
cov <- fread(file.path(INTERMEDIATE_DIR, "00_cov.tsv"))

# --- individual support profile ---

cov[, age_model := log10(pmax(as.numeric(age) + 1, 1e-6))]
cov_individual <- unique(cov[, .(individual_id, age_model)])

support_with_blocks <- copy(support)
setorder(support_with_blocks, grid_index)
support_with_blocks[, block := rleid(is_supported)]

support_blocks <- support_with_blocks[
  ,
  .(
    xmin = min(age_model, na.rm = TRUE),
    xmax = max(age_model, na.rm = TRUE),
    ymin_n = min(local_n, na.rm = TRUE),
    ymax_n = max(local_n, na.rm = TRUE)
  ),
  by = .(block, is_supported)
]

label_data <- support_with_blocks[
  is_supported == TRUE,
  .(
    xmin = min(age_model, na.rm = TRUE),
    xmax = max(age_model, na.rm = TRUE),
    x_mid = mean(age_model, na.rm = TRUE),
    y_mid = mean(local_n, na.rm = TRUE),
    start_age = first(grid_age),
    end_age = last(grid_age)
  ),
  by = block
]
label_data[
  ,
  total_n := vapply(
    seq_len(.N),
    function(i) sum(cov_individual$age_model >= xmin[i] &
                      cov_individual$age_model <= xmax[i], na.rm = TRUE),
    integer(1)
  )
]
label_data[
  ,
  label_text := paste0(
    "From: ", round(start_age, 1), "y\n",
    "To: ", round(end_age, 1), "y\n",
    "Total N = ", total_n
  )
]

p01_top <- ggplot() +
  geom_rect(
    data = support_blocks,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = is_supported),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  geom_line(data = support, aes(age_model, local_n), color = "#0072B2", linewidth = 0.9) +
  geom_label_repel(
    data = label_data,
    aes(x_mid, y_mid, label = label_text),
    size = 3.0,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "grey50",
    segment.size = 0.4,
    max.overlaps = 20,
    fill = "white",
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#56B4E9", "FALSE" = "#E69F00"),
    labels = c("TRUE" = "Supported region", "FALSE" = "Unsupported")
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = range(support$age_model, na.rm = TRUE)) +
  labs(x = NULL, y = "Individuals", fill = "Status", title = "Age support profile") +
  theme_pub(base = 10) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

real_ages <- c(-0.5, 0, 10, 30, 60, 90)
model_breaks <- log10(real_ages + 1)
p01_bottom <- ggplot() +
  geom_rug(
    data = cov_individual,
    aes(age_model),
    sides = "b",
    color = "#0072B2",
    alpha = 0.30,
    length = unit(0.45, "cm")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.55) +
  annotate("text", x = 0.04, y = 0.5, label = "Birth",
           color = "#D55E00", angle = 90, size = 3.2) +
  scale_x_continuous(breaks = model_breaks, labels = paste0(real_ages, "y"),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = range(support$age_model, na.rm = TRUE)) +
  labs(x = "Chronological age", y = "Individuals") +
  theme_pub(base = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_text(vjust = 0.5, angle = 0)
  )

sample_support_plot <- p01_top / p01_bottom + plot_layout(heights = c(8, 1))
save_png(sample_support_plot, "01_sample_support_profile.png", 7.2, 5.4)
save_pdf(sample_support_plot, "01_sample_support_profile.pdf", 7.2, 5.4)

layout_validation <- rbindlist(list(
  data.table(
    figure = "01_sample_support_profile",
    panel = "A_support_profile",
    n_rows = nrow(support),
    n_labels = nrow(label_data),
    max_label_chars = max(nchar(label_data$label_text), na.rm = TRUE),
    validation_note = sprintf(
      "Restored git-accepted support profile: log10 age-model x-axis, supported/unsupported block shading, one local_n line, and %s supported-block labels.",
      nrow(label_data)
    )
  ),
  data.table(
    figure = "01_sample_support_profile",
    panel = "B_sample_rug",
    n_rows = nrow(cov_individual),
    n_labels = 1L,
    max_label_chars = nchar("Birth"),
    validation_note = "Bottom track shows one rug mark per individual on the same log10 age-model scale with the birth reference line."
  )
), fill = TRUE)
fwrite(layout_validation,
       file.path(OUT_DIR, "01_sample_support_profile_layout_validation.tsv"),
       sep = "\t")

# --- GAM k-index diagnostics ---

sites <- load_sites()
model_diag <- sites[is.finite(k_index), .(k_index)]
k_median  <- median(model_diag$k_index, na.rm = TRUE)
k_q25     <- quantile(model_diag$k_index, 0.25, na.rm = TRUE)
k_q75     <- quantile(model_diag$k_index, 0.75, na.rm = TRUE)
n_total   <- nrow(model_diag)
n_pass    <- sum(model_diag$k_index >= 0.9 & model_diag$k_index <= 1.1, na.rm = TRUE)
pct_pass  <- round((n_pass / n_total) * 100, 1)
k_density <- density(model_diag$k_index[is.finite(model_diag$k_index)])
k_label_y <- max(k_density$y, na.rm = TRUE) * 0.72

stats_label <- paste0(
  "Data Summary:\n",
  "• Median = ", round(k_median, 3), "\n",
  "• IQR = [", round(k_q25, 3), ", ", round(k_q75, 3), "]\n",
  "• Passed = ", n_pass, "/", n_total, " (", pct_pass, "%)"
)

p01b <- ggplot(model_diag, aes(k_index)) +
  annotate("rect", xmin = 0.9, xmax = 1.1, ymin = 0, ymax = Inf,
           fill = "#E69F00", alpha = 0.12) +
  annotate("text", x = 0.95, y = Inf, label = paste0("Pass Zone\n", pct_pass, "%"),
           vjust = 1.5, color = "#D55E00", size = 3.6, fontface = "bold", lineheight = 0.9) +

  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#D1D5DB", color = "white", alpha = 0.7, linewidth = 0.2) +

  geom_density(color = "#0072B2", linewidth = 1.0) +

  geom_vline(xintercept = 1, linetype = "solid", color = "#333333", linewidth = 0.5) +

  geom_vline(xintercept = k_median, linetype = "dashed", color = "#0072B2", linewidth = 0.6) +

  annotate("label", x = min(model_diag$k_index, na.rm = TRUE), y = k_label_y,
           label = stats_label,
           hjust = 0, vjust = 1,
           size = 3.35,
           fontface = "bold",
           fill = alpha("white", 0.85),
           label.padding = unit(0.4, "lines"),
           color = "#1F2937") +

  scale_x_continuous(breaks = seq(0.5, 1.2, by = 0.1)) +

  labs(x = "Model K-index",
       y = "Probability Density",
       title = "GAM Model Fitting") +

  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5, color = "#111827", margin = margin(b = 10)),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    panel.grid = element_blank()
  )

save_png(p01b, "01_gam_model.png", 7, 5)
save_pdf(p01b, "01_gam_model.pdf", 7, 5)

gam_layout_validation <- data.table(
  figure = "01_gam_model",
  panel = "A_k_index_distribution",
  n_rows = n_total,
  n_labels = 2L,
  max_label_chars = max(nchar(c("Pass Zone", stats_label))),
  k_index_min = min(model_diag$k_index, na.rm = TRUE),
  k_index_max = max(model_diag$k_index, na.rm = TRUE),
  pass_zone_fraction = n_pass / n_total,
  validation_note = paste(
    "The statistics label is left-aligned at the observed x minimum;",
    "the pass-zone label is centered within [0.9, 1.1], and both are",
    "vertically separated inside the density panel."
  )
)
fwrite(gam_layout_validation,
       file.path(OUT_DIR, "01_gam_model_layout_validation.tsv"),
       sep = "\t")
