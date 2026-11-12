# Adjusted LOEUF constraint across APA direction classes, with pairwise contrasts.
suppressPackageStartupMessages({
  library(patchwork)
})
source("fig_common.R")

# --- config ---

GENE_AGG_FILE <- file.path(OUT_DIR, "05_gene_aggregation.tsv")
LOEUF_FILE <- file.path(OUT_DIR, "ed4c_loeuf.tsv")
LOEUF_MODEL_FILE <- file.path(OUT_DIR, "11_loeuf_adj_models.tsv")

direction_levels <- c("C1-spec", "C2-spec",
                      "C1-main", "C2-main", "Balanced")
direction_colors <- c("C1-spec" = "#56B4E9",
                      "C2-spec" = "#E69F00",
                      "C1-main" = "#0072B2",
                      "C2-main" = "#D55E00",
                      "Balanced" = "#7A7A7A")

group_band_colors <- c("Specific" = "#EAF7F2",
                       "Shared" = "#FBEAF3")

# --- loading data ---

gene_agg <- fread(GENE_AGG_FILE)
loeuf_raw <- fread(LOEUF_FILE, select = c("gene_symbol", "loeuf"))
loeuf_models <- fread(LOEUF_MODEL_FILE)

gene_agg[, direction_class := fcase(
  n_c1 > 0 & n_c2 == 0, "C1-spec",
  n_c2 > 0 & n_c1 == 0, "C2-spec",
  n_c1 > n_c2, "C1-main",
  n_c2 > n_c1, "C2-main",
  n_c1 == n_c2, "Balanced"
)]
gene_agg[, direction_class := factor(direction_class, levels = direction_levels)]

loeuf_raw <- unique(loeuf_raw[is.finite(loeuf)], by = "gene_symbol")
loeuf_df <- merge(gene_agg[, .(gene_symbol, direction_class)],
                  loeuf_raw, by = "gene_symbol")
loeuf_df[, direction_class := factor(direction_class, levels = direction_levels)]

# --- figure 09 data: LOEUF median and IQR ---

loeuf_summary <- loeuf_df[, .(
  n_genes = .N,
  median_loeuf = median(loeuf),
  q25_loeuf = quantile(loeuf, 0.25),
  q75_loeuf = quantile(loeuf, 0.75)
), by = direction_class]
loeuf_summary <- loeuf_summary[order(direction_class)]
loeuf_summary[, y := rev(seq_len(.N))]
loeuf_summary[, y_label := sprintf("%s (n=%s)", direction_class, comma(n_genes))]
loeuf_summary[, value_label := sprintf("%.2f", median_loeuf)]

group_bands <- data.table(
  band = c("Specific", "Shared"),
  ymin = c(3.5, 0.5),
  ymax = c(5.5, 3.5),
  fill = c(group_band_colors["Specific"], group_band_colors["Shared"])
)

pair_list <- combn(direction_levels, 2, simplify = FALSE)
loeuf_tests <- rbindlist(lapply(pair_list, function(pair) {
  x <- loeuf_df[direction_class == pair[1], loeuf]
  y <- loeuf_df[direction_class == pair[2], loeuf]
  data.table(
    comparison = paste(pair[1], "vs", pair[2]),
    n_a = length(x),
    median_a_raw_loeuf = median(x),
    n_b = length(y),
    median_b_raw_loeuf = median(y),
    pvalue = wilcox.test(x, y)$p.value
  )
}))
loeuf_tests[, fdr := p.adjust(pvalue, "BH")]
loeuf_tests[, delta_median := median_a_raw_loeuf - median_b_raw_loeuf]
loeuf_tests[, sig := fdr < FDR_THR]
loeuf_tests[, star := fcase(
  fdr < 0.001, "***",
  fdr < 0.01, "**",
  fdr < 0.05, "*",
  default = "ns"
)]
loeuf_tests[, contrast := factor(comparison, levels = rev(comparison[order(fdr)]))]
loeuf_tests[, star_x := delta_median + fifelse(delta_median >= 0, 0.012, -0.012)]
loeuf_tests[, star_hjust := fifelse(delta_median >= 0, 0, 1)]

kw_test <- kruskal.test(loeuf ~ direction_class, data = loeuf_df)
kw_label <- ifelse(kw_test$p.value < 0.001, "Kruskal-Wallis P < 0.001",
                   sprintf("Kruskal-Wallis P = %.3f", kw_test$p.value))
adj_row <- loeuf_models[model == "membership" & term == "membershipShared"]
adj_label <- if (nrow(adj_row) == 1L) {
  sprintf("Adjusted Shared beta = %.3f, P = %.2g", adj_row$estimate, adj_row$pvalue)
} else {
  "Adjusted Shared effect unavailable"
}

# --- figure 09A: LOEUF median with IQR ---

p_summary <- ggplot() +
  geom_rect(data = group_bands,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = band),
            alpha = 0.65, inherit.aes = FALSE) +
  geom_vline(xintercept = 1, color = "#8A8A8A", linetype = "dashed", linewidth = 0.35) +
  geom_segment(data = loeuf_summary,
               aes(x = q25_loeuf, xend = q75_loeuf, y = y, yend = y,
                   color = direction_class),
               linewidth = 2.2, lineend = "round") +
  geom_point(data = loeuf_summary,
             aes(x = median_loeuf, y = y, fill = direction_class),
             shape = 21, size = 3.0, color = "black", stroke = 0.35) +
  geom_text(data = loeuf_summary,
            aes(x = 1.125, y = y, label = value_label),
            hjust = 0, size = 3.35, color = "#111827", fontface = "bold") +
  geom_text(data = group_bands,
            aes(x = 0.38, y = ymax - 0.20, label = band),
            hjust = 0, size = 3.45, color = "#374151",
            fontface = "bold", inherit.aes = FALSE) +
  scale_y_continuous(breaks = loeuf_summary$y,
                     labels = loeuf_summary$y_label,
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_x_continuous(breaks = seq(0.4, 1.1, 0.1),
                     expand = expansion(mult = c(0.01, 0.12))) +
  scale_color_manual(values = direction_colors, guide = "none") +
  scale_fill_manual(values = c(group_band_colors, direction_colors), guide = "none") +
  coord_cartesian(xlim = c(0.34, 1.24), clip = "off") +
  labs(x = "LOEUF median with IQR", y = NULL,
       title = "Constraint by APA direction class",
       caption = paste(kw_label, adj_label, sep = "; ")) +
  theme_pub(base = 11) +
  theme(panel.grid.major.x = element_line(color = "#E5E7EB", linewidth = 0.35),
        axis.text.y = element_text(size = 10, face = "bold"),
        plot.caption = element_text(size = 8.8, hjust = 0, face = "bold"),
        plot.margin = margin(5.5, 12, 5.5, 5.5))

# --- figure 09B: pairwise median LOEUF contrasts ---

p_contrast <- ggplot(loeuf_tests, aes(delta_median, contrast, color = sig)) +
  geom_vline(xintercept = 0, color = "#8A8A8A", linetype = "dashed", linewidth = 0.32) +
  geom_segment(aes(x = 0, xend = delta_median, y = contrast, yend = contrast),
               color = "#C7C7C7", linewidth = 0.42) +
  geom_point(size = 2.25) +
  geom_text(aes(x = star_x, y = contrast, label = star, hjust = star_hjust),
            inherit.aes = FALSE, size = 3.15, fontface = "bold", color = "black") +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#666666"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.20, 0.24))) +
  labs(x = "Median LOEUF difference", y = NULL,
       title = "Pairwise contrasts",
       caption = "* FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001") +
  theme_pub(base = 11) +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12.5, face = "bold"),
        plot.caption = element_text(size = 8.8, hjust = 0, face = "bold"),
        plot.margin = margin(5.5, 8, 5.5, 5.5))

# --- assemble figure 09 ---

p09 <- p_summary | p_contrast +
  plot_layout(widths = c(1.30, 1.08)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_png(p09, "09_loeuf_constraint.png", 11.2, 4.8)
save_pdf(p09, "09_loeuf_constraint.pdf", 11.2, 4.8)

# --- layout validation ---

layout_validation <- data.table(
  figure = "09_loeuf_constraint",
  panel = c("A_LOEUF_median_IQR_lollipop", "B_pairwise_median_delta"),
  n_rows = c(nrow(loeuf_df), nrow(loeuf_tests)),
  n_categories_or_contrasts = c(uniqueN(loeuf_df$direction_class),
                                uniqueN(loeuf_tests$comparison)),
  max_label_chars = c(max(nchar(loeuf_summary$y_label)),
                      max(nchar(loeuf_tests$comparison))),
  validation_note = c(
    "Horizontal median + IQR lollipop by five spec/main direction-composition classes; no significance labels on the summary panel.",
    "Right-side panel shows all pairwise Wilcoxon median LOEUF differences with FDR stars."
  )
)
fwrite(layout_validation,
       file.path(OUT_DIR, "09_loeuf_constraint_layout_validation.tsv"), sep = "\t")

summary_annotation <- loeuf_summary[, .(
  direction_class,
  n_genes,
  median_loeuf,
  q25_loeuf,
  q75_loeuf
)]
fwrite(summary_annotation,
       file.path(OUT_DIR, "09_loeuf_constraint_sum.tsv"), sep = "\t")

contrast_annotation <- loeuf_tests[, .(
  comparison,
  n_a,
  median_a_raw_loeuf,
  n_b,
  median_b_raw_loeuf,
  delta_median,
  pvalue,
  fdr,
  star
)]
fwrite(contrast_annotation,
       file.path(OUT_DIR, "09_loeuf_constraint_contrasts.tsv"), sep = "\t")
