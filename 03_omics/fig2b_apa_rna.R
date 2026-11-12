#!/usr/bin/env Rscript
# Fig 2b: fraction of APA-RNA pairs positive and negative at FDR < 0.05, observed vs trans-gene null.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))

# --- config ---
omics_dir <- file.path(Sys.getenv("APA_ROOT"), "data", "03_omics")
figure_dir <- file.path(Sys.getenv("APA_ROOT"), "figures", "03_omics")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# --- helpers ---
g <- function(m) s[metric == m]

# --- data ---
s <- fread(file.path(omics_dir, "S1_apa_rna_cis_vs_trans_371.tsv"))
pos <- g("frac_pos_q05"); neg <- g("frac_neg_q05")
d <- rbindlist(list(
  data.table(sign = "Positive", src = "Observed", est = pos$observed, lo = NA_real_, hi = NA_real_),
  data.table(sign = "Positive", src = "Trans-gene null", est = pos$null_mean, lo = pos$null_q05, hi = pos$null_q95),
  data.table(sign = "Negative", src = "Observed", est = neg$observed, lo = NA_real_, hi = NA_real_),
  data.table(sign = "Negative", src = "Trans-gene null", est = neg$null_mean, lo = neg$null_q05, hi = neg$null_q95)))
d[, sign := factor(sign, levels = c("Positive", "Negative"))]
d[, src := factor(src, levels = c("Observed", "Trans-gene null"))]

# --- plot ---
p <- ggplot(d, aes(sign, est, colour = src, shape = src, group = src)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.3,
                position = position_dodge(0.45), na.rm = TRUE) +
  geom_point(aes(fill = src), size = 1.6, position = position_dodge(0.45)) +
  annotate("text", x = 2.4, y = 0.335, label = "pi[1]~0.67~vs~0.64~(null)", parse = TRUE,
           hjust = 1, size = 5.5 / .pt, colour = INK, family = FONT) +
  scale_colour_manual(values = c("Observed" = INK, "Trans-gene null" = NEUTRAL), name = NULL) +
  scale_shape_manual(values = c("Observed" = 21, "Trans-gene null" = 21), name = NULL) +
  scale_fill_manual(values = c("Observed" = INK, "Trans-gene null" = "white"), name = NULL) +
  scale_y_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05), name = "Fraction of APA-mRNA pairs") +
  scale_x_discrete(name = NULL) +
  theme_nature() + theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "fig2b_apa_rna_pairs.png"), p,
       width = 58, height = 48, units = "mm", dpi = 600, bg = "white")
