#!/usr/bin/env Rscript
# Fig 2d: APA-protein association beyond mRNA and age, 371 donors, observed vs trans-gene null.

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
mk <- function(m, lab, f) {
  r <- g(m)
  data.table(panel = lab, src = c("Observed", "Trans-gene null"), est = c(r$observed, r$null_mean),
             lo = c(NA, r$null_q05), hi = c(NA, r$null_q95), value_txt = c(f(r$observed), f(r$null_mean)))
}

# --- data ---
s <- fread(file.path(omics_dir, "S1_apa_prot_cis_vs_trans_371.tsv"))
d <- rbindlist(list(
  mk("n_q05", "Pairs (count)", function(x) if (abs(x - round(x)) < 1e-9) sprintf("%d", as.integer(round(x))) else sprintf("%.1f", x)),
  mk("pi1", "\\u03c0\\u2081", function(x) sprintf("%.3f", x)),
  mk("frac_neg", "Negative fraction", function(x) sprintf("%.1f%%", 100 * x))))
d[, panel := factor(panel, levels = c("Pairs (count)", "\\u03c0\\u2081", "Negative fraction"))]
d[, src := factor(src, levels = c("Observed", "Trans-gene null"))]
# observed 85 vs null 49.5; pi1 0.160 vs 0.118; 51.7% vs 49.4%

# --- plot ---
p <- ggplot(d, aes(src, est, colour = src, shape = src)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.3, na.rm = TRUE) +
  geom_point(size = 1.8) +
  geom_text(aes(label = value_txt), hjust = -0.25, vjust = -0.6, size = 5.5 / .pt, family = FONT, show.legend = FALSE) +
  facet_wrap(~panel, scales = "free_y", nrow = 1, strip.position = "top") +
  scale_colour_manual(values = c("Observed" = INK, "Trans-gene null" = NEUTRAL), guide = "none") +
  scale_shape_manual(values = c("Observed" = 16, "Trans-gene null" = 1), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)), name = "Value") +
  scale_x_discrete(name = NULL, labels = c("Observed" = "Obs.", "Trans-gene null" = "Null"), expand = expansion(add = c(0.4, 0.9))) +
  theme_nature() + theme(strip.text = element_text(size = 5.5))
ggsave(file.path(figure_dir, "fig2d_apa_protein_excess.png"), p,
       width = 88, height = 65, units = "mm", dpi = 600, bg = "white")
