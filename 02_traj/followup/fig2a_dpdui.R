#!/usr/bin/env Rscript
# Fig 2a: developmental dPDUI versus host-gene change per directional event, in residual-SD units
# Input is written by figures/producers/fig2a_prep_event_coupling.py
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrastr)
})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))

# --- load ---

d <- fread(file.path(Sys.getenv("APA_ROOT"), "data/02_traj/fig2a_event_coupling.tsv"))
d <- d[is.finite(apa_dev) & is.finite(ex_dev)]
d[, cl := factor(ifelse(cl == 1, 'C1', 'C2'), levels = c('C1', 'C2'))]
rr <- d[, .(n = .N, rho = round(cor(apa_dev, ex_dev, method = 'spearman'), 3)), by = cl]
lim <- ceiling(max(abs(c(d$apa_dev, d$ex_dev))))
lab <- rr[order(cl), .(cl, txt = paste0(as.character(cl), ' (n = ', format(n, big.mark = ','), ', rho = ',
                                 sprintf('%.3f', rho), ')'))]

# --- plot ---

p <- ggplot(d[sample(.N)], aes(apa_dev, ex_dev, colour = cl)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = NEUTRAL) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = NEUTRAL) +
  rasterise(geom_point(size = 0.35, alpha = 0.15, stroke = 0), dpi = 400) +
  annotate('text', x = -lim * 0.97, y = c(lim * 0.97, lim * 0.82), label = lab$txt, hjust = 0, vjust = 1,
           size = 5.5 / .pt, colour = c(C1_COLOR, C2_COLOR), family = FONT) +
  scale_colour_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR), guide = 'none') +
  scale_x_continuous(limits = c(-lim, lim), name = 'PDUI change (residual SD)') +
  scale_y_continuous(limits = c(-lim, lim), name = 'Host-gene change (residual SD)') +
  theme_nature()

# --- save ---

fig_dir <- file.path(FIG_DIR, "02_traj")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(fig_dir, "fig2a_dpdui_vs_drna.png"), p,
       width = 88, height = 65, units = "mm", dpi = 600, bg = "white")
