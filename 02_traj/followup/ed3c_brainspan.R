#!/usr/bin/env Rscript
# BrainSpan PFC concordance of C1/C2 events (external age-association replication).
# Reads the 12_external_c1c2 tables from data/02_traj/.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))

# --- load ---
data_dir <- file.path(Sys.getenv("APA_ROOT"), "data", "02_traj")
tp <- fread(file.path(data_dir, "12_external_c1c2_traj_plot.tsv"))
ev <- fread(file.path(data_dir, "12_external_c1c2_event_replication.tsv"))
ev[, cls := factor(ifelse(discovery_cluster_id == 1, "C1", "C2"), c("C1", "C2"))]
lab <- ev[, .(n = .N, nq = sum(external_age_q < 0.05, na.rm = TRUE)), by = cls]
lab[, txt := sprintf("%d events; q < 0.05 for %d", n, nq)]
lab[, `:=`(x = -Inf, y = Inf)]
tp[, cls := factor(ifelse(discovery_cluster_id == 1, "C1", "C2"), c("C1", "C2"))]

# --- plot ---
cols <- c(C1 = C1_COLOR, C2 = C2_COLOR)
p <- ggplot(tp, aes(age, colour = cls, fill = cls)) +
  geom_ribbon(aes(ymin = external_ci_low, ymax = external_ci_high), alpha = 0.18, colour = NA) +
  geom_line(aes(y = external_z), linewidth = 0.5) +
  geom_line(aes(y = discovery_z), linewidth = 0.3, linetype = "22", colour = INK) +
  geom_text(data = lab, aes(x = x, y = y, label = txt), inherit.aes = FALSE, hjust = -0.05, vjust = 1.6,
            size = 5.5 / .pt, colour = INK) +
  facet_wrap(~cls, nrow = 1) +
  scale_colour_manual(values = cols, guide = "none") + scale_fill_manual(values = cols, guide = "none") +
  labs(x = "Age (years, BrainSpan PFC)", y = "Standardized change (z)")

# --- save ---
fig_dir <- file.path(FIG_DIR, "02_traj")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(fig_dir, "ed3c_brainspan_concordance.png"), p,
       width = 88, height = 65, units = "mm", dpi = 600, bg = "white")