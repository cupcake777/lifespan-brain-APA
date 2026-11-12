# Cross-layer coupling (APA, transcript, gene) for the converged clusters,
# with a neutral class for |rho| <= 0.10.
suppressPackageStartupMessages({
  library(ggalluvial)
  library(patchwork)
})
source("fig_common.R")

# --- helpers ---

direction_call <- function(x) {
  fcase(x > 0.10, "pos", x < -0.10, "neg", default = "neutral")
}

# --- config ---

OUT_PNG <- "06_crosslayer.png"
OUT_PDF <- "06_crosslayer.pdf"
OUT_DATA <- file.path(OUT_DIR, "06_crosslayer_plot.tsv")
OUT_LAYOUT <- file.path(OUT_DIR, "06_crosslayer_layout_validation.tsv")

C_COLORS <- c("C1" = "#0072B2", "C2" = "#D55E00")
REL_COLORS <- c("APA-transcript" = "#0072B2",
                "Transcript-gene" = "#8C8C8C",
                "APA-gene" = "#D55E00")
REL_LEVELS <- c("APA-transcript", "Transcript-gene", "APA-gene")
DIR_COLORS <- c("pos" = "#0072B2",
                "neutral" = "#999999",
                "neg" = "#D55E00")

set.seed(7L)

# --- loading data ---

sync <- fread(file.path(OUT_DIR, "06_cross_layer_trend.tsv"))
sync[, cluster_label := paste0("C", cluster_id)]
sync[, cluster_label := factor(cluster_label, levels = c("C1", "C2"))]

# --- figure 06A data: cross-layer rho distributions ---

rho_long <- melt(sync,
                 id.vars = c("site_id", "gene_symbol", "cluster_label"),
                 measure.vars = c("apa_tx_spearman",
                                  "tx_gene_spearman",
                                  "apa_gene_spearman"),
                 variable.name = "relationship",
                 value.name = "rho")
rho_long <- rho_long[is.finite(rho)]
rho_long[, relationship := factor(
  relationship,
  levels = c("apa_tx_spearman", "tx_gene_spearman", "apa_gene_spearman"),
  labels = REL_LEVELS
)]

rho_summary <- rho_long[, .(
  n = .N,
  median_rho = median(rho),
  q25_rho = quantile(rho, 0.25),
  q75_rho = quantile(rho, 0.75)
), by = .(cluster_label, relationship)]

# --- figure 06B data: three-layer direction flow ---

rho_wide <- sync[
  is.finite(apa_tx_spearman) &
    is.finite(tx_gene_spearman) &
    is.finite(apa_gene_spearman)
]
rho_wide[, attenuation := apa_tx_spearman - apa_gene_spearman]

flow_dt <- rho_wide[, .(
  site_id,
  gene_symbol,
  cluster_label,
  apa_tx_dir = direction_call(apa_tx_spearman),
  tx_gene_dir = direction_call(tx_gene_spearman),
  apa_gene_dir = direction_call(apa_gene_spearman),
  attenuation
)]
flow_dt[, `:=`(
  apa_tx_dir = factor(apa_tx_dir, levels = c("pos", "neutral", "neg")),
  tx_gene_dir = factor(tx_gene_dir, levels = c("pos", "neutral", "neg")),
  apa_gene_dir = factor(apa_gene_dir, levels = c("pos", "neutral", "neg"))
)]

flow_summary <- flow_dt[, .(
  n = .N,
  median_attenuation = median(attenuation)
), by = .(cluster_label, apa_tx_dir, tx_gene_dir, apa_gene_dir)]
flow_summary[, fraction := n / sum(n), by = cluster_label]
flow_summary[, path_label := paste(apa_tx_dir, tx_gene_dir, apa_gene_dir, sep = " -> ")]

# --- figure 06C data: representative examples ---

example_sites <- rho_wide[
  apa_tx_spearman > 0.45 &
    tx_gene_spearman > 0.45 &
    apa_gene_spearman < 0.10 &
    !grepl("^ENSG", gene_symbol) &
    abs(apa_tx_spearman) < 0.99 &
    abs(tx_gene_spearman) < 0.99 &
    abs(apa_gene_spearman) < 0.99
]
example_sites[, attenuation := apa_tx_spearman - apa_gene_spearman]
example_sites[, example_score := attenuation + pmin(apa_tx_spearman, tx_gene_spearman)]
example_sites <- example_sites[order(-example_score)]
example_sites <- unique(example_sites, by = "gene_symbol")
example_sites <- example_sites[seq_len(min(.N, 6L))]
example_sites[, gene_label := factor(gene_symbol, levels = rev(gene_symbol))]

example_long <- melt(example_sites,
                     id.vars = c("site_id", "gene_symbol", "gene_label", "cluster_label"),
                     measure.vars = c("apa_tx_spearman",
                                      "tx_gene_spearman",
                                      "apa_gene_spearman"),
                     variable.name = "relationship",
                     value.name = "rho")
example_long[, relationship := factor(
  relationship,
  levels = c("apa_tx_spearman", "tx_gene_spearman", "apa_gene_spearman"),
  labels = REL_LEVELS
)]

# --- figure 06A: pairwise rho landscape ---

p_dist <- ggplot(rho_long, aes(relationship, rho, fill = cluster_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#555555", linewidth = 0.35) +
  geom_violin(position = position_dodge(width = 0.74),
              width = 0.72, alpha = 0.28, color = NA, scale = "width") +
  geom_boxplot(position = position_dodge(width = 0.74),
               width = 0.16, outlier.shape = NA, linewidth = 0.34,
               alpha = 0.78) +
  scale_fill_manual(values = C_COLORS, name = NULL) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = NULL, y = "Spearman rho", title = "Cross-layer correlation") +
  theme_pub(base = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "top",
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12.5, face = "bold"))

# --- figure 06B: cross-layer direction flow ---

p_flow <- ggplot(flow_summary,
                 aes(y = fraction,
                     axis1 = apa_tx_dir,
                     axis2 = tx_gene_dir,
                     axis3 = apa_gene_dir)) +
  geom_alluvium(aes(fill = apa_gene_dir),
                width = 0.14, alpha = 0.58, color = "white", linewidth = 0.08,
                discern = TRUE) +
  geom_stratum(width = 0.18, fill = "grey96", color = "grey45", linewidth = 0.24,
               discern = TRUE) +
  facet_wrap(~ cluster_label, nrow = 1) +
  scale_x_discrete(limits = REL_LEVELS,
                   labels = c("APA\nTx", "Tx\nGene", "APA\nGene"),
                   expand = c(0.09, 0.05)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_fill_manual(values = DIR_COLORS, name = "APA-gene dir") +
  labs(x = NULL, y = "Complete events",
       title = "Three-layer direction flow") +
  theme_pub(base = 11) +
  theme(strip.text = element_text(size = 10.5, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        panel.spacing.x = unit(0.85, "lines"),
        plot.title = element_text(size = 12.5, face = "bold"))

# --- figure 06C: representative cross-layer examples ---

p_examples <- ggplot(example_long, aes(rho, gene_label, fill = relationship)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#555555", linewidth = 0.32) +
  geom_col(position = position_dodge(width = 0.72),
           width = 0.58, color = "white", linewidth = 0.22) +
  scale_fill_manual(values = REL_COLORS, name = NULL) +
  coord_cartesian(xlim = c(-0.45, 0.85)) +
  labs(x = "Spearman rho", y = NULL, title = "Representative examples") +
  theme_pub(base = 11) +
  theme(legend.position = "top",
        legend.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10.2, face = "italic"),
        plot.title = element_text(size = 12.5, face = "bold"))

# --- assemble figure 06 ---

p06 <- (p_dist | p_flow | p_examples) +
  plot_layout(widths = c(1.00, 1.34, 1.05), guides = "keep")

save_png(p06, OUT_PNG, 13.8, 4.6)
ggsave(file.path(PANEL_DIR, OUT_PDF), p06, width = 13.8, height = 4.6, bg = "white")

# --- save figure 06 data ---

plot_data <- rbindlist(list(
  rho_summary[, .(panel = "A_rho_distribution",
                  cluster_label,
                  relationship,
                  label = NA_character_,
                  n,
                  x = NA_real_,
                  y = median_rho,
                  q25 = q25_rho,
                  q75 = q75_rho)],
  flow_summary[, .(panel = "B_direction_flow",
                 cluster_label,
                 relationship = "APA-transcript_to_Transcript-gene_to_APA-gene",
                 label = path_label,
                 n,
                 x = fraction,
                 y = median_attenuation,
                 q25 = NA_real_,
                 q75 = NA_real_)],
  example_long[, .(panel = "C_examples",
                   cluster_label,
                   relationship,
                   label = gene_symbol,
                   n = NA_integer_,
                   x = rho,
                   y = NA_real_,
                   q25 = NA_real_,
                   q75 = NA_real_)]
), fill = TRUE)
fwrite(plot_data, OUT_DATA, sep = "\t")

layout_validation <- data.table(
  figure = "06_crosslayer",
  panel = c("A_rho_distribution", "B_direction_flow", "C_examples"),
  n_rows = c(nrow(rho_long), nrow(flow_summary), nrow(example_long)),
  n_labels = c(uniqueN(rho_long$relationship),
               uniqueN(flow_summary$path_label),
               uniqueN(example_long$gene_symbol)),
  max_label_chars = c(max(nchar(as.character(rho_long$relationship))),
                      max(nchar(flow_summary$path_label)),
                      max(nchar(as.character(example_long$gene_symbol)))),
  validation_note = c(
    "No explanatory annotation layers; distribution panel uses violin plus boxplot.",
    "Panel B uses positive, neutral (|rho| <= 0.10), and negative states; unsupported ages do not enter rho.",
    "No explanatory annotation layers; representative examples use grouped horizontal rho bars."
  )
)
fwrite(layout_validation, OUT_LAYOUT, sep = "\t")
