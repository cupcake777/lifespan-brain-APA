#!/usr/bin/env Rscript
# Within-GWAS BH FDR, figures and top-gene tables for APA-QTL MAGMA gene-set enrichment.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")

# --- load packages ---

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
magma_dir <- file.path(apa_root, "data", "06_disease")
summary_path <- file.path(magma_dir, "magma_gsa_summary.tsv")
figure_dir <- file.path(apa_root, "figures", "06_disease")

out_fig1_pdf <- file.path(figure_dir, "magma_apa_qtl_enrichment_within_gwas_figure1.pdf")
out_fig1_png <- file.path(figure_dir, "magma_apa_qtl_enrichment_within_gwas_figure1.png")
out_fig2_pdf <- file.path(figure_dir, "magma_apa_qtl_enrichment_within_gwas_figure2_top_genes.pdf")
out_fig2_png <- file.path(figure_dir, "magma_apa_qtl_enrichment_within_gwas_figure2_top_genes.png")
out_hits <- file.path(magma_dir, "magma_apa_qtl_enrichment_within_gwas_significant.tsv")
out_notes <- file.path(magma_dir, "magma_apa_qtl_enrichment_within_gwas_notes.tsv")
out_top_genes <- file.path(magma_dir, "magma_apa_qtl_enrichment_within_gwas_top_genes.tsv")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

set_order <- c("prenatal", "postnatal", "interaction")
set_labels <- c(
  prenatal = "Prenatal",
  postnatal = "Postnatal",
  interaction = "Interaction"
)

set_colors <- c(
  prenatal = "#0072B2",
  postnatal = "#009E73",
  interaction = "#D55E00"
)

population_colors <- c(
  EAS = "#7A7A7A",
  EUR = "#222222"
)

# --- helpers ---

theme_magma <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_blank(),
      plot.caption = element_text(hjust = 0, color = "grey25", size = base_size - 1),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank()
    )
}

# --- load enrichment results ---

x <- fread(summary_path)
x[, p := as.numeric(p)]
x[, beta := as.numeric(beta)]
x[, beta_std := as.numeric(beta_std)]
x[, se := as.numeric(se)]
x[, ngenes := as.integer(ngenes)]
x[, set_name := factor(set_name, levels = set_order)]
x[, fdr_within_gwas := as.numeric(fdr_within_gwas)]
x[, fdr_within_gwas := p.adjust(p, method = "BH"), by = gwas_id]
x[, significant_within_gwas := fdr_within_gwas < 0.05]
if ("fdr_global" %in% names(x)) x[, fdr_global := NULL]
stopifnot(nrow(x) == 57L, all(x[, .N, by = gwas_id]$N == 3L),
          !anyDuplicated(x[, .(gwas_id, set_name)]))
x[, trait_id := gsub("^(EAS|EUR)_", "", gwas_id)]
x[, trait_axis := paste0(population, "_", trait_id)]
x[, trait_label := trait_id]
x[, neg_log10_p := -log10(p)]
x[, hit_label := ifelse(significant_within_gwas, trait_axis, "")]

fwrite(x[order(population, trait_id, set_name)], summary_path, sep = "\t")

trait_order <- x[, .(min_p = min(p)), by = .(population, trait_axis)][order(population, min_p)]$trait_axis
x[, trait_axis := factor(trait_axis, levels = trait_order)]

fwrite(
  x[significant_within_gwas == TRUE][order(p),
    .(gwas_id, population, trait_id, set_name, ngenes, beta, beta_std, se, p, fdr_within_gwas, source_gsa)
  ],
  out_hits,
  sep = "\t"
)

fwrite(
  x[, .(gwas_id, population, trait_id, set_name, p, fdr_within_gwas, significant_within_gwas)][order(population, trait_id, set_name)],
  out_notes,
  sep = "\t"
)

# --- figure data: significant hits and full screen ---

hits <- x[significant_within_gwas == TRUE]

threshold_y <- -log10(0.05)

heatmap_data <- copy(x)
heatmap_data[, p_label := fifelse(significant_within_gwas, "*", "")]
heatmap_data[, star_note := " * within-GWAS FDR < 0.05"]

# --- figure 2 data: example top MAGMA genes ---

qtl_gene_table <- fread(file.path(magma_dir, "apa_qtl_gene_sets.genes.tsv"))
qtl_gene_table[, gene_id := as.character(gene_id)]

top_gene_rows <- list()

for (i in seq_len(nrow(hits))) {
  hit <- hits[i]
  baseline_path <- file.path(magma_dir, "enrichment", hit$gwas_id, "baseline.genes.out")
  baseline_genes <- fread(baseline_path)
  baseline_genes[, GENE := as.character(GENE)]

  gene_block <- merge(
    baseline_genes,
    qtl_gene_table[set_name == as.character(hit$set_name), .(gene_id, symbol)],
    by.x = "GENE",
    by.y = "gene_id",
    all = FALSE
  )
  setnames(gene_block, c("GENE", "ZSTAT", "P"), c("gene_id", "zstat", "gene_p"))
  gene_block[, gwas_id := hit$gwas_id]
  gene_block[, trait_id := hit$trait_id]
  gene_block[, population := hit$population]
  gene_block[, set_name := as.character(hit$set_name)]
  gene_block[, hit_group := paste0(hit$trait_axis, " ", set_labels[as.character(hit$set_name)])]
  top_gene_rows[[i]] <- gene_block[zstat > 0][order(gene_p)][1:3]
}

top_genes <- rbindlist(top_gene_rows, use.names = TRUE)
top_genes[, gene_label := paste0(symbol, "  ")]
top_hit_levels <- hits[order(p), paste0(trait_axis, " ", set_labels[as.character(set_name)])]
top_genes[, hit_group := factor(hit_group, levels = top_hit_levels)]

fwrite(
  top_genes[, .(gwas_id, population, trait_id, set_name, gene_id, symbol, zstat, gene_p, CHR, START, STOP, NSNPS, NPARAM)],
  out_top_genes,
  sep = "\t"
)

# --- figure 1A: significant within-GWAS hits ---

p_hits <- ggplot(x, aes(x = beta_std, y = neg_log10_p)) +
  geom_hline(yintercept = threshold_y, linewidth = 0.35, linetype = "dashed", color = "grey45") +
  geom_point(aes(fill = set_name, color = population, size = ngenes, alpha = significant_within_gwas), shape = 21, stroke = 0.45) +
  scale_fill_manual(values = set_colors, labels = set_labels) +
  scale_color_manual(values = population_colors, guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.32), guide = "none") +
  scale_size_continuous(range = c(1.8, 6.2), breaks = c(250, 500, 1000, 1800)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.13))) +
  labs(
    title = "MAGMA APA-QTL enrichment",
    x = "Standardized MAGMA set beta",
    y = expression(-log[10](P))
  ) +
  theme_magma()

# --- figure 1B: all GWAS-by-set results ---

p_heatmap <- ggplot(heatmap_data, aes(x = trait_axis, y = set_name)) +
  geom_tile(aes(fill = neg_log10_p), color = "white", linewidth = 0.45) +
  geom_text(aes(label = p_label), size = 3.2, fontface = "bold", color = "black") +
  facet_grid(. ~ population, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = setNames(x$trait_label, x$trait_axis)) +
  scale_y_discrete(labels = set_labels) +
  scale_fill_viridis_c(option = "magma", begin = 0.1, end = 0.92, name = expression(-log[10](P))) +
  labs(x = "GWAS", y = NULL) +
  theme_magma() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    legend.position = "right",
    plot.title = element_blank()
  )

# --- figure 2: starred-set example genes ---

p_genes <- ggplot(top_genes, aes(x = zstat, y = reorder(symbol, zstat))) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(aes(fill = set_name), shape = 21, size = 2.3, stroke = 0.35, color = "black") +
  facet_wrap(~ hit_group, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = set_colors, labels = set_labels) +
  labs(
    x = "MAGMA gene Z",
    y = NULL,
    title = "Top genes in starred sets"
  ) +
  theme_magma(base_size = 8) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", hjust = 0),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0, size = 10)
  )

p_star_note <- ggplot(data.table(x = 0.02, y = 0.65, label = "* within-GWAS FDR < 0.05")) +
  geom_text(aes(x = x, y = y, label = label), hjust = 0, size = 3.2) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# --- save figures ---

figure1 <- p_hits / p_heatmap / p_star_note
figure1 <- figure1 + plot_layout(heights = c(1.45, 1, 0.12))

figure2 <- p_genes

ggsave(out_fig1_pdf, figure1, width = 7.3, height = 6.35, device = cairo_pdf)
ggsave(out_fig1_png, figure1, width = 7.3, height = 6.35, dpi = 450)
ggsave(out_fig2_pdf, figure2, width = 3.8, height = 5.45, device = cairo_pdf)
ggsave(out_fig2_png, figure2, width = 3.8, height = 5.45, dpi = 450)
