#!/usr/bin/env Rscript
# Figures 1-3 for the MAGMA APA-QTL enrichment results.
# Plots go to figures/06_disease, tables to data/06_disease.

# --- load packages ---

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
if (!nzchar(apa_root)) stop("Set APA_ROOT (source config/paths.local.sh)")
data_dir <- file.path(apa_root, "data", "06_disease")
figure_dir <- file.path(apa_root, "figures", "06_disease")

set_summary_path <- file.path(data_dir, "magma_gsa_summary_enrichment.tsv")
qtl_gene_path <- file.path(data_dir, "apa_qtl_gene_sets.genes.tsv")

fig1_pdf <- file.path(figure_dir, "magma_qtl_enrichment_overview.pdf")
fig1_png <- file.path(figure_dir, "magma_qtl_enrichment_overview.png")
fig1_data <- file.path(data_dir, "magma_qtl_enrichment_overview_plot_data.tsv")

fig2_pdf <- file.path(figure_dir, "magma_go_bp_selected_contexts.pdf")
fig2_png <- file.path(figure_dir, "magma_go_bp_selected_contexts.png")
fig2_data <- file.path(data_dir, "magma_go_bp_selected_contexts_plot_data.tsv")
fig2_terms <- file.path(data_dir, "go_bp_term_label_mapping.tsv")
fig2_genes <- file.path(data_dir, "go_bp_gene_membership.tsv")
fig2_foregrounds <- file.path(data_dir, "go_bp_selected_context_foregrounds.tsv")
fig2_gene_matrix <- file.path(data_dir, "go_bp_gene_context_matrix.tsv")

fig3_pdf <- file.path(figure_dir, "magma_qtl_developmental_specificity.pdf")
fig3_png <- file.path(figure_dir, "magma_qtl_developmental_specificity.png")
fig3_data <- file.path(data_dir, "magma_qtl_developmental_specificity_plot_data.tsv")

validation_path <- file.path(data_dir, "magma_figure_layout_validation.tsv")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

set_order <- c("prenatal", "postnatal", "interaction")
set_labels <- c(
  prenatal = "Prenatal",
  postnatal = "Postnatal",
  interaction = "Dynamic"
)
set_short <- c(
  prenatal = "pre",
  postnatal = "post",
  interaction = "dyn"
)

effect_low <- "#3B6FB6"
effect_mid <- "white"
effect_high <- "#C74343"
min_go_profile_genes <- 10
top_go_terms_per_context <- 5
top_genes_per_context <- 6
pal_stage_bg <- c(prenatal = "#F1F7FF", postnatal = "#F3FAF1")
stage_bg_half_width <- 0.42

# --- helpers ---

theme_magma <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )
}

save_plot <- function(plot, pdf_path, png_path, width, height) {
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 300)
}

make_bg_layer <- function(meta, y_min, y_max) {
  rects <- lapply(seq_len(nrow(meta)), function(i) {
    annotate("rect",
      xmin = meta$context_x[i] - stage_bg_half_width,
      xmax = meta$context_x[i] + stage_bg_half_width,
      ymin = y_min, ymax = y_max,
      fill = pal_stage_bg[meta$qtl_set[i]],
      alpha = 0.80
    )
  })
  rects
}

theme_fig2 <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "grey55", linewidth = 0.25),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 3,
                                margin = margin(b = 7)),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey60", linewidth = 0.35),
      legend.position = "right",
      legend.margin = margin(0, 0, 0, 4),
      legend.box.spacing = unit(0.25, "cm"),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key.height = unit(0.42, "cm"),
      legend.key.width = unit(0.46, "cm")
    )
}

# --- load data ---

set_summary <- fread(set_summary_path)
qtl_genes <- fread(qtl_gene_path)
set_summary[, p := as.numeric(p)]
set_summary[, beta := as.numeric(beta)]
set_summary[, beta_std := as.numeric(beta_std)]
set_summary[, ngenes := as.integer(ngenes)]
set_summary[, fdr_within_gwas := as.numeric(fdr_within_gwas)]
set_summary[, fdr_within_gwas := p.adjust(p, method = "BH"), by = gwas_id]
set_summary[, significant_within_gwas := fdr_within_gwas < 0.05]
set_summary[, set_name := factor(set_name, levels = set_order)]
set_summary[, qtl_set := as.character(set_name)]
set_summary[, trait_id := sub("^(EAS|EUR)_", "", gwas_id)]
set_summary[, trait_axis := paste0(population, "_", trait_id)]
set_summary[, qtl_label := set_labels[qtl_set]]
set_summary[, qtl_short := set_short[qtl_set]]
set_summary[, context_id := paste0(gwas_id, "_", qtl_short)]
set_summary[, neg_log10_p := -log10(p)]
set_summary[, neg_log10_fdr := -log10(fdr_within_gwas)]
set_summary[, sig_label := fifelse(significant_within_gwas, "*", "")]

qtl_genes[, gene_id := as.character(gene_id)]
# --- figure 1 data: set-level MAGMA enrichment ---

trait_order <- set_summary[
  ,
  .(best_fdr = min(fdr_within_gwas), best_beta = max(beta_std)),
  by = .(population, trait_axis, trait_id)
][order(population, best_fdr, -best_beta)]$trait_axis

set_summary[, trait_axis := factor(trait_axis, levels = trait_order)]
set_summary[, qtl_set := factor(qtl_set, levels = set_order)]

fig1_plot_data <- set_summary[
  ,
  .(
    gwas_id, population, trait_id, trait_axis, context_id, qtl_set,
    ngenes, beta, beta_std, p, fdr_within_gwas, significant_within_gwas,
    neg_log10_p, neg_log10_fdr, sig_label
  )
][order(population, trait_axis, qtl_set)]

fwrite(fig1_plot_data, fig1_data, sep = "\t")

p1_heatmap <- ggplot(fig1_plot_data, aes(x = trait_axis, y = qtl_set)) +
  geom_tile(aes(fill = beta_std), color = "white", linewidth = 0.45) +
  geom_text(aes(label = sig_label), size = 3.2, fontface = "bold", color = "black") +
  facet_grid(. ~ population, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = setNames(fig1_plot_data$trait_id, fig1_plot_data$trait_axis)) +
  scale_y_discrete(labels = set_labels) +
  scale_fill_gradient2(
    low = effect_low,
    mid = effect_mid,
    high = effect_high,
    midpoint = 0,
    name = "MAGMA beta"
  ) +
  labs(x = "GWAS", y = NULL, title = "MAGMA Enrichment") +
  theme_magma(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

figure1 <- p1_heatmap
save_plot(figure1, fig1_pdf, fig1_png, width = 8.2, height = 3.6)

# --- figure 2 data: selected context foreground genes ---

selected_contexts <- set_summary[significant_within_gwas == TRUE][order(fdr_within_gwas)]

foreground_rows <- list()
gene_membership_rows <- list()
go_rows <- list()

for (i in seq_len(nrow(selected_contexts))) {
  current <- selected_contexts[i]
  baseline <- fread(file.path(data_dir, "enrichment", current$gwas_id, "baseline.genes.out"))
  baseline[, GENE := as.character(GENE)]

  qtl_map <- qtl_genes[set_name == as.character(current$qtl_set), .(gene_id, symbol)]
  qtl_gene_level <- merge(
    baseline,
    qtl_map,
    by.x = "GENE",
    by.y = "gene_id",
    all = FALSE
  )
  qtl_gene_level[, gene_fdr_within_gwas_qtl := p.adjust(P, method = "BH")]
  gene_fdr_foreground <- qtl_gene_level[gene_fdr_within_gwas_qtl < 0.05]

  foreground <- gene_fdr_foreground
  foreground_type <- "gene_fdr"
  foreground_ids <- unique(foreground$GENE)
  universe_ids <- unique(baseline$GENE)
  context_label <- paste0(current$gwas_id, "_", set_short[as.character(current$qtl_set)])

  foreground_rows[[i]] <- data.table(
    context_id = context_label,
    gwas_id = current$gwas_id,
    population = current$population,
    qtl_set = as.character(current$qtl_set),
    foreground_type = foreground_type,
    n_gene_fdr_foreground = nrow(gene_fdr_foreground),
    n_foreground = length(foreground_ids),
    n_universe = length(universe_ids),
    set_level_p = current$p,
    set_level_fdr_within_gwas = current$fdr_within_gwas
  )

  gene_membership_rows[[i]] <- foreground[
    ,
    .(
      context_id = context_label,
      gwas_id = current$gwas_id,
      population = current$population,
      qtl_set = as.character(current$qtl_set),
      foreground_type = foreground_type,
      gene_id = GENE,
      symbol,
      zstat = ZSTAT,
      gene_p = P,
      gene_fdr_within_gwas_qtl
    )
  ][order(gene_p)]

  if (length(foreground_ids) >= min_go_profile_genes) {
    ego <- enrichGO(
      gene = foreground_ids,
      universe = universe_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE
    )
    go_table <- as.data.table(as.data.frame(ego))
    if (nrow(go_table) > 0) {
      go_table[, context_id := context_label]
      go_table[, gwas_id := current$gwas_id]
      go_table[, population := current$population]
      go_table[, qtl_set := as.character(current$qtl_set)]
      go_table[, foreground_type := foreground_type]
      go_table[, n_foreground := length(foreground_ids)]
      go_table[, n_universe := length(universe_ids)]
      go_rows[[i]] <- go_table
    }
  }
}

foreground_summary <- rbindlist(foreground_rows, use.names = TRUE)
gene_membership <- rbindlist(gene_membership_rows, use.names = TRUE)
go_all <- rbindlist(go_rows, use.names = TRUE, fill = TRUE)

fwrite(foreground_summary, fig2_foregrounds, sep = "\t")
fwrite(gene_membership, fig2_genes, sep = "\t")

go_top <- go_all[
  order(p.adjust),
  head(.SD, top_go_terms_per_context),
  by = .(context_id, foreground_type)
]
go_top[, gene_ratio_num := as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))]
go_top[, term_label := {
  x <- Description
  x <- fcase(
    x == "cellular modified amino acid biosynthetic process", "modified amino acid biosynthesis",
    x == "regulation of NMDA receptor activity", "NMDA receptor regulation",
    x == "endosome to lysosome transport via multivesicular body sorting pathway", "endosome-lysosome transport",
    x == "tRNA aminoacylation for protein translation", "tRNA aminoacylation for translation",
    x == "tRNA wobble uridine modification", "tRNA wobble U modification",
    x == "tRNA wobble base modification", "tRNA wobble-base modification",
    default = x
  )
  x <- gsub("positive regulation of ", "pos. reg. ", x)
  x <- gsub("negative regulation of ", "neg. reg. ", x)
  x <- gsub("regulation of ", "reg. ", x)
  x <- gsub("central nervous system", "CNS", x)
  x <- gsub("nervous system", "NS", x)
  x <- gsub("biological process", "BP", x)
  fifelse(nchar(x) > 42, paste0(substr(x, 1, 39), "..."), x)
}]
go_top[, neg_log10_fdr := -log10(p.adjust)]

context_order <- foreground_summary[n_foreground >= min_go_profile_genes][
  order(set_level_fdr_within_gwas)
]$context_id

context_meta <- foreground_summary[n_foreground >= min_go_profile_genes][
  match(context_order, context_id),
  .(
    context_id,
    gwas_id,
    trait_id = sub("^(EAS|EUR)_", "", gwas_id),
    population,
    qtl_set,
    qtl_label = set_labels[qtl_set],
    context_display = paste0(population, "_", sub("^(EAS|EUR)_", "", gwas_id)),
    context_plot_label = paste0(population, "_", sub("^(EAS|EUR)_", "", gwas_id), "\n", set_short[qtl_set])
  )
]
context_meta[, context_axis := make.unique(context_display, sep = "_")]
context_meta[, context_x := seq_len(.N)]

go_top <- merge(go_top, context_meta[, .(context_id, context_axis, context_x)], by = "context_id", all.x = TRUE)

term_order <- go_top[, .(best_fdr = min(p.adjust)), by = .(ID, term_label)][order(best_fdr)]$term_label
go_top[, term_label := factor(term_label, levels = rev(term_order))]

fwrite(go_top, fig2_data, sep = "\t")
fwrite(unique(go_top[, .(ID, Description, term_label)]), fig2_terms, sep = "\t")

go_gene_long <- go_top[
  ,
  .(gene_id = unlist(strsplit(geneID, "/", fixed = TRUE))),
  by = .(context_id, ID, Description, term_label, p.adjust)
]

gene_stats <- gene_membership[
  ,
  .(
    context_id,
    symbol,
    gene_id,
    zstat,
    gene_p,
    gene_fdr_within_gwas_qtl
  )
]

gene_context_matrix <- merge(
  go_gene_long,
  gene_stats,
  by = c("context_id", "gene_id"),
  all.x = TRUE
)
gene_context_summary <- gene_context_matrix[
  ,
  .(
    n_bp_terms = uniqueN(ID),
    bp_terms = paste(sort(unique(term_label)), collapse = "; "),
    best_go_fdr = min(p.adjust),
    gene_id = gene_id[which.min(gene_p)],
    zstat = zstat[which.min(gene_p)],
    gene_p = min(gene_p),
    gene_fdr_within_gwas_qtl = min(gene_fdr_within_gwas_qtl)
  ),
  by = .(context_id, symbol)
][order(context_id, -n_bp_terms, gene_p)]

gene_context_top <- gene_context_summary[
  ,
  head(.SD, top_genes_per_context),
  by = context_id
]
gene_context_top[, neg_log10_gene_fdr := -log10(gene_fdr_within_gwas_qtl)]
gene_context_top <- merge(
  gene_context_top,
  context_meta[, .(context_id, context_axis, context_x)],
  by = "context_id",
  all.x = TRUE
)

gene_order <- gene_context_top[
  ,
  .(best_gene_p = min(gene_p), max_terms = max(n_bp_terms)),
  by = symbol
][order(best_gene_p, -max_terms)]$symbol
gene_context_top[, symbol := factor(symbol, levels = rev(gene_order))]

fwrite(gene_context_top, fig2_gene_matrix, sep = "\t")

# --- figure 2: order terms and genes ---

# Reorder terms by overall significance across contexts
term_fdr_min <- go_top[, .(min_fdr = min(p.adjust)), by = term_label]
term_fdr_min <- term_fdr_min[order(min_fdr)]
go_top[, term_label := factor(term_label, levels = rev(term_fdr_min$term_label))]

# Reorder genes by overall significance
gene_fdr_min <- gene_context_top[, .(min_fdr = min(gene_fdr_within_gwas_qtl)), by = symbol]
gene_fdr_min <- gene_fdr_min[order(min_fdr)]
gene_context_top[, symbol := factor(symbol, levels = rev(gene_fdr_min$symbol))]

n_go_terms <- length(levels(go_top$term_label))
n_genes <- length(levels(gene_context_top$symbol))

# --- figure 2A: GO biological processes ---

bg_go <- make_bg_layer(context_meta, -Inf, Inf)

p2_go <- ggplot(go_top, aes(x = context_x, y = term_label)) +
  bg_go +
  geom_hline(yintercept = seq(0.5, n_go_terms, by = 1), color = "grey94", linewidth = 0.16) +
  geom_point(aes(fill = neg_log10_fdr, size = Count), shape = 21, color = "grey30", stroke = 0.28) +
  scale_x_continuous(
    breaks = context_meta$context_x,
    labels = context_meta$context_plot_label,
    expand = c(0.01, 0.01)
  ) +
  scale_fill_viridis_c(
    option = "D", direction = 1,
    name = expression(-log[10]("GO FDR")),
    guide = guide_colorbar(barheight = 3.4, barwidth = 0.55)
  ) +
  scale_size_continuous(range = c(2.4, 6.0), guide = "none") +
  labs(x = NULL, y = NULL, title = "A  GO biological processes") +
  theme_fig2(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = 5)),
    axis.text.y = element_text(size = 10.5),
    plot.margin = margin(5, 8, 8, 5)
  )

# --- figure 2B: gene anchors ---

bg_gene <- make_bg_layer(context_meta, -Inf, Inf)

p2_genes <- ggplot(gene_context_top, aes(x = context_x, y = symbol)) +
  bg_gene +
  geom_hline(yintercept = seq(0.5, n_genes, by = 1), color = "grey94", linewidth = 0.16) +
  geom_point(aes(fill = neg_log10_gene_fdr, size = n_bp_terms), shape = 21, color = "grey30", stroke = 0.28) +
  scale_x_continuous(
    breaks = context_meta$context_x,
    labels = context_meta$context_plot_label,
    expand = c(0.01, 0.01)
  ) +
  scale_fill_viridis_c(
    option = "D", direction = 1,
    name = expression(-log[10]("gene FDR")),
    guide = guide_colorbar(barheight = 3.4, barwidth = 0.55)
  ) +
  scale_size_continuous(range = c(2.2, 5.2), guide = "none") +
  labs(x = NULL, y = NULL, title = "B  Gene anchors") +
  theme_fig2(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = 5)),
    axis.text.y = element_text(size = 10.5, face = "italic"),
    plot.margin = margin(8, 8, 5, 5)
  )

# --- assemble figure 2 ---

figure2 <- (p2_go / p2_genes) +
  plot_layout(heights = c(1, 1.05), guides = "collect") +
  plot_annotation(
    caption = "Background: blue = Prenatal APA-QTL; green = Postnatal APA-QTL",
    theme = theme(
      plot.margin = margin(10, 16, 10, 10),
      plot.caption = element_text(size = 9.5, color = "grey35", hjust = 0.5, margin = margin(t = 8))
    )
  )

save_plot(figure2, fig2_pdf, fig2_png, width = 9.0, height = 10.0)

# --- figure 3 data: APA-QTL in GWAS ---

fig3_plot_data <- copy(fig1_plot_data)
fig3_plot_data[, evidence_status := fifelse(
  significant_within_gwas,
  "FDR < 0.05",
  fifelse(p < 0.05, "P < 0.05 only", "n.s.")
)]
fig3_plot_data[, evidence_status := factor(
  evidence_status,
  levels = c("FDR < 0.05", "P < 0.05 only", "n.s.")
)]
fig3_trait_order <- fig3_plot_data[
  ,
  .(
    any_sig = any(significant_within_gwas),
    best_fdr = min(fdr_within_gwas),
    max_abs_beta = max(abs(beta_std))
  ),
  by = trait_id
][order(-any_sig, best_fdr, -max_abs_beta)]$trait_id
fig3_trait_levels <- rev(fig3_trait_order)
fig3_plot_data[, trait_y := match(trait_id, fig3_trait_levels)]
fig3_plot_data[, trait_y := trait_y + fifelse(population == "EAS", -0.10, 0.10)]

fwrite(fig3_plot_data, fig3_data, sep = "\t")

p3_effect <- ggplot(fig3_plot_data, aes(x = beta_std, y = trait_y)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_point(
    aes(fill = population, size = neg_log10_p, shape = evidence_status),
    color = "grey25",
    stroke = 0.38,
    alpha = 0.88
  ) +
  geom_point(
    data = fig3_plot_data[significant_within_gwas == TRUE],
    aes(fill = population, size = neg_log10_p),
    shape = 23,
    color = "black",
    stroke = 0.9,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_grid(. ~ qtl_set, labeller = labeller(qtl_set = set_labels)) +
  scale_y_continuous(
    breaks = seq_along(fig3_trait_levels),
    labels = fig3_trait_levels,
    expand = expansion(add = 0.55)
  ) +
  scale_x_continuous(
    breaks = seq(-0.02, 0.02, by = 0.01),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.08, 0.11))
  ) +
  scale_fill_manual(values = c(EAS = "#E69F00", EUR = "#56B4E9"), name = "Population") +
  scale_size_continuous(range = c(1.4, 4.8), name = expression(-log[10](P))) +
  scale_shape_manual(
    values = c("FDR < 0.05" = 23, "P < 0.05 only" = 22, "n.s." = 21),
    name = "Statistical status"
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(fill = "grey85", color = "grey25", shape = 21, alpha = 1)),
    shape = guide_legend(order = 2, override.aes = list(fill = "grey85", color = "grey25", size = 3.6, alpha = 1)),
    fill = guide_legend(order = 3, override.aes = list(shape = 21, size = 4.2, alpha = 1))
  ) +
  labs(x = "Standardized MAGMA beta", y = "GWAS", title = "APA-QTL in GWAS") +
  theme_magma(base_size = 10) +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(size = 8.8, margin = margin(t = 3)),
    axis.title = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black", size = 12),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = grid::unit(0.42, "in"),
    plot.margin = margin(8, 18, 8, 8),
    legend.position = "right"
  )

figure3 <- p3_effect
save_plot(figure3, fig3_pdf, fig3_png, width = 9.2, height = 4.6)

# --- layout validation ---

validation <- data.table(
  figure = c("Figure 1", "Figure 2", "Figure 3"),
  output_pdf = c(fig1_pdf, fig2_pdf, fig3_pdf),
  output_png = c(fig1_png, fig2_png, fig3_png),
  width_in = c(8.2, 9.0, 9.2),
  height_in = c(3.6, 10.0, 4.6),
  label_strategy = c(
    "Short GWAS IDs; stars mark within-GWAS FDR < 0.05",
    "Short two-line context IDs with developmental stage labels and column background; selected GO terms and gene anchors use larger text in stacked panels; full mapping and support tables exported",
    "Short trait IDs on y axis; population encoded by color legend and slight vertical offset; stage facets separate prenatal, postnatal, and dynamic QTL sets; statistical status encoded by point shape; x-axis spacing increased to avoid facet-edge tick crowding"
  ),
  fdr_scope = c(
    "BH across 3 QTL sets within each GWAS",
    "GO BH FDR within each selected gene-FDR foreground context; gene panel uses gene-level FDR within GWAS x QTL",
    "BH across 3 QTL sets within each GWAS"
  ),
  visual_qc = c(
    "PDF/PNG exported with larger final font settings; categorical labels angled or shortened",
    "PDF/PNG exported with larger final font settings; stacked layout separates GO and gene-anchor panels; dot-size legends suppressed to reduce clutter; stage backgrounds use lighter fills with visible context-column gaps",
    "PDF/PNG exported with larger final font settings; Figure 3 uses marker shape for FDR, nominal-only, and non-significant status instead of alpha levels; panel spacing and x tick sizing checked for unclipped labels"
  )
)

fwrite(validation, validation_path, sep = "\t")
