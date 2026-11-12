#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")

# Gene-level MAGMA summary: BH FDR < 0.05 within each GWAS x QTL-set family.
# GO enrichment per QTL set, using that set's unique genes as background.

# --- load packages ---

library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
magma_dir <- file.path(apa_root, "data", "06_disease")
figure_dir <- file.path(apa_root, "figures", "06_disease")

out_gene_table <- file.path(magma_dir, "magma_gene_level_groups.tsv")
out_group_counts <- file.path(magma_dir, "magma_gene_level_gwas_qtl_counts.tsv")
out_stage_counts <- file.path(magma_dir, "magma_gene_level_stage_counts.tsv")
out_go_prenatal <- file.path(magma_dir, "magma_gene_level_GO_prenatal.tsv")
out_go_postnatal <- file.path(magma_dir, "magma_gene_level_GO_postnatal.tsv")
out_go_interaction <- file.path(magma_dir, "magma_gene_level_GO_interaction.tsv")
out_pdf <- file.path(figure_dir, "magma_gene_level_go_enrichment.pdf")
out_png <- file.path(figure_dir, "magma_gene_level_go_enrichment.png")
out_validation <- file.path(magma_dir, "magma_gene_level_go_enrichment_layout_validation.tsv")

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
go_min_gs_size <- 10
go_max_gs_size <- 500

# --- helpers ---

theme_gene <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black")
    )
}

go_table <- function(entrez_ids, background_ids) {
  ego <- enrichGO(
    gene = entrez_ids,
    universe = background_ids,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = go_min_gs_size,
    maxGSSize = go_max_gs_size,
    readable = TRUE
  )
  as.data.table(as.data.frame(ego))
}

# --- compute gene-level families per GWAS x QTL set ---

manifest <- fread(file.path(magma_dir, "gwas_manifest.tsv"))
qtl_gene_table <- fread(file.path(magma_dir, "apa_qtl_gene_sets.genes.tsv"))
qtl_gene_table[, gene_id := as.character(gene_id)]

all_gene_rows <- list()
group_counts <- list()

for (i in seq_len(nrow(manifest))) {
  gwas_id <- manifest$gwas_id[i]
  population <- manifest$ancestry[i]
  baseline <- fread(file.path(magma_dir, "enrichment", gwas_id, "baseline.genes.out"))
  baseline[, GENE := as.character(GENE)]

  for (qtl_set in set_order) {
    current_gwas <- gwas_id
    current_population <- population
    current_set <- qtl_set
    gene_map <- qtl_gene_table[set_name == qtl_set, .(gene_id, symbol)]
    group_genes <- merge(
      baseline[GENE %in% gene_map$gene_id],
      gene_map,
      by.x = "GENE",
      by.y = "gene_id",
      all = FALSE
    )
    group_genes[, gene_fdr_within_gwas_qtl := p.adjust(P, method = "BH")]
    group_genes[, significant_gene := gene_fdr_within_gwas_qtl < 0.05]
    group_genes[, gwas_id := current_gwas]
    group_genes[, population := current_population]
    group_genes[, qtl_set := current_set]

    all_gene_rows[[length(all_gene_rows) + 1]] <- group_genes[
      significant_gene == TRUE,
      .(gwas_id, population, qtl_set, gene_id = GENE, symbol, zstat = ZSTAT, gene_p = P, gene_fdr_within_gwas_qtl)
    ]

    group_counts[[length(group_counts) + 1]] <- data.table(
      gwas_id = gwas_id,
      population = population,
      qtl_set = qtl_set,
      n_tested_genes = nrow(group_genes),
      n_significant_genes = sum(group_genes$significant_gene)
    )
  }
}

sig_gene_table <- rbindlist(all_gene_rows, use.names = TRUE)
group_count_table <- rbindlist(group_counts, use.names = TRUE)

stage_gene_sets <- sig_gene_table[, .(n_unique_genes = uniqueN(gene_id)), by = qtl_set]
stage_gene_sets[, qtl_set := factor(qtl_set, levels = set_order)]

prenatal_genes <- unique(sig_gene_table[qtl_set == "prenatal", gene_id])
postnatal_genes <- unique(sig_gene_table[qtl_set == "postnatal", gene_id])
interaction_genes <- unique(sig_gene_table[qtl_set == "interaction", gene_id])

prenatal_background <- unique(qtl_gene_table[set_name == "prenatal", gene_id])
postnatal_background <- unique(qtl_gene_table[set_name == "postnatal", gene_id])
interaction_background <- unique(qtl_gene_table[set_name == "interaction", gene_id])

fwrite(sig_gene_table[order(gwas_id, qtl_set, gene_fdr_within_gwas_qtl)], out_gene_table, sep = "\t")
fwrite(group_count_table[order(population, gwas_id, qtl_set)], out_group_counts, sep = "\t")
fwrite(stage_gene_sets[order(qtl_set)], out_stage_counts, sep = "\t")

# --- run GO enrichment by QTL set ---

go_prenatal <- go_table(prenatal_genes, prenatal_background)
go_postnatal <- go_table(postnatal_genes, postnatal_background)
go_interaction <- go_table(interaction_genes, interaction_background)

fwrite(go_prenatal, out_go_prenatal, sep = "\t")
fwrite(go_postnatal, out_go_postnatal, sep = "\t")
fwrite(go_interaction, out_go_interaction, sep = "\t")

# --- figure 1A: significant genes per family ---

group_count_table[, qtl_set := factor(qtl_set, levels = set_order)]
group_count_table[, trait_id := gsub("^(EAS|EUR)_", "", gwas_id)]
trait_order <- group_count_table[, .(max_n = max(n_significant_genes)), by = .(population, gwas_id, trait_id)][order(population, -max_n)]$gwas_id
group_count_table[, gwas_id := factor(gwas_id, levels = trait_order)]
group_count_table[, label_color := ifelse(n_significant_genes >= 180, "black", "white")]

p_heatmap <- ggplot(group_count_table, aes(x = gwas_id, y = qtl_set, fill = n_significant_genes)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = n_significant_genes, color = label_color), size = 2.5) +
  facet_grid(. ~ population, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = setNames(group_count_table$trait_id, group_count_table$gwas_id)) +
  scale_y_discrete(labels = set_labels) +
  scale_color_identity() +
  scale_fill_viridis_c(option = "mako", name = "Genes") +
  labs(x = "GWAS", y = NULL, title = "A  Significant genes per GWAS x QTL set") +
  theme_gene(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank()
  )

# --- figure 1B: unique genes by QTL set ---

p_counts <- ggplot(stage_gene_sets, aes(x = qtl_set, y = n_unique_genes, fill = qtl_set)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.2) +
  geom_text(aes(label = n_unique_genes), vjust = -0.35, size = 3) +
  scale_x_discrete(labels = set_labels) +
  scale_fill_manual(values = set_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Unique significant genes", title = "B  Union across GWAS") +
  theme_gene(base_size = 8)

# --- figure 1C: GO enrichment by QTL set ---

go_tables <- list(
  "Prenatal genes" = go_prenatal,
  "Postnatal genes" = go_postnatal,
  "Interaction genes" = go_interaction
)
go_tables <- go_tables[vapply(go_tables, nrow, integer(1)) > 0]
go_panel_letters <- LETTERS[seq_along(go_tables) + 2]
go_panels <- list()
for (k in seq_along(go_tables)) {
  z <- head(go_tables[[k]][order(p.adjust)], 8)
  z[, gene_ratio_num := as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))]
  z[, term_label := fifelse(nchar(Description) > 58, paste0(substr(Description, 1, 55), "..."), Description)]

  go_panels[[k]] <- ggplot(z, aes(x = gene_ratio_num, y = reorder(term_label, -log10(p.adjust)))) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_viridis_c(option = "plasma", direction = -1, name = "Adj. P") +
    scale_size_continuous(range = c(2.2, 6.2), name = "Genes") +
    labs(x = "Gene ratio", y = NULL, title = paste(go_panel_letters[k], names(go_tables)[k])) +
    theme_gene(base_size = 8) +
    theme(axis.text.y = element_text(size = 7))
}

# --- save figure ---

figure <- wrap_plots(c(list(p_heatmap, p_counts), go_panels), ncol = 1) +
  plot_layout(heights = c(1.2, 0.55, rep(0.95, length(go_panels))))

figure_height <- 5.2 + 2.6 * length(go_panels)

ggsave(out_pdf, figure, width = 7.4, height = figure_height, device = cairo_pdf)
ggsave(out_png, figure, width = 7.4, height = figure_height, dpi = 450)

# --- layout validation ---

validation <- data.table(
  item = c(
    "gwas_qtl_gene_families",
    "families_with_significant_genes",
    "total_family_significant_gene_rows",
    "prenatal_unique_genes",
    "postnatal_unique_genes",
    "interaction_unique_genes",
    "prenatal_go_input_genes",
    "postnatal_go_input_genes",
    "interaction_go_input_genes",
    "prenatal_go_background_genes",
    "postnatal_go_background_genes",
    "interaction_go_background_genes",
    "enrichgo_min_gs_size",
    "enrichgo_max_gs_size",
    "prenatal_go_terms",
    "postnatal_go_terms",
    "interaction_go_terms",
    "go_panels_shown",
    "figure_pdf_exists",
    "figure_png_exists",
    "figure_width_in",
    "figure_height_in",
    "gene_level_definition"
  ),
  value = c(
    nrow(group_count_table),
    sum(group_count_table$n_significant_genes > 0),
    nrow(sig_gene_table),
    length(prenatal_genes),
    length(postnatal_genes),
    length(interaction_genes),
    length(prenatal_genes),
    length(postnatal_genes),
    length(interaction_genes),
    length(prenatal_background),
    length(postnatal_background),
    length(interaction_background),
    go_min_gs_size,
    go_max_gs_size,
    nrow(go_prenatal),
    nrow(go_postnatal),
    nrow(go_interaction),
    length(go_panels),
    file.exists(out_pdf),
    file.exists(out_png),
    7.4,
    figure_height,
    "Each GWAS x QTL set is corrected internally across genes in that QTL set using BH FDR < 0.05; GO uses the unique genes from the corresponding input APA-QTL set as background"
  )
)

fwrite(validation, out_validation, sep = "\t")
