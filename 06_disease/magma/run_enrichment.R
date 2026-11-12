#!/usr/bin/env Rscript
# GO and KEGG enrichment of EUR significant APA-QTL genes (genome-wide FDR < 0.05 per GWAS).
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")

# --- load packages ---

suppressMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
})

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
magma_dir <- file.path(apa_root, "data", "06_disease")
figure_dir <- file.path(apa_root, "figures", "06_disease")
dir.create(figure_dir, showWarnings = FALSE)

# --- load data ---

qtl_gene_table <- fread(file.path(magma_dir, "apa_qtl_gene_sets.genes.tsv"))
qtl_gene_table[, gene_id := as.character(gene_id)]

# --- extract significant genes per GWAS with genome-wide FDR correction ---

extract_sig_genes <- function(gwas_id, gene_set) {
  baseline <- fread(file.path(magma_dir, "enrichment", gwas_id, "baseline.genes.out"))
  baseline[, GENE := as.character(GENE)]
  baseline[, fdr := p.adjust(P, method = "BH")]
  
  genes_in_set <- qtl_gene_table[set_name == gene_set, gene_id]
  sig_genes <- baseline[GENE %in% genes_in_set & fdr < 0.05]
  
  # Convert Entrez to symbol
  gene_map <- qtl_gene_table[set_name == gene_set, .(gene_id, symbol)]
  sig_merged <- merge(sig_genes, gene_map, by.x = "GENE", by.y = "gene_id")
  
  return(sig_merged[, .(entrez = GENE, symbol, zstat = ZSTAT, p = P, fdr)])
}

# EUR_SCZ prenatal and postnatal
eur_scz_prenatal <- extract_sig_genes("EUR_SCZ", "prenatal")
eur_scz_postnatal <- extract_sig_genes("EUR_SCZ", "postnatal")

# EUR_BIP prenatal and postnatal
eur_bip_prenatal <- extract_sig_genes("EUR_BIP", "prenatal")
eur_bip_postnatal <- extract_sig_genes("EUR_BIP", "postnatal")

# --- combine prenatal and postnatal across GWAS ---

prenatal_genes <- unique(rbind(
  eur_scz_prenatal[, .(entrez, symbol)],
  eur_bip_prenatal[, .(entrez, symbol)]
))

postnatal_genes <- unique(rbind(
  eur_scz_postnatal[, .(entrez, symbol)],
  eur_bip_postnatal[, .(entrez, symbol)]
))

# --- run GO enrichment ---

run_go_enrichment <- function(gene_symbols, ont = "BP", label) {
  # Convert symbols to Entrez IDs
  gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_entrez) < 5) {
    return(NULL)
  }
  
  ego <- enrichGO(gene = gene_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  return(ego)
}

go_prenatal <- run_go_enrichment(prenatal_genes$symbol, "BP", "prenatal")
go_postnatal <- run_go_enrichment(postnatal_genes$symbol, "BP", "postnatal")

# --- run KEGG enrichment ---

run_kegg_enrichment <- function(gene_symbols, label) {
  gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_entrez) < 5) {
    return(NULL)
  }
  
  ekegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                      organism = 'hsa',
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
  
  return(ekegg)
}

kegg_prenatal <- run_kegg_enrichment(prenatal_genes$symbol, "prenatal")
kegg_postnatal <- run_kegg_enrichment(postnatal_genes$symbol, "postnatal")

# --- save results ---

save_enrichment_results <- function(result, name) {
  if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
    fwrite(as.data.frame(result),
           file.path(magma_dir, sprintf("enrichment_%s.tsv", name)),
           sep = "\t")
  }
}

save_enrichment_results(go_prenatal, "GO_BP_prenatal")
save_enrichment_results(go_postnatal, "GO_BP_postnatal")
save_enrichment_results(kegg_prenatal, "KEGG_prenatal")
save_enrichment_results(kegg_postnatal, "KEGG_postnatal")

# --- visualization ---

theme_enrichment <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text.y = element_text(color = "black", size = base_size),
      axis.text.x = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      legend.position = "right"
    )
}

# dot plot for top terms
plot_top_terms <- function(result, title, n = 15) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    return(NULL)
  }
  
  df <- as.data.frame(result)
  df <- head(df[order(df$p.adjust), ], n)
  
  # Scale gene ratio for plotting
  df$GeneRatio_num <- sapply(df$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  })
  
  p <- ggplot(df, aes(x = GeneRatio_num, y = reorder(Description, -p.adjust))) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_viridis_c(option = "plasma", direction = -1, name = "Adj. P") +
    scale_size_continuous(range = c(2, 6), name = "Gene Count") +
    labs(x = "Gene Ratio", y = NULL, title = title) +
    theme_enrichment()
  
  return(p)
}

# plots
p_go_pre <- plot_top_terms(go_prenatal, "GO Biological Process - Prenatal APA-QTL genes")
p_go_post <- plot_top_terms(go_postnatal, "GO Biological Process - Postnatal APA-QTL genes")
p_kegg_pre <- plot_top_terms(kegg_prenatal, "KEGG Pathways - Prenatal APA-QTL genes")
p_kegg_post <- plot_top_terms(kegg_postnatal, "KEGG Pathways - Postnatal APA-QTL genes")

# combine plots
if (!is.null(p_go_pre) && !is.null(p_go_post)) {
  p_go_combined <- p_go_pre / p_go_post +
    plot_annotation(
      title = "GO Enrichment: Prenatal vs Postnatal APA-QTL Genes",
      subtitle = "EUR GWAS, genome-wide FDR < 0.05",
      tag_levels = "A"
    )
  
  ggsave(file.path(figure_dir, "enrichment_GO_prenatal_vs_postnatal.pdf"),
         p_go_combined, width = 10, height = 12, device = cairo_pdf)
  ggsave(file.path(figure_dir, "enrichment_GO_prenatal_vs_postnatal.png"),
         p_go_combined, width = 10, height = 12, dpi = 300)
}

if (!is.null(p_kegg_pre) && !is.null(p_kegg_post)) {
  p_kegg_combined <- p_kegg_pre / p_kegg_post +
    plot_annotation(
      title = "KEGG Enrichment: Prenatal vs Postnatal APA-QTL Genes",
      subtitle = "EUR GWAS, genome-wide FDR < 0.05",
      tag_levels = "A"
    )
  
  ggsave(file.path(figure_dir, "enrichment_KEGG_prenatal_vs_postnatal.pdf"),
         p_kegg_combined, width = 10, height = 10, device = cairo_pdf)
  ggsave(file.path(figure_dir, "enrichment_KEGG_prenatal_vs_postnatal.png"),
         p_kegg_combined, width = 10, height = 10, dpi = 300)
}

# --- summary table ---

summary_dt <- data.table(
  group = c("Prenatal", "Postnatal"),
  n_genes = c(nrow(prenatal_genes), nrow(postnatal_genes)),
  go_terms = c(
    if (!is.null(go_prenatal)) nrow(as.data.frame(go_prenatal)) else 0,
    if (!is.null(go_postnatal)) nrow(as.data.frame(go_postnatal)) else 0
  ),
  kegg_terms = c(
    if (!is.null(kegg_prenatal)) nrow(as.data.frame(kegg_prenatal)) else 0,
    if (!is.null(kegg_postnatal)) nrow(as.data.frame(kegg_postnatal)) else 0
  )
)

fwrite(summary_dt, file.path(magma_dir, "enrichment_summary.tsv"), sep = "\t")
