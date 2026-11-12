# Over-representation of APA-linked protein candidates against the measured primary background.
# GO BP/MF, KEGG (local snapshot) and optional GenAge/CellAge sets; each custom TERM2GENE gets a background sentinel.

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# --- config ---

if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
mod <- file.path(Sys.getenv("APA_ROOT"), "data/03_omics")

minimum_gene_set_size <- 10L
maximum_gene_set_size <- 500L

# --- helpers ---

format_enrichment <- function(result, database) {
  table <- as.data.table(as.data.frame(result))
  if (nrow(table) == 0L) {
    return(data.table())
  }
  table[, database := database]
  table
}

# --- freeze candidate and measured background genes ---

universe <- fread(file.path(mod, "01_analysis_universe.tsv"))
candidates <- fread(file.path(mod, "03_integrated_gene_cands.tsv"))
gene_set_manifest <- fread(file.path(mod, "enrich_gene_set_manifest.tsv"))

background_symbols <- sort(unique(universe[
  primary_universe == TRUE & !is.na(gene_symbol), gene_symbol
]))
candidate_symbols <- sort(unique(candidates[
  apa_linked_protein_candidate == TRUE, gene_symbol
]))

background_mapping <- bitr(
  background_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
candidate_mapping <- if (length(candidate_symbols) == 0L) {
  data.frame(SYMBOL = character(), ENTREZID = character())
} else {
  bitr(
    candidate_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
}

background_entrez <- unique(background_mapping$ENTREZID)
candidate_entrez <- intersect(
  unique(candidate_mapping$ENTREZID),
  background_entrez
)

# --- GO BP and MF over-representation analysis ---

go_bp <- NULL
go_mf <- NULL
if (length(candidate_entrez) > 0L) {
  go_bp <- enrichGO(
    gene = candidate_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = minimum_gene_set_size,
    maxGSSize = maximum_gene_set_size,
    readable = FALSE
  )
  go_mf <- enrichGO(
    gene = candidate_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = minimum_gene_set_size,
    maxGSSize = maximum_gene_set_size,
    readable = FALSE
  )
}

# --- KEGG over-representation analysis from the local snapshot ---

kegg_link <- fread(
  file.path(mod, "kegg_hsa_gene_pathway.tsv"),
  header = FALSE,
  col.names = c("gene", "term")
)
kegg_link[, `:=`(
  gene = sub("^hsa:", "", gene),
  term = sub("^path:", "", term)
)]
kegg_name <- fread(
  file.path(mod, "kegg_hsa_pathways.tsv"),
  header = FALSE,
  col.names = c("term", "name")
)
kegg_name[, name := sub(" - Homo sapiens \\(human\\)$", "", name)]

kegg <- NULL
if (length(candidate_entrez) > 0L) {
  kegg <- enricher(
    gene = candidate_entrez,
    universe = background_entrez,
    TERM2GENE = kegg_link[, .(term, gene)],
    TERM2NAME = kegg_name,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = minimum_gene_set_size,
    maxGSSize = maximum_gene_set_size
  )
}

# --- optional GenAge and CellAge collections ---

aging_gene_set_path <- file.path(mod, "aging_gene_sets.tsv")
aging_results <- data.table()
aging_sets <- data.table(
  collection = character(), term = character(), gene_symbol = character()
)
if (file.exists(aging_gene_set_path) && length(candidate_symbols) > 0L) {
  aging_sets <- fread(aging_gene_set_path)
  aging_sets <- aging_sets[
    collection %chin% c("GenAge", "CellAge") & gene_symbol %chin% background_symbols
  ]
  aging_results <- rbindlist(lapply(
    unique(aging_sets$collection),
    function(collection_name) {
      collection_sets <- aging_sets[collection == collection_name]
      term2gene <- rbind(
        collection_sets[, .(term, gene_symbol)],
        data.table(
          term = "__measured_background__",
          gene_symbol = background_symbols
        )
      )
      result <- enricher(
        gene = candidate_symbols,
        universe = background_symbols,
        TERM2GENE = term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = minimum_gene_set_size,
        maxGSSize = maximum_gene_set_size
      )
      format_enrichment(result, collection_name)[
        ID != "__measured_background__"
      ]
    }
  ), fill = TRUE)
}

# --- combine results and record interpretation themes ---

enrichment <- rbindlist(list(
  format_enrichment(go_bp, "GO_BP"),
  format_enrichment(go_mf, "GO_MF"),
  format_enrichment(kegg, "KEGG"),
  aging_results
), fill = TRUE)

if (nrow(enrichment) == 0L) {
  enrichment <- data.table(
    database = character(), ID = character(), Description = character(),
    prespecified_theme = character(), GeneRatio = character(), BgRatio = character(),
    pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
    global_fdr = numeric(), geneID = character(), Count = integer()
  )
} else {
  enrichment[, global_fdr := p.adjust(pvalue, method = "BH")]
  enrichment[, prespecified_theme := fifelse(
    grepl("translat|ribosom|mrna processing|rna binding", Description,
          ignore.case = TRUE),
    "translation_regulation",
    fifelse(
      grepl("proteasom|ubiquitin|protein.*degrad|protein.*stabil|autophag|lysosom",
            Description, ignore.case = TRUE),
      "protein_stability_degradation",
      fifelse(
        grepl("stress|unfolded protein|heat shock|oxidative|dna damage",
              Description, ignore.case = TRUE),
        "stress_response",
        "other"
      )
    )
  )]
  setcolorder(enrichment, c(
    "database", "ID", "Description", "prespecified_theme",
    "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "global_fdr",
    "geneID", "Count"
  ))
}

mapping_qc <- data.table(
  metric = c(
    "primary_background_symbols", "primary_background_entrez_mapped",
    "candidate_symbols", "candidate_entrez_mapped",
    "candidate_entrez_in_background", "kegg_gene_sets",
    "kegg_background_genes", "aging_gene_set_source_available",
    "genage_background_genes", "genage_candidate_overlap",
    "cellage_background_genes", "cellage_candidate_overlap"
  ),
  value = c(
    length(background_symbols), length(background_entrez),
    length(candidate_symbols), uniqueN(candidate_mapping$ENTREZID),
    length(candidate_entrez), uniqueN(kegg_link$term),
    uniqueN(kegg_link[gene %chin% background_entrez, gene]),
    file.exists(aging_gene_set_path),
    uniqueN(aging_sets[collection == "GenAge", gene_symbol]),
    uniqueN(aging_sets[
      collection == "GenAge" & gene_symbol %chin% candidate_symbols,
      gene_symbol
    ]),
    uniqueN(aging_sets[collection == "CellAge", gene_symbol]),
    uniqueN(aging_sets[
      collection == "CellAge" & gene_symbol %chin% candidate_symbols,
      gene_symbol
    ])
  ),
  note = c(
    "All genes eligible for primary three-layer selection",
    "Mapped with org.Hs.eg.db SYMBOL to ENTREZID",
    "Frozen APA-linked protein candidate list",
    "Mapped with org.Hs.eg.db SYMBOL to ENTREZID",
    "Candidate Entrez IDs retained in the measured background",
    "Local KEGG human snapshot",
    "Measured background represented in local KEGG snapshot",
    paste(gene_set_manifest[collection %chin% c("GenAge", "CellAge"),
                            paste(collection, status, sep = "=")], collapse = ";"),
    "GenAge genes represented in the 5,258-gene measured background",
    "Frozen candidate genes in GenAge",
    "CellAge genes represented in the 5,258-gene measured background",
    "Frozen candidate genes in CellAge"
  )
)

# --- write enrichment outputs ---

fwrite(enrichment, file.path(mod, "05_functional_enrich.tsv"), sep = "\t", na = "NA")
fwrite(mapping_qc, file.path(mod, "05_enrich_mapping_qc.tsv"), sep = "\t", na = "NA")
