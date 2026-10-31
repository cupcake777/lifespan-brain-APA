# Gene-set enrichment (Fisher, one-sided) of the HC3-core G x Stage genes and their
# direction classes, BH within foreground and library. Writes functional_enrichment_* to data/05_gxs.
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# --- paths and settings ---
APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
D <- file.path(APA_ROOT, "data", "05_gxs")
REF <- file.path(APA_ROOT, "ref", "msigdb")
core <- fread(file.path(D, "interaction_hc3_core_groups.tsv"))
tested <- fread(file.path(D, "interaction_group_summary.tsv"), select = "gene_id")

gmt_files <- c(
  Hallmark = file.path(REF, "h.all.v2024.1.Hs.symbols.gmt"),
  Reactome = file.path(REF, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"),
  WikiPathways = file.path(REF, "c2.cp.wikipathways.v2024.1.Hs.symbols.gmt"),
  KEGG = file.path(REF, "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"),
  HPO = file.path(APA_ROOT, "ref", "hpo", "Human_Phenotype_Ontology.gmt"),
  AllenBrainHuman = file.path(APA_ROOT, "ref", "allen", "Allen_Brain_Atlas_10x_scRNA_2021.gmt"),
  GOBP = file.path(REF, "c5.go.bp.v2024.1.Hs.symbols.gmt"),
  GOMF = file.path(REF, "c5.go.mf.v2024.1.Hs.symbols.gmt"),
  GOCC = file.path(REF, "c5.go.cc.v2024.1.Hs.symbols.gmt")
)
min_term_size <- 10L
max_term_size <- 500L
fdr_cutoff <- 0.05

background <- unique(tested[!is.na(gene_id) & gene_id != "", gene_id])
foregrounds <- list(
  "All HC3 core" = unique(core$gene_id),
  Opposing = unique(core[point_estimate_class == "Opposing", gene_id]),
  Amplified = unique(core[point_estimate_class == "Amplified", gene_id]),
  Attenuated = unique(core[point_estimate_class == "Attenuated", gene_id])
)

# --- enrichment per library and foreground ---
result_list <- list()
for (library_name in names(gmt_files)) {
  lines <- readLines(gmt_files[[library_name]], warn = FALSE)
  lines <- lines[nzchar(lines)]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 3L]
  gene_sets <- data.table(
    library = library_name,
    term = vapply(fields, `[`, character(1), 1L),
    genes = lapply(fields, function(x) unique(x[-c(1L, 2L)]))
  )
  if (library_name == "AllenBrainHuman") gene_sets <- gene_sets[grepl("^Human ", term)]

  all_lib_genes <- unique(unlist(gene_sets$genes, use.names = FALSE))
  annotated_background <- intersect(background, all_lib_genes)
  if (length(annotated_background) == 0L) next

  gene_sets[, term_genes_bg := lapply(genes, intersect, y = annotated_background)]
  gene_sets[, term_size := lengths(term_genes_bg)]
  valid_sets <- gene_sets[term_size >= min_term_size & term_size <= max_term_size]
  if (nrow(valid_sets) == 0L) next

  for (foreground_name in names(foregrounds)) {
    foreground <- intersect(foregrounds[[foreground_name]], annotated_background)
    reference <- setdiff(annotated_background, foreground)
    len_fg <- length(foreground)
    len_ref <- length(reference)

    rows <- lapply(seq_len(nrow(valid_sets)), function(i) {
      term_genes <- valid_sets$term_genes_bg[[i]]
      driver_genes <- intersect(foreground, term_genes)
      a <- length(driver_genes)
      b <- len_fg - a
      c <- length(intersect(reference, term_genes))
      d <- len_ref - c
      enrichment_test <- fisher.test(matrix(c(a, b, c, d), nrow = 2L), alternative = "greater")
      term <- valid_sets$term[i]
      data.table(
        foreground = foreground_name,
        library = library_name,
        term = term,
        term_label = str_to_sentence(tolower(gsub("_", " ", gsub("^(GOBP_|GOMF_|GOCC_|KEGG_|REACTOME_)", "", term)))),
        foreground_genes = len_fg,
        annotated_background_genes = length(annotated_background),
        term_size = valid_sets$term_size[i],
        driver_count = a,
        driver_genes = paste(sort(driver_genes), collapse = "/"),
        odds_ratio = unname(enrichment_test$estimate),
        odds_ratio_ci95_low = unname(enrichment_test$conf.int[1L]),
        odds_ratio_ci95_high = unname(enrichment_test$conf.int[2L]),
        pvalue = enrichment_test$p.value
      )
    })
    result_list[[paste(foreground_name, library_name, sep = "|")]] <- rbindlist(rows, fill = TRUE)
  }
}

# --- BH within foreground x library (terms with >= 1 driver) ---
results <- rbindlist(result_list, fill = TRUE)
results[, bh_eligible := driver_count > 0L]
results[, fdr_within_library := NA_real_]
results[bh_eligible == TRUE, fdr_within_library := p.adjust(pvalue, method = "BH"), by = .(foreground, library)]
setorder(results, foreground, library, fdr_within_library, pvalue)
fwrite(results, file.path(D, "functional_enrichment_all_terms.tsv"), sep = "\t")

supported <- results[bh_eligible == TRUE & fdr_within_library < fdr_cutoff]
if (nrow(supported) > 0L) {
  supported[, driver_signature := vapply(strsplit(driver_genes, "/", fixed = TRUE),
                                         function(x) paste(sort(x), collapse = "/"), character(1))]
  supported[, redundant_terms_for_drivers := .N, by = .(foreground, driver_signature)]
  setorder(supported, foreground, fdr_within_library, pvalue, library, term)
}
fwrite(supported, file.path(D, "functional_enrichment_supported_terms.tsv"), sep = "\t")

summ <- results[, .(
  audited_terms = .N,
  bh_tested_terms = sum(bh_eligible),
  zero_overlap_terms = sum(!bh_eligible),
  annotated_foreground_genes = unique(foreground_genes),
  minimum_fdr = if (any(bh_eligible)) min(fdr_within_library, na.rm = TRUE) else NA_real_,
  supported_terms = sum(bh_eligible & fdr_within_library < fdr_cutoff, na.rm = TRUE)
), by = .(foreground, library)]
summ[, interpretation := fifelse(supported_terms > 0L,
                                 "Supported after BH correction within this library",
                                 "No term supported after BH correction within this library")]
setorder(summ, foreground, library)
fwrite(summ, file.path(D, "functional_enrichment_summary.tsv"), sep = "\t")

# --- conclusion text ---
supported_lines <- if (nrow(supported) > 0L) {
  supported[, sprintf(
    "- **%s / %s:** %s; OR %.2f (95%% CI %.2f-%.2f), one-sided Fisher p %.3g, within-library BH FDR %.3g; driver genes: %s",
    foreground, library, term_label, odds_ratio, odds_ratio_ci95_low,
    odds_ratio_ci95_high, pvalue, fdr_within_library, driver_genes)]
} else {
  "- No term passed BH correction within its source library"
}
writeLines(c(
  "# Functional enrichment conclusion", "",
  "This analysis tests the complete HC3-core gene set and the Opposing, Amplified, and Attenuated pattern gene sets against their corresponding annotated tested-gene backgrounds.",
  "Multiple testing is controlled separately within each foreground and source gene-set library using BH FDR among terms with at least one foreground overlap.",
  "Zero-overlap terms remain in the full audit table with bh_eligible=FALSE and FDR=NA, but are excluded from the BH matrix.",
  "", "## Supported results", "", supported_lines
), file.path(D, "functional_enrichment_conclusion.md"))
message(sprintf("Audited %d foreground-term rows; BH-tested %d nonzero-overlap rows; retained %d terms at within-library BH FDR < %.2f",
                nrow(results), results[bh_eligible == TRUE, .N], nrow(supported), fdr_cutoff))
