#!/usr/bin/env Rscript
# GO BP enrichment of APA-QTL genes and GxStage mechanism/consequence tests.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})


# --- helpers ---

map_gene_ids <- function(genes) {
  genes <- unique(na.omit(as.character(genes)))
  ensembl <- genes[grepl("^ENSG", genes)]
  symbols <- setdiff(genes, ensembl)
  mapped <- list()
  if (length(symbols) > 0L) {
    mapped[[length(mapped) + 1L]] <- suppressMessages(
      bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )
  }
  if (length(ensembl) > 0L) {
    mapped[[length(mapped) + 1L]] <- suppressMessages(
      bitr(sub("\\..*$", "", ensembl), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )
  }
  if (length(mapped) == 0L) return(character())
  unique(rbindlist(mapped, fill = TRUE)$ENTREZID)
}

build_tested_stage <- function(path, stage_name) {
  table <- fread(path, select = 1:4, showProgress = FALSE)
  setnames(table, names(table), c("phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_id"))
  parts <- tstrsplit(table$phenotype_id, "|", fixed = TRUE, keep = 1:4)
  if (length(parts) != 4L || any(lengths(parts) != nrow(table))) {
    stop(sprintf("Unable to parse phenotype IDs from %s", path))
  }
  table[, `:=`(
    transcript_id = parts[[1L]],
    gene_name = parts[[2L]],
    phenotype_strand = parts[[4L]],
    stage = stage_name,
    source = basename(path)
  )]
  table <- table[nzchar(phenotype_id) & nzchar(gene_name)]
  unique(table[, .(
    phenotype_id, transcript_id, gene_name, phenotype_chr,
    phenotype_start, phenotype_end, phenotype_strand, stage, source
  )], by = c("stage", "phenotype_id"))
}

run_go_contrast <- function(contrast, line, foreground, background, background_contract) {
  foreground <- unique(na.omit(as.character(foreground)))
  background <- unique(na.omit(as.character(background)))
  foreground_in_background <- intersect(foreground, background)
  foreground_ids <- intersect(map_gene_ids(foreground), map_gene_ids(background))
  background_ids <- map_gene_ids(background)
  if (length(foreground_ids) < minimum_mapped_foreground) {
    return(data.table(
      contrast = contrast,
      line = line,
      ontology = "GO_BP",
      term_id = NA_character_,
      term_description = NA_character_,
      foreground_input_n = length(foreground),
      foreground_in_background_n = length(foreground_in_background),
      foreground_missing_background_n = length(setdiff(foreground, background)),
      foreground_mapped_n = length(foreground_ids),
      background_mapped_n = length(background_ids),
      gene_ratio = NA_character_,
      background_ratio = NA_character_,
      p_value = NA_real_,
      p_adjust = NA_real_,
      gene_ids = NA_character_,
      test_status = "insufficient_mapped_foreground",
      background_contract = background_contract
    ))
  }
  fit <- enrichGO(
    gene = foreground_ids,
    universe = background_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  table <- as.data.table(as.data.frame(fit))
  if (nrow(table) == 0L) {
    return(data.table(
      contrast = contrast,
      line = line,
      ontology = "GO_BP",
      term_id = NA_character_,
      term_description = NA_character_,
      foreground_input_n = length(foreground),
      foreground_in_background_n = length(foreground_in_background),
      foreground_missing_background_n = length(setdiff(foreground, background)),
      foreground_mapped_n = length(foreground_ids),
      background_mapped_n = length(background_ids),
      gene_ratio = NA_character_,
      background_ratio = NA_character_,
      p_value = NA_real_,
      p_adjust = NA_real_,
      gene_ids = NA_character_,
      test_status = "no_testable_terms",
      background_contract = background_contract
    ))
  }
  table[, .(
    contrast = contrast,
    line = line,
    ontology = "GO_BP",
    term_id = ID,
    term_description = Description,
    foreground_input_n = length(foreground),
    foreground_in_background_n = length(foreground_in_background),
    foreground_missing_background_n = length(setdiff(foreground, background)),
    foreground_mapped_n = length(foreground_ids),
    background_mapped_n = length(background_ids),
    gene_ratio = GeneRatio,
    background_ratio = BgRatio,
    p_value = pvalue,
    p_adjust = p.adjust,
    gene_ids = geneID,
    test_status = "ok",
    background_contract = background_contract
  )]
}


# --- config ---

data_dir <- file.path(Sys.getenv("APA_ROOT"), "data/07_coloc")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
tested_universe_path <- file.path(data_dir, "tested_genes.tsv")
association_path <- file.path(data_dir, "01_significant_qtl_universe.tsv.gz")
mechanism_path <- file.path(data_dir, "02_signal_mechanisms.tsv.gz")
event_path <- file.path(data_dir, "annotations.tsv.gz")
consequence_path <- file.path(data_dir, "03_event_consequences.tsv.gz")
enrichment_output <- file.path(data_dir, "06_functional_enrichment.tsv.gz")
gxstage_output <- file.path(data_dir, "gxstage.tsv")
summary_output <- file.path(data_dir, "summary.tsv")
membership_output <- file.path(data_dir, "membership.tsv")
phenotype_paths <- c(
  prenatal = file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/phenotype.Prenatal.bed.gz"),
  postnatal = file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/phenotype.Postnatal.bed.gz")
)
analysis_version <- "framework_2026-07-25_v5"
minimum_mapped_foreground <- 10L


# --- loading data ---

associations <- fread(association_path)
mechanisms <- fread(mechanism_path)
events <- fread(event_path)
consequences <- fread(consequence_path)
valid_gxstage_classes <- c("Opposing", "Amplified", "Attenuated", "Uncertain")
gxstage_links <- unique(associations[
  trimws(as.character(gxstage_classification_status)) == "formal" &
    trimws(as.character(gxstage_class)) %chin% valid_gxstage_classes,
  .(
    interaction_group_id, signal_id, phenotype_id, gene_name,
    interaction_representative_phenotype_id, interaction_representative_variant_id,
    gxstage_class,
    gxstage_classification_status, gxstage_class_source,
    pattern_tau, pattern_delta, pattern_method
  )
])
if (anyNA(gxstage_links$interaction_group_id) || any(!nzchar(gxstage_links$interaction_group_id))) {
  stop("Formal GxStage records must have an interaction_group_id")
}
if (any(gxstage_links[, uniqueN(interaction_group_id), by = signal_id]$V1 > 1L)) {
  stop("A stage-specific QTL signal maps to multiple formal interaction groups")
}
gxstage_groups <- unique(gxstage_links[, .(
  interaction_group_id, phenotype_id, gene_name,
  interaction_representative_phenotype_id, interaction_representative_variant_id,
  gxstage_class, gxstage_classification_status, gxstage_class_source,
  pattern_tau, pattern_delta, pattern_method
)])
if (anyDuplicated(gxstage_groups$interaction_group_id)) {
  stop("Formal GxStage metadata are not unique by interaction_group_id")
}
tested <- rbindlist(Map(build_tested_stage, phenotype_paths, names(phenotype_paths)), use.names = TRUE)
setorder(tested, stage, phenotype_id)
fwrite(tested, tested_universe_path, sep = "\t")

significant_keys <- unique(associations[, .(stage, phenotype_id)])
tested_keys <- unique(tested[, .(stage, phenotype_id)])
missing_significant_keys <- significant_keys[!tested_keys, on = .(stage, phenotype_id)]
if (nrow(missing_significant_keys) > 0L) {
  stop(sprintf(
    "Full tested universe is missing %d significant stage-phenotype keys",
    nrow(missing_significant_keys)
  ))
}


# --- tested-universe GO enrichment ---

background_all <- unique(tested$gene_name)
significant_all <- unique(associations$gene_name)
high_pip_all <- unique(associations[finemap_tier == "high_confidence", gene_name])

go_results <- list(
  run_go_contrast(
    "all_significant_vs_tested", "line_a_primary", significant_all, background_all,
    "all_genes_in_full_prenatal_and_postnatal_apa_phenotype_beds"
  ),
  run_go_contrast(
    "prenatal_significant_vs_prenatal_tested",
    "line_a_primary",
    unique(associations[stage == "prenatal", gene_name]),
    unique(tested[stage == "prenatal", gene_name]),
    "all_genes_in_full_prenatal_apa_phenotype_bed"
  ),
  run_go_contrast(
    "postnatal_significant_vs_postnatal_tested",
    "line_a_primary",
    unique(associations[stage == "postnatal", gene_name]),
    unique(tested[stage == "postnatal", gene_name]),
    "all_genes_in_full_postnatal_apa_phenotype_bed"
  ),
  run_go_contrast(
    "high_confidence_finemapped_vs_significant",
    "fine_mapped_sensitivity",
    high_pip_all,
    significant_all,
    "all_genes_with_significant_apa_qtl"
  )
)

gxstage_genes <- unique(gxstage_groups[, .(gene_name, gxstage_class)])
go_results[[length(go_results) + 1L]] <- run_go_contrast(
  "flip_vs_primary_interaction_genes",
  "line_b_gxstage",
  gxstage_genes[gxstage_class == "Opposing", gene_name],
  gxstage_genes$gene_name,
  "all_genes_with_primary_upstream_gxstage_class"
)

go_table <- rbindlist(go_results, fill = TRUE)
go_table[, analysis_version := analysis_version]
fwrite(go_table, enrichment_output, sep = "\t")

membership_qc <- unique(go_table[, .(
  contrast, line, background_contract, foreground_input_n,
  foreground_in_background_n, foreground_missing_background_n,
  foreground_subset_pass = foreground_missing_background_n == 0L,
  analysis_version
)])
membership_qc <- rbind(
  membership_qc,
  data.table(
    contrast = "significant_stage_phenotypes_vs_full_tested_universe",
    line = "input_contract",
    background_contract = "all_stage_phenotypes_in_full_apa_phenotype_beds",
    foreground_input_n = nrow(significant_keys),
    foreground_in_background_n = nrow(significant_keys) - nrow(missing_significant_keys),
    foreground_missing_background_n = nrow(missing_significant_keys),
    foreground_subset_pass = nrow(missing_significant_keys) == 0L,
    analysis_version = analysis_version
  ),
  fill = TRUE
)
fwrite(membership_qc, membership_output, sep = "\t")


# --- GxStage interaction-group mechanism and consequence universe ---

gxstage_linked <- merge(
  gxstage_links,
  mechanisms[, .(
    interaction_group_id, signal_id,
    direct_pas_evidence,
    local_rbp_motif_change,
    local_clip_support,
    local_structure_change,
    gxstage_mechanism_variant_link
  )],
  by = c("interaction_group_id", "signal_id"),
  all.x = TRUE
)
gxstage_linked <- merge(
  gxstage_linked,
  events[, .(
    phenotype_id,
    dur_length,
    mirna_density_supported,
    rbp_density_supported,
    are_density_supported,
    gre_density_supported,
    ure_density_supported
  )],
  by = "phenotype_id",
  all.x = TRUE
)
mediation <- consequences[
  downstream_consequence_evidence_label == "mediation_supported",
  .(mediation_supported = TRUE),
  by = signal_id
]
gxstage_linked <- merge(gxstage_linked, mediation, by = "signal_id", all.x = TRUE)
gxstage_linked[is.na(mediation_supported), mediation_supported := FALSE]

variables <- c(
  "direct_pas_evidence", "local_rbp_motif_change", "local_clip_support",
  "local_structure_change", "mirna_density_supported", "rbp_density_supported",
  "are_density_supported", "gre_density_supported", "ure_density_supported",
  "mediation_supported"
)
for (variable in variables) {
  set(gxstage_linked, j = variable, value = fifelse(
    is.na(as.logical(gxstage_linked[[variable]])),
    FALSE,
    as.logical(gxstage_linked[[variable]])
  ))
}

gxstage_group <- gxstage_linked[, .(
  phenotype_id = phenotype_id[1L],
  gene_name = gene_name[1L],
  interaction_representative_phenotype_id = interaction_representative_phenotype_id[1L],
  interaction_representative_variant_id = interaction_representative_variant_id[1L],
  gxstage_class = gxstage_class[1L],
  gxstage_classification_status = gxstage_classification_status[1L],
  gxstage_class_source = gxstage_class_source[1L],
  pattern_tau = pattern_tau[1L],
  pattern_delta = pattern_delta[1L],
  pattern_method = pattern_method[1L],
  linked_stage_signal_n = uniqueN(signal_id),
  mechanism_variant_links = paste(sort(unique(na.omit(gxstage_mechanism_variant_link))), collapse = ";"),
  exact_variant_mechanism_n = sum(gxstage_mechanism_variant_link == "same_variant", na.rm = TRUE),
  direct_pas_evidence = any(direct_pas_evidence),
  local_rbp_motif_change = any(local_rbp_motif_change),
  local_clip_support = any(local_clip_support),
  local_structure_change = any(local_structure_change),
  mirna_density_supported = any(mirna_density_supported),
  rbp_density_supported = any(rbp_density_supported),
  are_density_supported = any(are_density_supported),
  gre_density_supported = any(gre_density_supported),
  ure_density_supported = any(ure_density_supported),
  mediation_supported = any(mediation_supported)
), by = interaction_group_id]
gxstage_group[, gxstage_binary_group := fifelse(gxstage_class == "Opposing", "Flip", "non_Flip")]


# --- GxStage mechanism and consequence tests ---

test_rows <- list()
for (variable in variables) {
  values <- as.logical(gxstage_group[[variable]])
  values[is.na(values)] <- FALSE
  gxstage_group[, value_for_test := values]
  counts <- gxstage_group[, .(
    n = .N,
    positive_n = sum(value_for_test)
  ), by = gxstage_binary_group]
  flip_n <- counts[gxstage_binary_group == "Flip", n]
  nonflip_n <- counts[gxstage_binary_group == "non_Flip", n]
  flip_positive <- counts[gxstage_binary_group == "Flip", positive_n]
  nonflip_positive <- counts[gxstage_binary_group == "non_Flip", positive_n]
  status <- if (length(flip_n) == 0L || length(nonflip_n) == 0L || min(flip_n, nonflip_n) < 10L) {
    "insufficient_n"
  } else if ((flip_positive + nonflip_positive) %in% c(0L, flip_n + nonflip_n)) {
    "no_outcome_variation"
  } else {
    "ok"
  }
  fit <- if (status == "ok") fisher.test(matrix(
    c(flip_positive, flip_n - flip_positive, nonflip_positive, nonflip_n - nonflip_positive),
    nrow = 2,
    byrow = TRUE
  )) else NULL
  test_rows[[length(test_rows) + 1L]] <- data.table(
    comparison = "Flip_vs_non_Flip",
    analysis_unit = "interaction_group",
    variable = variable,
    flip_n = ifelse(length(flip_n), flip_n, 0L),
    flip_positive_n = ifelse(length(flip_positive), flip_positive, 0L),
    nonflip_n = ifelse(length(nonflip_n), nonflip_n, 0L),
    nonflip_positive_n = ifelse(length(nonflip_positive), nonflip_positive, 0L),
    odds_ratio = if (is.null(fit)) NA_real_ else unname(fit$estimate),
    p_value = if (is.null(fit)) NA_real_ else fit$p.value,
    test_status = status,
    adjustment_status = "not_estimable_for_flip_n_below_10"
  )
}

gxstage_tests <- rbindlist(test_rows)
gxstage_tests[, p_adjust := p.adjust(p_value, method = "BH")]
gxstage_tests[, analysis_version := analysis_version]
fwrite(gxstage_tests, gxstage_output, sep = "\t")


# --- run summary ---

summary <- data.table(
  metric = c(
    "tested_background_genes", "significant_genes", "high_confidence_finemapped_genes",
    "tested_stage_phenotypes", "significant_stage_phenotypes_missing_from_tested",
    "gxstage_group_n", "flip_group_n", "nonflip_group_n",
    "primary_go_fdr_hits", "sensitivity_go_fdr_hits", "gxstage_fdr_hits"
  ),
  value = c(
    uniqueN(background_all), uniqueN(significant_all), uniqueN(high_pip_all),
    nrow(tested_keys), nrow(missing_significant_keys),
    nrow(gxstage_group), sum(gxstage_group$gxstage_binary_group == "Flip"),
    sum(gxstage_group$gxstage_binary_group == "non_Flip"),
    sum(go_table$line == "line_a_primary" & go_table$p_adjust < 0.05, na.rm = TRUE),
    sum(go_table$line == "fine_mapped_sensitivity" & go_table$p_adjust < 0.05, na.rm = TRUE),
    sum(gxstage_tests$p_adjust < 0.05, na.rm = TRUE)
  )
)
summary[, analysis_version := analysis_version]
fwrite(summary, summary_output, sep = "\t")
print(summary)
