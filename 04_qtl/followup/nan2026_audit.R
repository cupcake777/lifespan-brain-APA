#!/usr/bin/env Rscript
# Source audit for the Nan 2026 comparison: sensitivity models for the interaction and stage-specific apaQTL overlaps.
# Writes R4_nan2026_source_contract.tsv, R4_nan2026_overlap_sens.tsv and R4_nan2026_model_sens.tsv to data/04_qtl/.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# --- helpers ---
fit_sensitivity <- function(data, predictor, covariates, model_id, unit) {
  full_formula <- reformulate(c(sprintf("factor(%s)", predictor), covariates), response = "replicated")
  reduced_formula <- reformulate(covariates, response = "replicated")
  full_fit <- glm(full_formula, data = data, family = binomial())
  reduced_fit <- glm(reduced_formula, data = data, family = binomial())
  lrt <- anova(reduced_fit, full_fit, test = "Chisq")
  coefs <- as.data.table(summary(full_fit)$coefficients, keep.rownames = "term")
  setnames(coefs, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
           c("estimate", "std_error", "statistic", "p_value"))
  coefs[, `:=`(
    model_id = model_id,
    unit = unit,
    test = "Wald",
    n_observations = nrow(data),
    n_positive = sum(data$replicated),
    formula = paste(deparse(full_formula), collapse = " ")
  )]
  overall <- data.table(
    model_id = model_id,
    unit = unit,
    test = "likelihood-ratio",
    term = predictor,
    estimate = NA_real_,
    std_error = NA_real_,
    statistic = lrt$Deviance[2L],
    p_value = lrt$`Pr(>Chi)`[2L],
    n_observations = nrow(data),
    n_positive = sum(data$replicated),
    formula = paste(deparse(full_formula), collapse = " ")
  )
  rbind(overall, coefs, fill = TRUE)
}

overlap_counts <- function(data, grouping, match_column, analysis, contract) {
  data[, .(
    n_local = .N,
    n_matched = sum(get(match_column)),
    match_rate = mean(get(match_column))
  ), by = .(group = as.character(get(grouping)))][, `:=`(
    analysis = analysis,
    match_contract = contract,
    grouping_variable = grouping
  )]
}

read_stage_leads <- function(path, prefix) {
  x <- fread(path, select = c("group_id", "qval", "af", "pval_nominal"))
  setorder(x, group_id, pval_nominal)
  x <- x[, .SD[1L], by = group_id]
  setnames(x, c("qval", "af", "pval_nominal"),
           paste0(prefix, c("_qval", "_af", "_pval_nominal")))
  x
}

# --- config ---
APA_ROOT <- Sys.getenv("APA_ROOT")
if (!nzchar(APA_ROOT)) stop("Set APA_ROOT (source config/paths.sh)")
out_dir <- file.path(APA_ROOT, "data", "04_qtl")
workbook <- file.path(APA_ROOT, "ref", "nan2026", "41588_2026_2758_MOESM3_ESM.xlsx")
nominal_dir <- file.path(APA_ROOT, "ref", "nan2026")
gtf_path <- file.path(APA_ROOT, "ref", "gencode.v40.annotation.gtf.gz")
interaction_path <- file.path(APA_ROOT, "data", "05_gxs", "intx_hc3_core_groups.tsv")
prenatal_path <- file.path(out_dir, "prenatal.cis_qtl.txt.gz")
postnatal_path <- file.path(out_dir, "postnatal.cis_qtl.txt.gz")
cell_types <- c("Ast", "Exc", "Inh", "Mic", "Oli", "Opc")

# --- Nan 2026 source table ---
nan <- as.data.table(read_excel(workbook, sheet = "Table 11", skip = 1L))
setnames(nan, "Phenotype_id (transcript ID|gene ID|chromosome|strand)", "phenotype_id")
raw_sheet_rows <- nrow(nan)
nan <- nan[!is.na(SNP) & !is.na(phenotype_id) & !is.na(gene) & !is.na(celltype)]
excluded_sheet_rows <- raw_sheet_rows - nrow(nan)
nan[, transcript_id := tstrsplit(phenotype_id, "|", fixed = TRUE)[[1L]]]
nan_gene_keys <- unique(nan[, paste(gene, celltype, sep = "\t")])
nan_genes <- unique(nan$gene)
nan_transcripts <- unique(nan$transcript_id)

source_contract <- data.table(
  resource = c(
    "Supplementary Table 11 workbook",
    paste0(cell_types, ".3aQTL.sumstats.txt.gz"),
    "Author analysis code"
  ),
  source_url = c(
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-026-02758-w/MediaObjects/41588_2026_2758_MOESM3_ESM.xlsx",
    paste0("https://zenodo.org/api/records/20740673/files/", cell_types, ".3aQTL.sumstats.txt.gz/content"),
    "https://github.com/xiongxslab/NanJ_Carles_et_al_3aQTL_project"
  ),
  bytes = c(10180414, 610926198, 3990224351, 1530770464, 167127690, 1413646231, 484447745, NA_real_),
  local_path = c(workbook, file.path(nominal_dir, paste0(cell_types, ".3aQTL.sumstats.txt.gz")), "not copied; public repository inspected"),
  observed_contract = c(
    sprintf("Table 11: %d valid significant lead g-APA rows; %d non-data worksheet rows excluded; %d genes; %d transcripts; cell types=%s",
            nrow(nan), excluded_sheet_rows, uniqueN(nan$gene), uniqueN(nan$transcript_id), paste(sort(unique(nan$celltype)), collapse = ",")),
    rep("Full nominal variant-phenotype statistics required for exact pair-level replication", 6L),
    "Cross-cell-type producer matches phenotype_id + variant_id against full nominal statistics"
  )
)
source_contract[, sha256 := NA_character_]
source_contract[1L, sha256 := strsplit(system2("sha256sum", workbook, stdout = TRUE), " ", fixed = TRUE)[[1L]][1L]]
fwrite(source_contract, file.path(out_dir, "R4_nan2026_source_contract.tsv"), sep = "\t", na = "NA")

# --- interaction overlap sensitivities ---
interaction <- fread(interaction_path)
interaction[, transcript_id := tstrsplit(lead_phenotype_id, "|", fixed = TRUE)[[1L]]]
interaction[, `:=`(
  gene_any_cell = gene_id %chin% nan_genes,
  transcript_any_cell = transcript_id %chin% nan_transcripts
)]

interaction_overlap <- rbindlist(lapply(
  c("point_estimate_class", "stage_group_status", "strict_ci_direction_reversal"),
  function(grouping) rbind(
    overlap_counts(interaction, grouping, "gene_any_cell", "interaction HC3 core", "gene has a significant Nan g-APA in any cell type"),
    overlap_counts(interaction, grouping, "transcript_any_cell", "interaction HC3 core", "lead transcript has a significant Nan g-APA in any cell type")
  )
), fill = TRUE)

interaction_long <- interaction[, .(celltype = cell_types),
  by = .(group_id, gene_id, lead_af, point_estimate_class, stage_group_status, strict_ci_direction_reversal)]
interaction_long[, replicated := paste(gene_id, celltype, sep = "\t") %chin% nan_gene_keys]
interaction[, replicated := as.integer(gene_any_cell)]

interaction_models <- rbindlist(list(
  fit_sensitivity(interaction, "point_estimate_class", "lead_af", "interaction_point_class_gene_any_cell", "one HC3-core gene"),
  fit_sensitivity(interaction, "stage_group_status", "lead_af", "interaction_stage_status_gene_any_cell", "one HC3-core gene"),
  fit_sensitivity(interaction, "strict_ci_direction_reversal", "lead_af", "interaction_strict_reversal_gene_any_cell", "one HC3-core gene"),
  fit_sensitivity(interaction_long, "point_estimate_class", c("lead_af", "factor(celltype)"), "interaction_point_class_celltype_expanded", "gene x cell type; repeated-gene sensitivity"),
  fit_sensitivity(interaction_long, "stage_group_status", c("lead_af", "factor(celltype)"), "interaction_stage_status_celltype_expanded", "gene x cell type; repeated-gene sensitivity"),
  fit_sensitivity(interaction_long, "strict_ci_direction_reversal", c("lead_af", "factor(celltype)"), "interaction_strict_reversal_celltype_expanded", "gene x cell type; repeated-gene sensitivity")
), fill = TRUE)

# --- separate-stage apaQTL overlap sensitivities ---
gtf <- fread(cmd = sprintf("gzip -cd %s | grep -v '^#'", shQuote(gtf_path)), sep = "\t", header = FALSE,
             select = c(3L, 9L), quote = "", fill = TRUE, showProgress = FALSE)
gtf <- gtf[V3 == "gene"]
gtf[, gene_id := sub(".*gene_id \"([^\"]+)\".*", "\\1", V9)]
gtf[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
gtf[, gene_symbol := sub(".*gene_name \"([^\"]+)\".*", "\\1", V9)]
gene_map <- unique(gtf[, .(gene_id, gene_symbol)], by = "gene_id")

stage <- merge(read_stage_leads(prenatal_path, "prenatal"),
               read_stage_leads(postnatal_path, "postnatal"), by = "group_id", all = TRUE)
stage[, qtl_class := "Neither"]
stage[!is.na(prenatal_qval) & prenatal_qval <= 0.05 & (is.na(postnatal_qval) | postnatal_qval > 0.05), qtl_class := "Prenatal-specific"]
stage[(is.na(prenatal_qval) | prenatal_qval > 0.05) & !is.na(postnatal_qval) & postnatal_qval <= 0.05, qtl_class := "Postnatal-specific"]
stage[!is.na(prenatal_qval) & prenatal_qval <= 0.05 & !is.na(postnatal_qval) & postnatal_qval <= 0.05, qtl_class := "Shared"]
stage <- stage[qtl_class != "Neither"]
stage[, gene_symbol := fifelse(startsWith(group_id, "ENSG"), gene_map$gene_symbol[match(group_id, gene_map$gene_id)], group_id)]
stage[, maf := fcase(
  qtl_class == "Prenatal-specific", prenatal_af,
  qtl_class == "Postnatal-specific", postnatal_af,
  default = rowMeans(cbind(prenatal_af, postnatal_af), na.rm = TRUE)
)]
stage[, replicated := as.integer(gene_symbol %chin% nan_genes)]

stage_overlap <- stage[, .(
  n_local = .N,
  n_matched = sum(replicated),
  match_rate = mean(replicated)
), by = .(group = qtl_class)][, `:=`(
  analysis = "separate-stage apaQTL",
  match_contract = "gene has a significant Nan g-APA in any cell type",
  grouping_variable = "qtl_class"
)]

stage_long <- stage[, .(celltype = cell_types), by = .(group_id, gene_symbol, qtl_class, maf)]
stage_long[, replicated := paste(gene_symbol, celltype, sep = "\t") %chin% nan_gene_keys]
stage_models <- rbindlist(list(
  fit_sensitivity(stage, "qtl_class", "maf", "stage_qtl_class_gene_any_cell", "one local apaQTL gene"),
  fit_sensitivity(stage_long, "qtl_class", c("maf", "factor(celltype)"), "stage_qtl_class_celltype_expanded", "gene x cell type; repeated-gene sensitivity")
), fill = TRUE)

overlap <- rbind(interaction_overlap, stage_overlap, fill = TRUE)
setcolorder(overlap, c("analysis", "match_contract", "grouping_variable", "group", "n_local", "n_matched", "match_rate"))
fwrite(overlap, file.path(out_dir, "R4_nan2026_overlap_sens.tsv"), sep = "\t", na = "NA")

models <- rbind(interaction_models, stage_models, fill = TRUE)
models[, status := "sensitivity only: public lead table cannot reconstruct pair-level adjusted model"]
setcolorder(models, c("model_id", "unit", "status", "test", "term", "estimate", "std_error", "statistic", "p_value", "n_observations", "n_positive", "formula"))
fwrite(models, file.path(out_dir, "R4_nan2026_model_sens.tsv"), sep = "\t", na = "NA")

# --- input contract: Table 11 and interaction row counts ---
stopifnot(nrow(nan) == 6706L, excluded_sheet_rows == 436L,
          setequal(unique(nan$celltype), cell_types), nrow(interaction) == 193L)
cat(sprintf("Nan Table 11 rows=%d genes=%d transcripts=%d; interaction=%d; stage-QTL genes=%d\n",
            nrow(nan), uniqueN(nan$gene), uniqueN(nan$transcript_id), nrow(interaction), nrow(stage)))
