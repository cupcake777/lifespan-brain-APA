#!/usr/bin/env Rscript
# Builds matched APA, RNA and protein matrices and the sample manifest for the 12_protein analysis.

suppressPackageStartupMessages(library(data.table))

# --- config ---
apa_root <- Sys.getenv("APA_ROOT")
data_dir <- file.path(apa_root, "data")
omics_dir <- file.path(data_dir, "03_omics")
apa_path <- file.path(data_dir, "04_qtl", "apa_after_combat.tsv")
rna_path <- file.path(data_dir, "01_quant", "expr_after_combat.tsv")
rna_meta_path <- file.path(data_dir, "04_qtl", "CHB.all.meta")
protein_path <- file.path(omics_dir, "prot_modresults_covadj_adj.RData")
protein_meta_path <- file.path(omics_dir, "meta_all.txt")
output_dir <- omics_dir

# --- helpers ---
normalize_id <- function(x) gsub("[^A-Z0-9]", "", toupper(x))

# --- load sample identifiers and metadata ---
rna_header <- names(fread(rna_path, nrows = 0L, check.names = FALSE))
rna_ids <- rna_header[-1]
apa_header <- names(fread(apa_path, nrows = 0L, check.names = FALSE))
apa_ids <- apa_header[-1]

protein_objects <- new.env(parent = emptyenv())
load(protein_path, envir = protein_objects)
protein_ids <- colnames(protein_objects$quant_final)

rna_norm <- normalize_id(rna_ids)
protein_norm <- normalize_id(protein_ids)
protein_index <- match(rna_norm, protein_norm)
matched_rna_index <- which(!is.na(protein_index))

stopifnot(
  !anyDuplicated(rna_norm),
  !anyDuplicated(protein_norm),
  !anyDuplicated(apa_ids),
  length(matched_rna_index) > 0L,
  all(rna_ids[matched_rna_index] %in% apa_ids),
  identical(colnames(protein_objects$quant_final), rownames(protein_objects$covar))
)

rna_meta <- fread(rna_meta_path, data.table = FALSE, check.names = FALSE)
protein_meta <- fread(protein_meta_path, data.table = FALSE, check.names = FALSE)
rna_meta_index <- match(rna_norm[matched_rna_index], normalize_id(rna_meta$sample))
protein_meta_index <- match(
  protein_norm[protein_index[matched_rna_index]],
  normalize_id(protein_meta$ID)
)

sample_manifest <- data.frame(
  sample_id = rna_ids[matched_rna_index],
  canonical_id = rna_norm[matched_rna_index],
  rna_sample_id = rna_ids[matched_rna_index],
  protein_sample_id = protein_ids[protein_index[matched_rna_index]],
  match_rule = ifelse(
    rna_ids[matched_rna_index] == protein_ids[protein_index[matched_rna_index]],
    "exact",
    "case_or_separator"
  ),
  age_years = rna_meta$age[rna_meta_index],
  analysis_age_source = "RNA CHB.all.meta",
  sex = substr(toupper(rna_meta$sex[rna_meta_index]), 1, 1),
  life_stage = rna_meta$LifeStage[rna_meta_index],
  rna_batch = rna_meta$batch[rna_meta_index],
  rna_rin = rna_meta$RIN[rna_meta_index],
  protein_age = protein_meta$age[protein_meta_index],
  protein_sex = substr(toupper(protein_meta$sex[protein_meta_index]), 1, 1),
  protein_rin = protein_meta$RIN[protein_meta_index],
  protein_digest_batch = protein_meta$digest_batch[protein_meta_index],
  protein_disease_class = protein_meta$disease_class[protein_meta_index],
  stringsAsFactors = FALSE
)
sample_manifest$include_analysis <- TRUE
sample_manifest$exclusion_reason <- ""

analysis_samples <- sample_manifest$sample_id[sample_manifest$include_analysis]
analysis_protein_ids <- sample_manifest$protein_sample_id[sample_manifest$include_analysis]

stopifnot(
  !anyNA(rna_meta_index),
  !anyNA(protein_meta_index),
  all(sample_manifest$sex[sample_manifest$include_analysis] ==
        sample_manifest$protein_sex[sample_manifest$include_analysis])
)

# --- prepare matched APA matrix ---
apa_data <- fread(
  apa_path,
  select = c(apa_header[1], analysis_samples),
  check.names = FALSE,
  showProgress = TRUE
)
apa_feature_id <- apa_data[[1]]
apa_data[[1]] <- NULL
apa_matrix <- as.matrix(apa_data[, analysis_samples, with = FALSE])
rownames(apa_matrix) <- apa_feature_id

# --- prepare matched RNA matrix ---
rna_data <- fread(
  rna_path,
  select = c(rna_header[1], analysis_samples),
  check.names = FALSE,
  showProgress = TRUE
)
rna_gene_id <- rna_data[[1]]
rna_data[[1]] <- NULL
rna_matrix <- as.matrix(rna_data[, analysis_samples, with = FALSE])
rownames(rna_matrix) <- rna_gene_id

# --- prepare matched protein matrix and annotation ---
protein_annotation <- as.data.frame(protein_objects$data_rmOut_imp$inf)
protein_accession <- as.character(protein_annotation$Protein.accession)

stopifnot(
  nrow(protein_annotation) == nrow(protein_objects$quant_final),
  identical(protein_accession, rownames(protein_objects$data_rmOut_imp$log2_value)),
  !anyDuplicated(protein_accession)
)

protein_matrix <- protein_objects$quant_final[, analysis_protein_ids, drop = FALSE]
rownames(protein_matrix) <- protein_accession
colnames(protein_matrix) <- analysis_samples

stopifnot(
  identical(colnames(apa_matrix), analysis_samples),
  identical(colnames(rna_matrix), analysis_samples),
  identical(colnames(protein_matrix), analysis_samples)
)

# --- save matched inputs ---
sample_manifest_path <- file.path(output_dir, "matched_sample_manifest.tsv")
protein_annotation_path <- file.path(output_dir, "prot_feature_annotation.tsv")
apa_output_path <- file.path(output_dir, "apa_matched.rds")
rna_output_path <- file.path(output_dir, "rna_matched.rds")
protein_output_path <- file.path(output_dir, "protein_matched.rds")
input_manifest_path <- file.path(output_dir, "input_manifest.tsv")

fwrite(sample_manifest, sample_manifest_path, sep = "\t", na = "NA")
fwrite(protein_annotation, protein_annotation_path, sep = "\t", na = "NA")
saveRDS(apa_matrix, apa_output_path, compress = "gzip")
saveRDS(rna_matrix, rna_output_path, compress = "gzip")
saveRDS(protein_matrix, protein_output_path, compress = "gzip")

input_manifest <- data.frame(
  artifact = c(
    "matched_sample_manifest", "protein_feature_annotation",
    "apa_matched", "rna_matched", "protein_matched"
  ),
  path = c(
    sample_manifest_path, protein_annotation_path,
    apa_output_path, rna_output_path, protein_output_path
  ),
  source = c(
    paste(rna_meta_path, protein_meta_path, sep = ";"),
    protein_path,
    apa_path, rna_path, protein_path
  ),
  n_features = c(NA, nrow(protein_annotation), nrow(apa_matrix),
                 nrow(rna_matrix), nrow(protein_matrix)),
  n_samples = c(nrow(sample_manifest), NA, ncol(apa_matrix),
                ncol(rna_matrix), ncol(protein_matrix)),
  notes = c(
    sprintf("%d ID-level matches; RNA metadata is authoritative for analysis age",
            nrow(sample_manifest)),
    "Rows aligned to protein_matched.rds by Protein.accession",
    sprintf("Columns are %d analysis-ready RNA sample IDs in manifest order",
            length(analysis_samples)),
    sprintf("Columns are %d analysis-ready RNA sample IDs in manifest order",
            length(analysis_samples)),
    "Covariate-adjusted quant_final; columns renamed to RNA sample IDs"
  ),
  stringsAsFactors = FALSE
)
input_manifest$md5 <- unname(tools::md5sum(input_manifest$path))
if (file.exists(input_manifest_path)) {
  existing <- fread(input_manifest_path)
  source_rows <- existing[!artifact %in% input_manifest$artifact]
  input_manifest <- rbindlist(list(as.data.table(input_manifest), source_rows),
                              use.names = TRUE, fill = TRUE)
}
fwrite(input_manifest, input_manifest_path, sep = "\t", na = "NA")
