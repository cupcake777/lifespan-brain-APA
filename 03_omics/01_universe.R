# Analysis universe for the APA-RNA-protein integration.
# Maps locked APA sites to genes with RNA and protein evidence, then builds sample QC and the support profile.

suppressPackageStartupMessages(library(data.table))

# --- config ---

if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
mod <- file.path(Sys.getenv("APA_ROOT"), "data/03_omics")

grid_size <- 128L
inner_window <- 0.04
outer_window <- 0.08
inner_min <- 3L
outer_min <- 5L
min_unique_peptides <- 2L

# --- helpers ---

strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

append_reason <- function(current, condition, reason) {
  ifelse(
    condition,
    ifelse(current == "", reason, paste(current, reason, sep = ";")),
    current
  )
}

# --- load matched matrices and sample metadata ---

input_manifest <- fread(file.path(mod, "input_manifest.tsv"))

apa <- readRDS(file.path(mod, "apa_matched.rds"))
rna <- readRDS(file.path(mod, "rna_matched.rds"))
protein <- readRDS(file.path(mod, "protein_matched.rds"))
sample_manifest <- fread(file.path(mod, "matched_sample_manifest.tsv"))

sample_ids <- sample_manifest$sample_id
stopifnot(
  all(sample_manifest$include_analysis),
  !anyDuplicated(sample_manifest$canonical_id),
  identical(colnames(apa), sample_ids),
  identical(colnames(rna), sample_ids),
  identical(colnames(protein), sample_ids),
  !anyNA(apa),
  !anyNA(rna),
  !anyNA(protein)
)

rna_gene_base <- strip_version(rownames(rna))
stopifnot(!anyDuplicated(rna_gene_base))

# --- construct the locked APA discovery universe ---

discovery <- fread(file.path(Sys.getenv("APA_ROOT"), "data/02_traj/02_sites_clustered.tsv"))[
  directional_pass == TRUE & !is.na(cluster_id),
  .(
    site_id, transcript_id, apa_gene_field = gene_id, locus, strand,
    discovery_age_qvalue = age_qvalue,
    discovery_delta_logit = delta_logit,
    discovery_amplitude_logit = amplitude_logit,
    cluster_id, cluster_silhouette, cluster_consensus
  )
]
discovery[, `:=`(
  source_order = .I,
  transcript_base = strip_version(transcript_id)
)]

# GENCODE v40 transcript-to-gene table (missing from the current tree; see missing_inputs)
tx2gene <- fread(file.path(Sys.getenv("APA_ROOT"), "ref/gencode.v40.tx2gene.tsv"))[, .(
  transcript_base = strip_version(tx_name),
  gencode_transcript_id = tx_name,
  gene_symbol = gene_name,
  gencode_gene_id = gene_id,
  gencode_gene_base = strip_version(gene_id),
  gene_type,
  transcript_type
)]
tx2gene <- unique(tx2gene)
stopifnot(tx2gene[, uniqueN(gencode_gene_base), by = transcript_base][, max(V1)] == 1L)
tx2gene <- unique(tx2gene, by = "transcript_base")

universe <- merge(
  discovery,
  tx2gene,
  by = "transcript_base",
  all.x = TRUE,
  sort = FALSE
)
setorder(universe, source_order)

symbol_mapping <- unique(universe[
  !is.na(gene_symbol) & gene_symbol != "" &
    !is.na(gencode_gene_base) & gencode_gene_base != "",
  .(gene_symbol, gencode_gene_base)
])[, .(symbol_gene_id_count = .N), by = gene_symbol]
universe <- merge(universe, symbol_mapping, by = "gene_symbol", all.x = TRUE, sort = FALSE)
setorder(universe, source_order)

# --- add RNA and protein measurement evidence ---

protein_annotation <- fread(file.path(mod, "prot_feature_annotation.tsv"))
protein_annotation[, gene_symbol := trimws(Gene.name)]
protein_by_symbol <- protein_annotation[
  !is.na(gene_symbol) & gene_symbol != "",
  .(
    protein_feature_count = .N,
    protein_accessions = paste(Protein.accession, collapse = ";"),
    protein_accession = if (.N == 1L) Protein.accession[1] else NA_character_,
    protein_description = if (.N == 1L) Protein.description[1] else NA_character_,
    protein_coverage = if (.N == 1L) Coverage[1] else NA_real_,
    protein_peptides = if (.N == 1L) Peptides[1] else NA_integer_,
    protein_unique_peptides = if (.N == 1L) Unique.peptides[1] else NA_integer_
  ),
  by = gene_symbol
]

universe <- merge(universe, protein_by_symbol, by = "gene_symbol", all.x = TRUE, sort = FALSE)
setorder(universe, source_order)

universe[, `:=`(
  gencode_mapped = !is.na(gencode_gene_base) & gencode_gene_base != "",
  rna_available = gencode_gene_base %chin% rna_gene_base,
  protein_available = !is.na(protein_feature_count),
  gene_symbol_unambiguous = !is.na(symbol_gene_id_count) & symbol_gene_id_count == 1L,
  protein_mapping_unambiguous = !is.na(protein_feature_count) & protein_feature_count == 1L
)]
universe[, sensitivity_universe :=
  gencode_mapped & rna_available & protein_available &
    gene_symbol_unambiguous & protein_mapping_unambiguous]
universe[, primary_universe :=
  sensitivity_universe & !is.na(protein_unique_peptides) &
    protein_unique_peptides >= min_unique_peptides]

universe[, exclusion_reason := ""]
universe[, exclusion_reason := append_reason(
  exclusion_reason, !gencode_mapped, "no_gencode_mapping"
)]
universe[, exclusion_reason := append_reason(
  exclusion_reason, gencode_mapped & !rna_available, "no_rna_measurement"
)]
universe[, exclusion_reason := append_reason(
  exclusion_reason, gencode_mapped & !protein_available, "no_protein_measurement"
)]
universe[, exclusion_reason := append_reason(
  exclusion_reason, protein_available & !gene_symbol_unambiguous,
  "gene_symbol_maps_multiple_ensembl_genes"
)]
universe[, exclusion_reason := append_reason(
  exclusion_reason, protein_available & !protein_mapping_unambiguous,
  "multiple_protein_features"
)]
universe[, exclusion_reason := append_reason(
  exclusion_reason,
  sensitivity_universe & protein_unique_peptides < min_unique_peptides,
  "protein_unique_peptides_lt2"
)]
universe[exclusion_reason == "", exclusion_reason := "included_primary"]

universe[, source_order := NULL]
setcolorder(universe, c(
  "site_id", "transcript_id", "transcript_base", "apa_gene_field",
  "gencode_transcript_id", "gencode_gene_id", "gencode_gene_base",
  "gene_symbol", "gene_type", "transcript_type", "locus", "strand",
  "cluster_id", "discovery_age_qvalue", "discovery_delta_logit",
  "discovery_amplitude_logit", "cluster_silhouette", "cluster_consensus",
  "rna_available", "protein_available", "symbol_gene_id_count",
  "protein_feature_count", "protein_accessions", "protein_accession",
  "protein_description", "protein_coverage", "protein_peptides",
  "protein_unique_peptides", "gencode_mapped", "gene_symbol_unambiguous",
  "protein_mapping_unambiguous", "sensitivity_universe", "primary_universe",
  "exclusion_reason"
))

# --- matched-sample support profile ---

sample_manifest[, age_model := log10(pmax(age_years + 1, 1e-6))]
grid_model <- seq(
  min(sample_manifest$age_model),
  max(sample_manifest$age_model),
  length.out = grid_size
)
grid_age <- 10^grid_model - 1
age_range <- diff(range(sample_manifest$age_model))

inner_membership <- abs(outer(grid_model, sample_manifest$age_model, "-")) <=
  inner_window * age_range
outer_membership <- abs(outer(grid_model, sample_manifest$age_model, "-")) <=
  outer_window * age_range

support <- data.table(
  grid_index = seq_len(grid_size) - 1L,
  age_model = grid_model,
  grid_age = grid_age,
  local_n = rowSums(outer_membership),
  local_n_inner = rowSums(inner_membership)
)
support[, is_supported := local_n >= outer_min & local_n_inner >= inner_min]
support[, `:=`(
  inner_window = inner_window,
  outer_window = outer_window,
  inner_min = inner_min,
  outer_min = outer_min
)]

support[, support_block := rleid(is_supported)]
support_segments <- support[is_supported == TRUE, .(
  start_age = min(grid_age),
  end_age = max(grid_age),
  n_grid = .N
), by = support_block]
support[, support_block := NULL]

# --- sample and universe QC ---

protein_age_years <- fifelse(
  sample_manifest$protein_age < 0,
  sample_manifest$protein_age / 52,
  sample_manifest$protein_age
)

sample_qc <- rbindlist(list(
  data.table(metric = c(
    "n_samples", "n_unique_canonical_ids", "n_prenatal", "n_postnatal",
    "age_min_years", "age_median_years", "age_max_years",
    "rna_rin_min", "rna_rin_median", "rna_rin_max",
    "n_age_discordant_gt1y", "n_md5_verified_artifacts",
    "apa_features", "rna_features", "protein_features",
    "locked_discovery_sites", "sensitivity_universe_sites",
    "sensitivity_universe_genes", "primary_universe_sites",
    "primary_universe_genes", "n_supported_grid_points",
    "n_support_segments", "support_segment_boundaries"
  ), stratum = "all", value = as.character(c(
    nrow(sample_manifest), uniqueN(sample_manifest$canonical_id),
    sum(sample_manifest$age_years < 0), sum(sample_manifest$age_years >= 0),
    min(sample_manifest$age_years), median(sample_manifest$age_years),
    max(sample_manifest$age_years), min(sample_manifest$rna_rin),
    median(sample_manifest$rna_rin), max(sample_manifest$rna_rin),
    sum(abs(sample_manifest$age_years - protein_age_years) > 1),
    nrow(input_manifest), nrow(apa), nrow(rna), nrow(protein),
    nrow(universe), sum(universe$sensitivity_universe),
    uniqueN(universe[sensitivity_universe == TRUE, gencode_gene_base]),
    sum(universe$primary_universe),
    uniqueN(universe[primary_universe == TRUE, gencode_gene_base]),
    sum(support$is_supported), nrow(support_segments),
    paste(sprintf("%.3f..%.3f", support_segments$start_age,
                  support_segments$end_age), collapse = ";")
  ))),
  sample_manifest[, .(metric = "sex_count", stratum = sex,
                      value = as.character(.N)), by = sex][, sex := NULL],
  sample_manifest[, .(metric = "life_stage_count", stratum = life_stage,
                      value = as.character(.N)), by = life_stage][, life_stage := NULL],
  sample_manifest[, .(metric = "rna_batch_count", stratum = rna_batch,
                      value = as.character(.N)), by = rna_batch][, rna_batch := NULL],
  sample_manifest[, .(metric = "protein_digest_batch_count",
                      stratum = protein_digest_batch,
                      value = as.character(.N)),
                  by = protein_digest_batch][, protein_digest_batch := NULL],
  sample_manifest[, .(metric = "protein_disease_class_count",
                      stratum = as.character(protein_disease_class),
                      value = as.character(.N)),
                  by = protein_disease_class][, protein_disease_class := NULL]
), use.names = TRUE)

# --- write outputs ---

fwrite(universe, file.path(mod, "01_analysis_universe.tsv"), sep = "\t", na = "NA")
fwrite(sample_qc, file.path(mod, "01_sample_qc.tsv"), sep = "\t", na = "NA")
fwrite(support, file.path(mod, "01_support_profile.tsv"), sep = "\t", na = "NA")
