#!/usr/bin/env Rscript
# APA-QTL colocalization with stage-matched eQTL, GTEx v11 brain apaQTL and ROSMAP pQTL.
# Compares formal coloc pairs with published 3'aQTL-trait associations.
# Then runs the NDUFA6-SCZ and novel-pair follow-up.

suppressMessages({
  library(coloc)
  library(data.table)
  library(parallel)
  library(yaml)
  library(Rsamtools)
  library(GenomicRanges)
})

# --- helpers ---

.expand_env <- function(x) rapply(x, function(s) { for (v in regmatches(s, gregexpr("\\$\\{[A-Za-z_]+\\}", s))[[1]]) s <- sub(v, Sys.getenv(gsub("[${}]", "", v)), s, fixed = TRUE); s }, classes = "character", how = "replace")

source_identity <- function(path) {
  info <- file.info(path)
  list(
    path = normalizePath(path),
    bytes = unname(as.numeric(info$size)),
    mtime_unix = unname(as.numeric(info$mtime))
  )
}

atomic_fwrite <- function(table, path) {
  temporary_path <- paste0(path, ".tmp")
  fwrite(table, temporary_path, sep = "\t")
  file.rename(temporary_path, path)
}

empty_coloc <- function(candidates, status) {
  data.table(
    candidate_id = candidates$candidate_id,
    phenotype_id = NA_character_,
    stage = candidates$stage,
    eQTL_coloc = FALSE,
    eqtl_coloc_pp4 = NA_real_,
    eqtl_coloc_status = status,
    n_shared_snps = 0L,
    top_apa_snp = NA_character_,
    top_eqtl_snp = NA_character_
  )
}

variant_key <- function(ids) {
  ids <- sub("_b38$", "", ids)
  ids <- gsub("_", ":", ids, fixed = TRUE)
  parts <- tstrsplit(ids, ":", fixed = TRUE)
  alleles_1 <- toupper(parts[[3]])
  alleles_2 <- toupper(parts[[4]])
  paste(parts[[1]], parts[[2]], pmin(alleles_1, alleles_2), pmax(alleles_1, alleles_2), sep = ":")
}

run_coloc <- function(apa, external, apa_n, external_n, external_has_maf = TRUE, external_sdY = NULL) {
  apa <- apa[!is.na(variant_key) & !duplicated(variant_key)]
  external <- external[!is.na(variant_key) & !duplicated(variant_key)]
  common <- intersect(apa$variant_key, external$variant_key)
  if (length(common) < min_shared_snps) {
    return(list(status = "insufficient_shared_snps", n = length(common), n_before_filter = length(common), n_dropped_incomplete = 0L))
  }
  apa <- apa[match(common, variant_key)]
  external <- external[match(common, variant_key)]
  complete <- is.finite(apa$slope) & is.finite(apa$slope_se) & apa$slope_se > 0 &
    is.finite(apa$af) & apa$af > 0 & apa$af < 1 &
    is.finite(external$slope) & is.finite(external$slope_se) & external$slope_se > 0
  if (external_has_maf) {
    complete <- complete & is.finite(external$af) & external$af > 0 & external$af < 1
  }
  n_before_filter <- length(common)
  apa <- apa[complete]
  external <- external[complete]
  common <- common[complete]
  if (length(common) < min_shared_snps) {
    return(list(
      status = "insufficient_shared_snps_after_complete_case_filter",
      n = length(common), n_before_filter = n_before_filter,
      n_dropped_incomplete = n_before_filter - length(common)
    ))
  }
  dataset1 <- list(
    beta = apa$slope,
    varbeta = apa$slope_se^2,
    snp = common,
    type = "quant",
    N = apa_n,
    MAF = apa$af
  )
  dataset2 <- list(
    beta = external$slope,
    varbeta = external$slope_se^2,
    snp = common,
    type = "quant",
    N = external_n
  )
  if (external_has_maf) dataset2$MAF <- external$af
  if (!is.null(external_sdY)) dataset2$sdY <- external_sdY
  fit <- tryCatch({
    invisible(capture.output(result <- coloc.abf(
      dataset1 = dataset1,
      dataset2 = dataset2,
      p1 = as.numeric(params$p1),
      p2 = as.numeric(params$p2),
      p12 = as.numeric(params$p12)
    )))
    result
  },
    error = function(error) error
  )
  if (inherits(fit, "error")) {
    return(list(
      status = paste0("coloc_error:", conditionMessage(fit)), n = length(common),
      n_before_filter = n_before_filter,
      n_dropped_incomplete = n_before_filter - length(common)
    ))
  }
  list(
    status = "ok", n = length(common), summary = fit$summary,
    n_before_filter = n_before_filter,
    n_dropped_incomplete = n_before_filter - length(common)
  )
}

is_true <- function(value) {
  length(value) == 1L && !is.na(value) &&
    tolower(trimws(as.character(value))) %chin% c("true", "1", "yes")
}

is_present <- function(value) {
  length(value) == 1L && !is.na(value) &&
    !tolower(trimws(as.character(value))) %chin% c("", "na", "nan", "false", "none", "not_supported")
}

collapse_unique <- function(values) {
  paste(sort(unique(na.omit(as.character(values)))), collapse = "; ")
}

read_phenotype_rows <- function(path, phenotype) {
  command <- sprintf(
    "zcat %s | awk -F '\\t' -v target=%s 'NR==1 || $1==target'",
    shQuote(path), shQuote(phenotype)
  )
  fread(cmd = command)
}

canonical_variant <- function(x) {
  x <- sub("_b38$", "", x)
  gsub("_", ":", x, fixed = TRUE)
}

run_coloc_pair <- function(pair_name, first, second, first_type, second_type,
                           first_n, second_n, first_s = NULL, second_s = NULL,
                           p1_value = p1, p2_value = p2, p12_value = p12) {
  shared <- intersect(first$variant_id, second$variant_id)
  first <- first[match(shared, variant_id)]
  second <- second[match(shared, variant_id)]

  dataset1 <- list(
    beta = first$beta,
    varbeta = first$se^2,
    snp = shared,
    type = first_type,
    N = first_n,
    MAF = pmin(first$af, 1 - first$af)
  )
  dataset2 <- list(
    beta = second$beta,
    varbeta = second$se^2,
    snp = shared,
    type = second_type,
    N = second_n,
    MAF = pmin(second$af, 1 - second$af)
  )
  if (!is.null(first_s)) dataset1$s <- first_s
  if (!is.null(second_s)) dataset2$s <- second_s

  fit <- coloc.abf(
    dataset1, dataset2,
    p1 = p1_value, p2 = p2_value, p12 = p12_value
  )
  posterior <- as.data.table(fit$results)[order(-SNP.PP.H4)]
  posterior[, `:=`(
    pair = pair_name,
    posterior_rank = seq_len(.N),
    cumulative_conditional_SNP_PP_H4 = cumsum(SNP.PP.H4)
  )]
  cutoff_rank <- which(
    posterior$cumulative_conditional_SNP_PP_H4 >= credible_set_probability
  )[1]
  posterior[, in_conditional_95pct_H4_set := posterior_rank <= cutoff_rank]

  summary <- data.table(
    pair = pair_name,
    n_shared_snps = length(shared),
    PP.H0 = unname(fit$summary[["PP.H0.abf"]]),
    PP.H1 = unname(fit$summary[["PP.H1.abf"]]),
    PP.H2 = unname(fit$summary[["PP.H2.abf"]]),
    PP.H3 = unname(fit$summary[["PP.H3.abf"]]),
    PP.H4 = unname(fit$summary[["PP.H4.abf"]]),
    coloc_supported = unname(fit$summary[["PP.H4.abf"]]) >= pp4_threshold,
    conditional_H4_95pct_set_size = cutoff_rank,
    candidate_conditional_SNP_PP_H4 = posterior[snp == candidate_variant, SNP.PP.H4],
    candidate_in_conditional_95pct_H4_set = posterior[
      snp == candidate_variant, in_conditional_95pct_H4_set
    ],
    conditional_H4_set_interpretable = unname(fit$summary[["PP.H4.abf"]]) >= pp4_threshold
  )
  list(summary = summary, posterior = posterior)
}

source_contract_maf <- function(merged, columns) {
  available <- columns[columns %in% names(merged)]
  if (length(available) == 0L) {
    return(rep(maf_default, nrow(merged)))
  }
  maf <- pmin(as.numeric(merged[[available[[1L]]]]), 1 - as.numeric(merged[[available[[1L]]]]))
  maf[is.na(maf)] <- maf_na_replacement
  pmax(maf_epsilon, pmin(1 - maf_epsilon, maf))
}

run_local_gwas_coloc <- function(merged, gwas_type, gwas_n, gwas_prop,
                                 p1_value, p2_value, p12_value) {
  gwas_n_value <- if ("N" %in% names(merged)) {
    as.numeric(merged$N)
  } else {
    gwas_n
  }
  gwas_maf <- source_contract_maf(
    merged,
    c("EAF.gwas", "MAF.gwas", "af.gwas", "EAF", "MAF", "af")
  )
  qtl_maf <- source_contract_maf(merged, c("EAF.qtl", "MAF.qtl", "af.qtl"))
  dataset1 <- list(
    snp = merged$variant_id,
    beta = merged$BETA.gwas,
    varbeta = merged$SE.gwas^2,
    type = gwas_type,
    N = gwas_n_value,
    MAF = gwas_maf
  )
  if (gwas_type == "cc") dataset1$s <- gwas_prop

  dataset2 <- list(
    snp = merged$variant_id,
    beta = merged$BETA.qtl,
    varbeta = merged$SE.qtl^2,
    type = "quant",
    N = qtl_sample_sizes[[unique(merged$QTL_ID)]],
    MAF = qtl_maf
  )

  coloc.abf(
    dataset1, dataset2,
    p1 = p1_value, p2 = p2_value, p12 = p12_value
  )
}

# --- config ---

out_dir <- file.path(Sys.getenv("APA_ROOT"), "data/07_coloc")
config_path <- file.path(Sys.getenv("APA_ROOT"), "release/07_coloc/func/config.yml")
candidate_path <- file.path(out_dir, "02_signal_mechanisms.tsv.gz")
output_path <- file.path(out_dir, "internal_eqtl.tsv")
config <- .expand_env(read_yaml(config_path))
params <- config$parameters$internal_eqtl_coloc
analysis_version <- config$project$analysis_version
pipeline_run_id <- Sys.getenv("APA_FUNC_RUN_ID", "manual")
workers <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

nominal_paths <- list(
  prenatal = list(
    apa = config$inputs$prenatal_nominal_apaqtl$path,
    eqtl = config$inputs$prenatal_internal_eqtl_nominal$path,
    apa_sample = config$inputs$prenatal_qtl_samples$path,
    eqtl_sample = config$inputs$prenatal_internal_eqtl_samples$path
  ),
  postnatal = list(
    apa = config$inputs$postnatal_nominal_apaqtl$path,
    eqtl = config$inputs$postnatal_internal_eqtl_nominal$path,
    apa_sample = config$inputs$postnatal_qtl_samples$path,
    eqtl_sample = config$inputs$postnatal_internal_eqtl_samples$path
  )
)

# --- load selected signals ---

candidates <- unique(fread(candidate_path, select = c("candidate_id", "phenotype_id", "stage")))

# --- stage-matched APA/eQTL colocalization ---

rows <- list()
for (stage_name in names(nominal_paths)) {
  stage_candidates <- candidates[stage == stage_name]
  if (nrow(stage_candidates) == 0L) next
  paths <- nominal_paths[[stage_name]]
  apa_n <- nrow(fread(paths$apa_sample, header = FALSE))
  eqtl_n <- nrow(fread(paths$eqtl_sample, header = FALSE))
  apa <- fread(
    paths$apa,
    sep = "\t",
    select = c("phenotype_id", "variant_id", "af", "ma_samples", "slope", "slope_se")
  )
  apa <- apa[phenotype_id %chin% unique(stage_candidates$phenotype_id)]
  if (nrow(apa) == 0L) {
    rows[[length(rows) + 1L]] <- empty_coloc(stage_candidates, "no_apa_nominal_rows")
    next
  }
  eqtl_all <- fread(
    paths$eqtl,
    sep = "\t",
    select = c("phenotype_id", "variant_id", "af", "ma_samples", "slope", "slope_se")
  )
  setkey(eqtl_all, variant_id)
  stage_rows <- mclapply(
    seq_len(nrow(stage_candidates)),
    function(i) {
      candidate <- stage_candidates[i]
      apa_i <- apa[phenotype_id == candidate$phenotype_id]
      if (nrow(apa_i) == 0L) return(empty_coloc(candidate, "no_apa_nominal_rows"))

      shared_eqtl <- eqtl_all[J(apa_i$variant_id), nomatch = 0L]
      if (nrow(shared_eqtl) == 0L) return(empty_coloc(candidate, "no_eqtl_overlap"))

      rows <- list()
      for (egene_id in unique(shared_eqtl$phenotype_id)) {
        eqtl_i <- shared_eqtl[phenotype_id == egene_id]
        common <- intersect(apa_i$variant_id, eqtl_i$variant_id)
        if (length(common) < as.integer(params$min_shared_snps)) next

        apa_m <- apa_i[variant_id %chin% common][match(common, variant_id)]
        eqtl_m <- eqtl_i[variant_id %chin% common][match(common, variant_id)]
        invisible(capture.output(
          fit <- coloc.abf(
            dataset1 = list(
              beta = apa_m$slope, varbeta = apa_m$slope_se^2, snp = common,
              type = "quant", N = apa_n, MAF = apa_m$af
            ),
            dataset2 = list(
              beta = eqtl_m$slope, varbeta = eqtl_m$slope_se^2, snp = common,
              type = "quant", N = eqtl_n, MAF = eqtl_m$af
            ),
            p1 = as.numeric(params$p1),
            p2 = as.numeric(params$p2),
            p12 = as.numeric(params$p12)
          )
        ))
        summary <- fit$summary
        rows[[length(rows) + 1L]] <- data.table(
          candidate_id = candidate$candidate_id,
          phenotype_id = egene_id,
          stage = stage_name,
          eQTL_coloc = summary[["PP.H4.abf"]] >= as.numeric(params$pp4_threshold),
          eqtl_coloc_pp4 = summary[["PP.H4.abf"]],
          eqtl_coloc_status = "ok",
          n_shared_snps = length(common),
          top_apa_snp = apa_m$variant_id[which.max(abs(apa_m$slope / apa_m$slope_se))],
          top_eqtl_snp = eqtl_m$variant_id[which.max(abs(eqtl_m$slope / eqtl_m$slope_se))]
        )
      }
      if (length(rows)) rbindlist(rows, fill = TRUE) else empty_coloc(candidate, "insufficient_shared_snps")
    },
    mc.cores = workers,
    mc.preschedule = TRUE
  )
  rows[[length(rows) + 1L]] <- rbindlist(stage_rows, fill = TRUE)
}

# --- provenance and output ---

output <- rbindlist(rows, fill = TRUE)
provenance <- rbindlist(lapply(names(nominal_paths), function(stage_name) {
  paths <- nominal_paths[[stage_name]]
  apa <- source_identity(paths$apa)
  eqtl <- source_identity(paths$eqtl)
  data.table(
    stage = stage_name,
    apa_nominal_path = apa$path,
    apa_nominal_bytes = apa$bytes,
    apa_nominal_mtime_unix = apa$mtime_unix,
    eqtl_nominal_path = eqtl$path,
    eqtl_nominal_bytes = eqtl$bytes,
    eqtl_nominal_mtime_unix = eqtl$mtime_unix,
    apa_sample_n = nrow(fread(paths$apa_sample, header = FALSE)),
    eqtl_sample_n = nrow(fread(paths$eqtl_sample, header = FALSE))
  )
}))
output <- merge(output, provenance, by = "stage", all.x = TRUE, sort = FALSE)
output[, `:=`(
  analysis_version = analysis_version,
  pipeline_run_id = pipeline_run_id,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  config_path = normalizePath(config_path)
)]
setcolorder(output, c(
  "candidate_id", "phenotype_id", "stage", "eQTL_coloc", "eqtl_coloc_pp4", "eqtl_coloc_status",
  "n_shared_snps", "top_apa_snp", "top_eqtl_snp", "apa_nominal_path", "apa_nominal_bytes",
  "apa_nominal_mtime_unix", "eqtl_nominal_path", "eqtl_nominal_bytes", "eqtl_nominal_mtime_unix",
  "apa_sample_n", "eqtl_sample_n", "analysis_version", "pipeline_run_id", "generated_at_utc", "config_path"
))
atomic_fwrite(output, output_path)

# --- external QTL colocalization ---

inputs <- config$inputs
resources <- config$reference_resources
params <- config$parameters$external_qtl_coloc
external_cache_mode <- "force_rebuild"
all_output <- file.path(out_dir, "external_all.tsv")
summary_output <- file.path(out_dir, "external_summary.tsv")
source_audit_output <- file.path(out_dir, "external_sources.tsv")

pp4_threshold <- as.numeric(params$pp4_threshold)
min_shared_snps <- as.integer(params$min_shared_snps)

# --- candidate genes and local APA phenotypes ---

candidates <- fread(candidate_path)
candidates <- unique(candidates[, .(candidate_id, stage, phenotype_id, transcript_id, gene_name, variant_id)])
candidates <- candidates[stage %chin% c("prenatal", "postnatal")]
candidates[, local_event_coordinate := tstrsplit(phenotype_id, "|", fixed = TRUE)[[3]]]

gtf_command <- sprintf("zcat %s | awk -F '\\t' '$3==\"transcript\"'", shQuote(resources$gene_annotation_gtf$path))
gtf <- fread(cmd = gtf_command, sep = "\t", header = FALSE, select = 9, quote = "")
gtf[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', V9)]
gtf[, gene_id_base := sub("\\.[0-9]+$", "", gene_id)]
gtf[, transcript_id := sub('.*transcript_id "([^"]+)".*', '\\1', V9)]
transcript_map <- unique(gtf[, .(transcript_id, gene_id, gene_id_base)])
candidates <- merge(candidates, transcript_map, by = "transcript_id", all.x = TRUE)

group_paths <- list(
  prenatal = file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/phenotype.Prenatal.group"),
  postnatal = file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/phenotype.Postnatal.group")
)
for (stage_name in names(group_paths)) {
  group_map <- fread(
    group_paths[[stage_name]], sep = "\t", header = FALSE,
    col.names = c("phenotype_id", "local_apa_group_id")
  )
  group_map[, group_gene_name := tstrsplit(phenotype_id, "|", fixed = TRUE)[[2L]]]
  exact_group <- setNames(group_map$local_apa_group_id, group_map$phenotype_id)
  groups_by_gene <- group_map[, .(group_ids = list(unique(local_apa_group_id))), by = group_gene_name]
  groups_by_gene <- setNames(groups_by_gene$group_ids, groups_by_gene$group_gene_name)
  candidate_rows <- which(candidates$stage == stage_name)
  resolved_group <- unname(exact_group[candidates$phenotype_id[candidate_rows]])
  unresolved <- is.na(resolved_group)
  if (any(unresolved)) {
    resolved_group[unresolved] <- vapply(
      candidates$gene_name[candidate_rows][unresolved],
      function(gene_name) {
        matches <- groups_by_gene[[gene_name]]
        if (length(matches) == 1L) matches else NA_character_
      },
      character(1)
    )
  }
  exact_event_present <- candidates$phenotype_id[candidate_rows] %chin% group_map$phenotype_id
  candidates[candidate_rows, `:=`(
    local_apa_group_id = resolved_group,
    local_apa_phenotype_id = fifelse(!is.na(resolved_group), phenotype_id, NA_character_),
    local_apa_match_scope = fifelse(
      exact_event_present, "exact_event_in_01QTL_group",
      "same_gene_group_without_exact_event"
    )
  )]
}
if (anyNA(candidates$local_apa_group_id) || anyNA(candidates$local_apa_phenotype_id)) {
  stop("Every candidate must resolve to a gene-matched 01_QTL correction group")
}

apa_paths <- list(
  prenatal = inputs$prenatal_nominal_apaqtl$path,
  postnatal = inputs$postnatal_nominal_apaqtl$path
)
apa_sample_sizes <- c(
  prenatal = length(readLines(inputs$prenatal_qtl_samples$path)),
  postnatal = length(readLines(inputs$postnatal_qtl_samples$path))
)

apa_tables <- list()
for (stage_name in names(apa_paths)) {
  cache_path <- file.path(out_dir, paste0("local_apa_", stage_name, ".tsv.gz"))
  pattern_path <- paste0(cache_path, ".phenotypes")
  fwrite(
    data.table(phenotype_id = unique(na.omit(candidates[stage == stage_name]$local_apa_phenotype_id))),
    pattern_path, col.names = FALSE
  )
  awk <- "NR==FNR {wanted[$1]=1; next} FNR==1 || ($1 in wanted)"
  temporary_path <- paste0(cache_path, ".tmp")
  unlink(temporary_path)
  command <- sprintf(
    "zcat %s | awk -F '\\t' '%s' %s - | gzip -c > %s",
    shQuote(apa_paths[[stage_name]]), awk, shQuote(pattern_path), shQuote(temporary_path)
  )
  status <- system(command)
  if (status != 0) stop("APA subset extraction failed: ", apa_paths[[stage_name]])
  file.rename(temporary_path, cache_path)
  unlink(pattern_path)
  apa_tables[[stage_name]] <- fread(cache_path)
  apa_tables[[stage_name]][, variant_key := variant_key(variant_id)]
  setkeyv(apa_tables[[stage_name]], "phenotype_id")
}

# --- GTEx v11 apaQTL colocalization ---

gtex_sources <- list(
  list(name = "GTEx_v11_Brain_Cortex_apaQTL", resource = resources$gtex_v11_brain_cortex_apaqtl, qtl_type = "apaQTL", context = "adult_brain"),
  list(name = "GTEx_v11_Brain_BA9_apaQTL", resource = resources$gtex_v11_brain_ba9_apaqtl, qtl_type = "apaQTL", context = "adult_brain")
)

all_rows <- list()
source_audit <- list()
setDTthreads(1L)
for (source in gtex_sources) {
  cache_path <- file.path(out_dir, paste0(source$name, ".tsv.gz"))
  source_path <- source$resource$path
  covariate_path <- source$resource$covariates
  pattern_path <- paste0(cache_path, ".genes")
  fwrite(data.table(gene_id = unique(na.omit(candidates$gene_id_base))), pattern_path, col.names = FALSE)
  awk <- "NR==FNR {wanted[$1]=1; next} FNR==1 {print; next} {key=$1; sub(/_.*/,\"\",key); sub(/\\.[0-9]+$/,\"\",key); if(key in wanted) print}"
  temporary_path <- paste0(cache_path, ".tmp")
  unlink(temporary_path)
  command <- sprintf(
    "zcat %s | awk -F '\\t' '%s' %s - | gzip -c > %s",
    shQuote(source_path), awk, shQuote(pattern_path), shQuote(temporary_path)
  )
  status <- system(command)
  if (status != 0) stop("Reference subset extraction failed: ", source_path)
  file.rename(temporary_path, cache_path)
  unlink(pattern_path)
  external <- fread(cache_path)
  setnames(external, "phenotype_id", "external_phenotype_id")
  external[, gene_id_base := sub("\\.[0-9]+$", "", sub("_.*", "", external_phenotype_id))]
  external[, external_event_coordinate := sub("^[^_]+_", "", external_phenotype_id)]
  external[, variant_key := variant_key(variant_id)]
  setkeyv(external, "gene_id_base")
  external_n <- length(strsplit(readLines(covariate_path, n = 1), "\t", fixed = TRUE)[[1]]) - 1L

  source_audit[[length(source_audit) + 1]] <- data.table(
    source = source$name,
    qtl_type = source$qtl_type,
    source_context = source$context,
    sample_size = external_n,
    candidate_gene_rows = nrow(external),
    variant_mapping = "chr:position:unordered_allele_pair",
    phenotype_matching = "same_gene_GRCh38_apa_interval_overlap_or_same_gene_other_event",
    population_relation = "independent_external_population_not_matched_to_local_cohort",
    status = "tested"
  )

  source_rows <- mclapply(seq_len(nrow(candidates)), function(candidate_index) {
    candidate <- candidates[candidate_index]
    apa <- apa_tables[[candidate$stage]][.(candidate$local_apa_phenotype_id), nomatch = 0L]
    ext <- external[.(candidate$gene_id_base), nomatch = 0L]
    if (nrow(apa) == 0 || nrow(ext) == 0) {
      return(data.table(
        candidate_id = candidate$candidate_id, stage = candidate$stage,
        source = source$name, qtl_type = source$qtl_type, source_context = source$context,
        population_relation = "independent_external_population_not_matched_to_local_cohort",
        local_event_coordinate = candidate$local_event_coordinate,
        external_event_coordinate = NA_character_, event_match_scope = NA_character_,
        external_phenotype_id = NA_character_, n_shared_snps = 0L,
        PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_, PP.H3 = NA_real_, PP.H4 = NA_real_,
        coloc_supported = FALSE, external_validation_eligible = FALSE,
        evidence_role = if (nrow(apa) == 0) "no_local_apa_rows" else "no_external_gene_rows",
        coloc_status = if (nrow(apa) == 0) "no_local_apa_rows" else "no_external_gene_rows"
      ))
    }
    rbindlist(lapply(unique(ext$external_phenotype_id), function(external_id) {
      ext_one <- ext[external_phenotype_id == external_id]
      fit <- run_coloc(apa, ext_one, apa_sample_sizes[[candidate$stage]], external_n)
      external_coord <- ext_one$external_event_coordinate[[1L]]
      overlap <- FALSE
      if (!is.na(external_coord) && !is.na(candidate$local_event_coordinate)) {
        l <- tstrsplit(candidate$local_event_coordinate, "[:-]")
        e <- tstrsplit(external_coord, "[:-]")
        if (length(l) >= 3L && length(e) >= 3L) {
          lstart <- suppressWarnings(as.numeric(l[[2L]]))
          lend <- suppressWarnings(as.numeric(l[[3L]]))
          estart <- suppressWarnings(as.numeric(e[[2L]]))
          eend <- suppressWarnings(as.numeric(e[[3L]]))
          overlap <- !is.na(lstart) && !is.na(lend) && !is.na(estart) && !is.na(eend) &&
            identical(l[[1L]], e[[1L]]) && lstart <= eend && estart <= lend
        }
      }
      event_scope <- if (overlap) "same_gene_apa_overlap" else "same_gene_other_event"
      supported <- fit$status == "ok" && fit$summary[["PP.H4.abf"]] >= pp4_threshold
      eligible <- supported && event_scope == "same_gene_apa_overlap"
      coloc_summary <- if (fit$status == "ok") fit$summary else NULL
      data.table(
        candidate_id = candidate$candidate_id, stage = candidate$stage,
        source = source$name, qtl_type = source$qtl_type, source_context = source$context,
        population_relation = "independent_external_population_not_matched_to_local_cohort",
        local_event_coordinate = candidate$local_event_coordinate,
        external_event_coordinate = external_coord, event_match_scope = event_scope,
        external_phenotype_id = external_id, n_shared_snps = fit$n,
        PP.H0 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H0.abf"]],
        PP.H1 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H1.abf"]],
        PP.H2 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H2.abf"]],
        PP.H3 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H3.abf"]],
        PP.H4 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H4.abf"]],
        coloc_supported = supported,
        external_validation_eligible = eligible,
        evidence_role = if (event_scope == "same_gene_apa_overlap") {
          "secondary_cross_population_same_gene_apa_overlap_replication"
        } else {
          "secondary_cross_population_same_gene_other_apa_context"
        },
        coloc_status = fit$status
      )
    }), fill = TRUE)
  }, mc.cores = workers, mc.preschedule = TRUE)
  all_rows[[length(all_rows) + 1L]] <- rbindlist(source_rows, fill = TRUE)
  rm(external)
  gc()
}

# --- ROSMAP pQTL colocalization ---

rosmap <- fread(cmd = paste("zcat", shQuote(resources$rosmap_pqtl$path)))
rosmap[, gene_name := sub("^[^.]+\\.", "", Protein_GeneSymbol)]
rosmap[, external_phenotype_id := Protein_UniProt]
rosmap_chain <- resources$rosmap_pqtl$liftover_chain
rosmap_fasta <- resources$genome_fasta$path
rosmap_source_build <- resources$rosmap_pqtl$source_build
rosmap_target_build <- resources$rosmap_pqtl$target_build
if (!identical(rosmap_source_build, "GRCh37") || !identical(rosmap_target_build, "GRCh38")) {
  stop("ROSMAP source and target builds must be GRCh37 and GRCh38")
}

# GRCh37 to GRCh38 liftover with reference-allele validation
rosmap[, rosmap_row_id := .I]
liftover_bin <- Sys.which("liftOver")
input_path <- file.path(out_dir, "rosmap_grch37_input.bed")
mapped_path <- file.path(out_dir, "rosmap_grch38_mapped.bed")
unmapped_path <- file.path(out_dir, "rosmap_grch37_unmapped.bed")
lift_input <- rosmap[, .(
  chrom = paste0("chr", Chr), start = as.integer(BP) - 1L,
  end = as.integer(BP), rosmap_row_id
)]
fwrite(lift_input, input_path, sep = "\t", col.names = FALSE)
unlink(c(mapped_path, unmapped_path))
if (system2(liftover_bin, c(input_path, rosmap_chain, mapped_path, unmapped_path)) != 0L) {
  stop("ROSMAP liftover failed")
}
mapped <- if (file.exists(mapped_path) && file.info(mapped_path)$size > 0) {
  fread(mapped_path, header = FALSE, col.names = c("lifted_chr", "lifted_start0", "lifted_end", "rosmap_row_id"))
} else {
  data.table(lifted_chr = character(), lifted_start0 = integer(), lifted_end = integer(), rosmap_row_id = integer())
}
mapped[, rosmap_row_id := as.integer(rosmap_row_id)]
mapping_counts <- mapped[, .N, by = rosmap_row_id]
mapped <- mapping_counts[N == 1L][mapped, on = "rosmap_row_id", nomatch = 0L][
  , .(rosmap_row_id, lifted_chr, lifted_start0, lifted_end)
]
rosmap <- merge(rosmap, mapped, by = "rosmap_row_id", all.x = TRUE, sort = FALSE)
rosmap[, lifted_pos := as.integer(lifted_start0) + 1L]
rosmap[, liftover_status := fifelse(
  rosmap_row_id %in% mapping_counts[N > 1L]$rosmap_row_id, "ambiguous_liftover",
  fifelse(is.na(lifted_pos), "unmapped", "mapped")
)]

fasta <- FaFile(rosmap_fasta)
fasta_seqnames <- as.character(seqnames(scanFaIndex(fasta)))
rosmap[, lifted_reference := NA_character_]
in_fasta <- rosmap$liftover_status == "mapped" & rosmap$lifted_chr %chin% fasta_seqnames
if (any(in_fasta)) {
  open(fasta)
  gr <- GRanges(
    seqnames = rosmap$lifted_chr[in_fasta],
    ranges = IRanges::IRanges(start = rosmap$lifted_pos[in_fasta], width = 1L)
  )
  rosmap[in_fasta, lifted_reference := as.character(scanFa(fasta, param = gr))]
  close(fasta)
}

rosmap[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
direct_allele_match <- rosmap$lifted_reference == rosmap$A1 | rosmap$lifted_reference == rosmap$A2
reverse_A1 <- chartr("ACGT", "TGCA", rosmap$A1)
reverse_A2 <- chartr("ACGT", "TGCA", rosmap$A2)
reverse_allele_match <- rosmap$lifted_reference == reverse_A1 | rosmap$lifted_reference == reverse_A2
rosmap[, liftover_status := fifelse(
  liftover_status != "mapped", liftover_status,
  fifelse(nchar(A1) != 1L | nchar(A2) != 1L | A1 == A2, "non_snv_alleles",
    fifelse(is.na(lifted_reference), "target_contig_not_in_fasta",
      fifelse(direct_allele_match, "ok",
        fifelse(reverse_allele_match, "ok_reverse_complement", "reference_allele_mismatch")
      )
    )
  )
)]
rosmap[, variant_key := NA_character_]
valid <- rosmap$liftover_status %chin% c("ok", "ok_reverse_complement")
normalized_A1 <- ifelse(rosmap$liftover_status == "ok_reverse_complement", reverse_A1, rosmap$A1)
normalized_A2 <- ifelse(rosmap$liftover_status == "ok_reverse_complement", reverse_A2, rosmap$A2)
rosmap[valid, variant_key := paste0(
  lifted_chr, ":", lifted_pos, ":",
  pmin(normalized_A1[valid], normalized_A2[valid]), ":",
  pmax(normalized_A1[valid], normalized_A2[valid])
)]
rosmap_audit <- rosmap[, .(
  raw_rows = .N,
  lifted_b38_rows = sum(liftover_status %chin% c("ok", "ok_reverse_complement"), na.rm = TRUE),
  reverse_complement_rows = sum(liftover_status == "ok_reverse_complement", na.rm = TRUE),
  unmapped_rows = sum(liftover_status == "unmapped", na.rm = TRUE),
  ambiguous_liftover_rows = sum(liftover_status == "ambiguous_liftover", na.rm = TRUE),
  reference_allele_mismatch_rows = sum(liftover_status == "reference_allele_mismatch", na.rm = TRUE),
  other_liftover_excluded_rows = sum(!liftover_status %chin% c("ok", "ok_reverse_complement", "unmapped", "ambiguous_liftover", "reference_allele_mismatch"), na.rm = TRUE)
)]

setnames(rosmap, c("Beta", "SE"), c("slope", "slope_se"))
setkeyv(rosmap, "gene_name")
rosmap_n <- as.integer(resources$rosmap_pqtl$sample_size)
source_audit[[length(source_audit) + 1]] <- cbind(
  data.table(
    source = "ROSMAP_DLPFC_pQTL", qtl_type = "pQTL", source_context = "adult_DLPFC",
    sample_size = rosmap_n, candidate_gene_rows = rosmap[gene_name %chin% candidates$gene_name, .N],
    source_build = rosmap_source_build, target_build = rosmap_target_build,
    variant_mapping = "GRCh37_liftover_to_GRCh38:chr:position:unordered_allele_pair",
    phenotype_matching = "exact_same_gene_protein",
    population_relation = "independent_external_population_not_harmonized_to_local_cohort",
    status = "tested_after_GRCh37_to_GRCh38_liftover_with_reference_allele_validation"
  ),
  rosmap_audit
)

rosmap_rows <- mclapply(seq_len(nrow(candidates)), function(candidate_index) {
  candidate <- candidates[candidate_index]
  apa <- apa_tables[[candidate$stage]][.(candidate$local_apa_phenotype_id), nomatch = 0L]
  ext <- rosmap[.(candidate$gene_name), nomatch = 0L]
  if (nrow(apa) == 0 || nrow(ext) == 0) {
    return(data.table(
      candidate_id = candidate$candidate_id, stage = candidate$stage,
      source = "ROSMAP_DLPFC_pQTL", qtl_type = "pQTL", source_context = "adult_DLPFC",
      population_relation = "independent_external_population_not_harmonized_to_local_cohort",
      local_event_coordinate = candidate$local_event_coordinate,
      external_event_coordinate = NA_character_, event_match_scope = "same_gene_protein",
      external_phenotype_id = NA_character_, n_shared_snps = 0L,
      PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_, PP.H3 = NA_real_, PP.H4 = NA_real_,
      coloc_supported = FALSE, external_validation_eligible = FALSE,
      evidence_role = if (nrow(apa) == 0) "no_local_apa_rows" else "no_external_gene_rows",
      coloc_status = if (nrow(apa) == 0) "no_local_apa_rows" else "no_external_gene_rows"
    ))
  }
  rbindlist(lapply(unique(ext$external_phenotype_id), function(external_id) {
    fit <- run_coloc(
      apa, ext[external_phenotype_id == external_id],
      apa_sample_sizes[[candidate$stage]], rosmap_n,
      external_has_maf = FALSE, external_sdY = 1
    )
    coloc_summary <- if (fit$status == "ok") fit$summary else NULL
    supported <- fit$status == "ok" && fit$summary[["PP.H4.abf"]] >= pp4_threshold
    data.table(
      candidate_id = candidate$candidate_id, stage = candidate$stage,
      source = "ROSMAP_DLPFC_pQTL", qtl_type = "pQTL", source_context = "adult_DLPFC",
      population_relation = "independent_external_population_not_harmonized_to_local_cohort",
      local_event_coordinate = candidate$local_event_coordinate,
      external_event_coordinate = NA_character_, event_match_scope = "same_gene_protein",
      external_phenotype_id = external_id, n_shared_snps = fit$n,
      PP.H0 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H0.abf"]],
      PP.H1 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H1.abf"]],
      PP.H2 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H2.abf"]],
      PP.H3 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H3.abf"]],
      PP.H4 = if (is.null(coloc_summary)) NA_real_ else coloc_summary[["PP.H4.abf"]],
      coloc_supported = supported, external_validation_eligible = supported,
      evidence_role = "secondary_cross_population_same_gene_protein_support",
      coloc_status = if (fit$status == "ok") "ok_sdY_1_standardized_protein_assumption" else fit$status
    )
  }), fill = TRUE)
}, mc.cores = workers, mc.preschedule = TRUE)
all_rows[[length(all_rows) + 1L]] <- rbindlist(rosmap_rows, fill = TRUE)

# --- save coloc outputs ---

all_results <- rbindlist(all_rows, fill = TRUE)
all_results[, `:=`(analysis_version = analysis_version, config_path = config_path)]
setorder(all_results, candidate_id, source, -external_validation_eligible, -PP.H4, na.last = TRUE)
atomic_fwrite(all_results, all_output)

summary_results <- all_results[, .SD[1L], by = .(candidate_id, source)]
summary_results[, any_external_qtl_coloc := any(external_validation_eligible), by = candidate_id]
atomic_fwrite(summary_results, summary_output)
source_audit_table <- rbindlist(source_audit, fill = TRUE)
source_audit_table[, `:=`(
  external_cache_mode = external_cache_mode,
  analysis_version = analysis_version,
  config_path = config_path
)]
atomic_fwrite(source_audit_table, source_audit_output)

# --- published colocalization comparison ---

mechanism_path <- file.path(out_dir, "02_signal_mechanisms.tsv.gz")
consequence_path <- file.path(out_dir, "03_event_consequences.tsv.gz")
coloc_run <- file.path(out_dir, "coloc_run")
local_coloc_path <- file.path(coloc_run, "all_colocalization.csv")
published_path <- file.path(out_dir, "Traits_diseases_associated_3aQTLs.txt.gz")
published_output <- file.path(out_dir, "published.tsv")
published_match_output <- file.path(out_dir, "published_matches.tsv")
published_summary_output <- file.path(out_dir, "published_summary.json")
formal_pp4_threshold <- 0.75

trait_definitions <- list(
  SCZ = list("Schizophrenia", "\\bschizophrenia\\b"),
  BIP = list("Bipolar disorder", "\\bbipolar\\b"),
  INT = list("Intelligence", "\\bintelligen|\\bcognitive (?:ability|performance)|\\bgeneral cognitive ability|\\bverbal-numerical reasoning"),
  ADHD = list("Attention-deficit/hyperactivity disorder", "\\battention[- ]deficit|\\bADHD\\b"),
  MDD = list("Major depressive disorder", "\\bdepress|\\bmajor depressive\\b"),
  AD = list("Alzheimer disease", "\\balzheimer|\\bdementia\\b"),
  PD = list("Parkinson disease", "\\bparkinson\\b"),
  PTSD = list("Post-traumatic stress disorder", "\\bpost[- ]traumatic|\\bPTSD\\b"),
  ALS = list("Amyotrophic lateral sclerosis", "\\bamyotrophic|\\bALS\\b"),
  AIS = list("Ischemic stroke", "\\bstroke\\b|\\bischemic\\b")
)

mechanisms <- fread(mechanism_path)

local_coloc <- fread(local_coloc_path)
setnames(
  local_coloc,
  c("QTL_ID", "Phenotype", "GWAS_ID", "PP4", "Locus"),
  c("stage", "phenotype_id", "gwas_trait", "local_coloc_pp4", "local_coloc_locus_rsid")
)
local_coloc[, local_coloc_pp4 := as.numeric(local_coloc_pp4)]
setorder(local_coloc, -local_coloc_pp4)
local_coloc <- unique(local_coloc, by = c("stage", "phenotype_id", "gwas_trait"))
local_coloc <- local_coloc[local_coloc_pp4 >= formal_pp4_threshold]

pairs <- merge(
  mechanisms,
  local_coloc[, .(stage, phenotype_id, gwas_trait, local_coloc_pp4, local_coloc_locus_rsid, n_snps)],
  by = c("stage", "phenotype_id"),
  allow.cartesian = TRUE
)

consequences <- fread(consequence_path)
consequence_order <- c(
  mediation_supported = 0L,
  conditional_association_only = 1L,
  orthogonal_support_only = 2L,
  predicted_consequence_only = 3L,
  consequence_unresolved = 4L
)
consequences[, consequence_order := consequence_order[downstream_consequence_evidence_label]]
consequences[is.na(consequence_order), consequence_order := 99L]
consequences[, indirect_fdr := as.numeric(indirect_fdr)]
setorder(consequences, consequence_order, indirect_fdr, na.last = TRUE)
consequences <- unique(consequences, by = "candidate_id")[, .(
  candidate_id, downstream_consequence_evidence_label, outcome_type, outcome_id, indirect_fdr
)]
pairs <- merge(pairs, consequences, by = "candidate_id", all.x = TRUE)

external <- fread(all_output)
external[, supported := vapply(external_validation_eligible, is_true, logical(1))]
external[, `PP.H4` := as.numeric(`PP.H4`)]
external_counts <- external[, .(n_external_qtl_coloc = sum(supported)), by = candidate_id]
external_best <- copy(external)
setorder(external_best, candidate_id, -supported, -`PP.H4`, na.last = TRUE)
external_best <- unique(external_best, by = "candidate_id")
external_best <- external_best[, .(
  candidate_id,
  strongest_external_qtl_source = source,
  strongest_external_qtl_type = qtl_type,
  strongest_external_qtl_pp4 = `PP.H4`,
  strongest_external_event_match_scope = event_match_scope,
  strongest_external_evidence_role = evidence_role,
  strongest_external_population_relation = population_relation
)]
pairs <- merge(pairs, external_counts, by = "candidate_id", all.x = TRUE)
pairs <- merge(pairs, external_best, by = "candidate_id", all.x = TRUE)
pairs[is.na(n_external_qtl_coloc), n_external_qtl_coloc := 0L]

pairs[, `:=`(
  mechanism_class = snp_to_apa_evidence_label,
  functional_support_category = fifelse(is.na(local_annotation_status), "not_available", local_annotation_status),
  chain_completeness = "current_signal_formal_coloc_and_event_matched_external_contract",
  analysis_version = analysis_version,
  interpretation_limit = "Formal coloc supports a shared regional signal but does not by itself establish a causal variant or mechanism."
)]

functional_summaries <- lapply(seq_len(nrow(pairs)), function(index) {
  row <- pairs[index]
  evidence <- character()
  if (is_true(row$direct_pas_evidence)) evidence <- c(evidence, "direct_PAS_evidence")
  if (is_present(row$local_rbp_motif_change) || is_present(row$local_clip_support)) {
    evidence <- c(evidence, "local_RBP_support")
  }
  if (is_true(row$local_structure_change) || is_present(row$local_structure_delta_mfe)) {
    evidence <- c(evidence, "local_structure_support")
  }
  aparent <- suppressWarnings(as.numeric(row$APARENT_predicted_PDUI_delta))
  if (!is.na(aparent) && abs(aparent) > 0) {
    evidence <- c(evidence, sprintf("APARENT_delta=%.4g", aparent))
  }
  if (is_true(row$mechanism_star_eligible)) evidence <- c(evidence, "mechanism_star_eligible")
  if (as.integer(row$n_external_qtl_coloc) > 0L) {
    evidence <- c(evidence, sprintf(
      "external_%s=%s", row$strongest_external_qtl_type, row$strongest_external_qtl_source
    ))
  }
  consequence <- as.character(row$downstream_consequence_evidence_label)
  if (!is.na(consequence) && !consequence %chin% c("", "consequence_unresolved")) {
    evidence <- c(evidence, paste0("consequence=", consequence))
  }
  mechanism <- as.character(row$mechanism_class)
  status <- if (is_true(row$mechanism_star_eligible)) {
    "mechanism_supported"
  } else if (!is.na(mechanism) && !mechanism %chin% c("", "association_mechanism_unresolved")) {
    "partial_functional_context"
  } else if (length(evidence) > 0L) {
    "partial_functional_context"
  } else {
    "shared_locus_without_resolved_mechanism"
  }
  list(
    status = status,
    evidence = if (length(evidence) > 0L) paste(evidence, collapse = "; ") else "no_resolved_functional_evidence"
  )
})
pairs[, functional_explanation_status := vapply(
  functional_summaries,
  function(value) if (value$status == "mechanism_supported") "resolved" else "unresolved",
  character(1)
)]
pairs[, functional_evidence_summary := vapply(functional_summaries, `[[`, character(1), "evidence")]

published <- fread(published_path)
published[, published_gene := sub("^.*@", "", as.character(APA_gene))]
comparison_rows <- vector("list", nrow(pairs))
match_tables <- list()
for (index in seq_len(nrow(pairs))) {
  row <- pairs[index]
  trait_code <- sub("^[^_]+_", "", as.character(row$gwas_trait))
  definition <- trait_definitions[[trait_code]]
  if (is.null(definition)) definition <- list(trait_code, trait_code)
  gene_rows <- copy(published[published_gene == row$gene_name])
  gene_rows[, trait_match := grepl(definition[[2]], as.character(DISEASE_TRAIT), ignore.case = TRUE, perl = TRUE)]
  gene_rows[, exact_locus_match :=
    as.character(`3aQTL`) == row$local_coloc_locus_rsid |
    as.character(GWAS_SNP) == row$local_coloc_locus_rsid]
  trait_rows <- gene_rows[trait_match == TRUE]
  exact_rows <- trait_rows[exact_locus_match == TRUE]
  novelty_class <- if (nrow(exact_rows) > 0L) {
    "exact_locus_gene_trait_overlap"
  } else if (nrow(trait_rows) > 0L) {
    "gene_trait_overlap_different_locus"
  } else if (nrow(gene_rows) > 0L) {
    "gene_only_overlap_different_trait"
  } else {
    "not_reported_in_published_table"
  }
  comparison_rows[[index]] <- data.table(
    candidate_id = row$candidate_id,
    stage = row$stage,
    phenotype_id = row$phenotype_id,
    gwas_trait = row$gwas_trait,
    trait_name = definition[[1]],
    published_gene_present = nrow(gene_rows) > 0L,
    published_trait_match = nrow(trait_rows) > 0L,
    exact_locus_match = nrow(exact_rows) > 0L,
    n_published_gene_rows = nrow(gene_rows),
    n_published_gene_traits = uniqueN(na.omit(gene_rows$DISEASE_TRAIT)),
    n_published_gene_trait_rows = nrow(trait_rows),
    published_matching_traits = collapse_unique(trait_rows$DISEASE_TRAIT),
    published_matching_tissues = collapse_unique(trait_rows$Tissue),
    published_matching_rsids = collapse_unique(trait_rows[["3aQTL"]]),
    published_matching_pubmed_ids = collapse_unique(trait_rows$PUBMEDID),
    novelty_class = novelty_class
  )
  if (nrow(gene_rows) > 0L) {
    gene_rows[, `:=`(
      candidate_id = row$candidate_id,
      local_gwas_trait = row$gwas_trait,
      local_coloc_locus_rsid = row$local_coloc_locus_rsid
    )]
    setcolorder(gene_rows, c(
      "candidate_id", "local_gwas_trait", "local_coloc_locus_rsid",
      setdiff(names(gene_rows), c("candidate_id", "local_gwas_trait", "local_coloc_locus_rsid"))
    ))
    match_tables[[length(match_tables) + 1L]] <- gene_rows
  }
}

comparison <- rbindlist(comparison_rows, fill = TRUE)
comparison_keys <- c("candidate_id", "stage", "phenotype_id", "gwas_trait")
if (anyDuplicated(comparison, by = comparison_keys)) stop("Published comparison keys are not unique")
pairs <- merge(pairs, comparison, by = comparison_keys, all.x = TRUE)

output_columns <- c(
  "candidate_id", "gene_name", "stage", "phenotype_id", "variant_id", "gwas_trait", "trait_name",
  "local_coloc_pp4", "local_coloc_locus_rsid", "n_snps", "novelty_class",
  "published_gene_present", "published_trait_match", "exact_locus_match", "n_published_gene_rows",
  "n_published_gene_traits", "n_published_gene_trait_rows", "published_matching_traits",
  "published_matching_tissues", "published_matching_rsids", "published_matching_pubmed_ids",
  "mechanism_class", "snp_to_apa_evidence_label", "mechanism_star_eligible", "functional_support_category",
  "direct_pas_evidence", "APARENT_predicted_PDUI_delta", "local_rbp_motif_change", "local_clip_support",
  "local_structure_change", "local_structure_delta_mfe", "downstream_consequence_evidence_label",
  "outcome_type", "outcome_id", "indirect_fdr", "n_external_qtl_coloc", "strongest_external_qtl_source",
  "strongest_external_qtl_type", "strongest_external_qtl_pp4",
  "strongest_external_event_match_scope", "strongest_external_evidence_role",
  "strongest_external_population_relation", "chain_completeness", "analysis_version",
  "functional_explanation_status", "functional_evidence_summary", "interpretation_limit"
)

setorderv(pairs, c("novelty_class", "local_coloc_pp4"), c(1L, -1L))
fwrite(pairs[, ..output_columns], published_output, sep = "\t")

published_matches <- rbindlist(match_tables, fill = TRUE)
if (ncol(published_matches) == 0L) {
  published_matches <- data.table(
    candidate_id = character(),
    local_gwas_trait = character(),
    local_coloc_locus_rsid = character()
  )
}
fwrite(published_matches, published_match_output, sep = "\t")

novelty_counts <- as.list(table(pairs$novelty_class))
published_summary <- list(
  analysis_version = analysis_version,
  n_current_independent_signal_formal_coloc_pairs = nrow(pairs),
  formal_coloc_pp4_threshold = formal_pp4_threshold,
  novelty_class_counts = novelty_counts,
  n_pairs_with_eligible_external_qtl_support = sum(pairs$n_external_qtl_coloc > 0L),
  n_pairs_with_mechanism_star_eligibility = sum(vapply(pairs$mechanism_star_eligible, is_true, logical(1))),
  published_reference_scope = "single-SNP or LD-based 3'aQTL-GWAS association evidence; not formal coloc posterior",
  novelty_definition = "Novelty is relative only to the supplied published association table.",
  current_input_contract = c(
    "02_signal_mechanisms.tsv.gz",
    "03_event_consequences.tsv.gz",
    "coloc_run/all_colocalization.csv",
    "external_all.tsv"
  )
)
jsonlite::write_json(published_summary, published_summary_output, pretty = TRUE, auto_unbox = TRUE)

# --- NDUFA6-SCZ locus config ---

target_candidate_id <- "postnatal::ENST00000498737.8|NDUFA6|chr22:42085526-42086314|-::chr22:42070901:C:T"
candidate_variant <- "chr22:42070901:C:T"
local_phenotype <- "ENST00000498737.8|NDUFA6|chr22:42085526-42086314|-"
external_phenotype <- "ENSG00000184983.11_chr22:42085526-42086314"
gwas_id <- "EUR_SCZ"
gwas_locus <- "rs2743467"
gwas_n <- 130644
gwas_case_fraction <- 0.4086
pp4_threshold <- 0.75
credible_set_probability <- 0.95
analysis_version <- "framework_2026-07-25_v5"
p1 <- 1e-4
p2 <- 1e-4
p12 <- 1e-5

# --- load local APA-SCZ locus ---

rds_path <- file.path(coloc_run, "rds", "EUR_SCZ_postnatal_chr22_42548880_rs2743467_NDUFA6_coloc.rds")
locus <- readRDS(rds_path)
merged <- as.data.table(locus$merged_data)

local_scz_local <- merged[, .(
  variant_id,
  beta = BETA.qtl,
  se = SE.qtl,
  af
)]
local_scz_gwas <- merged[, .(
  variant_id,
  beta = BETA.gwas,
  se = SE.gwas,
  af = EAF
)]

# --- load full local and GTEx BA9 APA signals ---

local_full <- read_phenotype_rows(file.path(out_dir, "local_apa_postnatal.tsv.gz"), local_phenotype)
local_full[, variant_id := canonical_variant(variant_id)]
local_full <- local_full[, .(
  variant_id,
  beta = slope,
  se = slope_se,
  af,
  p = pval_nominal
)]

ba9_full <- read_phenotype_rows(file.path(out_dir, "GTEx_v11_Brain_BA9_apaQTL.tsv.gz"), external_phenotype)
ba9_full[, variant_id := canonical_variant(variant_id)]
ba9_full <- ba9_full[, .(
  variant_id,
  beta = slope,
  se = slope_se,
  af,
  p = pval_nominal
)]

postnatal_n <- length(readLines(file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/apa.sample.Postnatal.final")))
ba9_covariates <- paste0(Sys.getenv("APA_HPC_HOME"), "/share_group_folder/.ref/GTEx_v11/GTEx_Analysis_v11_eQTL/Brain_Frontal_Cortex_BA9.v11.covariates.txt")
ba9_n <- length(strsplit(readLines(ba9_covariates, n = 1), "\t", fixed = TRUE)[[1]]) - 1L

# --- pairwise coloc triangulation ---

local_scz <- run_coloc_pair(
  "local_postnatal_APA_vs_EUR_SCZ",
  local_scz_local, local_scz_gwas,
  "quant", "cc", postnatal_n, gwas_n,
  second_s = gwas_case_fraction
)
local_ba9 <- run_coloc_pair(
  "local_postnatal_APA_vs_GTEx_BA9_APA",
  local_full, ba9_full,
  "quant", "quant", postnatal_n, ba9_n
)

scz_ba9_gwas <- local_scz_gwas[variant_id %chin% ba9_full$variant_id]
scz_ba9 <- run_coloc_pair(
  "EUR_SCZ_vs_GTEx_BA9_APA",
  scz_ba9_gwas, ba9_full,
  "cc", "quant", gwas_n, ba9_n,
  first_s = gwas_case_fraction
)

pair_summary <- rbindlist(list(
  local_scz$summary,
  local_ba9$summary,
  scz_ba9$summary
))

existing_external <- fread(all_output)[
  candidate_id == target_candidate_id &
    source == "GTEx_v11_Brain_BA9_apaQTL" &
    external_phenotype_id == external_phenotype &
    event_match_scope == "same_gene_apa_overlap" &
    external_validation_eligible == TRUE
]
existing_external_pp4 <- if (nrow(existing_external) == 0L) {
  NA_real_
} else {
  max(existing_external$PP.H4, na.rm = TRUE)
}
pair_summary[, previous_reported_PP4 := NA_real_]
pair_summary[
  pair == "local_postnatal_APA_vs_EUR_SCZ",
  previous_reported_PP4 := locus$pp4
]
pair_summary[
  pair == "local_postnatal_APA_vs_GTEx_BA9_APA",
  previous_reported_PP4 := existing_external_pp4
]
pair_summary[, absolute_recompute_delta := abs(PP.H4 - previous_reported_PP4)]
stopifnot(all(pair_summary[!is.na(previous_reported_PP4), absolute_recompute_delta] < 0.001))

# --- prior sensitivity ---

prior_grid <- data.table(
  setting = c("default", "higher_first_trait", "higher_second_trait", "higher_both", "conservative"),
  p1 = c(1e-4, 1e-3, 1e-4, 1e-3, 1e-5),
  p2 = c(1e-4, 1e-4, 1e-3, 1e-3, 1e-5),
  p12 = c(1e-5, 1e-5, 1e-5, 1e-5, 1e-6)
)

prior_sensitivity <- rbindlist(lapply(seq_len(nrow(prior_grid)), function(index) {
  prior <- prior_grid[index]
  fits <- suppressWarnings(list(
    run_coloc_pair(
      "local_postnatal_APA_vs_EUR_SCZ",
      local_scz_local, local_scz_gwas,
      "quant", "cc", postnatal_n, gwas_n,
      second_s = gwas_case_fraction,
      p1_value = prior$p1, p2_value = prior$p2, p12_value = prior$p12
    )$summary,
    run_coloc_pair(
      "local_postnatal_APA_vs_GTEx_BA9_APA",
      local_full, ba9_full,
      "quant", "quant", postnatal_n, ba9_n,
      p1_value = prior$p1, p2_value = prior$p2, p12_value = prior$p12
    )$summary,
    run_coloc_pair(
      "EUR_SCZ_vs_GTEx_BA9_APA",
      scz_ba9_gwas, ba9_full,
      "cc", "quant", gwas_n, ba9_n,
      first_s = gwas_case_fraction,
      p1_value = prior$p1, p2_value = prior$p2, p12_value = prior$p12
    )$summary
  ))
  result <- rbindlist(fits)
  result[, `:=`(
    setting = prior$setting,
    p1 = prior$p1,
    p2 = prior$p2,
    p12 = prior$p12
  )]
  result
}))
prior_sensitivity <- prior_sensitivity[, .(
  setting, pair, p1, p2, p12, n_shared_snps, PP.H4, coloc_supported
)]

# --- shared-variant posterior evidence ---

posterior <- rbindlist(list(
  local_scz$posterior,
  local_ba9$posterior,
  scz_ba9$posterior
), fill = TRUE)
posterior <- posterior[
  in_conditional_95pct_H4_set == TRUE | snp == candidate_variant,
  .(pair, variant_id = snp, SNP.PP.H4, posterior_rank,
    cumulative_conditional_SNP_PP_H4, in_conditional_95pct_H4_set)
]

local_effects <- copy(local_full)
setnames(
  local_effects,
  c("beta", "se", "af", "p"),
  c("beta_local_APA", "se_local_APA", "af_local_APA", "p_local_APA")
)
gwas_effects <- copy(local_scz_gwas)
setnames(
  gwas_effects,
  c("beta", "se", "af"),
  c("beta_EUR_SCZ", "se_EUR_SCZ", "af_EUR_SCZ")
)
ba9_effects <- copy(ba9_full)
setnames(
  ba9_effects,
  c("beta", "se", "af", "p"),
  c("beta_GTEx_BA9_APA", "se_GTEx_BA9_APA", "af_GTEx_BA9_APA", "p_GTEx_BA9_APA")
)
effect_evidence <- Reduce(
  function(first, second) merge(first, second, by = "variant_id"),
  list(local_effects, gwas_effects, ba9_effects)
)
rsid_map <- unique(merged[, .(variant_id, current_b38_rsid = rsid.gwas)])
effect_evidence <- merge(effect_evidence, rsid_map, by = "variant_id", all.x = TRUE)
posterior <- merge(posterior, effect_evidence, by = "variant_id", all.x = TRUE)
setorder(posterior, pair, posterior_rank)

# --- stage specificity and locus competitors ---

all_coloc <- fread(local_coloc_path)
locus_competitors <- all_coloc[GWAS_ID == gwas_id & Locus == gwas_locus]
locus_competitors[, gene_name := tstrsplit(Phenotype, "|", fixed = TRUE)[[2]]]
setorder(locus_competitors, -PP4)
locus_competitors <- locus_competitors[, .(
  GWAS_ID, QTL_ID, Locus, gene_name, Phenotype, PP4, n_snps,
  passes_pp4_threshold = PP4 >= pp4_threshold
)]

prenatal_ndufa6_values <- locus_competitors[
  QTL_ID == "prenatal" & gene_name == "NDUFA6", PP4
]
prenatal_ndufa6_pp4 <- if (length(prenatal_ndufa6_values) == 0L) {
  NA_real_
} else {
  max(prenatal_ndufa6_values)
}
postnatal_competitors <- locus_competitors[
  QTL_ID == "postnatal" & passes_pp4_threshold == TRUE
]

# --- integrated interpretation ---

candidate <- fread(file.path(out_dir, "02_signal_mechanisms.tsv.gz"))[
  candidate_id == target_candidate_id
]
stopifnot(nrow(candidate) == 1L)
published <- fread(file.path(out_dir, "published.tsv"))[
  candidate_id == target_candidate_id &
    gwas_trait == gwas_id &
    stage == "postnatal" &
    phenotype_id == local_phenotype
]
stopifnot(nrow(published) == 1L)
candidate_rsid <- merged[variant_id == candidate_variant, unique(rsid.gwas)]
candidate_effects <- effect_evidence[variant_id == candidate_variant]

pair_ranges <- prior_sensitivity[, .(
  prior_sensitivity_PP4_min = min(PP.H4),
  prior_sensitivity_PP4_max = max(PP.H4),
  prior_settings_supported = sum(coloc_supported)
), by = pair]

summary <- data.table(
  candidate_id = target_candidate_id,
  gene_name = "NDUFA6",
  stage = "postnatal",
  candidate_variant_b38 = candidate_variant,
  candidate_current_b38_rsid = candidate_rsid,
  historical_gwas_hg19_locus_label = gwas_locus,
  rsid_roles = paste0(
    gwas_locus, "=historical_hg19_GWAS_clump_label;",
    candidate_rsid, "=current_GRCh38_candidate_variant"
  ),
  local_finemap_PIP = candidate$representative_variant_pip,
  local_APA_vs_SCZ_PP4 = pair_summary[pair == "local_postnatal_APA_vs_EUR_SCZ", PP.H4],
  local_APA_vs_GTEx_BA9_APA_PP4 = pair_summary[pair == "local_postnatal_APA_vs_GTEx_BA9_APA", PP.H4],
  SCZ_vs_GTEx_BA9_APA_PP4 = pair_summary[pair == "EUR_SCZ_vs_GTEx_BA9_APA", PP.H4],
  all_three_pairwise_coloc_supported = all(pair_summary$coloc_supported),
  local_APA_beta_ALT_T = candidate_effects$beta_local_APA,
  GTEx_BA9_APA_beta_ALT_T = candidate_effects$beta_GTEx_BA9_APA,
  EUR_SCZ_beta_ALT_T = candidate_effects$beta_EUR_SCZ,
  candidate_effect_directions_concordant = length(unique(sign(c(
    candidate_effects$beta_local_APA,
    candidate_effects$beta_GTEx_BA9_APA,
    candidate_effects$beta_EUR_SCZ
  )))) == 1L,
  prenatal_NDUFA6_APA_vs_SCZ_PP4 = prenatal_ndufa6_pp4,
  postnatal_specific_coloc_support = pair_summary[
    pair == "local_postnatal_APA_vs_EUR_SCZ", coloc_supported
  ] & (is.na(prenatal_ndufa6_pp4) | prenatal_ndufa6_pp4 < pp4_threshold),
  n_postnatal_APA_genes_PP4_ge_0_75_at_locus = nrow(postnatal_competitors),
  postnatal_competing_APA_genes = paste(postnatal_competitors$gene_name, collapse = ";"),
  mechanism_class = published$mechanism_class,
  functional_support_category = published$functional_support_category,
  mechanism_resolved = published$functional_explanation_status == "resolved",
  interpretation = if (
    pair_summary[pair == "local_postnatal_APA_vs_EUR_SCZ", coloc_supported] &&
    pair_summary[pair == "local_postnatal_APA_vs_GTEx_BA9_APA", coloc_supported] &&
    !pair_summary[pair == "EUR_SCZ_vs_GTEx_BA9_APA", coloc_supported]
  ) {
    "postnatal_local_APA_SCZ_coloc_and_BA9_APA_replication_without_direct_SCZ_BA9_APA_coloc"
  } else {
    "partial_pairwise_coloc_support_gene_and_molecular_mechanism_not_resolved"
  }
)
summary[, `:=`(
  local_APA_vs_SCZ_prior_PP4_min = pair_ranges[
    pair == "local_postnatal_APA_vs_EUR_SCZ", prior_sensitivity_PP4_min
  ],
  local_APA_vs_SCZ_prior_PP4_max = pair_ranges[
    pair == "local_postnatal_APA_vs_EUR_SCZ", prior_sensitivity_PP4_max
  ],
  local_APA_vs_SCZ_prior_settings_supported = pair_ranges[
    pair == "local_postnatal_APA_vs_EUR_SCZ", prior_settings_supported
  ],
  local_APA_vs_GTEx_BA9_APA_prior_PP4_min = pair_ranges[
    pair == "local_postnatal_APA_vs_GTEx_BA9_APA", prior_sensitivity_PP4_min
  ],
  local_APA_vs_GTEx_BA9_APA_prior_PP4_max = pair_ranges[
    pair == "local_postnatal_APA_vs_GTEx_BA9_APA", prior_sensitivity_PP4_max
  ],
  local_APA_vs_GTEx_BA9_APA_prior_settings_supported = pair_ranges[
    pair == "local_postnatal_APA_vs_GTEx_BA9_APA", prior_settings_supported
  ],
  SCZ_vs_GTEx_BA9_APA_prior_PP4_min = pair_ranges[
    pair == "EUR_SCZ_vs_GTEx_BA9_APA", prior_sensitivity_PP4_min
  ],
  SCZ_vs_GTEx_BA9_APA_prior_PP4_max = pair_ranges[
    pair == "EUR_SCZ_vs_GTEx_BA9_APA", prior_sensitivity_PP4_max
  ],
  SCZ_vs_GTEx_BA9_APA_prior_settings_supported = pair_ranges[
    pair == "EUR_SCZ_vs_GTEx_BA9_APA", prior_settings_supported
  ]
)]
summary[, default_pairwise_support_prior_robust :=
  local_APA_vs_SCZ_prior_settings_supported == nrow(prior_grid) &
  local_APA_vs_GTEx_BA9_APA_prior_settings_supported == nrow(prior_grid) &
  SCZ_vs_GTEx_BA9_APA_prior_settings_supported == nrow(prior_grid)]
summary[, interpretation := paste0(interpretation, "_prior_sensitive")]
summary[, analysis_version := analysis_version]
pair_summary[, analysis_version := analysis_version]
prior_sensitivity[, analysis_version := analysis_version]
posterior[, analysis_version := analysis_version]
locus_competitors[, analysis_version := analysis_version]

# --- save outputs ---

fwrite(summary, file.path(out_dir, "case_summary.tsv"), sep = "\t")
fwrite(pair_summary, file.path(out_dir, "case_pairs.tsv"), sep = "\t")
fwrite(prior_sensitivity, file.path(out_dir, "case_sensitivity.tsv"), sep = "\t")
fwrite(posterior, file.path(out_dir, "case_posteriors.tsv"), sep = "\t")
fwrite(locus_competitors, file.path(out_dir, "case_competitors.tsv"), sep = "\t")

# --- broader locus config ---

prior_grid <- data.table(
  setting = c("default", "higher_GWAS", "higher_APA", "higher_both", "conservative"),
  p1 = c(1e-4, 1e-3, 1e-4, 1e-3, 1e-5),
  p2 = c(1e-4, 1e-4, 1e-3, 1e-3, 1e-5),
  p12 = c(1e-5, 1e-5, 1e-5, 1e-5, 1e-6)
)
qtl_sample_sizes <- c(prenatal = 189, postnatal = 508)
maf_default <- 0.1
maf_na_replacement <- 0.05
maf_epsilon <- 1e-6

# --- load prioritized novel pairs ---

published <- fread(file.path(out_dir, "published.tsv"))
targets <- published[novelty_class == "not_reported_in_published_table"]
candidate_table <- fread(file.path(out_dir, "02_signal_mechanisms.tsv.gz"))[, .(
  candidate_id,
  local_finemap_PIP = representative_variant_pip
)]
targets <- merge(targets, candidate_table, by = "candidate_id", all.x = TRUE)
all_coloc <- fread(local_coloc_path)
gwas_config <- read_yaml(file.path(Sys.getenv("APA_ROOT"), "release/07_coloc/easycoloc/config/gwas.yml"))$datasets
rds_dir <- file.path(coloc_run, "rds")
rds_files <- list.files(rds_dir, pattern = "_coloc\\.rds$", full.names = TRUE)

# --- local APA-GWAS robustness and variant evidence ---

summary_rows <- list()
prior_rows <- list()
posterior_rows <- list()
competitor_rows <- list()

for (target_index in seq_len(nrow(targets))) {
  target <- targets[target_index]
  pattern <- paste0(
    "^", target$gwas_trait, "_", target$stage, "_.*_",
    target$local_coloc_locus_rsid, "_", target$gene_name, "_coloc\\.rds$"
  )
  rds_path <- rds_files[grepl(pattern, basename(rds_files))]
  if (length(rds_path) != 1L) {
    summary_rows[[target_index]] <- data.table(
      candidate_id = target$candidate_id,
      gene_name = target$gene_name,
      stage = target$stage,
      gwas_trait = target$gwas_trait,
      phenotype_id = target$phenotype_id,
      candidate_variant_b38 = target$variant_id,
      historical_gwas_locus_label = target$local_coloc_locus_rsid,
      local_finemap_PIP = target$local_finemap_PIP,
      previous_reported_PP4 = target$local_coloc_pp4,
      recomputed_default_PP4 = NA_real_,
      absolute_recompute_delta = NA_real_,
      recompute_classification_matches = NA,
      prior_PP4_min = NA_real_,
      prior_PP4_max = NA_real_,
      prior_settings_supported = NA_integer_,
      prior_settings_tested = nrow(prior_grid),
      prior_robust_local_coloc = FALSE,
      candidate_variant_in_coloc_dataset = NA,
      candidate_conditional_SNP_PP_H4 = NA_real_,
      candidate_in_conditional_95pct_H4_set = NA,
      APA_beta_ALT = NA_real_,
      GWAS_beta_ALT = NA_real_,
      APA_and_GWAS_effect_direction_concordant = NA,
      opposite_stage_PP4 = NA_real_,
      stage_specific_default_support = NA,
      n_same_stage_genes_PP4_ge_0_75_at_locus = NA_integer_,
      same_stage_competing_genes = NA_character_,
      gene_specific_at_locus = NA,
      n_external_qtl_coloc = target$n_external_qtl_coloc,
      strongest_external_qtl_source = target$strongest_external_qtl_source,
      strongest_external_qtl_type = target$strongest_external_qtl_type,
      strongest_external_qtl_PP4 = target$strongest_external_qtl_pp4,
      mechanism_class = target$mechanism_class,
      mechanism_resolved = target$functional_explanation_status == "resolved",
      evidence_class = if (length(rds_path) == 0L) {
        "missing_locus_data"
      } else {
        "ambiguous_locus_data"
      },
      rds_status = if (length(rds_path) == 0L) "missing" else "ambiguous"
    )
    next
  }

  locus <- readRDS(rds_path)
  merged <- as.data.table(locus$merged_data)
  merged[, QTL_ID := target$stage]

  gwas <- Filter(
    function(dataset) identical(dataset$id, target$gwas_trait),
    gwas_config
  )[[1]]
  gwas_prop <- if (gwas$type == "cc") as.numeric(gwas$prop) else NULL

  sensitivity <- rbindlist(lapply(seq_len(nrow(prior_grid)), function(prior_index) {
    prior <- prior_grid[prior_index]
    fit <- suppressWarnings(run_local_gwas_coloc(
      merged, gwas$type, as.numeric(gwas$sample_size_n), gwas_prop,
      prior$p1, prior$p2, prior$p12
    ))
    data.table(
      candidate_id = target$candidate_id,
      gene_name = target$gene_name,
      stage = target$stage,
      gwas_trait = target$gwas_trait,
      setting = prior$setting,
      p1 = prior$p1,
      p2 = prior$p2,
      p12 = prior$p12,
      n_shared_snps = nrow(merged),
      PP.H4 = unname(fit$summary[["PP.H4.abf"]]),
      coloc_supported = unname(fit$summary[["PP.H4.abf"]]) >= pp4_threshold
    )
  }))
  prior_rows[[target_index]] <- sensitivity

  default_fit <- suppressWarnings(run_local_gwas_coloc(
    merged, gwas$type, as.numeric(gwas$sample_size_n), gwas_prop,
    prior_grid[setting == "default", p1],
    prior_grid[setting == "default", p2],
    prior_grid[setting == "default", p12]
  ))
  default_posterior <- as.data.table(default_fit$results)[order(-SNP.PP.H4)]
  default_posterior[, `:=`(
    posterior_rank = seq_len(.N),
    cumulative_conditional_SNP_PP_H4 = cumsum(SNP.PP.H4)
  )]
  cutoff_rank <- which(
    default_posterior$cumulative_conditional_SNP_PP_H4 >= credible_set_probability
  )[1]
  default_posterior[, in_conditional_95pct_H4_set := posterior_rank <= cutoff_rank]
  candidate_posterior <- default_posterior[snp == target$variant_id]
  candidate_variant_in_coloc_dataset <- nrow(candidate_posterior) == 1L

  posterior_output <- default_posterior[
    in_conditional_95pct_H4_set == TRUE | snp == target$variant_id,
    .(
      candidate_id = target$candidate_id,
      gene_name = target$gene_name,
      stage = target$stage,
      gwas_trait = target$gwas_trait,
      variant_id = snp,
      SNP.PP.H4,
      posterior_rank,
      cumulative_conditional_SNP_PP_H4,
      in_conditional_95pct_H4_set
    )
  ]
  variant_effects <- merged[, .(
    variant_id,
    current_b38_rsid = rsid.gwas,
    beta_APA = BETA.qtl,
    se_APA = SE.qtl,
    beta_GWAS = BETA.gwas,
    se_GWAS = SE.gwas
  )]
  posterior_rows[[target_index]] <- merge(
    posterior_output, variant_effects,
    by = "variant_id", all.x = TRUE
  )

  competitors <- all_coloc[
    GWAS_ID == target$gwas_trait &
    QTL_ID == target$stage &
    Locus == target$local_coloc_locus_rsid
  ]
  competitors[, gene_name := tstrsplit(Phenotype, "|", fixed = TRUE)[[2]]]
  setorder(competitors, -PP4)
  competitor_rows[[target_index]] <- competitors[, .(
    candidate_id = target$candidate_id,
    GWAS_ID, QTL_ID, Locus, gene_name, Phenotype, PP4, n_snps,
    passes_pp4_threshold = PP4 >= pp4_threshold
  )]

  opposite_stage <- if (target$stage == "prenatal") "postnatal" else "prenatal"
  opposite_stage_rows <- all_coloc[
    GWAS_ID == target$gwas_trait &
    QTL_ID == opposite_stage &
    Locus == target$local_coloc_locus_rsid &
    Phenotype == target$phenotype_id
  ]
  opposite_stage_pp4 <- if (nrow(opposite_stage_rows) == 0) {
    NA_real_
  } else {
    max(opposite_stage_rows$PP4)
  }

  candidate_effect <- merged[variant_id == target$variant_id]
  if (nrow(candidate_effect) == 0) {
    candidate_effect <- data.table(
      rsid.gwas = NA_character_,
      BETA.qtl = NA_real_,
      BETA.gwas = NA_real_
    )
  }
  supported_competitors <- competitors[PP4 >= pp4_threshold]
  prior_supported <- sum(sensitivity$coloc_supported)

  summary_rows[[target_index]] <- data.table(
    candidate_id = target$candidate_id,
    gene_name = target$gene_name,
    stage = target$stage,
    gwas_trait = target$gwas_trait,
    phenotype_id = target$phenotype_id,
    candidate_variant_b38 = target$variant_id,
    candidate_current_b38_rsid = candidate_effect$rsid.gwas,
    historical_gwas_locus_label = target$local_coloc_locus_rsid,
    local_finemap_PIP = target$local_finemap_PIP,
    previous_reported_PP4 = target$local_coloc_pp4,
    recomputed_default_PP4 = sensitivity[setting == "default", PP.H4],
    absolute_recompute_delta = abs(
      sensitivity[setting == "default", PP.H4] - target$local_coloc_pp4
    ),
    recompute_classification_matches =
      (sensitivity[setting == "default", PP.H4] >= pp4_threshold) ==
      (target$local_coloc_pp4 >= pp4_threshold),
    prior_PP4_min = min(sensitivity$PP.H4),
    prior_PP4_max = max(sensitivity$PP.H4),
    prior_settings_supported = prior_supported,
    prior_settings_tested = nrow(prior_grid),
    prior_robust_local_coloc = prior_supported == nrow(prior_grid),
    candidate_variant_in_coloc_dataset = candidate_variant_in_coloc_dataset,
    candidate_conditional_SNP_PP_H4 = if (candidate_variant_in_coloc_dataset) {
      candidate_posterior$SNP.PP.H4
    } else {
      NA_real_
    },
    candidate_in_conditional_95pct_H4_set = if (candidate_variant_in_coloc_dataset) {
      candidate_posterior$in_conditional_95pct_H4_set
    } else {
      NA
    },
    APA_beta_ALT = candidate_effect$BETA.qtl,
    GWAS_beta_ALT = candidate_effect$BETA.gwas,
    APA_and_GWAS_effect_direction_concordant = sign(candidate_effect$BETA.qtl) == sign(candidate_effect$BETA.gwas),
    opposite_stage_PP4 = opposite_stage_pp4,
    stage_specific_default_support = target$local_coloc_pp4 >= pp4_threshold &
      (is.na(opposite_stage_pp4) | opposite_stage_pp4 < pp4_threshold),
    n_same_stage_genes_PP4_ge_0_75_at_locus = nrow(supported_competitors),
    same_stage_competing_genes = paste(supported_competitors$gene_name, collapse = ";"),
    gene_specific_at_locus = nrow(supported_competitors) == 1L,
    n_external_qtl_coloc = target$n_external_qtl_coloc,
    strongest_external_qtl_source = target$strongest_external_qtl_source,
    strongest_external_qtl_type = target$strongest_external_qtl_type,
    strongest_external_qtl_PP4 = target$strongest_external_qtl_pp4,
    mechanism_class = target$mechanism_class,
    mechanism_resolved = target$functional_explanation_status == "resolved",
    evidence_class = if (!candidate_variant_in_coloc_dataset) {
      "regional_coloc_candidate_variant_absent_from_shared_dataset"
    } else if (target$n_external_qtl_coloc > 0) {
      "external_QTL_replication_local_coloc_prior_sensitive"
    } else if (!candidate_posterior$in_conditional_95pct_H4_set) {
      "regional_coloc_candidate_variant_outside_conditional_95pct_H4_set"
    } else if (prior_supported == nrow(prior_grid)) {
      "prior_robust_local_coloc_without_external_QTL_replication"
    } else {
      "prior_sensitive_local_coloc_without_external_QTL_replication"
    },
    rds_status = "available"
  )
}

# --- prioritize follow-up candidates ---

summary <- rbindlist(summary_rows, fill = TRUE)
classification_mismatches <- summary[
  !is.na(recompute_classification_matches) & !recompute_classification_matches,
  .(
    candidate_id, gene_name, stage, gwas_trait,
    previous_reported_PP4, recomputed_default_PP4, absolute_recompute_delta
  )
]
if (nrow(classification_mismatches) > 0L) {
  stop(sprintf("Current-v5 coloc classification mismatch for %d available targets", nrow(classification_mismatches)))
}
setorder(
  summary,
  -n_external_qtl_coloc,
  -candidate_in_conditional_95pct_H4_set,
  -prior_robust_local_coloc,
  -candidate_variant_in_coloc_dataset,
  -gene_specific_at_locus,
  -prior_settings_supported,
  -recomputed_default_PP4,
  na.last = TRUE
)
summary[, rank := seq_len(.N)]
setcolorder(summary, c("rank", setdiff(names(summary), "rank")))

prior_sensitivity <- rbindlist(prior_rows)
posterior <- rbindlist(posterior_rows, fill = TRUE)
competitors <- rbindlist(competitor_rows, fill = TRUE)

# --- save outputs ---

summary[, analysis_version := analysis_version]
prior_sensitivity[, analysis_version := analysis_version]
posterior[, analysis_version := analysis_version]
competitors[, analysis_version := analysis_version]

fwrite(summary, file.path(out_dir, "novel_summary.tsv"), sep = "\t")
fwrite(prior_sensitivity, file.path(out_dir, "novel_sensitivity.tsv"), sep = "\t")
fwrite(posterior, file.path(out_dir, "novel_posteriors.tsv"), sep = "\t")
fwrite(competitors, file.path(out_dir, "novel_competitors.tsv"), sep = "\t")

print(summary[, .(
  rank, gene_name, stage, gwas_trait,
  recomputed_default_PP4, prior_settings_supported,
  n_external_qtl_coloc, gene_specific_at_locus,
  mechanism_resolved, evidence_class
)])
