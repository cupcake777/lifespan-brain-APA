#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(glue)
})

source("src/utils_config.R")

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list(dataset = NULL, force = FALSE, threads = NULL)
  idx <- 1L
  while (idx <= length(args)) {
    arg <- args[[idx]]
    if (arg %in% c("--global", "--gwas", "--qtl", "--dataset", "--threads")) {
      if (idx == length(args)) {
        stop("Missing value for ", arg, call. = FALSE)
      }
      parsed[[substring(arg, 3L)]] <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--force") {
      parsed$force <- TRUE
      idx <- idx + 1L
    } else if (arg %in% c("--help", "-h")) {
      cat("Usage: Rscript tools/harmonize_gwas.R [--global PATH] [--gwas PATH] [--qtl PATH] [--dataset ID] [--threads N] [--force]\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  parsed
}

opts <- parse_args(args)
if (!is.null(opts$global)) Sys.setenv(EASYCOLOC_GLOBAL_CONFIG = opts$global)
if (!is.null(opts$gwas)) Sys.setenv(EASYCOLOC_GWAS_CONFIG = opts$gwas)
if (!is.null(opts$qtl)) Sys.setenv(EASYCOLOC_QTL_CONFIG = opts$qtl)

cfg_bundle <- easycoloc_read_configs()
cfg_global <- cfg_bundle$global
cfg_gwas <- cfg_bundle$gwas
cfg_qtl <- cfg_bundle$qtl

utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

qtl_build <- if (is.null(cfg_qtl$qtl_info$build)) "hg38" else cfg_qtl$qtl_info$build
n_threads <- if (!is.null(opts$threads)) as.integer(opts$threads) else if (!is.null(cfg_global$n_cores)) as.integer(cfg_global$n_cores) else 1L
n_threads <- max(1L, n_threads)

datasets <- cfg_gwas$datasets
if (!is.null(opts$dataset)) {
  datasets <- Filter(function(ds) identical(as.character(ds$id), as.character(opts$dataset)), datasets)
  if (length(datasets) == 0) {
    stop(glue("No GWAS dataset found for --dataset {opts$dataset}"), call. = FALSE)
  }
}

if (length(datasets) == 0) {
  stop("No GWAS datasets configured", call. = FALSE)
}

if (is.null(cfg_global$harmonize_dir) || !nzchar(cfg_global$harmonize_dir)) {
  stop("harmonize_dir is not configured", call. = FALSE)
}
dir.create(cfg_global$harmonize_dir, recursive = TRUE, showWarnings = FALSE)

message(glue("[HARMONIZE] Target build: {qtl_build}"))
message(glue("[HARMONIZE] Output dir: {cfg_global$harmonize_dir}"))
message(glue("[HARMONIZE] Threads per dataset: {n_threads}"))

rows <- list()
for (ds in datasets) {
  source_build <- if (is.null(ds$build)) "19" else gsub("hg", "", ds$build)
  target_build <- gsub("hg", "", qtl_build)
  harmony_settings <- if (!is.null(cfg_global$harmonization_settings)) cfg_global$harmonization_settings else list()
  compress_output <- if (!is.null(harmony_settings$compress_output)) isTRUE(harmony_settings$compress_output) else TRUE
  output_file <- easycoloc_harmonized_cache_path(
    save_dir = cfg_global$harmonize_dir,
    dataset_id = ds$id,
    source_build = source_build,
    target_build = target_build,
    compress_output = compress_output
  )
  cache_candidates <- easycoloc_harmonized_cache_candidates(
    save_dir = cfg_global$harmonize_dir,
    dataset_id = ds$id,
    source_build = source_build,
    target_build = target_build,
    compress_output = compress_output
  )
  cache_exists_before <- any(file.exists(cache_candidates))
  message(glue("[HARMONIZE] Start {ds$id}"))
  harm_dt <- prepare_gwas_harmony(ds, qtl_build = qtl_build, n_threads = n_threads, prefer_cache = !isTRUE(opts$force))
  existing_cache <- cache_candidates[file.exists(cache_candidates)]
  actual_output_file <- if (length(existing_cache) > 0L) existing_cache[[1L]] else output_file
  rows[[length(rows) + 1L]] <- data.table(
    GWAS_ID = ds$id,
    input_file = ds$file,
    output_file = actual_output_file,
    rows = nrow(harm_dt),
    columns = paste(names(harm_dt), collapse = ","),
    reused_cache = cache_exists_before && !isTRUE(opts$force)
  )
  message(glue("[HARMONIZE] Done {ds$id}: {nrow(harm_dt)} rows -> {actual_output_file}"))
}

summary_dt <- data.table::rbindlist(rows, fill = TRUE)
summary_file <- file.path(cfg_global$harmonize_dir, "harmonized_gwas_manifest.tsv")
data.table::fwrite(summary_dt, summary_file, sep = "\t", quote = FALSE, na = "NA")
message(glue("[HARMONIZE] Manifest: {summary_file}"))
