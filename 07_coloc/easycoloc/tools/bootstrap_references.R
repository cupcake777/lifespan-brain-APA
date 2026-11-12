#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
})

source("src/utils_config.R")
source("src/utils_download.R")
source("src/utils_reference.R")
source("src/utils_bootstrap.R")

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("Usage:\n")
  cat("  Rscript tools/bootstrap_references.R SOURCE_YAML [--global PATH] [--gwas PATH] [--qtl PATH] [--rewrite-config] [--force]\n")
  cat("  Rscript tools/bootstrap_references.R --demo DEST_DIR [--run] [--force]\n")
  cat("  Rscript tools/bootstrap_references.R --setup-1kg DEST_DIR [--build hg19|hg38] [--pop EAS] [--chromosomes 1-22] [--base-url URL] [--panel-url URL] [--force]\n")
  cat("  Rscript tools/bootstrap_references.R --fetch-gtex-meta DEST_DIR [--gtex-url URL] [--gtex-eqtl-dir DIR] [--gtex-sqtl-dir DIR] [--qtl-config-out PATH] [--gtex-qtl-type eqtl|sqtl] [--force]\n")
  quit(status = 0)
}

parse_args <- function(args) {
  parsed <- list(
    source_yaml = NULL,
    mode = "materialize",
    dest_dir = NULL,
    rewrite_config = FALSE,
    force = FALSE,
    run_demo = FALSE,
    build = "hg19",
    pop = "ALL",
    chromosomes = "1-22",
    base_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp",
    panel_url = NULL,
    gtex_url = "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
    gtex_eqtl_dir = NULL,
    gtex_sqtl_dir = NULL,
    qtl_config_out = NULL,
    gtex_qtl_type = "eqtl"
  )
  idx <- 1L
  while (idx <= length(args)) {
    arg <- args[[idx]]
    if (arg %in% c("--global", "--gwas", "--qtl")) {
      if (idx == length(args)) stop("Missing value for ", arg, call. = FALSE)
      env_name <- paste0("EASYCOLOC_", toupper(substring(arg, 3L)), "_CONFIG")
      do.call(Sys.setenv, as.list(setNames(args[[idx + 1L]], env_name)))
      idx <- idx + 2L
    } else if (arg == "--rewrite-config") {
      parsed$rewrite_config <- TRUE
      idx <- idx + 1L
    } else if (arg == "--force") {
      parsed$force <- TRUE
      idx <- idx + 1L
    } else if (arg == "--run") {
      parsed$run_demo <- TRUE
      idx <- idx + 1L
    } else if (arg == "--demo") {
      parsed$mode <- "demo"
      if (idx == length(args)) stop("Missing destination directory for --demo", call. = FALSE)
      parsed$dest_dir <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--setup-1kg") {
      parsed$mode <- "setup_1kg"
      if (idx == length(args)) stop("Missing destination directory for --setup-1kg", call. = FALSE)
      parsed$dest_dir <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--fetch-gtex-meta") {
      parsed$mode <- "fetch_gtex_meta"
      if (idx == length(args)) stop("Missing destination directory for --fetch-gtex-meta", call. = FALSE)
      parsed$dest_dir <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--pop") {
      parsed$pop <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--build") {
      parsed$build <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--chromosomes") {
      parsed$chromosomes <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--base-url") {
      parsed$base_url <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--panel-url") {
      parsed$panel_url <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--gtex-url") {
      parsed$gtex_url <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--gtex-eqtl-dir") {
      parsed$gtex_eqtl_dir <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--gtex-sqtl-dir") {
      parsed$gtex_sqtl_dir <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--qtl-config-out") {
      parsed$qtl_config_out <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg == "--gtex-qtl-type") {
      parsed$gtex_qtl_type <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg %in% c("--help", "-h")) {
      usage()
    } else if (startsWith(arg, "--")) {
      stop("Unknown argument: ", arg, call. = FALSE)
    } else if (parsed$mode != "materialize" && is.null(parsed$dest_dir)) {
      parsed$dest_dir <- arg
      idx <- idx + 1L
    } else if (is.null(parsed$source_yaml)) {
      parsed$source_yaml <- arg
      idx <- idx + 1L
    } else {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }
  }
  if (identical(parsed$mode, "materialize") && is.null(parsed$source_yaml)) {
    stop("SOURCE_YAML is required", call. = FALSE)
  }
  if (!identical(parsed$mode, "materialize") && is.null(parsed$dest_dir)) {
    stop("Destination directory is required for mode: ", parsed$mode, call. = FALSE)
  }
  parsed
}

parsed <- parse_args(args)
cfg_bundle <- NULL
if (identical(parsed$mode, "materialize") || isTRUE(parsed$rewrite_config)) {
  cfg_bundle <- easycoloc_read_configs()
}

if (identical(parsed$mode, "materialize")) {
  source_yaml <- normalizePath(parsed$source_yaml, mustWork = FALSE)
  source_cfg <- yaml::read_yaml(source_yaml)

  manifest <- easycoloc_bootstrap_apply(
    cfg_bundle = cfg_bundle,
    source_cfg = source_cfg,
    rewrite_config = parsed$rewrite_config,
    force = parsed$force
  )

  manifest_path <- file.path(cfg_bundle$project_root, "reference_bootstrap_manifest.tsv")
  json_path <- file.path(cfg_bundle$project_root, "reference_bootstrap_manifest.json")
  fwrite(manifest, manifest_path, sep = "\t")
  write_json(
    list(
      project_root = cfg_bundle$project_root,
      global_config = cfg_bundle$paths$global,
      source_yaml = source_yaml,
      rewrite_config = parsed$rewrite_config,
      force = parsed$force,
      manifest = manifest
    ),
    json_path,
    pretty = TRUE,
    auto_unbox = TRUE
  )

  cat("[BOOTSTRAP] source_yaml:", source_yaml, "\n")
  cat("[BOOTSTRAP] project_root:", cfg_bundle$project_root, "\n")
  cat("[BOOTSTRAP] rewrite_config:", parsed$rewrite_config, "\n")
  cat("[BOOTSTRAP] wrote:", manifest_path, "\n")
  cat("[BOOTSTRAP] wrote:", json_path, "\n")

  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i]
    cat(
      sprintf(
        "[BOOTSTRAP] %-16s action=%-8s mode=%-7s target=%s\n",
        row$resource,
        row$action,
        row$mode,
        row$target_rel
      )
    )
  }
} else if (identical(parsed$mode, "demo")) {
  manifest <- easycoloc_bootstrap_demo(
    dest_dir = parsed$dest_dir,
    force = parsed$force,
    run_pipeline = parsed$run_demo
  )
  print(manifest)
} else if (identical(parsed$mode, "setup_1kg")) {
  manifest <- easycoloc_setup_1kg(
    dest_dir = parsed$dest_dir,
    pop = parsed$pop,
    build = parsed$build,
    chromosomes = parsed$chromosomes,
    base_url = parsed$base_url,
    panel_url = parsed$panel_url,
    force = parsed$force,
    rewrite_config = parsed$rewrite_config,
    cfg_bundle = cfg_bundle
  )
  print(manifest)
} else if (identical(parsed$mode, "fetch_gtex_meta")) {
  if (is.null(parsed$qtl_config_out) && !is.null(cfg_bundle)) {
    parsed$qtl_config_out <- file.path(dirname(cfg_bundle$paths$qtl), "qtl_gtex_generated.yml")
  }
  manifest <- easycoloc_fetch_gtex_meta(
    dest_dir = parsed$dest_dir,
    url = parsed$gtex_url,
    force = parsed$force,
    eqtl_dir = parsed$gtex_eqtl_dir,
    sqtl_dir = parsed$gtex_sqtl_dir,
    config_output_file = parsed$qtl_config_out,
    config_qtl_type = parsed$gtex_qtl_type
  )
  print(manifest)
}
