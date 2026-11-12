#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(yaml)
})

source("src/utils_config.R")

args <- commandArgs(trailingOnly = TRUE)

parse_doctor_args <- function(args) {
  parsed <- list()
  idx <- 1L
  while (idx <= length(args)) {
    arg <- args[[idx]]
    if (arg %in% c("--global", "--gwas", "--qtl")) {
      if (idx == length(args)) {
        stop("Missing value for ", arg, call. = FALSE)
      }
      parsed[[substring(arg, 3L)]] <- args[[idx + 1L]]
      idx <- idx + 2L
    } else if (arg %in% c("--help", "-h")) {
      cat("Usage: Rscript tools/doctor_easycoloc.R [--global PATH] [--gwas PATH] [--qtl PATH]\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  parsed
}

doctor_args <- parse_doctor_args(args)
if (!is.null(doctor_args$global)) Sys.setenv(EASYCOLOC_GLOBAL_CONFIG = doctor_args$global)
if (!is.null(doctor_args$gwas)) Sys.setenv(EASYCOLOC_GWAS_CONFIG = doctor_args$gwas)
if (!is.null(doctor_args$qtl)) Sys.setenv(EASYCOLOC_QTL_CONFIG = doctor_args$qtl)

cfg_bundle <- easycoloc_read_configs()

results <- data.table(
  status = character(),
  category = character(),
  item = character(),
  detail = character()
)

record_check <- function(status, category, item, detail) {
  results <<- rbind(
    results,
    data.table(
      status = status,
      category = category,
      item = item,
      detail = as.character(detail)
    ),
    fill = TRUE
  )
}

check_scalar_path <- function(path_value, category, item, required = TRUE, allow_missing = FALSE) {
  if (is.list(path_value)) {
    if (length(path_value) == 0) {
      if (required) {
        record_check("FAIL", category, item, "path is missing")
      } else {
        record_check("WARN", category, item, "path not configured")
      }
      return(invisible(FALSE))
    }
    ok <- TRUE
    for (path_name in names(path_value)) {
      child_item <- if (nzchar(path_name)) paste0(item, "_", path_name) else item
      child_ok <- check_scalar_path(
        path_value[[path_name]],
        category,
        child_item,
        required = required,
        allow_missing = allow_missing
      )
      ok <- ok && isTRUE(child_ok)
    }
    return(invisible(ok))
  }
  if (length(path_value) > 1) {
    ok <- TRUE
    for (idx in seq_along(path_value)) {
      child_item <- paste0(item, "_", idx)
      child_ok <- check_scalar_path(
        path_value[[idx]],
        category,
        child_item,
        required = required,
        allow_missing = allow_missing
      )
      ok <- ok && isTRUE(child_ok)
    }
    return(invisible(ok))
  }
  if (is.null(path_value) || length(path_value) == 0 || is.na(path_value) || !nzchar(path_value)) {
    if (required) {
      record_check("FAIL", category, item, "path is missing")
    } else {
      record_check("WARN", category, item, "path not configured")
    }
    return(invisible(FALSE))
  }
  if (allow_missing || file.exists(path_value) || dir.exists(path_value)) {
    record_check("OK", category, item, path_value)
    return(invisible(TRUE))
  }
  if (required) {
    record_check("FAIL", category, item, paste("missing:", path_value))
  } else {
    record_check("WARN", category, item, paste("missing:", path_value))
  }
  invisible(FALSE)
}

check_plink_prefix <- function(prefix_path) {
  if (is.null(prefix_path) || is.na(prefix_path) || !nzchar(prefix_path)) {
    record_check("WARN", "reference", "plink_hg38", "not configured")
    return(invisible(FALSE))
  }
  suffixes <- c(".bed", ".bim", ".fam")
  missing <- suffixes[!file.exists(paste0(prefix_path, suffixes))]
  if (length(missing) == 0) {
    record_check("OK", "reference", "plink_hg38", prefix_path)
    return(invisible(TRUE))
  }
  record_check(
    "WARN",
    "reference",
    "plink_hg38",
    paste0("missing PLINK files for prefix ", prefix_path, ": ", paste(missing, collapse = ", "))
  )
  invisible(FALSE)
}

check_recomb_prefix <- function(recom_value, item = "recombination_map", required = FALSE) {
  if (is.list(recom_value)) {
    if (length(recom_value) == 0) {
      if (required) {
        record_check("FAIL", "reference", item, "path is missing")
      } else {
        record_check("WARN", "reference", item, "path not configured")
      }
      return(invisible(FALSE))
    }
    ok <- TRUE
    for (path_name in names(recom_value)) {
      child_item <- if (nzchar(path_name)) paste0(item, "_", path_name) else item
      child_ok <- check_recomb_prefix(
        recom_value[[path_name]],
        item = child_item,
        required = required
      )
      ok <- ok && isTRUE(child_ok)
    }
    return(invisible(ok))
  }

  if (is.null(recom_value) || length(recom_value) == 0 || is.na(recom_value) || !nzchar(recom_value)) {
    if (required) {
      record_check("FAIL", "reference", item, "path is missing")
    } else {
      record_check("WARN", "reference", item, "path not configured")
    }
    return(invisible(FALSE))
  }

  chroms <- as.character(seq_len(22))
  expected_files <- c(
    paste0(recom_value, "_recombination_map_hapmap_format_hg38_chr_", chroms, ".txt"),
    paste0(recom_value, "_recombination_map_hg38_chr_", chroms, ".bed")
  )
  if (any(file.exists(expected_files))) {
    record_check("OK", "reference", item, recom_value)
    return(invisible(TRUE))
  }

  detail <- paste0("missing recombination map files for prefix ", recom_value)
  if (required) {
    record_check("FAIL", "reference", item, detail)
  } else {
    record_check("WARN", "reference", item, detail)
  }
  invisible(FALSE)
}

check_required_names <- function(x, required_names, category, item) {
  missing <- setdiff(required_names, names(x))
  if (length(missing) == 0) {
    record_check("OK", category, item, "required keys present")
    return(invisible(TRUE))
  }
  record_check("FAIL", category, item, paste("missing keys:", paste(missing, collapse = ", ")))
  invisible(FALSE)
}

required_packages <- c(
  "data.table", "dplyr", "yaml", "glue", "jsonlite",
  "coloc", "susieR", "ggplot2"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) == 0) {
  record_check("OK", "runtime", "R_packages", "core packages available")
} else {
  record_check("WARN", "runtime", "R_packages", paste("missing:", paste(missing_packages, collapse = ", ")))
}

for (tool_name in c("Rscript", "plink", "python", "tabix", "liftOver")) {
  tool_path <- Sys.which(tool_name)
  if (nzchar(tool_path)) {
    record_check("OK", "runtime", paste0("tool_", tool_name), tool_path)
  } else {
    record_check("WARN", "runtime", paste0("tool_", tool_name), "not found on PATH")
  }
}

record_check("OK", "config", "global_yaml", cfg_bundle$paths$global)
record_check("OK", "config", "gwas_yaml", cfg_bundle$paths$gwas)
record_check("OK", "config", "qtl_yaml", cfg_bundle$paths$qtl)

cfg_global <- cfg_bundle$global
cfg_gwas <- cfg_bundle$gwas
cfg_qtl <- cfg_bundle$qtl

check_required_names(
  cfg_global,
  c("output_dir", "temp_dir", "clump", "coloc_settings", "harmonization_settings"),
  "config",
  "global_structure"
)

if (!is.null(cfg_global$output_dir) && nzchar(cfg_global$output_dir)) {
  output_parent <- dirname(cfg_global$output_dir)
  if (dir.exists(output_parent)) {
    record_check("OK", "config", "output_parent", output_parent)
  } else {
    record_check("WARN", "config", "output_parent", paste("parent directory missing:", output_parent))
  }
}

check_plink_prefix(cfg_global$plink_hg38)
check_scalar_path(cfg_global$plink_keep, "reference", "plink_keep", required = FALSE)
check_recomb_prefix(cfg_global$recom, "recombination_map", required = FALSE)
check_scalar_path(cfg_global$gene_anno, "reference", "gene_annotation", required = FALSE)
check_scalar_path(cfg_global$ref_genome_hg19, "reference", "ref_genome_hg19", required = FALSE)
check_scalar_path(cfg_global$ref_genome_hg38, "reference", "ref_genome_hg38", required = FALSE)
check_scalar_path(cfg_global$`1kg_af`, "reference", "1kg_af", required = FALSE)
check_scalar_path(cfg_global$dbsnp_hg19, "reference", "dbsnp_hg19", required = FALSE)
check_scalar_path(cfg_global$dbsnp_hg38, "reference", "dbsnp_hg38", required = FALSE)
check_scalar_path(cfg_global$harmonize_dir, "reference", "harmonize_dir", required = FALSE, allow_missing = TRUE)
check_scalar_path(cfg_global$harmonization_settings$liftover_chain, "reference", "liftover_chain", required = FALSE)
if (!is.null(cfg_global$hash_table_dir) && nzchar(cfg_global$hash_table_dir)) {
  check_scalar_path(cfg_global$hash_table_dir, "reference", "hash_table_dir", required = FALSE)
} else {
  record_check("WARN", "reference", "hash_table_dir", "not configured; SNP lookup falls back to direct/allelic matching")
}

if (is.null(cfg_gwas$datasets) || length(cfg_gwas$datasets) == 0) {
  record_check("FAIL", "gwas", "datasets", "no GWAS datasets configured")
} else {
  record_check("OK", "gwas", "datasets", paste(length(cfg_gwas$datasets), "dataset(s) configured"))
  required_dataset_fields <- c("id", "file", "build", "pop", "type", "columns")
  required_column_map <- c("snp", "chrom", "pos", "a1", "a2", "pval", "beta", "se")
  for (idx in seq_along(cfg_gwas$datasets)) {
    ds <- cfg_gwas$datasets[[idx]]
    label <- if (!is.null(ds$id)) ds$id else paste0("dataset_", idx)
    missing_fields <- required_dataset_fields[!vapply(required_dataset_fields, function(field) !is.null(ds[[field]]), logical(1))]
    if (length(missing_fields) > 0) {
      record_check("FAIL", "gwas", label, paste("missing fields:", paste(missing_fields, collapse = ", ")))
      next
    }
    check_scalar_path(ds$file, "gwas", paste0(label, "_file"), required = TRUE)
    missing_columns <- required_column_map[!required_column_map %in% names(ds$columns)]
    if (length(missing_columns) == 0) {
      record_check("OK", "gwas", paste0(label, "_columns"), "required column map present")
    } else {
      record_check("FAIL", "gwas", paste0(label, "_columns"), paste("missing:", paste(missing_columns, collapse = ", ")))
    }
  }
}

check_required_names(cfg_qtl, c("qtl_info", "QTL_all_header", "QTL_cols"), "qtl", "qtl_structure")
if (!is.null(cfg_qtl$qtl_info)) {
  check_scalar_path(cfg_qtl$qtl_info$file, "qtl", "qtl_summary_file", required = TRUE)
  if (!is.null(cfg_qtl$qtl_info$columns)) {
    missing_qtl_meta <- c("id", "sample_size", "all_filename", "sig_filename")
    missing_qtl_meta <- missing_qtl_meta[!missing_qtl_meta %in% names(cfg_qtl$qtl_info$columns)]
    if (length(missing_qtl_meta) == 0) {
      record_check("OK", "qtl", "qtl_metadata_columns", "required QTL metadata mappings present")
    } else {
      record_check("FAIL", "qtl", "qtl_metadata_columns", paste("missing:", paste(missing_qtl_meta, collapse = ", ")))
    }
  } else {
    record_check("FAIL", "qtl", "qtl_metadata_columns", "qtl_info.columns missing")
  }
}
required_qtl_cols <- c("phenotype", "chrom", "pos", "pval", "beta", "se")
missing_qtl_cols <- required_qtl_cols[!required_qtl_cols %in% names(cfg_qtl$QTL_cols)]
if (length(missing_qtl_cols) == 0) {
  record_check("OK", "qtl", "QTL_cols", "required mappings present")
} else {
  record_check("FAIL", "qtl", "QTL_cols", paste("missing:", paste(missing_qtl_cols, collapse = ", ")))
}

selected_qtl_ids <- easycoloc_selected_ids("EASYCOLOC_QTL_ID")
if (length(selected_qtl_ids) > 0L && file.exists(cfg_qtl$qtl_info$file)) {
  qtl_meta <- tryCatch(
    {
      meta <- data.table::fread(cfg_qtl$qtl_info$file)
      meta <- easycoloc_resolve_qtl_meta_paths(meta, cfg_qtl$qtl_info$file, cfg_qtl)
      easycoloc_filter_qtl_meta(meta, cfg_qtl, selected_ids = selected_qtl_ids)
    },
    error = function(e) {
      record_check("FAIL", "qtl", "selected_datasets", e$message)
      NULL
    }
  )
  if (!is.null(qtl_meta)) {
    record_check("OK", "qtl", "selected_datasets", paste(selected_qtl_ids, collapse = ", "))
    all_col <- cfg_qtl$qtl_info$columns$all_filename
    sig_col <- cfg_qtl$qtl_info$columns$sig_filename
    prefilter_sig_pairs <- isTRUE(cfg_global$analysis$prefilter_sig_pairs)
    for (row_idx in seq_len(nrow(qtl_meta))) {
      qtl_id <- as.character(qtl_meta[[cfg_qtl$qtl_info$columns$id]][row_idx])
      check_scalar_path(qtl_meta[[all_col]][row_idx], "qtl", paste0(qtl_id, "_all_pairs"), required = TRUE)
      check_scalar_path(
        qtl_meta[[sig_col]][row_idx],
        "qtl",
        paste0(qtl_id, "_sig_pairs"),
        required = prefilter_sig_pairs
      )
    }
  }
}

status_levels <- c("FAIL", "WARN", "OK")
results[, status := factor(status, levels = status_levels)]
setorder(results, status, category, item)
results[, status := as.character(status)]

for (row_idx in seq_len(nrow(results))) {
  row <- results[row_idx]
  tag <- switch(row$status,
    FAIL = "[FAIL]",
    WARN = "[WARN]",
    OK = "[OK]"
  )
  cat(tag, sprintf("%-10s", row$category), sprintf("%-30s", row$item), row$detail, "\n")
}

summary_counts <- results[, .N, by = status]
cat(
  "\n[DOCTOR] summary:",
  paste(sprintf("%s=%s", summary_counts$status, summary_counts$N), collapse = " "),
  "\n"
)

quit(status = if (any(results$status == "FAIL")) 1 else 0)
