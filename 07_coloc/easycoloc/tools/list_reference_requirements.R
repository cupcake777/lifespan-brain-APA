#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_config.R")
source("src/utils_reference.R")

args <- commandArgs(trailingOnly = TRUE)

include_qtl_files <- FALSE
parse_args <- function(args) {
  idx <- 1L
  while (idx <= length(args)) {
    arg <- args[[idx]]
    if (arg %in% c("--global", "--gwas", "--qtl")) {
      if (idx == length(args)) {
        stop("Missing value for ", arg, call. = FALSE)
      }
      env_name <- paste0("EASYCOLOC_", toupper(substring(arg, 3L)), "_CONFIG")
      do.call(Sys.setenv, as.list(setNames(args[[idx + 1L]], env_name)))
      idx <- idx + 2L
    } else if (arg == "--include-qtl-files") {
      include_qtl_files <<- TRUE
      idx <- idx + 1L
    } else if (arg %in% c("--help", "-h")) {
      cat("Usage: Rscript tools/list_reference_requirements.R [--global PATH] [--gwas PATH] [--qtl PATH] [--include-qtl-files]\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
}

parse_args(args)
cfg_bundle <- easycoloc_read_configs()
ref_dt <- easycoloc_collect_reference_requirements(cfg_bundle, include_qtl_files = include_qtl_files)
setorder(ref_dt, -required, stage, name)

cat("[REF] Configs\n")
cat("[REF] global:", cfg_bundle$paths$global, "\n")
cat("[REF] gwas:", cfg_bundle$paths$gwas, "\n")
cat("[REF] qtl:", cfg_bundle$paths$qtl, "\n\n")

for (row_idx in seq_len(nrow(ref_dt))) {
  row <- ref_dt[row_idx]
  req_label <- if (isTRUE(row$required)) "required" else "optional"
  size_label <- if (!is.na(row$size_human) && nzchar(row$size_human)) row$size_human else "NA"
  cat(
    sprintf(
      "[REF] %-9s %-16s %-28s status=%-18s size=%-10s format=%s\n",
      req_label,
      row$stage,
      row$name,
      row$status,
      size_label,
      row$format
    )
  )
  cat("      path:", ifelse(is.na(row$path), "NA", row$path), "\n")
  cat("      use :", row$used_for, "\n")
  if (!is.na(row$note) && nzchar(row$note)) {
    cat("      note:", row$note, "\n")
  }
}

summary_dt <- ref_dt[, .(
  total = .N,
  present = sum(status == "present", na.rm = TRUE),
  missing = sum(status %in% c("missing", "missing_path", "not_configured") & required, na.rm = TRUE)
) , by = .(required, stage)]

cat("\n[REF] summary\n")
print(summary_dt)
