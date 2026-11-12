#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
source("src/utils_config.R")

cfg_bundle <- if (length(args) >= 1 && nzchar(args[1])) NULL else easycoloc_try_read_configs()
output_dir <- easycoloc_resolve_output_dir_arg(args, cfg_bundle = cfg_bundle, required = TRUE)

if (!dir.exists(output_dir)) {
  stop("Output directory does not exist: ", output_dir, call. = FALSE)
}

classify_output_file <- function(rel_path) {
  if (grepl("^abf/.*_locus_results\\.csv$", rel_path)) return("abf_locus_results")
  if (grepl("^rds/.*\\.rds$", rel_path)) return("rds_bundle")
  if (grepl("^plots/.*\\.pdf$", rel_path)) return("plot_pdf")
  if (grepl("^plots/.*\\.png$", rel_path)) return("plot_png")
  if (grepl("^susie/.*_susie\\.csv$", rel_path)) return("susie_summary")
  if (grepl("^run_easycoloc_full.*\\.log$", rel_path)) return("run_log")
  if (grepl("all_pairs\\.txt\\.gz(\\.tbi)?$", rel_path)) return("qtl_tabix_input")
  if (grepl("signif_pairs\\.txt\\.gz(\\.tbi)?$", rel_path)) return("qtl_significant_input")
  if (grepl("\\.R$", rel_path)) return("aux_script")
  "other"
}

all_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]

manifest_dt <- data.table(
  path = all_files
)
manifest_dt[, rel_path := sub(paste0("^", normalizePath(output_dir, mustWork = FALSE), "/?"), "", path)]
manifest_dt[, size_bytes := file.info(path)$size]
manifest_dt[, mtime := format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S")]
manifest_dt[, extension := tools::file_ext(rel_path)]
manifest_dt[, category := vapply(rel_path, classify_output_file, character(1))]
manifest_dt[, size_mb := round(size_bytes / 1024^2, 3)]
setorder(manifest_dt, category, rel_path)

manifest_file <- file.path(output_dir, "output_manifest.tsv")
fwrite(
  manifest_dt[, .(rel_path, category, extension, size_bytes, size_mb, mtime)],
  manifest_file,
  sep = "\t"
)

category_summary <- manifest_dt[
  ,
  .(
    file_count = .N,
    total_size_bytes = sum(size_bytes, na.rm = TRUE),
    total_size_mb = round(sum(size_bytes, na.rm = TRUE) / 1024^2, 3)
  ),
  by = category
][order(-file_count, category)]

summary_payload <- list(
  output_dir = output_dir,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  total_files = nrow(manifest_dt),
  total_size_bytes = sum(manifest_dt$size_bytes, na.rm = TRUE),
  total_size_mb = round(sum(manifest_dt$size_bytes, na.rm = TRUE) / 1024^2, 3),
  category_summary = category_summary
)

summary_file <- file.path(output_dir, "output_manifest_summary.json")
write_json(summary_payload, summary_file, pretty = TRUE, auto_unbox = TRUE)

cat("[MANIFEST] wrote:", manifest_file, "\n")
cat("[MANIFEST] wrote:", summary_file, "\n")
cat("[MANIFEST] total files:", nrow(manifest_dt), "\n")
