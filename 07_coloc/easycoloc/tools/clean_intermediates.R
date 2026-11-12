#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || args[1] %in% c("--help", "-h")) {
  cat("Usage: Rscript tools/clean_intermediates.R RESULTS_OR_PROJECT_DIR [--apply]\n")
  quit(status = 0)
}

root <- normalizePath(path.expand(args[1]), mustWork = FALSE)
apply_cleanup <- "--apply" %in% args
if (basename(root) == "results") {
  project_dir <- normalizePath(file.path(root, ".."), mustWork = FALSE)
  results_dir <- root
} else {
  project_dir <- root
  results_dir <- file.path(project_dir, "results")
}

collect_existing <- function(paths, category, reason) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) {
    return(data.table())
  }
  data.table(
    path = normalizePath(paths, mustWork = FALSE),
    category = category,
    reason = reason,
    size_bytes = file.info(paths)$size
  )
}

runtime_dir <- file.path(results_dir, "runtime")
candidates <- rbindlist(
  list(
    collect_existing(
      c(
        file.path(results_dir, "deduplicated_colocalization_results.csv"),
        list.files(results_dir, pattern = "^significant_unique_trait_phenotype_PP4_.*\\.csv$", full.names = TRUE)
      ),
      "stale_summary",
      "Legacy PP4-max deduplicated output; no longer part of EasyColoc summary semantics."
    ),
    collect_existing(
      list.files(project_dir, pattern = "^(coloc|easycoloc)_.*\\.(out|err)$", full.names = TRUE),
      "job_log",
      "Scheduler stdout/stderr logs; keep only if debugging a failed run."
    ),
    collect_existing(
      list.files(file.path(project_dir, "temp"), pattern = "^(clump_output_|clump_input_|easycoloc_.*_keep_)", full.names = TRUE),
      "plink_temp",
      "Successful PLINK clump temporary files; failed-run diagnostics should be reviewed before deletion."
    ),
    collect_existing(
      file.path(runtime_dir, c("active_run.json", "monitor_snapshot.json")),
      "runtime_snapshot",
      "Ephemeral monitor state; not needed after run completion."
    )
  ),
  fill = TRUE
)

if (nrow(candidates) == 0) {
  cat("[CLEAN] no cleanable intermediate files found\n")
  quit(status = 0)
}

setorder(candidates, category, path)
candidates[, size_mb := round(size_bytes / 1024^2, 3)]
print(candidates[, .(category, size_mb, path, reason)])
cat("[CLEAN] candidates:", nrow(candidates), "files,", round(sum(candidates$size_bytes, na.rm = TRUE) / 1024^2, 3), "MB\n")

if (isTRUE(apply_cleanup)) {
  unlink(candidates$path)
  cat("[CLEAN] deleted candidates\n")
} else {
  cat("[CLEAN] dry run only; rerun with --apply to delete these files\n")
}
