#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || args[1] %in% c("--help", "-h")) {
  cat("Usage: Rscript tools/finalize_run.R RESULTS_DIR\n")
  quit(status = 0)
}

source("src/utils_config.R")

cfg_bundle <- if (length(args) >= 1 && nzchar(args[1])) NULL else easycoloc_try_read_configs()
output_dir <- easycoloc_resolve_output_dir_arg(args, cfg_bundle = cfg_bundle, required = TRUE)
if (length(args) >= 1 && nzchar(args[1])) {
  output_parent <- normalizePath(file.path(output_dir, ".."), mustWork = FALSE)
  sibling_config <- file.path(output_parent, "config")
  if (dir.exists(sibling_config) &&
    file.exists(file.path(sibling_config, "global.yml")) &&
    file.exists(file.path(sibling_config, "gwas.yml")) &&
    file.exists(file.path(sibling_config, "qtl.yml"))) {
    Sys.setenv(
      EASYCOLOC_GLOBAL_CONFIG = file.path(sibling_config, "global.yml"),
      EASYCOLOC_GWAS_CONFIG = file.path(sibling_config, "gwas.yml"),
      EASYCOLOC_QTL_CONFIG = file.path(sibling_config, "qtl.yml")
    )
    cfg_bundle <- easycoloc_try_read_configs()
  } else {
    cfg_bundle <- easycoloc_try_read_configs()
  }
} else {
  cfg_bundle <- easycoloc_try_read_configs()
}
cfg_global <- if (!is.null(cfg_bundle)) cfg_bundle$global else NULL

utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

pp4_threshold <- if (!is.null(cfg_global$coloc_settings$pp4_threshold)) {
  as.numeric(cfg_global$coloc_settings$pp4_threshold)
} else if (!is.null(cfg_global$sig_threshold)) {
  as.numeric(cfg_global$sig_threshold)
} else {
  0.8
}

message(glue("[FINALIZE] output_dir: {output_dir}"))
runtime_dir <- file.path(output_dir, "runtime")
shard_runtime_root <- file.path(runtime_dir, "shards")
if (dir.exists(shard_runtime_root)) {
  shard_task_files <- list.files(shard_runtime_root, pattern = "^task_state\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(shard_task_files) > 0) {
    task_state <- rbindlist(lapply(shard_task_files, function(path) {
      dt <- tryCatch(fread(path), error = function(e) data.table())
      if (nrow(dt) > 0) {
        dt[, shard_runtime := basename(dirname(path))]
      }
      dt
    }), fill = TRUE)
    if (nrow(task_state) > 0) {
      if ("task_id" %in% names(task_state) && "updated_at" %in% names(task_state)) {
        setorder(task_state, task_id, updated_at)
        task_state <- task_state[!duplicated(task_id, fromLast = TRUE)]
      }
      dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
      fwrite(task_state, file.path(runtime_dir, "task_state.tsv"), sep = "\t", na = "NA")
      message(glue("[FINALIZE] merged shard task states: {nrow(task_state)} task(s)"))
    }
  }

  shard_event_files <- list.files(shard_runtime_root, pattern = "^event_log\\.ndjson$", recursive = TRUE, full.names = TRUE)
  if (length(shard_event_files) > 0) {
    dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
    event_lines <- unlist(lapply(shard_event_files, function(path) {
      lines <- readLines(path, warn = FALSE)
      lines[nzchar(lines)]
    }), use.names = FALSE)
    writeLines(event_lines, file.path(runtime_dir, "event_log.ndjson"), useBytes = TRUE)
    message(glue("[FINALIZE] merged shard event logs: {length(event_lines)} event(s)"))
  }
}

merge_all_results(
  output_dir = output_dir,
  pp4_threshold = pp4_threshold,
  merge_susie = TRUE,
  save_summary = TRUE
)

report_file <- file.path(output_dir, "coloc_report.html")
if (file.exists(file.path("src", "utils_report.R"))) {
  source(file.path("src", "utils_report.R"))
}
if (exists("generate_html_report", mode = "function")) {
  generate_html_report(
    results_dir = output_dir,
    output_file = report_file,
    project_name = "EasyColoc Analysis"
  )
  message(glue("[FINALIZE] report: {report_file}"))
}

manifest_script <- file.path("tools", "build_output_manifest.R")
if (file.exists(manifest_script)) {
  status <- system2("Rscript", c(manifest_script, output_dir), stdout = TRUE, stderr = TRUE)
  cat(paste(status, collapse = "\n"), "\n")
}

check_script <- file.path("tools", "check_run_completion.R")
if (file.exists(check_script)) {
  status <- system2("Rscript", c(check_script, output_dir), stdout = TRUE, stderr = TRUE)
  cat(paste(status, collapse = "\n"), "\n")
  check_code <- attr(status, "status")
  if (!is.null(check_code) && !identical(as.integer(check_code), 0L)) {
    warning("[FINALIZE] completion check reported non-zero status; inspect runtime state for shard-level failures", call. = FALSE)
  }
}

message("[FINALIZE] complete")
