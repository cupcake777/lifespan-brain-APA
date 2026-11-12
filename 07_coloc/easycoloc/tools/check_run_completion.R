#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
source("src/utils_config.R")

cfg_bundle <- if (length(args) >= 1 && nzchar(args[1])) NULL else easycoloc_try_read_configs()
output_dir <- easycoloc_resolve_output_dir_arg(args, cfg_bundle = cfg_bundle, required = TRUE)
requested_log <- NULL

if (length(args) >= 2 && nzchar(args[2])) {
  requested_log <- normalizePath(args[2], mustWork = FALSE)
}

latest_log <- NULL
if (!is.null(requested_log) && file.exists(requested_log)) {
  latest_log <- requested_log
} else if (dir.exists(output_dir)) {
  timestamped_logs <- list.files(
    output_dir,
    pattern = "^run_easycoloc_full_[0-9]{8}_[0-9]{4,6}\\.log$",
    full.names = TRUE
  )
  fallback_logs <- list.files(
    output_dir,
    pattern = "^run_easycoloc_full.*\\.log$",
    full.names = TRUE
  )

  log_candidates <- if (length(timestamped_logs) > 0) timestamped_logs else fallback_logs
  if (length(log_candidates) > 0) {
    log_info <- file.info(log_candidates)
    latest_log <- rownames(log_info)[which.max(log_info$mtime)]
  }
}

log_lines <- character(0)
if (!is.null(latest_log) && file.exists(latest_log)) {
  log_lines <- readLines(latest_log, warn = FALSE)
}

runtime_dir <- file.path(output_dir, "runtime")
heartbeat_file <- file.path(runtime_dir, "heartbeat.json")
task_state_file <- file.path(runtime_dir, "task_state.tsv")
event_log_file <- file.path(runtime_dir, "event_log.ndjson")

read_json_safely <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

read_ndjson_safely <- function(path) {
  if (!file.exists(path)) return(data.table())
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(data.table())

  parsed <- lapply(lines, function(line) {
    tryCatch(
      {
        record <- as.list(fromJSON(line))
        record[lengths(record) > 0L]
      },
      error = function(e) NULL
    )
  })
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0) return(data.table())
  rbindlist(parsed, fill = TRUE)
}

heartbeat <- read_json_safely(heartbeat_file)
task_state_dt <- if (file.exists(task_state_file)) {
  tryCatch(as.data.table(fread(task_state_file)), error = function(e) data.table())
} else {
  data.table()
}
event_log_dt <- read_ndjson_safely(event_log_file)

current_event_log_dt <- event_log_dt
current_pid <- if (!is.null(heartbeat) && !is.null(heartbeat$pid)) {
  suppressWarnings(as.integer(heartbeat$pid))
} else {
  NA_integer_
}
if (!is.na(current_pid) && nrow(event_log_dt) > 0 && "pid" %in% names(event_log_dt)) {
  event_pid <- suppressWarnings(as.integer(event_log_dt$pid))
  current_event_log_dt <- event_log_dt[event_pid == current_pid]
}

has_complete_marker <- any(grepl("EasyColoc Analysis Complete!", log_lines, fixed = TRUE))
has_summary_marker <- any(grepl("[SUM] Merge complete!", log_lines, fixed = TRUE))
has_report_marker <- any(grepl("[REPORT] Interactive report saved:", log_lines, fixed = TRUE))
has_failed_tasks <- nrow(task_state_dt) > 0 &&
  "status" %in% names(task_state_dt) &&
  any(task_state_dt$status %in% c("failed"), na.rm = TRUE)
has_runtime_failures <- nrow(current_event_log_dt) > 0 &&
  "stage" %in% names(current_event_log_dt) &&
  any(current_event_log_dt$stage %in% c("gwas_failed", "summary_failed", "report_failed"), na.rm = TRUE)
heartbeat_failed <- !is.null(heartbeat) &&
  !is.null(heartbeat$stage) &&
  heartbeat$stage %in% c("pipeline_failed", "summary_failed", "report_failed")
heartbeat_complete <- !is.null(heartbeat) &&
  !is.null(heartbeat$stage) &&
  identical(as.character(heartbeat$stage), "pipeline_complete")
event_summary_complete <- nrow(current_event_log_dt) > 0 &&
  "stage" %in% names(current_event_log_dt) &&
  any(current_event_log_dt$stage == "summary_complete", na.rm = TRUE)
event_report_complete <- nrow(current_event_log_dt) > 0 &&
  "stage" %in% names(current_event_log_dt) &&
  any(current_event_log_dt$stage == "report_complete", na.rm = TRUE)

abf_count <- if (dir.exists(file.path(output_dir, "abf"))) {
  length(list.files(file.path(output_dir, "abf"), pattern = "_locus_results\\.csv$", full.names = TRUE))
} else 0L

rds_count <- if (dir.exists(file.path(output_dir, "rds"))) {
  length(list.files(file.path(output_dir, "rds"), pattern = "\\.rds$", full.names = TRUE))
} else 0L

plot_count <- if (dir.exists(file.path(output_dir, "plots"))) {
  length(list.files(file.path(output_dir, "plots"), pattern = "\\.(pdf|png)$", full.names = TRUE))
} else 0L

susie_count <- if (dir.exists(file.path(output_dir, "susie"))) {
  length(list.files(file.path(output_dir, "susie"), pattern = "_susie\\.csv$", full.names = TRUE))
} else 0L

manifest_file <- file.path(output_dir, "output_manifest.tsv")
manifest_summary_file <- file.path(output_dir, "output_manifest_summary.json")
merged_results_file <- file.path(output_dir, "all_colocalization_results.csv")
report_file <- file.path(output_dir, "coloc_report.html")

status <- if (has_failed_tasks || has_runtime_failures || heartbeat_failed) {
  "FAILED"
} else if (has_complete_marker && has_summary_marker && has_report_marker) {
  "COMPLETE"
} else if (heartbeat_complete && event_summary_complete && event_report_complete) {
  "COMPLETE"
} else if (heartbeat_complete && abf_count > 0) {
  "COMPLETE"
} else if (file.exists(merged_results_file) &&
  file.exists(report_file) &&
  file.exists(manifest_file) &&
  abf_count > 0) {
  "COMPLETE"
} else if (length(log_lines) > 0) {
  "INCOMPLETE"
} else {
  "NOT_STARTED"
}

cat("[CHECK] status:", status, "\n")
cat("[CHECK] output_dir:", output_dir, "\n")
cat("[CHECK] latest_log:", ifelse(is.null(latest_log), "NA", latest_log), "\n")
cat("[CHECK] abf:", abf_count, "rds:", rds_count, "plots:", plot_count, "susie:", susie_count, "\n")
cat("[CHECK] manifest:", ifelse(file.exists(manifest_file), manifest_file, "missing"), "\n")
cat("[CHECK] manifest_summary:", ifelse(file.exists(manifest_summary_file), manifest_summary_file, "missing"), "\n")
cat("[CHECK] heartbeat_stage:", ifelse(is.null(heartbeat$stage), "NA", heartbeat$stage), "\n")
cat("[CHECK] runtime_failures:", has_runtime_failures, "\n")
cat("[CHECK] failed_tasks:", has_failed_tasks, "\n")

if (length(log_lines) > 0) {
  cat("[CHECK] last_log_line:", tail(log_lines, 1), "\n")
}

quit(status = if (identical(status, "COMPLETE")) 0 else 1)
