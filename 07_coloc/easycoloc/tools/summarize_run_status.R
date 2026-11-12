#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
source("src/utils_config.R")
cfg_bundle <- if (length(args) >= 1 && nzchar(args[1])) NULL else easycoloc_try_read_configs()
cfg_gwas <- if (!is.null(cfg_bundle)) cfg_bundle$gwas else NULL

output_dir <- easycoloc_resolve_output_dir_arg(args, cfg_bundle = cfg_bundle, required = TRUE)

timestamped_logs <- if (dir.exists(output_dir)) {
  list.files(output_dir, pattern = "^run_easycoloc_(full|bg)_[0-9]{8}_[0-9]{4,6}\\.log$", full.names = TRUE)
} else {
  character(0)
}

latest_log <- NULL
if (length(timestamped_logs) > 0) {
  info <- file.info(timestamped_logs)
  # Filter out empty log files if there are non-empty ones available
  non_empty_logs <- timestamped_logs[info$size > 0]
  if (length(non_empty_logs) > 0) {
    info <- info[non_empty_logs, , drop = FALSE]
    latest_log <- rownames(info)[which.max(info$mtime)]
  } else {
    latest_log <- rownames(info)[which.max(info$mtime)]
  }
}

runtime_dir <- file.path(output_dir, "runtime")
active_run_file <- file.path(runtime_dir, "active_run.json")
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
    tryCatch(as.list(fromJSON(line)), error = function(e) NULL)
  })
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0) return(data.table())
  rbindlist(parsed, fill = TRUE)
}

active_run <- read_json_safely(active_run_file)
heartbeat <- read_json_safely(heartbeat_file)
if (!is.null(active_run) && !is.null(active_run$log_file) && file.exists(active_run$log_file)) {
  latest_log <- active_run$log_file
}

task_status_dt <- if (file.exists(task_state_file)) {
  tryCatch(as.data.table(fread(task_state_file)), error = function(e) data.table())
} else {
  data.table()
}

runtime_event_dt <- read_ndjson_safely(event_log_file)

log_lines <- if (!is.null(latest_log) && file.exists(latest_log)) readLines(latest_log, warn = FALSE) else character(0)
last_started_gwas <- {
  hits <- grep("^Processing GWAS Dataset:", log_lines, value = TRUE)
  if (length(hits) > 0) sub("^Processing GWAS Dataset: ", "", tail(hits, 1)) else NA_character_
}

infer_gwas_ids <- function(cfg_gwas = NULL, task_status_dt = data.table(), log_lines = character(0)) {
  ids <- character(0)

  if (!is.null(cfg_gwas$datasets) && length(cfg_gwas$datasets) > 0) {
    ids <- c(ids, vapply(cfg_gwas$datasets, `[[`, character(1), "id"))
  }

  if (nrow(task_status_dt) > 0 && "gwas_id" %in% names(task_status_dt)) {
    ids <- c(ids, as.character(task_status_dt$gwas_id))
  }

  if (length(log_lines) > 0) {
    ids <- c(ids, sub("^Processing GWAS Dataset: ", "", grep("^Processing GWAS Dataset:", log_lines, value = TRUE)))
  }

  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  sort(ids)
}

count_files_for_ids <- function(subdir, pattern_suffix, gwas_ids) {
  dir_path <- file.path(output_dir, subdir)
  if (!dir.exists(dir_path) || length(gwas_ids) == 0) {
    return(setNames(integer(length(gwas_ids)), gwas_ids))
  }

  files <- list.files(dir_path, full.names = FALSE)
  setNames(
    vapply(gwas_ids, function(id) {
      escaped_id <- easycoloc_regex_escape(id)
      sum(grepl(paste0("^", escaped_id, ".*", pattern_suffix, "$"), files))
    }, integer(1)),
    gwas_ids
  )
}

gwas_ids <- infer_gwas_ids(cfg_gwas = cfg_gwas, task_status_dt = task_status_dt, log_lines = log_lines)

if (length(gwas_ids) > 0) {
  abf_counts <- count_files_for_ids("abf", "_locus_results\\.csv", gwas_ids)
  rds_counts <- count_files_for_ids("rds", "\\.rds", gwas_ids)
  plot_counts <- count_files_for_ids("plots", "\\.(pdf|png)", gwas_ids)
  susie_counts <- count_files_for_ids("susie", "_susie\\.csv", gwas_ids)

  summary_dt <- data.table(
    gwas_id = gwas_ids,
    abf = as.integer(abf_counts[gwas_ids]),
    rds = as.integer(rds_counts[gwas_ids]),
    plots = as.integer(plot_counts[gwas_ids]),
    susie = as.integer(susie_counts[gwas_ids])
  )
  summary_dt[, started := gwas_id %in% gsub("^.*: ", "", grep("^Processing GWAS Dataset:", log_lines, value = TRUE))]
  is_terminal_heartbeat <- !is.null(heartbeat) &&
    !is.null(heartbeat$stage) &&
    heartbeat$stage %in% c("pipeline_complete", "pipeline_failed", "summary_failed", "report_failed")
  summary_dt[, current := !is_terminal_heartbeat & !is.na(last_started_gwas) & gwas_id == last_started_gwas]
} else {
  summary_dt <- data.table(
    gwas_id = character(),
    abf = integer(),
    rds = integer(),
    plots = integer(),
    susie = integer(),
    started = logical(),
    current = logical()
  )
}

txt_file <- file.path(output_dir, "run_status_summary.txt")
json_file <- file.path(output_dir, "run_status_summary.json")

status_text <- c(
  paste0("output_dir: ", output_dir),
  paste0("latest_log: ", ifelse(is.null(latest_log), "NA", latest_log)),
  paste0("last_started_gwas: ", ifelse(is.na(last_started_gwas), "NA", last_started_gwas)),
  paste0("active_run_pid: ", ifelse(is.null(active_run$pid), "NA", active_run$pid)),
  paste0("heartbeat_stage: ", ifelse(is.null(heartbeat$stage), "NA", heartbeat$stage)),
  paste0("heartbeat_time: ", ifelse(is.null(heartbeat$timestamp), "NA", heartbeat$timestamp)),
  "",
  capture.output(print(summary_dt)),
  "",
  if (nrow(task_status_dt) > 0) {
    c("task_status_counts:", capture.output(print(task_status_dt[, .N, by = status][order(status)])))
  } else {
    "task_status_counts: none"
  }
)

status_payload <- list(
  output_dir = output_dir,
  latest_log = latest_log,
  last_started_gwas = last_started_gwas,
  active_run = active_run,
  heartbeat = heartbeat,
  summary = summary_dt,
  task_status_counts = if (nrow(task_status_dt) > 0) task_status_dt[, .N, by = status][order(status)] else data.table(),
  runtime_event_counts = if (nrow(runtime_event_dt) > 0 && "stage" %in% names(runtime_event_dt)) {
    runtime_event_dt[, .N, by = stage][order(stage)]
  } else {
    data.table()
  }
)

if (easycoloc_is_writable_dir(output_dir)) {
  writeLines(status_text, txt_file)
  write_json(status_payload, json_file, pretty = TRUE, auto_unbox = TRUE)
  cat("[STATUS] wrote:", txt_file, "\n")
  cat("[STATUS] wrote:", json_file, "\n")
} else {
  cat("[STATUS] skipped writing summary files because output_dir is not writable\n")
}

if (nrow(summary_dt) == 0) {
  cat("[STATUS] no GWAS identifiers inferred from config, runtime state, or logs\n")
}
print(summary_dt)
