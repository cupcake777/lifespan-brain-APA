#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
source("src/utils_config.R")

snapshot_file <- NULL
write_snapshot <- TRUE
output_arg <- NULL

idx <- 1L
while (idx <= length(args)) {
  arg <- args[[idx]]
  if (identical(arg, "--snapshot-file")) {
    if (idx == length(args)) {
      stop("Missing value for --snapshot-file", call. = FALSE)
    }
    snapshot_file <- normalizePath(path.expand(args[[idx + 1L]]), mustWork = FALSE)
    idx <- idx + 2L
  } else if (identical(arg, "--no-write-snapshot")) {
    write_snapshot <- FALSE
    idx <- idx + 1L
  } else if (identical(arg, "--help") || identical(arg, "-h")) {
    cat("Usage: Rscript tools/monitor_easycoloc.R [--no-write-snapshot] [--snapshot-file PATH] [OUTPUT_DIR]\n")
    quit(status = 0)
  } else if (startsWith(arg, "--")) {
    stop("Unknown argument: ", arg, call. = FALSE)
  } else if (is.null(output_arg)) {
    output_arg <- arg
    idx <- idx + 1L
  } else {
    stop("Unexpected positional argument: ", arg, call. = FALSE)
  }
}

cfg_bundle <- if (is.null(output_arg)) easycoloc_try_read_configs() else NULL
output_dir <- easycoloc_resolve_output_dir_arg(
  args = if (is.null(output_arg)) character(0) else output_arg,
  cfg_bundle = cfg_bundle,
  required = TRUE
)

runtime_dir <- file.path(output_dir, "runtime")
active_run_file <- file.path(runtime_dir, "active_run.json")
heartbeat_file <- file.path(runtime_dir, "heartbeat.json")
task_state_file <- file.path(runtime_dir, "task_state.tsv")
default_snapshot_file <- file.path(runtime_dir, "monitor_snapshot.json")

read_json_safely <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

pid_is_alive <- function(pid) {
  pid_num <- suppressWarnings(as.integer(pid))
  if (is.na(pid_num) || pid_num < 1L) return(FALSE)
  status <- suppressWarnings(system2("kill", c("-0", as.character(pid_num)), stdout = FALSE, stderr = FALSE))
  identical(status, 0L)
}

parse_runtime_timestamp <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(NA)
  tz_name <- Sys.getenv("TZ", unset = "UTC")
  as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = tz_name)
}

active_run <- read_json_safely(active_run_file)
heartbeat <- read_json_safely(heartbeat_file)
task_state <- if (file.exists(task_state_file)) {
  tryCatch(as.data.table(fread(task_state_file)), error = function(e) data.table())
} else {
  data.table()
}

heartbeat_ts <- if (!is.null(heartbeat) && !is.null(heartbeat$timestamp)) {
  parse_runtime_timestamp(heartbeat$timestamp)
} else {
  as.POSIXct(NA)
}
heartbeat_recent <- !is.na(heartbeat_ts) &&
  as.numeric(difftime(Sys.time(), heartbeat_ts, units = "secs")) <= 600
heartbeat_matches_active_run <- !is.null(active_run) &&
  !is.null(active_run$pid) &&
  !is.null(heartbeat) &&
  !is.null(heartbeat$pid) &&
  identical(as.integer(active_run$pid), as.integer(heartbeat$pid))
active_run_pid_alive <- !is.null(active_run) &&
  !is.null(active_run$pid) &&
  pid_is_alive(active_run$pid)
heartbeat_is_terminal <- !is.null(heartbeat) &&
  !is.null(heartbeat$stage) &&
  heartbeat$stage %in% c("pipeline_complete", "pipeline_failed", "summary_failed", "report_failed")
active_run_state <- if (is.null(active_run)) {
  "missing"
} else if (isTRUE(heartbeat_is_terminal) && isTRUE(heartbeat_matches_active_run)) {
  "finished"
} else if (isTRUE(active_run_pid_alive) || (isTRUE(heartbeat_recent) && isTRUE(heartbeat_matches_active_run))) {
  "running"
} else {
  "stale"
}

log_candidates <- character(0)
if (dir.exists(output_dir)) {
  log_candidates <- list.files(
    output_dir,
    pattern = "^run_easycoloc_full.*\\.log$",
    full.names = TRUE
  )
}

active_log <- if (isTRUE(active_run_pid_alive) &&
                  !is.null(active_run) &&
                  !is.null(active_run$log_file) &&
                  file.exists(active_run$log_file)) {
  active_run$log_file
} else if (length(log_candidates) > 0) {
  info <- file.info(log_candidates)
  info$path <- rownames(info)
  info <- info[order(info$mtime, info$size, decreasing = TRUE), ]
  non_empty <- info$path[!is.na(info$size) & info$size > 0]
  if (length(non_empty) > 0) non_empty[[1]] else info$path[[1]]
} else {
  NA_character_
}

last_log_line <- NA_character_
if (!is.na(active_log) && file.exists(active_log)) {
  log_lines <- readLines(active_log, warn = FALSE)
  if (length(log_lines) > 0) last_log_line <- tail(log_lines, 1)
}

count_outputs <- function(subdir, pattern) {
  dir_path <- file.path(output_dir, subdir)
  if (!dir.exists(dir_path)) return(0L)
  length(list.files(dir_path, pattern = pattern, full.names = TRUE))
}

task_counts <- if (nrow(task_state) > 0 && "status" %in% names(task_state)) {
  as.list(table(task_state$status, useNA = "ifany"))
} else {
  list()
}

snapshot <- list(
  output_dir = output_dir,
  active_log = active_log,
  active_run = active_run,
  active_run_state = active_run_state,
  active_run_pid_alive = active_run_pid_alive,
  heartbeat_recent = heartbeat_recent,
  heartbeat = heartbeat,
  last_log_line = last_log_line,
  outputs = list(
    abf = count_outputs("abf", "_locus_results\\.csv$"),
    rds = count_outputs("rds", "\\.rds$"),
    plots = count_outputs("plots", "\\.(pdf|png)$"),
    susie = count_outputs("susie", "_susie\\.csv$")
  ),
  task_counts = task_counts,
  task_rows = nrow(task_state)
)

if (isTRUE(write_snapshot)) {
  if (is.null(snapshot_file)) {
    snapshot_file <- default_snapshot_file
  }

  snapshot_dir <- dirname(snapshot_file)
  if (dir.exists(snapshot_dir) && easycoloc_is_writable_dir(snapshot_dir)) {
    write_json(snapshot, snapshot_file, auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    warning("Skipping monitor snapshot write because directory is not writable: ", snapshot_dir, call. = FALSE)
  }
}

cat("[MONITOR] output_dir:", output_dir, "\n")
cat("[MONITOR] active_log:", ifelse(is.na(active_log), "NA", active_log), "\n")
if (!is.null(active_run)) {
  cat("[MONITOR] active_run_state:", active_run_state, "\n")
  cat("[MONITOR] active_run_pid:", ifelse(is.null(active_run$pid), "NA", active_run$pid), "\n")
  cat("[MONITOR] active_run_label:", ifelse(is.null(active_run$run_label), "NA", active_run$run_label), "\n")
}
if (!is.null(heartbeat)) {
  cat("[MONITOR] heartbeat_stage:", ifelse(is.null(heartbeat$stage), "NA", heartbeat$stage), "\n")
  cat("[MONITOR] heartbeat_time:", ifelse(is.null(heartbeat$timestamp), "NA", heartbeat$timestamp), "\n")
  cat("[MONITOR] heartbeat_recent:", heartbeat_recent, "\n")
}
cat("[MONITOR] outputs:",
    sprintf("abf=%s rds=%s plots=%s susie=%s",
            snapshot$outputs$abf,
            snapshot$outputs$rds,
            snapshot$outputs$plots,
            snapshot$outputs$susie),
    "\n")
if (length(task_counts) > 0) {
  cat("[MONITOR] task_status:",
      paste(sprintf("%s=%s", names(task_counts), unlist(task_counts)), collapse = " "),
      "\n")
}
cat("[MONITOR] last_log_line:", ifelse(is.na(last_log_line), "NA", last_log_line), "\n")
if (isTRUE(write_snapshot) && !is.null(snapshot_file)) {
  cat("[MONITOR] snapshot_file:", snapshot_file, "\n")
}
