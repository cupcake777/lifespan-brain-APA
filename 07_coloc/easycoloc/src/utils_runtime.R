suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(glue)
  library(tools)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

.runtime_tracker <- new.env(parent = emptyenv())
.runtime_state_flush_counter <- 0L

runtime_is_enabled <- function() {
  isTRUE(get0("enabled", envir = .runtime_tracker, ifnotfound = FALSE))
}

runtime_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

runtime_scalar_chr <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }
  x <- as.character(x[[1]])
  if (!nzchar(x) || identical(x, "NA")) {
    return(NA_character_)
  }
  x
}

runtime_scalar_num <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x[[1]]))
}

compute_runtime_config_fingerprint <- function(files) {
  files <- unique(files[file.exists(files)])
  if (length(files) == 0) {
    return(NA_character_)
  }

  hashes <- tryCatch(tools::md5sum(files), error = function(e) NULL)
  if (is.null(hashes) || length(hashes) == 0) {
    return(NA_character_)
  }

  paste(sprintf("%s=%s", basename(names(hashes)), unname(hashes)), collapse = ";")
}

runtime_empty_task_state <- function() {
  data.table(
    task_id = character(),
    status = character(),
    attempts = integer(),
    gwas_id = character(),
    locus_id = character(),
    qtl_id = character(),
    phenotype = character(),
    identification_method = character(),
    result_file = character(),
    plot_file = character(),
    locuscompare_plot_file = character(),
    locuscompare_json_file = character(),
    rds_file = character(),
    pp4 = numeric(),
    n_snps = numeric(),
    started_at = character(),
    completed_at = character(),
    updated_at = character(),
    config_fingerprint = character(),
    error_message = character()
  )
}

runtime_paths <- function() {
  get0("paths", envir = .runtime_tracker, ifnotfound = list())
}

runtime_task_state <- function() {
  get0("task_state", envir = .runtime_tracker, ifnotfound = runtime_empty_task_state())
}

runtime_task_index <- function() {
  get0("task_index", envir = .runtime_tracker, ifnotfound = new.env(parent = emptyenv()))
}

rebuild_runtime_task_index <- function(state = runtime_task_state()) {
  index_env <- new.env(parent = emptyenv())
  if (nrow(state) > 0 && "task_id" %in% names(state)) {
    for (idx in seq_len(nrow(state))) {
      task_id <- runtime_scalar_chr(state$task_id[idx])
      if (!is.na(task_id)) {
        assign(task_id, idx, envir = index_env)
      }
    }
  }
  assign("task_index", index_env, envir = .runtime_tracker)
  invisible(index_env)
}

persist_runtime_task_state <- function() {
  if (!runtime_is_enabled()) {
    return(invisible(FALSE))
  }

  paths <- runtime_paths()
  state <- copy(runtime_task_state())
  tmp_file <- paste0(paths$task_state, ".tmp")
  fwrite(state, tmp_file, sep = "\t", na = "NA")
  file.rename(tmp_file, paths$task_state)
  invisible(TRUE)
}

runtime_should_flush_task_state <- function(interval = 1L) {
  interval <- suppressWarnings(as.integer(interval))
  if (is.na(interval) || interval <= 1L) {
    return(TRUE)
  }
  .runtime_state_flush_counter <<- .runtime_state_flush_counter + 1L
  .runtime_state_flush_counter %% interval == 0L
}

initialize_runtime_tracker <- function(output_dir,
                                       enabled = TRUE,
                                       config_fingerprint = NA_character_) {
  runtime_dir <- file.path(output_dir, "runtime")
  paths <- list(
    runtime_dir = runtime_dir,
    active_run = file.path(runtime_dir, "active_run.json"),
    heartbeat = file.path(runtime_dir, "heartbeat.json"),
    task_state = file.path(runtime_dir, "task_state.tsv"),
    event_log = file.path(runtime_dir, "event_log.ndjson"),
    monitor_snapshot = file.path(runtime_dir, "monitor_snapshot.json")
  )

  assign("enabled", isTRUE(enabled), envir = .runtime_tracker)
  assign("config_fingerprint", runtime_scalar_chr(config_fingerprint), envir = .runtime_tracker)
  assign("paths", paths, envir = .runtime_tracker)
  .runtime_state_flush_counter <<- 0L

  if (!isTRUE(enabled)) {
    assign("task_state", runtime_empty_task_state(), envir = .runtime_tracker)
    rebuild_runtime_task_index(runtime_empty_task_state())
    return(invisible(FALSE))
  }

  dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)

  state <- if (file.exists(paths$task_state)) {
    tryCatch(fread(paths$task_state), error = function(e) runtime_empty_task_state())
  } else {
    runtime_empty_task_state()
  }

  state <- as.data.table(state)
  missing_cols <- setdiff(names(runtime_empty_task_state()), names(state))
  for (col in missing_cols) {
    state[[col]] <- runtime_empty_task_state()[[col]]
  }
  state <- state[, names(runtime_empty_task_state()), with = FALSE]

  assign("task_state", state, envir = .runtime_tracker)
  rebuild_runtime_task_index(state)
  invisible(TRUE)
}

register_active_run <- function(log_file = NA_character_,
                                run_label = NA_character_,
                                parent_pid = NA_integer_,
                                command = "Rscript run_coloc.R") {
  if (!runtime_is_enabled()) {
    return(invisible(FALSE))
  }

  payload <- list(
    timestamp = runtime_timestamp(),
    pid = Sys.getpid(),
    parent_pid = runtime_scalar_num(parent_pid),
    run_label = runtime_scalar_chr(run_label),
    log_file = runtime_scalar_chr(log_file),
    command = runtime_scalar_chr(command),
    config_fingerprint = runtime_scalar_chr(get0("config_fingerprint", envir = .runtime_tracker, ifnotfound = NA_character_))
  )

  write_json(payload, runtime_paths()$active_run, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(TRUE)
}

clear_active_run <- function() {
  if (!runtime_is_enabled()) {
    return(invisible(FALSE))
  }

  active_run_file <- runtime_paths()$active_run
  if (!is.null(active_run_file) && file.exists(active_run_file)) {
    unlink(active_run_file)
  }
  invisible(TRUE)
}

append_runtime_event <- function(level = "INFO",
                                 stage,
                                 message_text = NA_character_,
                                 gwas_id = NA_character_,
                                 locus_id = NA_character_,
                                 qtl_id = NA_character_,
                                 phenotype = NA_character_,
                                 task_id = NA_character_,
                                 extras = list()) {
  if (!runtime_is_enabled()) {
    return(invisible(FALSE))
  }

  payload <- list(
    timestamp = runtime_timestamp(),
    pid = Sys.getpid(),
    level = runtime_scalar_chr(level),
    stage = runtime_scalar_chr(stage),
    message = runtime_scalar_chr(message_text),
    gwas_id = runtime_scalar_chr(gwas_id),
    locus_id = runtime_scalar_chr(locus_id),
    qtl_id = runtime_scalar_chr(qtl_id),
    phenotype = runtime_scalar_chr(phenotype),
    task_id = runtime_scalar_chr(task_id),
    extras = extras
  )

  cat(
    toJSON(payload, auto_unbox = TRUE, null = "null"),
    "\n",
    file = runtime_paths()$event_log,
    append = TRUE
  )
  invisible(TRUE)
}

write_runtime_heartbeat <- function(stage,
                                    message_text = NA_character_,
                                    counters = list(),
                                    extras = list()) {
  if (!runtime_is_enabled()) {
    return(invisible(FALSE))
  }

  payload <- list(
    timestamp = runtime_timestamp(),
    pid = Sys.getpid(),
    stage = runtime_scalar_chr(stage),
    message = runtime_scalar_chr(message_text),
    counters = counters,
    extras = extras,
    config_fingerprint = runtime_scalar_chr(get0("config_fingerprint", envir = .runtime_tracker, ifnotfound = NA_character_))
  )

  write_json(payload, runtime_paths()$heartbeat, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(TRUE)
}

build_task_id <- function(gwas_id, locus_id, qtl_id, phenotype) {
  pieces <- c(gwas_id, locus_id, qtl_id, phenotype)
  pieces <- vapply(
    pieces,
    function(x) {
      x <- runtime_scalar_chr(x)
      if (is.na(x)) {
        return("NA")
      }
      gsub("[^[:alnum:]_.:-]+", "_", x)
    },
    character(1)
  )
  paste(pieces, collapse = "__")
}

runtime_task_completed <- function(task_id, config_fingerprint = NA_character_) {
  if (!runtime_is_enabled()) {
    return(FALSE)
  }

  task_id <- runtime_scalar_chr(task_id)
  if (is.na(task_id)) {
    return(FALSE)
  }

  idx <- get0(task_id, envir = runtime_task_index(), ifnotfound = NA_integer_)
  if (is.na(idx)) {
    return(FALSE)
  }
  state <- runtime_task_state()
  if (idx < 1L || idx > nrow(state)) {
    return(FALSE)
  }
  hit <- state[idx]

  if (!is.na(config_fingerprint) &&
    "config_fingerprint" %in% names(hit) &&
    !all(is.na(hit$config_fingerprint)) &&
    !any(hit$config_fingerprint == config_fingerprint, na.rm = TRUE)) {
    return(FALSE)
  }

  any(hit$status == "completed", na.rm = TRUE)
}

runtime_update_task <- function(task_id,
                                status = NULL,
                                gwas_id = NULL,
                                locus_id = NULL,
                                qtl_id = NULL,
                                phenotype = NULL,
                                identification_method = NULL,
                                result_file = NULL,
                                plot_file = NULL,
                                locuscompare_plot_file = NULL,
                                locuscompare_json_file = NULL,
                                rds_file = NULL,
                                pp4 = NULL,
                                n_snps = NULL,
                                error_message = NULL,
                                persist = TRUE) {
  if (!runtime_is_enabled()) {
    return(invisible(FALSE))
  }

  task_id <- runtime_scalar_chr(task_id)
  if (is.na(task_id)) {
    return(invisible(FALSE))
  }

  state <- runtime_task_state()
  idx <- get0(task_id, envir = runtime_task_index(), ifnotfound = NA_integer_)
  if (is.na(idx)) {
    new_row <- runtime_empty_task_state()[, lapply(.SD, function(col) {
      if (is.integer(col)) {
        return(NA_integer_)
      }
      if (is.numeric(col)) {
        return(NA_real_)
      }
      NA_character_
    })]
    state <- rbind(state, new_row, fill = TRUE)
    idx <- nrow(state)
    state$task_id[idx] <- task_id
    state$attempts[idx] <- 0L
    assign(task_id, idx, envir = runtime_task_index())
  }

  now <- runtime_timestamp()
  new_status <- runtime_scalar_chr(status)

  if (!is.na(new_status)) {
    state$status[idx] <- new_status
    if (identical(new_status, "running")) {
      current_attempts <- suppressWarnings(as.integer(state$attempts[idx]))
      if (is.na(current_attempts)) current_attempts <- 0L
      state$attempts[idx] <- current_attempts + 1L
      if (is.na(state$started_at[idx]) || !nzchar(state$started_at[idx])) {
        state$started_at[idx] <- now
      }
    }
    if (new_status %in% c("completed", "failed", "skipped")) {
      state$completed_at[idx] <- now
    }
  }

  updates <- list(
    gwas_id = runtime_scalar_chr(gwas_id),
    locus_id = runtime_scalar_chr(locus_id),
    qtl_id = runtime_scalar_chr(qtl_id),
    phenotype = runtime_scalar_chr(phenotype),
    identification_method = runtime_scalar_chr(identification_method),
    result_file = runtime_scalar_chr(result_file),
    plot_file = runtime_scalar_chr(plot_file),
    locuscompare_plot_file = runtime_scalar_chr(locuscompare_plot_file),
    locuscompare_json_file = runtime_scalar_chr(locuscompare_json_file),
    rds_file = runtime_scalar_chr(rds_file),
    pp4 = runtime_scalar_num(pp4),
    n_snps = runtime_scalar_num(n_snps),
    error_message = runtime_scalar_chr(error_message)
  )

  for (nm in names(updates)) {
    value <- updates[[nm]]
    if (!is.na(value)) state[[nm]][idx] <- value
  }

  state$updated_at[idx] <- now
  state$config_fingerprint[idx] <- runtime_scalar_chr(
    get0("config_fingerprint", envir = .runtime_tracker, ifnotfound = NA_character_)
  )

  assign("task_state", state, envir = .runtime_tracker)
  if (isTRUE(persist)) {
    persist_runtime_task_state()
  }
  invisible(TRUE)
}

runtime_get_task_record <- function(task_id) {
  if (!runtime_is_enabled()) {
    return(NULL)
  }

  task_id <- runtime_scalar_chr(task_id)
  if (is.na(task_id)) {
    return(NULL)
  }

  idx <- get0(task_id, envir = runtime_task_index(), ifnotfound = NA_integer_)
  if (is.na(idx)) {
    return(NULL)
  }
  state <- runtime_task_state()
  if (idx < 1L || idx > nrow(state)) {
    return(NULL)
  }
  state[idx]
}

.gene_label_cache <- new.env(parent = emptyenv())

resolve_gene_label <- function(phenotype_label) {
  phenotype_label <- runtime_scalar_chr(phenotype_label)
  if (is.na(phenotype_label)) {
    return("NA")
  }

  cache_key <- phenotype_label
  if (exists(cache_key, envir = .gene_label_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gene_label_cache))
  }

  resolved <- phenotype_label

  parsed <- tryCatch(parse_geneID(phenotype_label), error = function(e) NULL)
  if (!is.null(parsed) && !is.null(parsed$geneSymbol) &&
    !is.na(parsed$geneSymbol) && nzchar(parsed$geneSymbol)) {
    resolved <- parsed$geneSymbol
  } else {
    base <- sub("\\..*$", "", phenotype_label)
    if (grepl("^ENSG[0-9]{11}$", base) && requireNamespace("clusterProfiler", quietly = TRUE)) {
      converted <- tryCatch(
        suppressMessages(
          clusterProfiler::bitr(base, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
        ),
        error = function(e) NULL
      )
      if (!is.null(converted) && nrow(converted) > 0 && "SYMBOL" %in% names(converted)) {
        resolved <- converted$SYMBOL[1]
      }
    }
  }

  resolved <- sanitize_filename(resolved)
  if (!nzchar(resolved)) resolved <- "NA"
  assign(cache_key, resolved, envir = .gene_label_cache)
  resolved
}
