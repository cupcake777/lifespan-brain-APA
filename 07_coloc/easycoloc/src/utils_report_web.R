suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

resolve_current_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else as.character(ofile)
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]

  if (length(frame_files) == 0) {
    return(getwd())
  }

  dirname(normalizePath(frame_files[[length(frame_files)]], mustWork = FALSE))
}

utils_report_web_dir <- resolve_current_script_dir()
source(file.path(utils_report_web_dir, "utils_report.R"))

read_json_if_present <- function(path) {
  if (!file.exists(path)) {
    return(list(missing = TRUE, path = basename(path)))
  }

  parsed <- tryCatch(
    read_json(path, simplifyVector = FALSE),
    error = function(e) list(invalid_json = TRUE, path = basename(path), error = e$message)
  )
  parsed
}

report_clean_path <- function(path) {
  path <- as.character(path)
  if (length(path) == 0 || is.na(path[1]) || !nzchar(path[1])) {
    return("")
  }
  path <- gsub("\\\\", "/", path[1])
  path <- sub("^\\./", "", path)
  path
}

report_existing_rel_path <- function(results_dir, rel_path) {
  rel_path <- report_clean_path(rel_path)
  if (!nzchar(rel_path)) {
    return("")
  }
  if (file.exists(file.path(results_dir, rel_path))) {
    return(rel_path)
  }
  if (file.exists(rel_path)) {
    return(rel_path)
  }
  ""
}

report_match_asset <- function(row, assets, suffix_regex) {
  if (length(assets) == 0) {
    return("")
  }

  candidates <- assets[grepl(suffix_regex, assets, ignore.case = TRUE)]
  if (length(candidates) == 0) {
    return("")
  }

  fields <- c(
    report_value(row, "GWAS_ID", default = ""),
    report_value(row, "QTL_ID", default = ""),
    report_value(row, "Locus", default = ""),
    report_value(row, "Gene", fallback = "Phenotype", default = "")
  )
  tokens <- unique(vapply(fields, function(x) {
    token <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
    token <- gsub("^_+|_+$", "", token)
    token
  }, character(1)))
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) > 0) {
    scored <- vapply(candidates, function(name) {
      sum(vapply(tokens, function(token) grepl(token, name, fixed = TRUE), logical(1)))
    }, integer(1))
    best <- candidates[which.max(scored)]
    if (max(scored) >= 2L) {
      return(best)
    }
  }
  ""
}

normalize_report_results <- function(results_dt, results_dir = NULL, asset_index = NULL) {
  if (is.null(results_dt) || nrow(results_dt) == 0) {
    return(list())
  }

  rows <- split(results_dt, seq_len(nrow(results_dt)))
  out <- vector("list", length(rows))
  has_plot_file_col <- "plot_file" %in% names(results_dt)
  has_locuscompare_plot_file_col <- "locuscompare_plot_file" %in% names(results_dt)
  has_locuscompare_json_file_col <- "locuscompare_json_file" %in% names(results_dt)
  has_any_explicit_asset_col <- any(c(
    has_plot_file_col,
    has_locuscompare_plot_file_col,
    has_locuscompare_json_file_col
  ))
  plot_assets <- if (!is.null(asset_index) && !is.null(asset_index$plots)) {
    vapply(asset_index$plots, `[[`, "", "rel_path")
  } else {
    character(0)
  }
  json_assets <- if (!is.null(asset_index) && !is.null(asset_index$plot_data)) {
    vapply(asset_index$plot_data, `[[`, "", "rel_path")
  } else {
    character(0)
  }

  for (i in seq_along(rows)) {
    row <- rows[[i]]
    row_pp4 <- suppressWarnings(as.numeric(report_value(row, "PP4", default = NA_real_)))
    row_pp4 <- if (length(row_pp4) == 0) NA_real_ else row_pp4[[1]]
    legacy_asset_lookup <- !has_any_explicit_asset_col && (is.na(row_pp4) || row_pp4 >= 0.7)

    raw_plot_file <- report_clean_path(report_value(row, "plot_file", default = ""))
    raw_locuscompare_plot_file <- report_clean_path(report_value(row, "locuscompare_plot_file", default = ""))
    raw_locuscompare_json_file <- report_clean_path(report_value(row, "locuscompare_json_file", default = ""))

    plot_file <- report_existing_rel_path(results_dir, raw_plot_file)
    locuscompare_plot_file <- report_existing_rel_path(results_dir, raw_locuscompare_plot_file)
    locuscompare_json_file <- report_existing_rel_path(results_dir, raw_locuscompare_json_file)

    if (!nzchar(plot_file) && (nzchar(raw_plot_file) || (!has_plot_file_col && legacy_asset_lookup))) {
      matched <- report_match_asset(row, plot_assets, "_coloc\\.(png|pdf)$")
      plot_file <- report_existing_rel_path(results_dir, matched)
    }
    if (!nzchar(locuscompare_plot_file) && (nzchar(raw_locuscompare_plot_file) || (!has_locuscompare_plot_file_col && legacy_asset_lookup))) {
      matched <- report_match_asset(row, plot_assets, "_locuscompare\\.(png|pdf)$")
      locuscompare_plot_file <- report_existing_rel_path(results_dir, matched)
    }
    if (!nzchar(locuscompare_json_file) && (nzchar(raw_locuscompare_json_file) || (!has_locuscompare_json_file_col && legacy_asset_lookup))) {
      matched <- report_match_asset(row, json_assets, "_locuscompare\\.json$")
      locuscompare_json_file <- report_existing_rel_path(results_dir, matched)
    }

    out[[i]] <- list(
      gwas_id = as.character(report_value(row, "GWAS_ID", default = "")),
      qtl_id = as.character(report_value(row, "QTL_ID", default = "")),
      locus = as.character(report_value(row, "Locus", default = "")),
      phenotype = as.character(report_value(row, "Gene", fallback = "Phenotype", default = "")),
      pp4 = as.numeric(report_value(row, "PP4", default = NA_real_)),
      n_snps = as.integer(report_value(row, "n_snps", default = NA_integer_)),
      plot_file = plot_file,
      locuscompare_plot_file = locuscompare_plot_file,
      locuscompare_json_file = locuscompare_json_file,
      source_file = as.character(report_value(row, "source_file", default = ""))
    )
  }

  out
}

report_expected_plot_rows <- function(results_dt) {
  if (is.null(results_dt) || nrow(results_dt) == 0) {
    return(logical())
  }

  expected <- rep(FALSE, nrow(results_dt))
  asset_cols <- intersect(
    c("plot_file", "locuscompare_plot_file", "locuscompare_json_file"),
    names(results_dt)
  )

  if (length(asset_cols) > 0) {
    for (col in asset_cols) {
      expected <- expected | nzchar(vapply(results_dt[[col]], report_clean_path, character(1)))
    }
    return(expected)
  }

  if ("PP4" %in% names(results_dt)) {
    pp4 <- suppressWarnings(as.numeric(results_dt$PP4))
    return(!is.na(pp4) & pp4 >= 0.7)
  }

  rep(TRUE, nrow(results_dt))
}

index_report_assets <- function(results_dir) {
  plot_dir <- file.path(results_dir, "plots")
  plot_files <- if (dir.exists(plot_dir)) list.files(plot_dir, full.names = FALSE) else character(0)
  plots <- lapply(plot_files, function(name) {
    list(
      name = name,
      rel_path = file.path("plots", name)
    )
  })

  plot_data_dir <- file.path(results_dir, "plot_data")
  plot_data_files <- if (dir.exists(plot_data_dir)) list.files(plot_data_dir, full.names = FALSE) else character(0)
  plot_data <- lapply(plot_data_files, function(name) {
    list(
      name = name,
      rel_path = file.path("plot_data", name)
    )
  })

  top_level_entries <- list.files(results_dir, full.names = FALSE, recursive = FALSE)
  top_level_files <- top_level_entries[!dir.exists(file.path(results_dir, top_level_entries))]
  download_candidates <- top_level_files[grepl("\\.(csv|tsv|txt|json)$", top_level_files, ignore.case = TRUE)]
  download_files <- sort(setdiff(download_candidates, c("heartbeat.json", "monitor_snapshot.json")))

  downloads <- lapply(download_files, function(name) {
    list(
      name = name,
      rel_path = name
    )
  })

  list(
    plots = plots,
    plot_data = plot_data,
    downloads = downloads
  )
}

build_report_web_payload <- function(results_dir, project_name = NULL) {
  results_dt <- read_report_results(results_dir)
  if (is.null(results_dt) || nrow(results_dt) == 0) {
    stop("No ABF result files found under ", file.path(results_dir, "abf"), call. = FALSE)
  }

  project_name_resolved <- if (is.null(project_name) || !nzchar(project_name)) {
    basename(normalizePath(results_dir, mustWork = FALSE))
  } else {
    project_name
  }

  assets <- index_report_assets(results_dir)
  rows <- normalize_report_results(results_dt, results_dir = results_dir, asset_index = assets)
  expected_plot <- report_expected_plot_rows(results_dt)
  missing_plot_n <- sum(vapply(rows, function(row) {
    !nzchar(row$plot_file) && !nzchar(row$locuscompare_plot_file)
  }, logical(1)) & expected_plot)
  warnings <- if (missing_plot_n > 0) {
    list(sprintf("No coloc/locuscompare plot asset matched %d result row(s).", missing_plot_n))
  } else {
    list()
  }

  list(
    meta = list(
      project_name = project_name_resolved,
      results_dir = normalizePath(results_dir, mustWork = FALSE),
      generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      report_version = "report-web-v1"
    ),
    summary = get_summary_stats(results_dt),
    results = rows,
    assets = assets,
    runtime = list(
      heartbeat = read_json_if_present(file.path(results_dir, "heartbeat.json")),
      monitor_snapshot = read_json_if_present(file.path(results_dir, "monitor_snapshot.json"))
    ),
    warnings = warnings
  )
}

write_report_web_payload <- function(results_dir, project_name = NULL) {
  output_dir <- file.path(results_dir, "report_web")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  payload <- build_report_web_payload(results_dir, project_name = project_name)
  out_file <- file.path(output_dir, "report-data.json")
  write_json(payload, out_file, auto_unbox = TRUE, pretty = TRUE, null = "null")

  out_file
}
