easycoloc_get_config_path <- function(envvar, default_path) {
  value <- Sys.getenv(envvar, unset = default_path)
  if (!nzchar(value)) {
    default_path
  } else {
    path.expand(value)
  }
}

easycoloc_resolve_project_root <- function(cfg_global, global_config_path) {
  env_root <- Sys.getenv("EASYCOLOC_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(path.expand(env_root), mustWork = FALSE))
  }
  if (!is.null(cfg_global$project_root) && nzchar(cfg_global$project_root)) {
    base_dir <- easycoloc_path_dirname(global_config_path)
    expanded_root <- path.expand(cfg_global$project_root)
    if (grepl("^(/|[A-Za-z]:[/\\\\])", expanded_root)) {
      return(normalizePath(expanded_root, mustWork = FALSE))
    }
    return(normalizePath(file.path(base_dir, expanded_root), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

easycoloc_config_paths <- function() {
  list(
    global = easycoloc_get_config_path("EASYCOLOC_GLOBAL_CONFIG", "config/global.yml"),
    gwas = easycoloc_get_config_path("EASYCOLOC_GWAS_CONFIG", "config/gwas.yml"),
    qtl = easycoloc_get_config_path("EASYCOLOC_QTL_CONFIG", "config/qtl.yml")
  )
}

easycoloc_path_dirname <- function(path) {
  normalizePath(dirname(path), mustWork = FALSE)
}

easycoloc_expand_env_vars <- function(path_value) {
  if (is.null(path_value) || length(path_value) == 0) {
    return(path_value)
  }
  if (!is.character(path_value)) {
    return(path_value)
  }
  expanded <- vapply(path_value, function(single_path) {
    if (is.na(single_path) || !nzchar(single_path)) {
      return(single_path)
    }
    out <- single_path
    matches <- gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}|\\$[A-Za-z_][A-Za-z0-9_]*", out, perl = TRUE)
    tokens <- unique(regmatches(out, matches)[[1]])
    if (length(tokens) == 0 || identical(tokens, character(0))) {
      return(out)
    }
    for (token in tokens) {
      var <- sub("^\\$\\{?([^}]+)\\}?$", "\\1", token, perl = TRUE)
      value <- Sys.getenv(var, unset = token)
      out <- gsub(token, value, out, fixed = TRUE)
    }
    out
  }, character(1))
  unname(expanded)
}

easycoloc_resolve_path <- function(path_value, base_dir) {
  if (is.null(path_value) || length(path_value) == 0) {
    return(path_value)
  }
  if (!is.character(path_value)) {
    return(path_value)
  }
  resolved <- vapply(path_value, function(single_path) {
    if (is.na(single_path) || !nzchar(single_path)) {
      return(single_path)
    }
    expanded <- path.expand(easycoloc_expand_env_vars(single_path))
    if (grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) {
      return(expanded)
    }
    normalizePath(file.path(base_dir, expanded), mustWork = FALSE)
  }, character(1))
  unname(resolved)
}

easycoloc_resolve_named_paths <- function(path_list, base_dir) {
  if (is.null(path_list) || length(path_list) == 0) {
    return(path_list)
  }
  if (!is.list(path_list)) {
    return(easycoloc_resolve_path(path_list, base_dir))
  }
  for (nm in names(path_list)) {
    path_list[[nm]] <- easycoloc_resolve_path(path_list[[nm]], base_dir)
  }
  path_list
}

easycoloc_resolve_global_paths <- function(cfg_global, base_dir) {
  path_fields <- c(
    "output_dir", "temp_dir", "plink_hg19", "plink_hg38", "hash_table_dir",
    "gene_anno", "ref_genome_hg19", "ref_genome_hg38",
    "1kg_af", "dbsnp_hg19", "dbsnp_hg38", "harmonize_dir"
  )
  for (field in path_fields) {
    if (!is.null(cfg_global[[field]])) {
      cfg_global[[field]] <- easycoloc_resolve_path(cfg_global[[field]], base_dir)
    }
  }
  if (!is.null(cfg_global$plink_keep)) {
    cfg_global$plink_keep <- easycoloc_resolve_named_paths(cfg_global$plink_keep, base_dir)
  }
  if (!is.null(cfg_global$recom)) {
    cfg_global$recom <- easycoloc_resolve_named_paths(cfg_global$recom, base_dir)
  }
  if (!is.null(cfg_global$harmonization_settings$liftover_chain)) {
    cfg_global$harmonization_settings$liftover_chain <- easycoloc_resolve_path(
      cfg_global$harmonization_settings$liftover_chain,
      base_dir
    )
  }
  if (!is.null(cfg_global$runtime$disk_ld_cache_dir) &&
    is.character(cfg_global$runtime$disk_ld_cache_dir) &&
    nzchar(cfg_global$runtime$disk_ld_cache_dir)) {
    cfg_global$runtime$disk_ld_cache_dir <- easycoloc_resolve_path(
      cfg_global$runtime$disk_ld_cache_dir,
      base_dir
    )
  }
  cfg_global
}

easycoloc_resolve_gwas_paths <- function(cfg_gwas, base_dir) {
  if (!is.null(cfg_gwas$datasets) && length(cfg_gwas$datasets) > 0) {
    for (idx in seq_along(cfg_gwas$datasets)) {
      if (!is.null(cfg_gwas$datasets[[idx]]$file)) {
        cfg_gwas$datasets[[idx]]$file <- easycoloc_resolve_path(
          cfg_gwas$datasets[[idx]]$file,
          base_dir
        )
      }
    }
  }
  cfg_gwas
}

easycoloc_selected_ids <- function(envvar) {
  raw_value <- Sys.getenv(envvar, unset = "")
  if (!nzchar(raw_value)) {
    return(character())
  }
  selected <- trimws(strsplit(raw_value, ",", fixed = TRUE)[[1]])
  unique(selected[nzchar(selected)])
}

easycoloc_filter_gwas_datasets <- function(cfg_gwas,
                                           selected_ids = easycoloc_selected_ids("EASYCOLOC_GWAS_ID")) {
  if (length(selected_ids) == 0L) {
    return(cfg_gwas)
  }
  datasets <- cfg_gwas$datasets
  available_ids <- if (is.null(datasets)) {
    character()
  } else {
    vapply(datasets, function(dataset) {
      if (is.null(dataset$id)) NA_character_ else as.character(dataset$id)
    }, character(1))
  }
  missing_ids <- setdiff(selected_ids, available_ids)
  if (length(missing_ids) > 0L) {
    stop(
      sprintf(
        "Unknown GWAS ID(s): %s. Available IDs: %s",
        paste(missing_ids, collapse = ", "),
        if (length(available_ids) == 0L) "<none>" else paste(available_ids, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  cfg_gwas$datasets <- datasets[available_ids %in% selected_ids]
  cfg_gwas
}

easycoloc_filter_qtl_meta <- function(qtl_meta,
                                      cfg_qtl,
                                      selected_ids = easycoloc_selected_ids("EASYCOLOC_QTL_ID")) {
  if (length(selected_ids) == 0L) {
    return(qtl_meta)
  }
  id_col <- cfg_qtl$qtl_info$columns$id
  if (is.null(id_col) || !nzchar(id_col) || !id_col %in% names(qtl_meta)) {
    stop("Cannot select QTL data: qtl_info.columns.id is absent from the QTL metadata table", call. = FALSE)
  }
  available_ids <- unique(as.character(qtl_meta[[id_col]]))
  available_ids <- available_ids[!is.na(available_ids) & nzchar(available_ids)]
  missing_ids <- setdiff(selected_ids, available_ids)
  if (length(missing_ids) > 0L) {
    stop(
      sprintf(
        "Unknown QTL ID(s): %s. Available IDs: %s",
        paste(missing_ids, collapse = ", "),
        if (length(available_ids) == 0L) "<none>" else paste(available_ids, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  qtl_meta[as.character(qtl_meta[[id_col]]) %in% selected_ids, , drop = FALSE]
}

easycoloc_validate_gwas_column_maps <- function(cfg_gwas, config_path) {
  datasets <- cfg_gwas$datasets
  if (is.null(datasets) || length(datasets) == 0L) {
    return(invisible(TRUE))
  }

  for (dataset in datasets) {
    column_names <- names(dataset$columns)
    if (!is.null(column_names) && "FALSE" %in% column_names) {
      dataset_id <- if (is.null(dataset$id)) "<unknown>" else dataset$id
      stop(
        sprintf(
          "Invalid GWAS column map for %s in %s: YAML parsed an unquoted 'n' key as FALSE; write \"n\": instead",
          dataset_id,
          config_path
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

easycoloc_resolve_qtl_paths <- function(cfg_qtl, base_dir) {
  if (!is.null(cfg_qtl$qtl_info$file)) {
    cfg_qtl$qtl_info$file <- easycoloc_resolve_path(cfg_qtl$qtl_info$file, base_dir)
  }
  cfg_qtl
}

easycoloc_resolve_qtl_meta_paths <- function(qtl_meta, qtl_summary_path, qtl_cfg) {
  if (is.null(qtl_meta) || nrow(qtl_meta) == 0) {
    return(qtl_meta)
  }
  base_dir <- easycoloc_path_dirname(qtl_summary_path)
  path_cols <- c(
    qtl_cfg$qtl_info$columns$all_filename,
    qtl_cfg$qtl_info$columns$sig_filename
  )
  path_cols <- unique(path_cols[!is.na(path_cols) & nzchar(path_cols)])
  for (col in path_cols) {
    if (!is.null(col) && col %in% names(qtl_meta)) {
      qtl_meta[[col]] <- easycoloc_resolve_path(as.character(qtl_meta[[col]]), base_dir)
    }
  }
  qtl_meta
}

easycoloc_read_configs <- function() {
  paths <- easycoloc_config_paths()
  cfg_global <- yaml::read_yaml(paths$global)
  project_root <- easycoloc_resolve_project_root(cfg_global, paths$global)
  cfg_gwas <- yaml::read_yaml(paths$gwas)
  cfg_qtl <- yaml::read_yaml(paths$qtl)
  easycoloc_validate_gwas_column_maps(cfg_gwas, paths$gwas)
  cfg_gwas <- easycoloc_filter_gwas_datasets(cfg_gwas)
  output_override <- Sys.getenv("EASYCOLOC_OUTPUT_DIR", unset = "")
  if (nzchar(output_override)) {
    cfg_global$output_dir <- output_override
  }
  list(
    paths = paths,
    project_root = project_root,
    global = easycoloc_resolve_global_paths(cfg_global, project_root),
    gwas = easycoloc_resolve_gwas_paths(cfg_gwas, project_root),
    qtl = easycoloc_resolve_qtl_paths(cfg_qtl, project_root)
  )
}

easycoloc_try_read_configs <- function() {
  tryCatch(
    easycoloc_read_configs(),
    error = function(e) NULL
  )
}

easycoloc_resolve_output_dir_arg <- function(args = character(),
                                             arg_index = 1L,
                                             cfg_bundle = NULL,
                                             required = FALSE) {
  if (length(args) >= arg_index && nzchar(args[[arg_index]])) {
    return(normalizePath(path.expand(args[[arg_index]]), mustWork = FALSE))
  }

  if (is.null(cfg_bundle)) {
    cfg_bundle <- easycoloc_try_read_configs()
  }

  if (!is.null(cfg_bundle) &&
    !is.null(cfg_bundle$global$output_dir) &&
    nzchar(cfg_bundle$global$output_dir)) {
    return(normalizePath(cfg_bundle$global$output_dir, mustWork = FALSE))
  }

  if (isTRUE(required)) {
    stop("Output directory must be provided explicitly or resolvable from config/global.yml", call. = FALSE)
  }

  NA_character_
}

easycoloc_is_writable_dir <- function(path) {
  dir.exists(path) && file.access(path, 2) == 0
}

easycoloc_regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}
