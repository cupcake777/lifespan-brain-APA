# ------------------------------------------------------------------------------
# src/utils_plot.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(GenomicRanges)
  library(dplyr)
  library(data.table)
  library(glue)
  library(ggrepel)
  library(grid)
  library(viridis) # Perceptually uniform, colorblind-friendly color scales
})

resolve_col <- function(df, preferred, alternatives = NULL) {
  if (preferred %in% names(df)) {
    return(preferred)
  }
  if (!is.null(alternatives)) {
    for (alt in alternatives) {
      if (alt %in% names(df)) {
        return(alt)
      }
    }
  }
  all_cols <- names(df)
  match_idx <- which(tolower(all_cols) == tolower(preferred))
  if (length(match_idx) > 0) {
    return(all_cols[match_idx[1]])
  }

  if (!is.null(alternatives)) {
    for (alt in alternatives) {
      match_idx <- which(tolower(all_cols) == tolower(alt))
      if (length(match_idx) > 0) {
        return(all_cols[match_idx[1]])
      }
    }
  }
  return(NULL)
}

safe_df_convert <- function(df) {
  if (is.null(df)) {
    return(NULL)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    return(NULL)
  }
  return(df)
}

safe_numeric <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) {
    return(default)
  }
  suppressWarnings(as.numeric(as.character(x)))
}

safe_character <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) {
    return(default)
  }
  as.character(x)
}

normalize_chrom_label <- function(chrom) {
  chrom_clean <- gsub("^chr", "", safe_character(chrom), ignore.case = TRUE)
  paste0("chr", chrom_clean)
}

format_plot_label <- function(label, fallback = "Gene") {
  label_chr <- safe_character(label, default = fallback)
  if (is.na(label_chr) || !nzchar(label_chr)) {
    return(fallback)
  }

  if (grepl("|", label_chr, fixed = TRUE)) {
    label_parts <- strsplit(label_chr, "|", fixed = TRUE)[[1]]
    if (length(label_parts) >= 2 && nzchar(label_parts[2])) {
      return(label_parts[2])
    }
  }

  label_chr
}

format_qtl_label <- function(qtl_type) {
  qtl_chr <- safe_character(qtl_type, default = "")
  if (is.na(qtl_chr) || !nzchar(qtl_chr)) {
    return("")
  }
  gsub("_", " ", qtl_chr, fixed = TRUE)
}

.ld_reference_cache <- new.env(parent = emptyenv())

get_ld_plink_timeout <- function(default = 120L) {
  timeout <- getOption("easycoloc.ld_plink_timeout", default)
  if (is.null(timeout) || length(timeout) == 0 || is.na(timeout)) {
    return(default)
  }
  timeout <- as.integer(timeout[1])
  if (!is.finite(timeout) || timeout <= 0L) {
    return(default)
  }
  timeout
}

run_ld_plink <- function(plink_bin, args, timeout = get_ld_plink_timeout(), label = "LD") {
  cmd_text <- paste(c(shQuote(plink_bin), shQuote(args)), collapse = " ")
  message(glue("[LD] Running PLINK {label} with timeout {timeout}s: {cmd_text}"))
  status <- tryCatch(
    system2(
      plink_bin,
      args = args,
      stdout = FALSE,
      stderr = FALSE,
      timeout = timeout
    ),
    warning = function(w) {
      message(glue("[LD] PLINK {label} warning: {conditionMessage(w)}"))
      124L
    },
    error = function(e) {
      message(glue("[LD] PLINK {label} failed to start/run: {conditionMessage(e)}"))
      1L
    }
  )
  status <- as.integer(status)
  if (!identical(status, 0L)) {
    message(glue("[LD] PLINK {label} exited non-zero or timed out (status={status})"))
  }
  status
}

load_ld_reference_mapping <- function(bfile, variants = NULL) {
  mapping_file <- paste0(bfile, "_to_rsid.tsv")
  variants <- unique(as.character(variants))
  cache_key <- normalizePath(mapping_file, mustWork = FALSE)

  if (length(variants) == 0 && exists(cache_key, envir = .ld_reference_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .ld_reference_cache, inherits = FALSE))
  }

  mapping_dt <- NULL
  if (file.exists(mapping_file)) {
    mapping_dt <- if (length(variants) > 0) {
      variant_file <- tempfile(pattern = "ld_variants_", fileext = ".txt")
      writeLines(variants, variant_file)
      on.exit(unlink(variant_file), add = TRUE)

      awk_cmd <- sprintf(
        "awk -F '\\t' 'NR==FNR {wanted[$1]=1; next} FNR==1 || (($1 in wanted) || ($2 in wanted))' %s %s",
        shQuote(variant_file),
        shQuote(mapping_file)
      )

      tryCatch(
        data.table::fread(cmd = awk_cmd, select = c("unique_id", "rsid")),
        error = function(e) {
          message(glue("[LD] Failed to subset mapping file {mapping_file}: {e$message}"))
          NULL
        }
      )
    } else {
      tryCatch(
        data.table::fread(mapping_file, select = c("unique_id", "rsid")),
        error = function(e) {
          message(glue("[LD] Failed to read mapping file {mapping_file}: {e$message}"))
          NULL
        }
      )
    }

    if (!is.null(mapping_dt) && nrow(mapping_dt) > 0) {
      mapping_dt <- unique(mapping_dt[!is.na(unique_id) & nzchar(unique_id)])
      mapping_dt[, unique_id := as.character(unique_id)]
      mapping_dt[, rsid := as.character(rsid)]
      data.table::setkey(mapping_dt, unique_id)
      data.table::setindex(mapping_dt, rsid)
    }
  }

  if (length(variants) == 0) {
    assign(cache_key, mapping_dt, envir = .ld_reference_cache)
  }
  mapping_dt
}

resolve_ld_variant_ids <- function(variants, bfile) {
  variants <- unique(as.character(variants))
  mapping_dt <- load_ld_reference_mapping(bfile, variants = variants)
  mapping_file <- paste0(bfile, "_to_rsid.tsv")

  if (is.null(mapping_dt) || nrow(mapping_dt) == 0) {
    return(list(
      query_ids = variants,
      unique_to_rsid = setNames(character(), character())
    ))
  }

  unique_matches <- mapping_dt[J(variants), nomatch = 0L, .(unique_id, rsid)]
  unique_hits <- unique(unique_matches$unique_id)

  unresolved_variants <- variants[!variants %chin% unique_hits]
  rsid_matches <- if (length(unresolved_variants) > 0L) {
    mapping_dt[J(unresolved_variants), on = "rsid", nomatch = 0L, .(unique_id, rsid)]
  } else {
    mapping_dt[0, .(unique_id, rsid)]
  }

  resolved_rsids <- unique(rsid_matches$rsid)
  unmatched_variants <- unresolved_variants[!unresolved_variants %chin% resolved_rsids]
  query_ids <- unique(c(unique_hits, rsid_matches$unique_id, unmatched_variants))

  unique_to_rsid_dt <- unique(rbind(
    unique_matches[, .(unique_id, rsid)],
    rsid_matches[, .(unique_id, rsid)],
    fill = TRUE
  ))

  resolved_count <- data.table::uniqueN(resolved_rsids)
  if (resolved_count > 0) {
    message(glue("[LD] Resolved {resolved_count}/{length(variants)} rsID(s) via {basename(mapping_file)}"))
  }

  list(
    query_ids = query_ids,
    unique_to_rsid = setNames(unique_to_rsid_dt$rsid, unique_to_rsid_dt$unique_id)
  )
}

map_ld_result_to_input_ids <- function(ld_result, unique_to_rsid) {
  if (is.null(ld_result) || nrow(ld_result) == 0 || length(unique_to_rsid) == 0) {
    return(ld_result)
  }

  ld_dt <- data.table::as.data.table(ld_result)
  for (col in c("SNP_A", "SNP_B")) {
    mapped <- unique_to_rsid[ld_dt[[col]]]
    ld_dt[[col]] <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, ld_dt[[col]])
  }

  ld_dt <- ld_dt[
    ,
    .(R = R[which.max(abs(R))]),
    by = .(SNP_A, SNP_B)
  ]

  as.data.frame(ld_dt, stringsAsFactors = FALSE)
}

ld_matrix_to_pairs <- function(ld_res) {
  if (is.null(ld_res) || is.null(ld_res$R)) {
    return(NULL)
  }
  ld_mat <- as.matrix(ld_res$R)
  if (nrow(ld_mat) == 0 || ncol(ld_mat) == 0) {
    return(NULL)
  }
  ld_result <- as.data.frame(as.table(ld_mat), stringsAsFactors = FALSE)
  colnames(ld_result) <- c("SNP_A", "SNP_B", "R")
  ld_result <- ld_result[ld_result$SNP_A != ld_result$SNP_B, ]
  if (nrow(ld_result) == 0) {
    return(NULL)
  }
  ld_result
}

normalize_plot_ld_pairs <- function(ld_df) {
  if (is.null(ld_df)) return(NULL)
  ld_result <- if (is.list(ld_df) && !is.data.frame(ld_df) && !is.null(ld_df$R)) {
    ld_matrix_to_pairs(ld_df)
  } else {
    as.data.frame(ld_df, stringsAsFactors = FALSE)
  }
  if (is.null(ld_result) || !all(c("SNP_A", "SNP_B", "R") %in% names(ld_result))) {
    return(NULL)
  }
  ld_result
}

parse_feature_region <- function(label) {
  label_chr <- safe_character(label, default = "")
  if (is.na(label_chr) || !nzchar(label_chr)) {
    return(NULL)
  }

  region_match <- regexec("(chr[[:alnum:]]+):([0-9,]+)-([0-9,]+)", label_chr)
  region_parts <- regmatches(label_chr, region_match)[[1]]
  if (length(region_parts) != 4) {
    return(NULL)
  }

  list(
    chrom = region_parts[2],
    start = suppressWarnings(as.numeric(gsub(",", "", region_parts[3]))),
    end = suppressWarnings(as.numeric(gsub(",", "", region_parts[4])))
  )
}

clamp_unit_interval <- function(x) {
  pmin(pmax(x, 0), 1)
}

estimate_internal_legend_box <- function(legend_nrow = 1, legend_labels = NULL,
                                         plot_width = 10, plot_height = 8) {
  label_chars <- if (!is.null(legend_labels) && length(legend_labels) > 0) {
    suppressWarnings(max(nchar(as.character(legend_labels)), na.rm = TRUE))
  } else {
    8
  }
  if (!is.finite(label_chars)) label_chars <- 8

  width_frac <- 0.18 + 0.008 * max(0, label_chars - 6) + 0.01 * max(0, legend_nrow - 5)
  height_frac <- 0.11 + 0.045 * max(1, legend_nrow)

  if (is.finite(plot_width) && plot_width > 0) {
    width_frac <- width_frac * (6.8 / plot_width)^0.12
  }
  if (is.finite(plot_height) && plot_height > 0) {
    height_frac <- height_frac * (5.2 / plot_height)^0.10
  }

  width_frac <- min(0.34, max(0.18, width_frac))
  height_frac <- min(0.42, max(0.18, height_frac))

  list(
    width = width_frac,
    height = height_frac,
    margin_x = min(0.04, max(0.02, width_frac * 0.16)),
    margin_y = min(0.04, max(0.02, height_frac * 0.12))
  )
}

choose_external_legend_layout <- function(legend_nrow = 1,
                                          plot_width = 10, plot_height = 8) {
  use_bottom <- (is.finite(plot_width) && plot_width < 6.2) ||
    (is.finite(plot_height) && plot_height < 4.8) ||
    legend_nrow >= 7

  if (use_bottom) {
    return(list(
      position = "bottom",
      justification = "center",
      direction = "horizontal",
      guide_ncol = min(3, max(2, ceiling(legend_nrow / 2))),
      guide_nrow = max(1, ceiling(legend_nrow / min(3, max(2, ceiling(legend_nrow / 2))))),
      box_margin = margin(t = 3, r = 0, b = 0, l = 0),
      plot_margin = margin(t = 6, r = 6, b = 4, l = 6, unit = "pt")
    ))
  }

  list(
    position = "right",
    justification = "center",
    direction = "vertical",
    guide_ncol = 1,
    guide_nrow = legend_nrow,
    box_margin = margin(t = 0, r = 0, b = 0, l = 4),
    plot_margin = margin(t = 6, r = 2, b = 2, l = 6, unit = "pt")
  )
}

choose_internal_legend_position <- function(points_df, label_df = NULL,
                                            x_col = "x", y_col = "y",
                                            xlim = NULL, ylim = NULL,
                                            legend_nrow = 1, legend_labels = NULL,
                                            plot_width = 10, plot_height = 8) {
  default_layout <- list(
    corner = "top_left",
    position = c(0.02, 0.98),
    justification = c(0, 1),
    score = NA_real_
  )

  points_df <- safe_df_convert(points_df)
  if (is.null(points_df) || !all(c(x_col, y_col) %in% names(points_df))) {
    return(default_layout)
  }

  x_vals <- safe_numeric(points_df[[x_col]])
  y_vals <- safe_numeric(points_df[[y_col]])
  keep <- !is.na(x_vals) & !is.na(y_vals)
  if (!any(keep)) {
    return(default_layout)
  }

  x_vals <- x_vals[keep]
  y_vals <- y_vals[keep]
  point_df <- points_df[keep, , drop = FALSE]

  x_limits <- if (!is.null(xlim) && length(xlim) == 2 && all(is.finite(xlim))) {
    as.numeric(xlim)
  } else {
    range(x_vals, na.rm = TRUE)
  }
  y_limits <- if (!is.null(ylim) && length(ylim) == 2 && all(is.finite(ylim))) {
    as.numeric(ylim)
  } else {
    range(y_vals, na.rm = TRUE)
  }

  if (!all(is.finite(x_limits)) || diff(x_limits) <= 0 || !all(is.finite(y_limits)) || diff(y_limits) <= 0) {
    return(default_layout)
  }

  x_scaled <- clamp_unit_interval((x_vals - x_limits[1]) / diff(x_limits))
  y_scaled <- clamp_unit_interval((y_vals - y_limits[1]) / diff(y_limits))

  point_weight <- 1 + 0.65 * y_scaled
  if ("is_lead" %in% names(point_df)) {
    point_weight <- point_weight + ifelse(!is.na(point_df$is_lead) & point_df$is_lead, 1.8, 0)
  }
  if ("is_credible" %in% names(point_df)) {
    point_weight <- point_weight + ifelse(!is.na(point_df$is_credible) & point_df$is_credible, 0.9, 0)
  }

  legend_box <- estimate_internal_legend_box(
    legend_nrow = legend_nrow,
    legend_labels = legend_labels,
    plot_width = plot_width,
    plot_height = plot_height
  )

  candidates <- list(
    list(corner = "top_left", position = c(legend_box$margin_x, 1 - legend_box$margin_y), justification = c(0, 1)),
    list(corner = "top_right", position = c(1 - legend_box$margin_x, 1 - legend_box$margin_y), justification = c(1, 1)),
    list(corner = "bottom_right", position = c(1 - legend_box$margin_x, legend_box$margin_y), justification = c(1, 0)),
    list(corner = "bottom_left", position = c(legend_box$margin_x, legend_box$margin_y), justification = c(0, 0))
  )

  label_scaled <- NULL
  if (!is.null(label_df) && is.data.frame(label_df) && all(c("x", "y") %in% names(label_df))) {
    label_x <- safe_numeric(label_df$x)
    label_y <- safe_numeric(label_df$y)
    label_keep <- !is.na(label_x) & !is.na(label_y)
    if (any(label_keep)) {
      label_scaled <- data.frame(
        x = clamp_unit_interval((label_x[label_keep] - x_limits[1]) / diff(x_limits)),
        y = clamp_unit_interval((label_y[label_keep] - y_limits[1]) / diff(y_limits))
      )
    }
  }

  score_candidate <- function(candidate) {
    pos_x <- candidate$position[1]
    pos_y <- candidate$position[2]
    just_x <- candidate$justification[1]
    just_y <- candidate$justification[2]

    if (just_x == 0) {
      x_lo <- pos_x
      x_hi <- pos_x + legend_box$width
    } else {
      x_lo <- pos_x - legend_box$width
      x_hi <- pos_x
    }

    if (just_y == 0) {
      y_lo <- pos_y
      y_hi <- pos_y + legend_box$height
    } else {
      y_lo <- pos_y - legend_box$height
      y_hi <- pos_y
    }

    x_lo <- clamp_unit_interval(x_lo)
    x_hi <- clamp_unit_interval(x_hi)
    y_lo <- clamp_unit_interval(y_lo)
    y_hi <- clamp_unit_interval(y_hi)

    point_hit <- x_scaled >= x_lo & x_scaled <= x_hi & y_scaled >= y_lo & y_scaled <= y_hi
    point_score <- sum(point_weight[point_hit], na.rm = TRUE)

    label_score <- 0
    if (!is.null(label_scaled) && nrow(label_scaled) > 0) {
      label_hit <- label_scaled$x >= x_lo & label_scaled$x <= x_hi &
        label_scaled$y >= y_lo & label_scaled$y <= y_hi
      label_score <- 4 * sum(label_hit, na.rm = TRUE)
    }

    point_score + label_score
  }

  scores <- vapply(candidates, score_candidate, numeric(1))
  best_idx <- which.min(scores)[1]
  best_candidate <- candidates[[best_idx]]
  best_candidate$score <- scores[[best_idx]]
  best_candidate
}

resolve_plot_window <- function(leadSNP_DF, pos_col, chr_col,
                                lead_snp_val = NULL,
                                plot_window_bp = 200000,
                                phenotype_info = NULL) {
  positions <- safe_numeric(leadSNP_DF[[pos_col]])
  positions <- positions[!is.na(positions) & positions > 0]
  if (length(positions) == 0) {
    return(NULL)
  }

  if (!is.null(lead_snp_val)) {
    snp_col <- resolve_col(leadSNP_DF, "rsid", c("snp", "SNP", "marker"))
    if (!is.null(snp_col)) {
      lead_row <- leadSNP_DF[leadSNP_DF[[snp_col]] == lead_snp_val, ]
      if (nrow(lead_row) > 0) {
        center_pos <- safe_numeric(lead_row[[pos_col]][1])
        min_pos <- center_pos - plot_window_bp
        max_pos <- center_pos + plot_window_bp
      } else {
        min_pos <- min(positions)
        max_pos <- max(positions)
      }
    } else {
      min_pos <- min(positions)
      max_pos <- max(positions)
    }
  } else {
    min_pos <- min(positions)
    max_pos <- max(positions)
  }

  chrom_num <- safe_character(leadSNP_DF[[chr_col]])[1]
  chrom_label <- normalize_chrom_label(chrom_num)
  base_min_pos <- min_pos
  base_max_pos <- max_pos

  feature_region <- parse_feature_region(phenotype_info)
  if (!is.null(feature_region) &&
    identical(normalize_chrom_label(feature_region$chrom), chrom_label) &&
    !is.na(feature_region$start) && !is.na(feature_region$end)) {
    min_pos <- min(min_pos, feature_region$start)
    max_pos <- max(max_pos, feature_region$end)
  }

  list(
    min_pos = min_pos,
    max_pos = max_pos,
    base_min_pos = base_min_pos,
    base_max_pos = base_max_pos,
    chrom_num = chrom_num,
    chrom_label = chrom_label,
    feature_region = feature_region,
    expanded = (min_pos != base_min_pos || max_pos != base_max_pos)
  )
}

placeholder_plot <- function(label, detail = NULL) {
  label_text <- if (!is.null(detail) && nzchar(detail)) {
    paste(label, detail, sep = "\n")
  } else {
    label
  }

  ggplot(data.frame(x = 0.5, y = 0.5, label = label_text), aes(x = x, y = y)) +
    geom_text(aes(label = label), size = 4, lineheight = 1.1) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

ld_extract <- function(variants, bfile, plink_bin) {
  variants <- unique(as.character(variants))
  if (length(variants) == 0) {
    return(NULL)
  }

  fn <- tempfile()
  on.exit(unlink(paste0(fn, "*")), add = TRUE)

  tryCatch(
    {
      resolved_ids <- resolve_ld_variant_ids(variants, bfile)
      query_ids <- resolved_ids$query_ids
      if (length(query_ids) == 0) {
        message("[LD] No variants remained after ID resolution")
        return(NULL)
      }

      data.table::fwrite(list(query_ids), fn, col.names = FALSE)

      if (!file.exists(paste0(bfile, ".bed"))) {
        message("[LD] PLINK bfile not found: ", bfile)
        return(NULL)
      }

      plink_args <- c(
        "--bfile", bfile,
        "--extract", fn,
        "--allow-extra-chr",
        "--r", "square",
        "--make-just-bim",
        "--out", fn
      )
      sys_out <- run_ld_plink(plink_bin, plink_args, label = glue("plot LD ({length(query_ids)} variants)"))
      if (!identical(sys_out, 0L)) {
        return(NULL)
      }

      ld_file <- paste0(fn, ".ld")
      bim_file <- paste0(fn, ".bim")

      if (!file.exists(ld_file)) {
        message("[LD] LD file not created. Check PLINK output.")
        return(NULL)
      }

      ld_mat <- as.matrix(data.table::fread(ld_file))
      if (!file.exists(bim_file)) {
        message("[LD] BIM file not found")
        return(NULL)
      }
      bim <- data.table::fread(bim_file, col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))

      colnames(ld_mat) <- bim$SNP
      rownames(ld_mat) <- bim$SNP

      ld_result <- as.data.frame(as.table(ld_mat))
      colnames(ld_result) <- c("SNP_A", "SNP_B", "R")
      ld_result <- ld_result[ld_result$SNP_A != ld_result$SNP_B, ]
      ld_result <- map_ld_result_to_input_ids(ld_result, resolved_ids$unique_to_rsid)

      n_snps_in_bim <- nrow(bim)
      message(glue("[LD] Extracted LD for {length(unique(c(ld_result$SNP_A, ld_result$SNP_B)))}/{length(variants)} SNPs (PLINK queried {length(query_ids)} IDs; reference has {n_snps_in_bim} SNPs in region)"))

      res_inv <- data.frame(SNP_A = ld_result$SNP_B, SNP_B = ld_result$SNP_A, R = ld_result$R)

      return(dplyr::bind_rows(ld_result, res_inv) %>% dplyr::distinct())
    },
    error = function(e) {
      message(glue("[LD] Extraction failed: {e$message}"))
      return(NULL)
    }
  )
}

LD_plot <- function(df, ld_df = NULL, lead_snps = NULL,
                    bfile = NULL, plink_bin = "plink",
                    plot_title = NULL, plot_subtitle = NULL,
                    region_recomb = NULL, xlim = NULL,
                    show_lead_line = FALSE,
                    colocalized_snp = NULL,
                    credible_set = NULL,
                    gwas_display_threshold = NULL,
                    plot_width = 10, plot_height = 8) {
  df <- safe_df_convert(df)
  if (is.null(df)) {
    return(placeholder_plot("No Input Data"))
  }

  col_snp <- resolve_col(df, "rsid", c("SNPID.gwas", "snp", "SNP", "SNPID", "marker"))
  col_chr <- resolve_col(df, "chromosome", c("CHR", "chr", "CHR.qtl", "CHR.gwas"))
  col_pos <- resolve_col(df, "position", c("POS", "pos", "POS.qtl", "POS.gwas", "BP"))
  col_p <- resolve_col(df, "p_value", c("P", "pval", "P.qtl", "P.gwas"))

  if (is.null(col_snp) || is.null(col_chr) || is.null(col_pos) || is.null(col_p)) {
    return(placeholder_plot("Missing Columns"))
  }

  plot_df <- data.frame(
    rsid = safe_character(df[[col_snp]]),
    chromosome = safe_character(df[[col_chr]]),
    position = safe_numeric(df[[col_pos]]),
    p_value = safe_numeric(df[[col_p]]),
    stringsAsFactors = FALSE
  )
  col_gwas_p <- resolve_col(df, "P.gwas", c("p_gwas", "gwas_p", "pval_gwas"))
  if (!is.null(col_gwas_p)) {
    plot_df$gwas_p <- safe_numeric(df[[col_gwas_p]])
  } else {
    plot_df$gwas_p <- NA_real_
  }
  rsid_sample <- head(plot_df$rsid[!is.na(plot_df$rsid)], 3)
  is_rsid_format <- all(grepl("^rs", rsid_sample))
  message(glue("[LD_plot] SNP column: {col_snp}, Sample: {paste(rsid_sample, collapse=', ')}, Is rsID format: {is_rsid_format}"))

  plot_df <- plot_df[
    !is.na(plot_df$rsid) & plot_df$rsid != "" &
      !is.na(plot_df$position) &
      !is.na(plot_df$p_value) & plot_df$p_value > 0,
  ]

  if (nrow(plot_df) == 0) {
    return(placeholder_plot("No Valid SNPs"))
  }

  if (!is.null(xlim) && length(xlim) == 2) {
    plot_df <- plot_df[plot_df$position >= xlim[1] & plot_df$position <= xlim[2], ]
    if (nrow(plot_df) == 0) {
      return(placeholder_plot("No SNPs in Window"))
    }
  }

  message(glue("[LD_plot] Plotting {nrow(plot_df)} SNPs"))
  if (!is.null(lead_snps) && length(lead_snps) > 0) {
    lead_snp <- lead_snps[1]
    if (!lead_snp %in% plot_df$rsid) {
      lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
    }
  } else {
    lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
  }
  lead_pos <- plot_df$position[plot_df$rsid == lead_snp][1]

  message(glue("[LD_plot] Lead SNP: {lead_snp}"))
  message(glue("[LD_plot] Lead SNP in data: {lead_snp %in% plot_df$rsid}"))
  if (!(lead_snp %in% plot_df$rsid)) {
    close_matches <- grep(lead_snp, plot_df$rsid, value = TRUE)
    message(glue("[LD] Lead SNP '{lead_snp}' not found. Close matches: {paste(head(close_matches, 5), collapse=', ')}"))
    trimmed_match <- grep(gsub("^\\s+|\\s+$", "", lead_snp), plot_df$rsid, value = TRUE)
    if (length(trimmed_match) > 0) {
      message(glue("[LD] Found after trimming whitespace"))
    }
  }
  sample_snps <- head(unique(plot_df$rsid), 5)
  message(glue("[LD_plot] Sample SNPs in data: {paste(sample_snps, collapse=', ')}"))

  plot_df$r2 <- NA_real_
  if (!is.null(ld_df)) {
    ld_result <- normalize_plot_ld_pairs(ld_df)
    if (!is.null(ld_result)) {
      message(glue("[LD_plot] Using supplied LD table/matrix with {nrow(ld_result)} pairs"))
    } else {
      message("[LD_plot] Supplied ld_df is not a usable LD table/matrix")
      ld_result <- NULL
    }
  } else if (!is.null(bfile) && file.exists(paste0(bfile, ".bed"))) {
    snps_for_ld <- unique(c(lead_snp, plot_df$rsid))
    message(glue("[LD] Extracting LD for {length(snps_for_ld)} SNPs including lead SNP..."))

    ld_result <- ld_extract(snps_for_ld, bfile, plink_bin)
  } else {
    message(glue("[LD_plot] PLINK bfile not found: {bfile}"))
    ld_result <- NULL
  }

  if (!is.null(ld_result)) {
      lead_in_ld <- any(ld_result$SNP_A == lead_snp | ld_result$SNP_B == lead_snp)
      message(glue("[LD] Lead SNP in LD results: {lead_in_ld}"))

      if (lead_in_ld) {
        ld_to_lead <- ld_result[ld_result$SNP_A == lead_snp, ]
        message(glue("[LD_plot] LD pairs with lead SNP: {nrow(ld_to_lead)}"))
      } else {
        message(glue("[LD] WARNING: Lead SNP {lead_snp} not in PLINK reference!"))
        snps_in_ld <- unique(c(ld_result$SNP_A, ld_result$SNP_B))
        plot_df$in_ld <- plot_df$rsid %in% snps_in_ld

        if (any(plot_df$in_ld)) {
          snps_in_ld_df <- plot_df[plot_df$in_ld, ]
          if (nrow(snps_in_ld_df) > 0) {
            alt_lead <- snps_in_ld_df$rsid[which.min(snps_in_ld_df$p_value)]
            message(glue("[LD] Using {alt_lead} as alternative lead SNP for LD visualization"))
            ld_to_lead <- ld_result[ld_result$SNP_A == alt_lead, ]
            if (nrow(ld_to_lead) > 0) {
              ld_to_lead$r2_calc <- abs(ld_to_lead$R)^2
              plot_df <- merge(
                plot_df,
                ld_to_lead[, c("SNP_B", "r2_calc")],
                by.x = "rsid",
                by.y = "SNP_B",
                all.x = TRUE
              )

              plot_df$r2 <- ifelse(!is.na(plot_df$r2_calc), plot_df$r2_calc, plot_df$r2)
              plot_df$r2_calc <- NULL
              plot_df$r2[plot_df$rsid == alt_lead] <- 1.0
              lead_snp <- alt_lead
            }
          }
        }
        ld_to_lead <- data.frame()
      }

      if (nrow(ld_to_lead) > 0) {
        ld_to_lead$r2_calc <- abs(ld_to_lead$R)^2

        plot_df <- merge(
          plot_df,
          ld_to_lead[, c("SNP_B", "r2_calc")],
          by.x = "rsid",
          by.y = "SNP_B",
          all.x = TRUE
        )

        plot_df$r2 <- ifelse(!is.na(plot_df$r2_calc), plot_df$r2_calc, plot_df$r2)
        plot_df$r2_calc <- NULL
      }
  } else {
    message("[LD_plot] LD extraction returned NULL - check PLINK bfile and SNP IDs")
  }
  plot_df$r2[plot_df$rsid == lead_snp] <- 1.0
  plot_df$r2[is.na(plot_df$r2)] <- 0.05

  high_ld_count <- sum(plot_df$r2 > 0.2 & plot_df$rsid != lead_snp)
  message(glue("[LD_plot] Found {high_ld_count} proxy SNPs with r² > 0.2"))

  plot_df$is_lead <- plot_df$rsid == lead_snp
  # Mark credible set SNPs (excluding lead SNP — it has its own styling)
  cs_snps <- character(0)
  cs_pp <- numeric(0)
  if (!is.null(credible_set) && is.data.frame(credible_set) && nrow(credible_set) > 0) {
    cs_snps <- as.character(credible_set$snp)
    cs_pp <- credible_set$SNP.PP.H4
  }

  # Create is_credible column based on credible_set
  plot_df$is_credible <- FALSE
  if (length(cs_snps) > 0 && nrow(plot_df) > 0) {
    # Match SNPs, ensuring proper handling of NA values
    rsid_clean <- as.character(plot_df$rsid)
    cs_snps_clean <- as.character(cs_snps)
    credible_match <- rsid_clean %in% cs_snps_clean
    # Handle NA rsids by setting them to FALSE
    credible_match[is.na(credible_match)] <- FALSE
    # Ensure length matches
    if (length(credible_match) == nrow(plot_df)) {
      plot_df$is_credible <- credible_match
    } else {
      message("[LD_plot] Warning: Length mismatch in credible_set assignment")
      plot_df$is_credible <- rep(FALSE, nrow(plot_df))
    }
  }

  if (!is.null(gwas_display_threshold) && !all(is.na(plot_df$gwas_p))) {
    keep_mask <- (!is.na(plot_df$gwas_p) & plot_df$gwas_p <= gwas_display_threshold) |
      plot_df$is_lead | plot_df$is_credible
    kept_n <- sum(keep_mask, na.rm = TRUE)
    if (kept_n > 0 && kept_n < nrow(plot_df)) {
      message(glue("[LD_plot] Filtering displayed SNPs by GWAS threshold {signif(gwas_display_threshold, 3)}: {kept_n}/{nrow(plot_df)} retained"))
      plot_df <- plot_df[keep_mask, , drop = FALSE]
    }
  }

  # ==========================================================================
  # LD Color Scale: Viridis (perceptually uniform, colorblind-friendly)
  # ==========================================================================
  # Replaces previous RdYlBu scheme which had poor yellow distinguishability
  # and red-green confusion issues for colorblind users (~8% of males).
  # viridis(n=5, option="D") produces: dark purple -> blue -> green -> yellow
  # Option "magma" (sequential) or "cividis" (optimized for colorblind) available.
  vir_colors <- viridis(5, option = "D")
  plot_df$color_code <- cut(
    plot_df$r2,
    breaks = c(-0.01, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = vir_colors,
    include.lowest = TRUE
  )
  plot_df$color_code <- as.character(plot_df$color_code)
  plot_df$color_code[plot_df$is_lead] <- "#7F3C8D"
  plot_df$color_code[plot_df$is_credible & !plot_df$is_lead] <- "#E67E22"

  maxlogP <- max(-log10(plot_df$p_value), na.rm = TRUE)
  if (is.infinite(maxlogP) || maxlogP < 1) maxlogP <- 10
  minlogP <- min(-log10(plot_df$p_value), na.rm = TRUE)
  if (is.infinite(minlogP) || is.na(minlogP)) minlogP <- 0
  # Y-axis lower bound: don't always start from 0.
  # For strong signals (maxlogP > 7), trim the non-informative bottom by starting
  # from -log10(0.1) = 1. SNPs with p > 0.1 are pure noise in a colocalization
  # context and trimming them makes the significant region more prominent.
  if (maxlogP > 5) {
    ymin_auto <- 1
  } else {
    ymin_auto <- 0
  }
  p <- ggplot()

  if (!is.null(region_recomb) && is.data.frame(region_recomb) && nrow(region_recomb) > 0) {
    r_pos_col <- resolve_col(region_recomb, "Position(bp)", c("position", "pos", "bp"))
    r_rate_col <- resolve_col(region_recomb, "Rate(cM/Mb)", c("rate", "recomb", "cM/Mb"))

    if (!is.null(r_pos_col) && !is.null(r_rate_col)) {
      region_recomb$pos_mapped <- safe_numeric(region_recomb[[r_pos_col]])
      region_recomb$rate_mapped <- safe_numeric(region_recomb[[r_rate_col]])

      region_recomb <- region_recomb[!is.na(region_recomb$pos_mapped) & !is.na(region_recomb$rate_mapped), ]

      if (!is.null(xlim) && length(xlim) == 2) {
        region_recomb <- region_recomb[
          region_recomb$pos_mapped >= xlim[1] &
            region_recomb$pos_mapped <= xlim[2],
        ]
      }

      if (nrow(region_recomb) > 0) {
        recomb_max <- max(region_recomb$rate_mapped, na.rm = TRUE)
        if (!is.infinite(recomb_max) && recomb_max > 0) {
          # Scale recombination rate to fit within the visible y-range
          y_range <- maxlogP * 0.85 - ymin_auto
          scale_factor <- y_range / recomb_max
          region_recomb$scaled_rate <- region_recomb$rate_mapped * scale_factor + ymin_auto

          p <- p +
            geom_ribbon(
              data = region_recomb,
              aes(x = pos_mapped / 1e6, ymin = ymin_auto, ymax = pmax(scaled_rate, ymin_auto)),
              fill = "#DCEAF7", alpha = 0.65, inherit.aes = FALSE
            ) +
            geom_line(
              data = region_recomb,
              aes(x = pos_mapped / 1e6, y = scaled_rate),
              color = "#2E6F95", linewidth = 0.75, alpha = 0.9, inherit.aes = FALSE
            )
        }
      }
    }
  }
  if (show_lead_line && !is.na(lead_pos)) {
    p <- p +
      geom_vline(
        xintercept = lead_pos / 1e6,
        linetype = "dotted",
        color = "#7F3C8D",
        linewidth = 0.8,
        alpha = 0.6
      )
  }

  # --- Plot points in 3 layers: background, credible set, lead SNP ---
  df_bg <- plot_df[!plot_df$is_lead & !plot_df$is_credible, ]
  df_cs <- plot_df[plot_df$is_credible & !plot_df$is_lead, ]
  df_lead <- plot_df[plot_df$is_lead, ]

  # Layer 1: Background SNPs (circles)
  if (nrow(df_bg) > 0) {
    p <- p + geom_point(
      data = df_bg,
      aes(x = position / 1e6, y = -log10(p_value), fill = color_code),
      shape = 21, color = "#FFFFFF",
      stroke = 0.18, size = 2.8, alpha = 0.66
    )
  }
  # Layer 2: Credible set SNPs (circles, orange, slightly larger)
  # Use aes(fill=color_code) so the fill value goes through scale_fill_identity legend
  if (nrow(df_cs) > 0) {
    p <- p + geom_point(
      data = df_cs,
      aes(x = position / 1e6, y = -log10(p_value), fill = color_code),
      shape = 21, color = "#D35400",
      stroke = 0.72, size = 3.3, alpha = 0.9
    )
  }
  # Layer 3: Lead SNP (diamond, purple)
  # Use aes(fill=color_code) so the fill value goes through scale_fill_identity legend
  if (nrow(df_lead) > 0) {
    p <- p + geom_point(
      data = df_lead,
      aes(x = position / 1e6, y = -log10(p_value), fill = color_code),
      shape = 23, color = "#5B2C6F",
      stroke = 0.82, size = 3.8, alpha = 0.96
    )
  }

  p <- p + geom_hline(
    yintercept = -log10(5e-8),
    linetype = "dashed", color = "#6B7280", linewidth = 0.6
  )

  # --- Label lead SNP + credible set SNPs using ggrepel ---
  label_df <- data.frame(
    x = numeric(0), y = numeric(0), label = character(0),
    color = character(0), is_lead = logical(0),
    stringsAsFactors = FALSE
  )

  # Add lead SNP label
  if (nrow(df_lead) > 0) {
    lead_y <- -log10(df_lead$p_value[1])
    label_df <- rbind(label_df, data.frame(
      x = df_lead$position[1] / 1e6, y = lead_y,
      label = lead_snp, color = "#7F3C8D", is_lead = TRUE,
      stringsAsFactors = FALSE
    ))
  }

  # Add credible set labels (top 5 by PP.H4)
  if (nrow(df_cs) > 0 && length(cs_snps) > 0) {
    max_cs_labels <- min(2, nrow(df_cs))
    cs_pp_for_sort <- cs_pp[match(df_cs$rsid, cs_snps)]
    cs_order <- order(-cs_pp_for_sort, na.last = TRUE)
    df_cs_top <- df_cs[cs_order[seq_len(max_cs_labels)], ]
    for (i in seq_len(nrow(df_cs_top))) {
      cs_rsid <- df_cs_top$rsid[i]
      pp_val <- cs_pp[match(cs_rsid, cs_snps)]
      pp_label <- cs_rsid
      label_df <- rbind(label_df, data.frame(
        x = df_cs_top$position[i] / 1e6,
        y = -log10(df_cs_top$p_value[i]),
        label = pp_label, color = "#D35400", is_lead = FALSE,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Use geom_text_repel for all labels with segments
  if (nrow(label_df) > 0) {
    p <- p + geom_text_repel(
      data = label_df,
      aes(x = x, y = y, label = label),
      color = label_df$color,
      size = ifelse(label_df$is_lead, 2.9, 2.45),
      fontface = "bold",
      direction = "both",
      nudge_y = (maxlogP - minlogP) * 0.035,
      segment.color = "#A0AAB4",
      segment.size = 0.22,
      segment.linetype = "dotted",
      min.segment.length = 0.1,
      box.padding = 0.42,
      point.padding = 0.32,
      max.overlaps = 8,
      force = 2.4,
      force_pull = 0.7,
      max.time = 3,
      max.iter = 30000,
      inherit.aes = FALSE
    )
  }

  # --- Build legend with lead SNP (diamond) + credible set (orange) + LD colors ---
  # Use the same viridis colors as defined above for data coloring
  ld_colors_all <- vir_colors
  ld_labels_all <- c("0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "< 0.2")
  # Only keep LD bins that actually have SNPs in the data
  present_ld <- ld_colors_all %in% unique(plot_df$color_code)
  ld_colors <- ld_colors_all[present_ld]
  ld_labels <- ld_labels_all[present_ld]

  # Build legend breaks and labels
  # Note: override.aes shape must be SCALAR (not vector) for scale_fill_identity
  legend_breaks <- c("#7F3C8D", "#E67E22", ld_colors)
  legend_labels <- c("Lead SNP", "95% CS", ld_labels)

  # Remove credible set if not present
  if (nrow(df_cs) == 0) {
    legend_breaks <- c("#7F3C8D", ld_colors)
    legend_labels <- c("Lead SNP", ld_labels)
  }

  legend_ncol <- 1
  legend_nrow <- length(legend_breaks)
  ylim <- c(ymin_auto, ceiling(maxlogP * 1.1))
  legend_layout <- choose_external_legend_layout(
    legend_nrow = legend_nrow,
    plot_width = plot_width,
    plot_height = plot_height
  )
  message(glue("[LD_plot] Legend placement selected: {legend_layout$position}"))

  p <- p +
    scale_fill_identity(
      name = "LD (r2)",
      labels = legend_labels,
      breaks = legend_breaks,
      drop = FALSE,
      guide = guide_legend(
        override.aes = list(shape = 21, size = 3.1, alpha = 1, stroke = 0.4),
        title.position = "top",
        title.hjust = 0,
        nrow = legend_layout$guide_nrow,
        ncol = legend_layout$guide_ncol,
        byrow = TRUE
      )
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = NULL,
      y = expression(-log[10] ~ italic(P)),
      caption = glue(
        "lead {lead_snp} | proxies {high_ld_count} | cs95 {nrow(df_cs)}"
      )
    ) +
    theme_classic(base_size = 9.2, base_family = "sans") +
    theme(
      legend.position = legend_layout$position,
      legend.justification = legend_layout$justification,
      legend.direction = legend_layout$direction,
      legend.title = element_text(size = 7.8, face = "bold"),
      legend.text = element_text(size = 7.1),
      legend.key.size = unit(0.5, "lines"),
      legend.key.height = unit(0.55, "lines"),
      legend.background = element_rect(fill = "#FFFFFFF2", color = "#D8E1EA", linewidth = 0.5),
      legend.margin = margin(t = 4, r = 4, b = 4, l = 4),
      legend.box.margin = legend_layout$box_margin,
      legend.box.spacing = unit(0, "pt"),
      axis.line = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks = element_line(color = "#243447", linewidth = 0.4),
      axis.ticks.length = unit(0.08, "cm"),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8.2, color = "#203040"),
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = 9.1, face = "bold", margin = margin(r = 4)),
      plot.title = element_text(size = 11.5, face = "bold", color = "#12233B", hjust = 0, margin = margin(b = 1)),
      plot.subtitle = element_text(size = 8.1, color = "#4A5A6A", hjust = 0, margin = margin(b = 4)),
      plot.caption = element_text(size = 7.1, color = "#5B6775", hjust = 0, margin = margin(t = 3)),
      plot.margin = legend_layout$plot_margin,
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#FBFCFD", color = NA),
      panel.border = element_rect(fill = NA, color = "#D8E1EA", linewidth = 0.45),
      panel.grid.major.y = element_line(color = "#E7EDF3", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (!is.null(xlim) && length(xlim) == 2) {
    p <- p + scale_x_continuous(limits = xlim / 1e6, expand = expansion(mult = c(0.01, 0.01)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)))
  }
  if (!is.null(region_recomb) && exists("scale_factor")) {
    p <- p + scale_y_continuous(
      name = expression(-log[10] ~ italic(P)),
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = sec_axis(
        ~ (. - ymin_auto) / scale_factor,
        name = "Recombination rate (cM/Mb)"
      )
    ) +
      theme(
        axis.title.y.right = element_text(size = 8.7, face = "bold", color = "#2E6F95", margin = margin(l = 4)),
        axis.text.y.right = element_text(size = 7.8, color = "#2E6F95")
      )
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  }

  p <- p + coord_cartesian(ylim = ylim)

  return(p)
}

genetrack <- function(chrom_str, xlim_bp, gtf_path = NULL, show_strand_legend = FALSE) {
  chrom_str <- normalize_chrom_label(chrom_str)

  if (is.null(xlim_bp) || length(xlim_bp) != 2) {
    xlim_bp <- c(0, 1e6)
  }

  BPStart <- xlim_bp[1]
  BPStop <- xlim_bp[2]

  gr <- GRanges(seqnames = chrom_str, ranges = IRanges(start = BPStart, end = BPStop))
  tx_df <- NULL
  if (!is.null(gtf_path) && file.exists(gtf_path)) {
    tryCatch(
      {
        subset_gtf <- rtracklayer::import(gtf_path, which = gr, feature.type = "exon")
        if (length(subset_gtf) > 0) {
          df_gtf <- as.data.frame(subset_gtf)

          if ("gene_name" %in% names(df_gtf)) {
            df_gtf$symbol <- df_gtf$gene_name
          } else if ("gene_id" %in% names(df_gtf)) {
            df_gtf$symbol <- df_gtf$gene_id
          } else {
            df_gtf$symbol <- "Gene"
          }

          if ("strand" %in% names(df_gtf)) {
            df_gtf$gene_strand <- as.character(df_gtf$strand)
          } else {
            df_gtf$gene_strand <- "+"
          }

          df_gtf$gene_strand[is.na(df_gtf$gene_strand) | df_gtf$gene_strand == "*"] <- "+"

          if ("gene_type" %in% names(df_gtf)) {
            df_gtf <- df_gtf[df_gtf$gene_type == "protein_coding", ]
          }

          if (nrow(df_gtf) > 0) {
            tx_summary <- df_gtf %>%
              dplyr::group_by(symbol) %>%
              dplyr::summarize(
                start = min(start, na.rm = TRUE),
                end = max(end, na.rm = TRUE),
                strand = dplyr::first(gene_strand),
                .groups = "drop"
              )
            tx_df <- as.data.frame(tx_summary)
          }
        }
      },
      error = function(e) {
        message(glue("[Gene Track] GTF error: {e$message}"))
      }
    )
  }
  if (is.null(tx_df) || nrow(tx_df) == 0) {
    tryCatch(
      {
        required_pkgs <- c(
          "TxDb.Hsapiens.UCSC.hg38.knownGene",
          "GenomicFeatures",
          "AnnotationDbi",
          "org.Hs.eg.db"
        )
        missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
        if (length(missing_pkgs) > 0) {
          stop(glue("optional gene-track packages unavailable: {paste(missing_pkgs, collapse=', ')}"))
        }
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        subset_tx <- GenomicFeatures::transcriptsByOverlaps(txdb, gr)

        if (length(subset_tx) > 0) {
          gene_ids <- subset_tx$gene_id
          gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]

          symbols <- NULL
          if (length(gene_ids) > 0) {
            symbols <- suppressMessages(
              AnnotationDbi::mapIds(
                org.Hs.eg.db::org.Hs.eg.db,
                keys = gene_ids,
                column = "SYMBOL",
                keytype = "ENTREZID",
                multiVals = "first"
              )
            )
          }

          tx_df <- data.frame(
            gene_id = subset_tx$gene_id,
            tx_name = subset_tx$tx_name,
            start = start(subset_tx),
            end = end(subset_tx),
            strand = as.character(strand(subset_tx)),
            stringsAsFactors = FALSE
          )

          if (!is.null(symbols) && length(symbols) > 0) {
            symbol_df <- data.frame(
              gene_id = names(symbols),
              symbol = as.character(symbols),
              stringsAsFactors = FALSE
            )
            tx_df <- merge(tx_df, symbol_df, by = "gene_id", all.x = TRUE)
            tx_df$symbol <- ifelse(is.na(tx_df$symbol), tx_df$tx_name, tx_df$symbol)
          } else {
            tx_df$symbol <- tx_df$tx_name
          }

          tx_df <- tx_df[!grepl("^(LOC|LINC|MIR|SNOR|RN7|RNU)", tx_df$symbol), ]

          if (nrow(tx_df) > 0) {
            tx_df <- tx_df %>%
              dplyr::group_by(symbol) %>%
              dplyr::summarize(
                start = min(start, na.rm = TRUE),
                end = max(end, na.rm = TRUE),
                strand = dplyr::first(strand),
                .groups = "drop"
              ) %>%
              as.data.frame()
          }
        }
      },
      error = function(e) {
        message(glue("[Gene Track] TxDb error: {e$message}"))
      }
    )
  }

  if (is.null(tx_df) || nrow(tx_df) == 0) {
    return(
      ggplot() +
        theme_void() +
        geom_text(aes(x = mean(xlim_bp) / 1e6, y = 0.5, label = "No genes in region"),
          size = 3, color = "grey50", fontface = "italic"
        ) +
        scale_x_continuous(
          limits = xlim_bp / 1e6,
          expand = c(0.01, 0.01),
          name = paste0("Position on ", chrom_str, " (Mb)")
        ) +
        theme(
          plot.margin = margin(t = 2, r = 5, b = 8, l = 8, unit = "pt"),
          axis.text.x = element_text(size = 8.2, color = "black"),
          axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 5)),
          axis.ticks.x = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length.x = unit(0.12, "cm"),
          axis.line.x = element_line(color = "black", linewidth = 0.6)
        )
    )
  }

  message(glue("[Gene Track] Drawing {nrow(tx_df)} genes"))

  tx_df <- tx_df[order(tx_df$start), ]

  tx_df$strand[is.na(tx_df$strand) | tx_df$strand == "*" | tx_df$strand == ""] <- "+"
  tx_df$strand <- ifelse(tx_df$strand %in% c("+", "-"), tx_df$strand, "+")
  tx_df$y <- 1
  if (nrow(tx_df) > 1) {
    overlap_margin <- (BPStop - BPStart) * 0.03
    for (i in 2:nrow(tx_df)) {
      # Check overlap against ALL previously placed genes, not just the previous one
      candidate_y <- 1
      placed <- TRUE
      for (level in 1:10) { # support up to 10 stacking levels
        conflict <- FALSE
        for (j in 1:(i - 1)) {
          if (tx_df$y[j] == level && tx_df$start[i] < (tx_df$end[j] + overlap_margin)) {
            conflict <- TRUE
            break
          }
        }
        if (!conflict) {
          candidate_y <- level
          break
        }
        candidate_y <- level + 1
      }
      tx_df$y[i] <- candidate_y
    }
  }
  max_y_level <- max(tx_df$y)

  tx_df$color <- ifelse(tx_df$strand == "+", "#D73027", "#4575B4")
  lane_df <- data.frame(
    ymin = seq(0.5, max_y_level - 0.5, by = 1),
    ymax = seq(1.5, max_y_level + 0.5, by = 1),
    fill = rep(c("#F8FAFC", "#F1F5F9"), length.out = max_y_level),
    stringsAsFactors = FALSE
  )
  p <- ggplot(tx_df)
  p <- p + geom_rect(
    data = lane_df,
    aes(xmin = xlim_bp[1] / 1e6, xmax = xlim_bp[2] / 1e6, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = lane_df$fill,
    alpha = 1
  )
  p <- p + geom_segment(
    aes(x = start / 1e6, xend = end / 1e6, y = y, yend = y, color = strand),
    linewidth = 1.8,
    lineend = "round"
  )
  arrow_length <- (BPStop - BPStart) * 0.015

  fwd_genes <- tx_df[tx_df$strand == "+", ]
  if (nrow(fwd_genes) > 0) {
    p <- p + geom_segment(
      data = fwd_genes,
      aes(x = (end - arrow_length) / 1e6, xend = end / 1e6, y = y, yend = y),
      arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
      color = "#D73027",
      linewidth = 0.95
    )
  }

  rev_genes <- tx_df[tx_df$strand == "-", ]
  if (nrow(rev_genes) > 0) {
    p <- p + geom_segment(
      data = rev_genes,
      aes(x = (start + arrow_length) / 1e6, xend = start / 1e6, y = y, yend = y),
      arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
      color = "#4575B4",
      linewidth = 0.95
    )
  }
  p <- p + geom_text_repel(
    aes(x = (start + end) / 2e6, y = y, label = symbol, color = strand),
    size = 2.2,
    fontface = "bold.italic",
    direction = "y",
    nudge_y = 0.18,
    segment.color = "#A0AAB4",
    segment.size = 0.2,
    segment.linetype = "dotted",
    min.segment.length = 0.2,
    box.padding = 0.12,
    point.padding = 0.1,
    max.overlaps = 18,
    force = 1.1,
    show.legend = FALSE
  )
  p <- p +
    scale_color_manual(
      values = c("+" = "#D73027", "-" = "#4575B4"),
      labels = c("+" = "Forward", "-" = "Reverse"),
      name = "Strand"
    ) +
    scale_x_continuous(
      limits = xlim_bp / 1e6,
      expand = expansion(mult = c(0.01, 0.01)),
      name = paste0("Genomic position on ", chrom_str, " (Mb)")
    ) +
    scale_y_continuous(limits = c(0.3, max_y_level + 1.2), expand = c(0, 0)) +
    theme_void() +
    theme(
      plot.margin = margin(t = 2, r = 6, b = 4, l = 6, unit = "pt"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "#D8E1EA", linewidth = 0.45),
      axis.text.x = element_text(size = 7.8, color = "#203040"),
      axis.title.x = element_text(size = 8.6, face = "bold", margin = margin(t = 4)),
      axis.ticks.x = element_line(color = "#243447", linewidth = 0.4),
      axis.ticks.length.x = unit(0.08, "cm"),
      axis.line.x = element_blank(),
      legend.position = if (show_strand_legend) "bottom" else "none",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7.2)
    )

  return(p)
}

plot_qtl_association <- function(qtl_all_chrom, qtl_all_pvalue, leadSNP_DF,
                                 ld_df = NULL, gtf_path = NULL, region_recomb = NULL,
                                 recomb_path = NULL,
                                 show_lead_line = FALSE,
                                 qtl_type = NULL, phenotype_info = NULL,
                                 significance_threshold = 5e-8,
                                 significance_label = NULL,
                                 title_phenotype_field = "gene",
                                 plot_width = 10, plot_height = 8,
                                 credible_set = NULL,
                                 plot_window_bp = 200000) {
  message("[plot_qtl_association] Starting...")

  leadSNP_DF <- safe_df_convert(leadSNP_DF)
  if (is.null(leadSNP_DF)) {
    return(NULL)
  }

  lead_snp_val <- if (exists("lead_SNP", envir = .GlobalEnv)) get("lead_SNP", envir = .GlobalEnv) else NULL
  gene_sym <- if (exists("geneSymbol", envir = .GlobalEnv)) get("geneSymbol", envir = .GlobalEnv) else "Gene"
  plink_bfile_val <- if (exists("plink_bfile", envir = .GlobalEnv)) get("plink_bfile", envir = .GlobalEnv) else NULL

  pos_col <- resolve_col(leadSNP_DF, "POS.qtl", c("POS", "pos", "BP", "POS.gwas", "position"))
  chr_col <- resolve_col(leadSNP_DF, qtl_all_chrom, c("CHR", "CHR.qtl", "chr", "chromosome"))

  if (is.null(pos_col) || is.null(chr_col)) {
    message("[ERROR] Cannot resolve position/chr columns")
    return(NULL)
  }

  window_info <- resolve_plot_window(
    leadSNP_DF = leadSNP_DF,
    pos_col = pos_col,
    chr_col = chr_col,
    lead_snp_val = lead_snp_val,
    plot_window_bp = plot_window_bp,
    phenotype_info = phenotype_info
  )
  if (is.null(window_info)) {
    message("[ERROR] Cannot resolve plotting window")
    return(NULL)
  }

  min_pos <- window_info$min_pos
  max_pos <- window_info$max_pos
  chrom_num <- window_info$chrom_num
  chrom_label <- window_info$chrom_label

  # Build title: format 'QTLtype_Phenotypegene' (e.g., 'stage3_apa_ARL3')
  display_label <- format_plot_label(gene_sym, fallback = format_plot_label(phenotype_info))

  qtl_label <- format_qtl_label(qtl_type)

  if (!is.null(qtl_type) && !is.null(phenotype_info)) {
    plot_title_str <- paste0(qtl_label, " / ", display_label)
  } else if (!is.null(qtl_type)) {
    plot_title_str <- paste0(qtl_label, " / ", display_label)
  } else {
    plot_title_str <- paste(lead_snp_val, "-", display_label, "QTL")
  }

  plot_subtitle_str <- paste0(
    chrom_label, ":",
    format(as.integer(min_pos), big.mark = ","),
    "-",
    format(as.integer(max_pos), big.mark = ",")
  )

  if (isTRUE(window_info$expanded)) {
    plot_subtitle_str <- paste0(
      chrom_label, ":",
      format(as.integer(min_pos), big.mark = ","),
      "-",
      format(as.integer(max_pos), big.mark = ",")
    )
    message(glue("[plot_qtl_association] Expanded plot window to cover phenotype transcript: {plot_subtitle_str}"))
  }

  recomb_needs_reload <- is.null(region_recomb) || !is.data.frame(region_recomb) || nrow(region_recomb) == 0
  if (!recomb_needs_reload && !is.null(recomb_path)) {
    recomb_pos_col <- resolve_col(region_recomb, "Position(bp)", c("position", "pos", "bp"))
    if (!is.null(recomb_pos_col)) {
      recomb_pos <- safe_numeric(region_recomb[[recomb_pos_col]])
      recomb_pos <- recomb_pos[!is.na(recomb_pos)]
      if (length(recomb_pos) > 0) {
        recomb_needs_reload <- min(recomb_pos) > min_pos || max(recomb_pos) < max_pos
      } else {
        recomb_needs_reload <- TRUE
      }
    } else {
      recomb_needs_reload <- TRUE
    }
  }

  if (recomb_needs_reload && !is.null(recomb_path)) {
    region_recomb <- tryCatch(
      {
        load_recomb_map(chrom_label, min_pos, max_pos, recomb_path)
      },
      error = function(e) {
        message(glue("[plot_qtl_association] Recombination reload failed: {e$message}"))
        region_recomb
      }
    )
  }

  p_assoc <- tryCatch(
    {
      LD_plot(
        df = leadSNP_DF,
        ld_df = ld_df,
        lead_snps = lead_snp_val,
        bfile = plink_bfile_val,
        plink_bin = "plink",
        plot_title = plot_title_str,
        plot_subtitle = plot_subtitle_str,
        region_recomb = region_recomb,
        xlim = c(min_pos, max_pos),
        show_lead_line = show_lead_line,
        colocalized_snp = lead_snp_val,
        credible_set = credible_set,
        gwas_display_threshold = significance_threshold,
        plot_width = plot_width,
        plot_height = plot_height
      )
    },
    error = function(e) {
      message(glue("[ERROR] LD plot failed: {e$message}"))
      placeholder_plot("LD Plot Error", e$message)
    }
  )

  p_track <- tryCatch(
    {
      genetrack(chrom_num, c(min_pos, max_pos), gtf_path = gtf_path)
    },
    error = function(e) {
      message(glue("[ERROR] Gene track failed: {e$message}"))
      placeholder_plot("Gene Track Error", e$message)
    }
  )

  # Combine LD plot and gene track
  tryCatch(
    {
      ggarrange(p_assoc, p_track,
        ncol = 1,
        heights = c(3.8, 0.95),
        align = "v"
      )
    },
    error = function(e) {
      message(glue("[ERROR] ggarrange failed: {e$message}"))
      # Return just the LD plot if arrangement fails
      message("[WARN] Returning LD plot only (no gene track)")
      p_assoc
    }
  )
}
