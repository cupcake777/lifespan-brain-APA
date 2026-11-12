suppressPackageStartupMessages({
  library(data.table)
  library(glue)
})

easycoloc_file_size_bytes <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(NA_real_)
  }
  info <- tryCatch(file.info(path), error = function(e) NULL)
  if (is.null(info) || nrow(info) == 0 || isTRUE(info$isdir[1])) {
    return(NA_real_)
  }
  as.numeric(info$size[1])
}

easycoloc_dir_size_bytes <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !dir.exists(path)) {
    return(NA_real_)
  }
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files)]
  if (length(files) == 0) {
    return(0)
  }
  sum(file.info(files)$size, na.rm = TRUE)
}

easycoloc_plink_prefix_size_bytes <- function(prefix_path) {
  if (is.null(prefix_path) || is.na(prefix_path) || !nzchar(prefix_path)) {
    return(NA_real_)
  }
  suffixes <- c(".bed", ".bim", ".fam")
  files <- paste0(prefix_path, suffixes)
  if (!all(file.exists(files))) {
    return(NA_real_)
  }
  sum(file.info(files)$size, na.rm = TRUE)
}

easycoloc_format_size <- function(size_bytes) {
  if (is.na(size_bytes)) {
    return(NA_character_)
  }
  if (size_bytes < 1024) {
    return(paste0(size_bytes, " B"))
  }
  units <- c("B", "KB", "MB", "GB", "TB")
  size <- size_bytes
  unit_idx <- 1L
  while (size >= 1024 && unit_idx < length(units)) {
    size <- size / 1024
    unit_idx <- unit_idx + 1L
  }
  sprintf("%.2f %s", size, units[unit_idx])
}

easycoloc_plink_prefix_status <- function(prefix_path) {
  if (is.null(prefix_path) || is.na(prefix_path) || !nzchar(prefix_path)) {
    return("missing_path")
  }
  suffixes <- c(".bed", ".bim", ".fam")
  files <- paste0(prefix_path, suffixes)
  missing <- suffixes[!file.exists(files)]
  if (length(missing) == 0) {
    return("present")
  }
  paste0("missing:", paste(missing, collapse = ","))
}

easycoloc_required_builds <- function(cfg_gwas) {
  if (is.null(cfg_gwas$datasets) || length(cfg_gwas$datasets) == 0) {
    return(character(0))
  }
  unique(vapply(
    cfg_gwas$datasets,
    function(ds) if (!is.null(ds$build)) as.character(ds$build) else "hg38",
    character(1)
  ))
}

easycoloc_collect_reference_requirements <- function(cfg_bundle, include_qtl_files = FALSE) {
  cfg_global <- cfg_bundle$global
  cfg_gwas <- cfg_bundle$gwas
  cfg_qtl <- cfg_bundle$qtl
  builds <- easycoloc_required_builds(cfg_gwas)
  requires_hg19 <- "hg19" %in% builds
  requires_hg38 <- length(builds) == 0 || "hg38" %in% builds || isTRUE(any(gsub("hg", "", builds) == "38"))
  qtl_build <- if (!is.null(cfg_qtl$qtl_info$build)) as.character(cfg_qtl$qtl_info$build) else "hg38"

  add_row <- function(rows, name, required, path, path_type, format, used_for, stage, note = "") {
    if (is.list(path) && !is.data.frame(path)) {
      path_names <- names(path)
      if (is.null(path_names)) path_names <- as.character(seq_along(path))
      for (idx in seq_along(path)) {
        rows <- add_row(
          rows = rows,
          name = paste0(name, "_", path_names[[idx]]),
          required = required,
          path = path[[idx]],
          path_type = path_type,
          format = format,
          used_for = used_for,
          stage = stage,
          note = note
        )
      }
      return(rows)
    }
    if (length(path) > 1) {
      path_names <- names(path)
      if (is.null(path_names)) path_names <- as.character(seq_along(path))
      for (idx in seq_along(path)) {
        rows <- add_row(
          rows = rows,
          name = paste0(name, "_", path_names[[idx]]),
          required = required,
          path = path[[idx]],
          path_type = path_type,
          format = format,
          used_for = used_for,
          stage = stage,
          note = note
        )
      }
      return(rows)
    }

    exists_flag <- FALSE
    status <- "missing"
    size_bytes <- NA_real_

    if (!is.null(path) && !is.na(path) && nzchar(path)) {
      if (identical(path_type, "plink_prefix")) {
        status <- easycoloc_plink_prefix_status(path)
        exists_flag <- identical(status, "present")
        size_bytes <- easycoloc_plink_prefix_size_bytes(path)
      } else if (identical(path_type, "directory")) {
        exists_flag <- dir.exists(path)
        status <- if (exists_flag) "present" else "missing"
        size_bytes <- easycoloc_dir_size_bytes(path)
      } else {
        exists_flag <- file.exists(path)
        status <- if (exists_flag) "present" else "missing"
        size_bytes <- easycoloc_file_size_bytes(path)
      }
    } else {
      status <- "not_configured"
    }

    rbind(
      rows,
      data.table(
        name = as.character(name),
        required = required,
        stage = as.character(stage),
        used_for = as.character(used_for),
        format = as.character(format),
        path_type = as.character(path_type),
        path = if (is.null(path)) NA_character_ else as.character(path),
        status = as.character(status),
        exists = exists_flag,
        size_bytes = size_bytes,
        size_human = easycoloc_format_size(size_bytes),
        note = as.character(note)
      ),
      fill = TRUE
    )
  }

  rows <- data.table()
  if (identical(qtl_build, "hg19")) {
    rows <- add_row(
      rows,
      name = "plink_hg19",
      required = TRUE,
      path = cfg_global$plink_hg19,
      path_type = "plink_prefix",
      format = ".bed/.bim/.fam prefix",
      used_for = "locus clumping + SuSiE LD matrix for hg19 QTL build",
      stage = "core"
    )
  } else {
    rows <- add_row(
      rows,
      name = "plink_hg38",
      required = TRUE,
      path = cfg_global$plink_hg38,
      path_type = "plink_prefix",
      format = ".bed/.bim/.fam prefix",
      used_for = "locus clumping + SuSiE LD matrix for hg38 QTL build",
      stage = "core"
    )
  }
  rows <- add_row(
    rows,
    name = "plink_keep",
    required = FALSE,
    path = cfg_global$plink_keep,
    path_type = "file",
    format = "PLINK --keep file, recommended two columns: FID IID",
    used_for = "population-specific LD subset",
    stage = "optional",
    note = "If missing, EasyColoc falls back to all samples in the PLINK reference. One-column IID files are accepted and converted at runtime."
  )
  rows <- add_row(
    rows,
    name = "qtl_summary",
    required = TRUE,
    path = cfg_qtl$qtl_info$file,
    path_type = "file",
    format = "CSV metadata table",
    used_for = "locate QTL tabix resources",
    stage = "core"
  )
  rows <- add_row(
    rows,
    name = "gene_anno",
    required = FALSE,
    path = cfg_global$gene_anno,
    path_type = "file",
    format = "GTF/GTF.GZ",
    used_for = "gene track plotting",
    stage = "plotting"
  )
  rows <- add_row(
    rows,
    name = "recombination_map",
    required = FALSE,
    path = cfg_global$recom,
    path_type = "file",
    format = "prefix for chr-specific txt/bed maps",
    used_for = "recombination ribbon in locus plots",
    stage = "plotting"
  )
  rows <- add_row(
    rows,
    name = "hash_table_dir",
    required = FALSE,
    path = cfg_global$hash_table_dir,
    path_type = "directory",
    format = "chr_*_hash_table.rds directory",
    used_for = "rsID to chr:pos rescue matching",
    stage = "matching",
    note = "Optional but useful when rsID direct overlap is poor."
  )
  rows <- add_row(
    rows,
    name = "harmonize_dir",
    required = FALSE,
    path = cfg_global$harmonize_dir,
    path_type = "directory",
    format = "working directory",
    used_for = "native harmonized GWAS caches",
    stage = "harmonization"
  )

  if (requires_hg19) {
    rows <- add_row(
      rows,
      name = "ref_genome_hg19",
      required = TRUE,
      path = cfg_global$ref_genome_hg19,
      path_type = "file",
      format = "FASTA/FASTA.GZ",
      used_for = "hg19 GWAS harmonization",
      stage = "harmonization"
    )
    rows <- add_row(
      rows,
      name = "dbsnp_hg19",
      required = TRUE,
      path = cfg_global$dbsnp_hg19,
      path_type = "file",
      format = "dbSNP VCF/BED.GZ",
      used_for = "hg19 rsID normalization",
      stage = "harmonization"
    )
    rows <- add_row(
      rows,
      name = "liftover_chain",
      required = FALSE,
      path = cfg_global$harmonization_settings$liftover_chain,
      path_type = "file",
      format = "chain.gz",
      used_for = "fallback hg19 to hg38 coordinate conversion",
      stage = "harmonization"
    )
  }

  if (requires_hg38) {
    rows <- add_row(
      rows,
      name = "ref_genome_hg38",
      required = TRUE,
      path = cfg_global$ref_genome_hg38,
      path_type = "file",
      format = "FASTA/FASTA.GZ",
      used_for = "hg38 GWAS harmonization",
      stage = "harmonization"
    )
    rows <- add_row(
      rows,
      name = "dbsnp_hg38",
      required = TRUE,
      path = cfg_global$dbsnp_hg38,
      path_type = "file",
      format = "dbSNP VCF/BED.GZ",
      used_for = "hg38 rsID normalization",
      stage = "harmonization"
    )
  }

  if (!is.null(cfg_gwas$datasets) && length(cfg_gwas$datasets) > 0) {
    pop_build_dt <- unique(data.table(
      build = vapply(cfg_gwas$datasets, function(ds) if (!is.null(ds$build)) as.character(ds$build) else "hg38", character(1)),
      pop = vapply(cfg_gwas$datasets, function(ds) if (!is.null(ds$pop)) as.character(ds$pop) else "EUR", character(1))
    ))
    for (idx in seq_len(nrow(pop_build_dt))) {
      ref_vcf <- file.path(
        cfg_global[["1kg_af"]],
        glue("1KG_hg{gsub('hg', '', pop_build_dt$build[idx])}_{pop_build_dt$pop[idx]}_AF.tsv.gz")
      )
      rows <- add_row(
        rows,
        name = glue("1kg_af_{pop_build_dt$build[idx]}_{pop_build_dt$pop[idx]}"),
        required = TRUE,
        path = ref_vcf,
        path_type = "file",
        format = "TSV.GZ allele frequency table",
        used_for = "native harmonization prior/reference AF",
        stage = "harmonization"
      )
    }
  }

  if (isTRUE(include_qtl_files) && file.exists(cfg_qtl$qtl_info$file)) {
    qtl_meta <- tryCatch(fread(cfg_qtl$qtl_info$file), error = function(e) NULL)
    if (!is.null(qtl_meta) && nrow(qtl_meta) > 0) {
      qtl_meta <- easycoloc_resolve_qtl_meta_paths(qtl_meta, cfg_qtl$qtl_info$file, cfg_qtl)
      id_col <- cfg_qtl$qtl_info$columns$id
      all_col <- cfg_qtl$qtl_info$columns$all_filename
      sig_col <- cfg_qtl$qtl_info$columns$sig_filename
      for (i in seq_len(nrow(qtl_meta))) {
        qtl_id <- as.character(qtl_meta[[id_col]][i])
        rows <- add_row(
          rows,
          name = glue("qtl_allpairs_{qtl_id}"),
          required = TRUE,
          path = as.character(qtl_meta[[all_col]][i]),
          path_type = "file",
          format = "bgzip tabix table",
          used_for = "QTL regional query",
          stage = "qtl_data"
        )
        if (!is.null(sig_col) && sig_col %in% names(qtl_meta)) {
          rows <- add_row(
            rows,
            name = glue("qtl_sigpairs_{qtl_id}"),
            required = FALSE,
            path = as.character(qtl_meta[[sig_col]][i]),
            path_type = "file",
            format = "bgzip tabix table",
            used_for = "fast phenotype prefilter",
            stage = "qtl_data",
            note = "Optional but strongly recommended for speed."
          )
        }
      }
    }
  }

  rows
}
