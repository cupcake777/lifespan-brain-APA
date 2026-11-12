# ------------------------------------------------------------------------------
# src/utils_format.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glue)
  library(jsonlite)
  library(rlang)
})

`%notin%` <- Negate(`%in%`)
if (file.exists("src/utils_hash.R")) {
  source("src/utils_hash.R")
} else {
  warning("Hash conversion script not exist!.")
}

.ld_cache <- new.env(parent = emptyenv())
.tabix_cache <- new.env(parent = emptyenv())

easycoloc_prepare_plink_keep_file <- function(keep_file,
                                              bfile = NULL,
                                              temp_dir = tempdir(),
                                              label = "keep",
                                              verbose = TRUE) {
  if (is.null(keep_file) || !nzchar(keep_file) || !file.exists(keep_file)) {
    return(NULL)
  }
  keep_dt <- tryCatch(
    data.table::fread(keep_file, header = FALSE, fill = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(keep_dt) || nrow(keep_dt) == 0) {
    stop(glue("[PLINK] Keep file is empty or unreadable: {keep_file}"), call. = FALSE)
  }
  if (ncol(keep_dt) >= 2) {
    return(keep_file)
  }

  iid <- unique(trimws(as.character(keep_dt[[1]])))
  iid <- iid[!is.na(iid) & nzchar(iid)]
  if (length(iid) == 0) {
    stop(glue("[PLINK] One-column keep file has no sample IDs: {keep_file}"), call. = FALSE)
  }

  keep_out <- data.table::data.table(FID = "0", IID = iid)
  fam_file <- if (!is.null(bfile) && nzchar(bfile)) paste0(bfile, ".fam") else NULL
  if (!is.null(fam_file) && file.exists(fam_file)) {
    fam_dt <- data.table::fread(
      fam_file,
      header = FALSE,
      select = c(1, 2),
      col.names = c("FID", "IID"),
      showProgress = FALSE
    )
    fam_dt[, IID := as.character(IID)]
    matched <- fam_dt[IID %in% iid, .(FID = as.character(FID), IID)]
    missing_iid <- setdiff(iid, matched$IID)
    if (length(missing_iid) > 0 && isTRUE(verbose)) {
      message(glue("[PLINK] Keep file {basename(keep_file)}: {length(missing_iid)} IDs not found in {basename(fam_file)} and will be omitted"))
    }
    keep_out <- unique(matched)
    if (nrow(keep_out) == 0) {
      stop(glue("[PLINK] No keep-file samples overlap PLINK .fam: {keep_file}"), call. = FALSE)
    }
  } else if (isTRUE(verbose)) {
    message(glue("[PLINK] Converting one-column keep file {basename(keep_file)} to FID/IID format with FID=0"))
  }

  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  safe_label <- gsub("[^A-Za-z0-9_.-]+", "_", label)
  out_file <- tempfile(pattern = paste0("easycoloc_", safe_label, "_keep_"), tmpdir = temp_dir, fileext = ".sample")
  data.table::fwrite(keep_out, out_file, sep = "\t", col.names = FALSE, quote = FALSE)
  if (isTRUE(verbose)) {
    message(glue("[PLINK] Converted one-column keep file for {label}: {basename(keep_file)} -> {basename(out_file)} ({nrow(keep_out)} samples)"))
  }
  out_file
}

ld_cache_group_key <- function(bfile, keep_file) {
  paste(
    normalizePath(bfile, mustWork = FALSE),
    if (is.null(keep_file) || !nzchar(keep_file)) "NO_KEEP" else normalizePath(keep_file, mustWork = FALSE),
    sep = "::"
  )
}

ld_cache_entry_key <- function(variants) {
  paste(sort(unique(as.character(variants))), collapse = "|")
}

ld_disk_cache_dir <- function() {
  cache_dir <- getOption("easycoloc.disk_ld_cache_dir", NULL)
  if (is.null(cache_dir) || length(cache_dir) == 0 || is.na(cache_dir[1]) || !nzchar(cache_dir[1])) {
    return(NULL)
  }
  normalizePath(path.expand(as.character(cache_dir[1])), mustWork = FALSE)
}

ld_reference_signature <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path)) {
    return(paste(path, "missing", sep = "::"))
  }
  info <- file.info(path)
  paste(path, info$size, as.numeric(info$mtime), sep = "::")
}

ld_disk_cache_key <- function(variants, bfile, keep_file = NULL) {
  variants <- sort(unique(as.character(variants)))
  ref_paths <- c(paste0(bfile, ".bed"), paste0(bfile, ".bim"), paste0(bfile, ".fam"))
  ref_sig <- paste(vapply(ref_paths, ld_reference_signature, character(1)), collapse = "||")
  keep_sig <- if (is.null(keep_file) || !nzchar(keep_file)) {
    "NO_KEEP"
  } else {
    ld_reference_signature(keep_file)
  }
  digest_src <- paste(
    normalizePath(bfile, mustWork = FALSE),
    keep_sig,
    ref_sig,
    paste(variants, collapse = "|"),
    sep = "\n"
  )
  key_file <- tempfile("easycoloc_ld_key_")
  on.exit(unlink(key_file), add = TRUE)
  writeLines(digest_src, key_file, useBytes = TRUE)
  unname(tools::md5sum(key_file))
}

ld_disk_cache_path <- function(variants, bfile, keep_file = NULL) {
  cache_dir <- ld_disk_cache_dir()
  if (is.null(cache_dir)) {
    return(NULL)
  }
  file.path(cache_dir, paste0(ld_disk_cache_key(variants, bfile, keep_file), ".rds"))
}

ld_disk_cache_read <- function(variants, bfile, keep_file = NULL) {
  cache_path <- ld_disk_cache_path(variants, bfile, keep_file)
  if (is.null(cache_path) || !file.exists(cache_path)) {
    return(NULL)
  }
  cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
  if (is.null(cached) || is.null(cached$R)) {
    return(NULL)
  }
  ld_cache_store(variants, bfile, keep_file = keep_file, ld_res = cached)
  cached
}

ld_disk_cache_write <- function(variants, bfile, keep_file = NULL, ld_res) {
  cache_path <- ld_disk_cache_path(variants, bfile, keep_file)
  if (is.null(cache_path) || is.null(ld_res) || is.null(ld_res$R)) {
    return(invisible(FALSE))
  }
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(cache_path, ".tmp.", Sys.getpid())
  ok <- tryCatch(
    {
      saveRDS(ld_res, tmp_path)
      file.rename(tmp_path, cache_path)
    },
    error = function(e) FALSE
  )
  if (file.exists(tmp_path)) {
    unlink(tmp_path)
  }
  invisible(isTRUE(ok))
}

ld_cache_lookup <- function(variants, bfile, keep_file = NULL) {
  req_snps <- sort(unique(as.character(variants)))
  if (length(req_snps) == 0) {
    return(NULL)
  }
  group_key <- ld_cache_group_key(bfile, keep_file)
  if (!exists(group_key, envir = .ld_cache, inherits = FALSE)) {
    return(NULL)
  }
  group_entries <- get(group_key, envir = .ld_cache, inherits = FALSE)
  exact_key <- ld_cache_entry_key(req_snps)
  if (!is.null(group_entries[[exact_key]])) {
    entry <- group_entries[[exact_key]]
    available <- req_snps[req_snps %in% rownames(entry$R) & req_snps %in% colnames(entry$R)]
    if (length(available) == 0) {
      return(NULL)
    }
    snp_info <- entry$snp_info[match(available, entry$snp_info$SNP), ]
    return(list(
      R = entry$R[available, available, drop = FALSE],
      snp_info = snp_info
    ))
  }

  candidate_keys <- names(group_entries)
  if (length(candidate_keys) == 0) {
    return(NULL)
  }
  best_entry <- NULL
  best_size <- Inf
  for (entry_key in candidate_keys) {
    entry <- group_entries[[entry_key]]
    if (is.null(entry) || is.null(entry$R)) next
    if (length(entry$snps) < length(req_snps)) next
    if (all(req_snps %in% entry$snps) && length(entry$snps) < best_size) {
      best_entry <- entry
      best_size <- length(entry$snps)
    }
  }
  if (is.null(best_entry)) {
    return(NULL)
  }
  available <- req_snps[
    req_snps %in% best_entry$snps &
      req_snps %in% rownames(best_entry$R) &
      req_snps %in% colnames(best_entry$R)
  ]
  if (length(available) == 0) {
    return(NULL)
  }
  snp_info <- best_entry$snp_info[match(available, best_entry$snp_info$SNP), ]
  list(
    R = best_entry$R[available, available, drop = FALSE],
    snp_info = snp_info
  )
}

ld_cache_store <- function(variants, bfile, keep_file = NULL, ld_res, max_entries = 12L) {
  if (is.null(ld_res) || is.null(ld_res$R)) {
    return(invisible(FALSE))
  }
  mat_snps <- intersect(rownames(ld_res$R), colnames(ld_res$R))
  snps <- sort(unique(as.character(mat_snps)))
  if (length(snps) == 0) {
    return(invisible(FALSE))
  }
  ld_res$R <- ld_res$R[snps, snps, drop = FALSE]
  if (!is.null(ld_res$snp_info) && "SNP" %in% names(ld_res$snp_info)) {
    ld_res$snp_info <- ld_res$snp_info[match(snps, ld_res$snp_info$SNP), ]
  }
  group_key <- ld_cache_group_key(bfile, keep_file)
  group_entries <- if (exists(group_key, envir = .ld_cache, inherits = FALSE)) {
    get(group_key, envir = .ld_cache, inherits = FALSE)
  } else {
    list()
  }
  entry_key <- ld_cache_entry_key(snps)
  group_entries[[entry_key]] <- list(
    snps = snps,
    R = ld_res$R,
    snp_info = ld_res$snp_info,
    cached_at = Sys.time()
  )
  if (length(group_entries) > max_entries) {
    ord <- order(vapply(group_entries, function(x) as.numeric(x$cached_at), numeric(1)), decreasing = TRUE)
    group_entries <- group_entries[ord[seq_len(max_entries)]]
  }
  assign(group_key, group_entries, envir = .ld_cache)
  invisible(TRUE)
}

tabix_cache_key <- function(tabix_file, chrom, start, end) {
  paste(normalizePath(tabix_file, mustWork = FALSE), as.character(chrom), as.integer(start), as.integer(end), sep = "::")
}

# =============================================================================
# APA Gene ID Parser: Handle alternative polyadenylation (APA) transcript format
# =============================================================================
# Parses gene IDs in format: "ENST001|GENE1|..." (GTEx APA format)
# or standard ENSEMBL IDs: "ENSG00000123456.7"
# Returns list with geneSymbol, is_apa flag, and optional transcript_id
# =============================================================================
parse_geneID <- function(geneID, org_db = "org.Hs.eg.db") {
  geneID <- as.character(geneID)

  if (is.na(geneID) || geneID == "") {
    return(list(
      geneSymbol = NA_character_,
      original_id = geneID,
      is_apa = FALSE,
      is_ensembl = FALSE,
      transcript_id = NULL,
      ensembl_id = NULL
    ))
  }

  result <- list(
    geneSymbol = geneID,
    original_id = geneID,
    is_apa = FALSE,
    is_ensembl = FALSE,
    transcript_id = NULL,
    ensembl_id = NULL
  )

  if (grepl("\\|", geneID)) {
    parts <- strsplit(geneID, "\\|")[[1]]
    if (length(parts) >= 2) {
      result$is_apa <- TRUE
      result$geneSymbol <- parts[2]
      result$transcript_id <- parts[1]
      message(sprintf("[parse_geneID] APA format detected: %s -> %s", geneID, result$geneSymbol))
      return(result)
    }
  }

  if (grepl("^ENSG[0-9]{11}", geneID)) {
    result$is_ensembl <- TRUE
    geneID_noDOT <- gsub("\\..*$", "", geneID)
    result$ensembl_id <- geneID_noDOT

    tryCatch(
      {
        symbol <- AnnotationDbi::mapIds(
          get(org_db, envir = asNamespace(org_db)),
          keys = geneID_noDOT,
          column = "SYMBOL",
          keytype = "ENTREZID",
          multiVals = "first"
        )
        if (!is.na(symbol)) {
          result$geneSymbol <- symbol
        }
      },
      error = function(e) {
        tryCatch(
          {
            symbol <- AnnotationDbi::mapIds(
              get(org_db, envir = asNamespace(org_db)),
              keys = geneID_noDOT,
              column = "SYMBOL",
              keytype = "ENSEMBL",
              multiVals = "first"
            )
            if (!is.na(symbol)) {
              result$geneSymbol <- symbol
            }
          },
          error = function(e2) {
            message(sprintf("[parse_geneID] ENSEMBL ID could not be converted: %s", geneID_noDOT))
          }
        )
      }
    )
  }

  if (result$is_ensembl && result$geneSymbol == geneID) {
    message(sprintf("[parse_geneID] Using ENSEMBL ID as symbol: %s", result$geneSymbol))
  }

  return(result)
}

add_variant_id <- function(df, chr_col = "CHR", pos_col = "POS", allele1_col = "NEA", allele2_col = "EA") {
  req_cols <- c(chr_col, pos_col, allele1_col, allele2_col)
  if (!all(req_cols %in% names(df))) {
    stop(glue("Missing columns for ID creation: {paste(setdiff(req_cols, names(df)), collapse=', ')}"))
  }
  # Remove "chr" prefix if already present to avoid duplication
  chr_clean <- gsub("^chr", "", df[[chr_col]])
  v_id <- paste0("chr", chr_clean, ":", df[[pos_col]], ":", df[[allele1_col]], ":", df[[allele2_col]])
  v_id_flip <- paste0("chr", chr_clean, ":", df[[pos_col]], ":", df[[allele2_col]], ":", df[[allele1_col]])
  if (inherits(df, "data.table")) {
    df[, variant_id := v_id]
    df[, variant_id_flip := v_id_flip]
  } else {
    df$variant_id <- v_id
    df$variant_id_flip <- v_id_flip
  }

  return(df)
}

format_sumstats <- function(df, type, col_map, case_control = FALSE, sample_size_n = NULL) {
  df_std <- as.data.table(copy(df))

  target_names <- c(
    snp = "SNPID",
    chrom = "CHR",
    pos = "POS",
    a1 = "EA",
    a2 = "NEA",
    beta = "BETA",
    se = "SE",
    pval = "P",
    n = "N",
    maf = "EAF",
    af = "EAF"
  )

  for (key in names(col_map)) {
    if (key %in% names(target_names) && col_map[[key]] %in% names(df_std)) {
      setnames(df_std, old = col_map[[key]], new = target_names[[key]])
    }
  }

  numeric_cols <- c("BETA", "SE", "P", "POS", "N", "EAF")
  for (col in numeric_cols) {
    if (col %in% names(df_std)) df_std[[col]] <- as.numeric(df_std[[col]])
  }
  if (case_control && "BETA" %notin% names(df_std) && "OR" %in% names(df_std)) {
    df_std[, BETA := log(OR)]
  }
  if (("N" %notin% names(df_std) || all(is.na(df_std$N))) &&
      !is.null(sample_size_n) &&
      !is.na(suppressWarnings(as.numeric(sample_size_n)))) {
    df_std <- data.table::setalloccol(df_std)
    df_std[, N := as.numeric(sample_size_n)]
  }

  return(df_std)
}

easycoloc_standardize_harmonized_gwas <- function(dt, sample_size_n = NULL) {
  res <- as.data.table(copy(dt))
  rename_first <- function(candidates, target) {
    hit <- intersect(candidates, names(res))[1]
    if (!is.na(hit) && nzchar(hit) && hit != target) {
      if (target %in% names(res)) {
        res[, (hit) := NULL]
      } else {
        data.table::setnames(res, hit, target)
      }
    }
  }

  original_id_col <- intersect(c("SNPID", "rsID", "rsid", "SNP", "ID"), names(res))[1]
  if (!is.na(original_id_col) && nzchar(original_id_col) && original_id_col != "SNPID") {
    if ("SNPID" %in% names(res)) {
      res[, (original_id_col) := NULL]
    } else {
      data.table::setnames(res, original_id_col, "SNPID")
    }
  }
  rename_first(c("CHR", "CHROM", "chrom", "chromosome"), "CHR")
  rename_first(c("POS", "BP", "pos", "base_pair_location"), "POS")
  rename_first(c("EA", "A1", "effect_allele", "EFFECT_ALLELE"), "EA")
  rename_first(c("NEA", "A2", "other_allele", "OTHER_ALLELE"), "NEA")
  rename_first(c("EAF", "AF", "MAF", "effect_allele_frequency"), "EAF")
  rename_first(c("BETA", "beta"), "BETA")
  rename_first(c("SE", "standard_error", "se"), "SE")
  rename_first(c("P", "PVAL", "p_value", "pval"), "P")
  rename_first(c("N", "n", "TotalN", "N_analyzed"), "N")

  if ("CHR" %in% names(res)) {
    res[, CHR := gsub("^chr", "", as.character(CHR), ignore.case = TRUE)]
  }
  if ("POS" %in% names(res)) {
    res[, POS := as.integer(POS)]
  }
  for (allele_col in intersect(c("EA", "NEA"), names(res))) {
    res[, (allele_col) := toupper(as.character(get(allele_col)))]
  }
  for (num_col in intersect(c("EAF", "BETA", "SE", "P", "N"), names(res))) {
    res[, (num_col) := as.numeric(get(num_col))]
  }
  if (("N" %notin% names(res) || all(is.na(res$N))) &&
      !is.null(sample_size_n) &&
      !is.na(suppressWarnings(as.numeric(sample_size_n)))) {
    res[, N := as.numeric(sample_size_n)]
  }
  if (all(c("CHR", "POS", "NEA", "EA") %in% names(res))) {
    res[, variant_id := paste0("chr", gsub("^chr", "", CHR, ignore.case = TRUE), ":", POS, ":", NEA, ":", EA)]
  }
  if ("rsID" %notin% names(res)) {
    res[, rsID := NA_character_]
  }
  if ("SNPID" %in% names(res)) {
    snpid_chr <- as.character(res$SNPID)
    rsid_from_snpid <- is.na(res$rsID) & grepl("^rs[0-9]+$", snpid_chr, ignore.case = TRUE)
    res[rsid_from_snpid, rsID := snpid_chr[rsid_from_snpid]]
    if ("variant_id" %in% names(res)) {
      res[is.na(rsID) & snpid_chr == variant_id, rsID := NA_character_]
    }
  }
  if ("SNPID" %notin% names(res)) {
    res[, SNPID := if ("rsID" %in% names(res)) rsID else NA_character_]
  }
  if ("variant_id" %in% names(res)) {
    res[is.na(SNPID) | !nzchar(as.character(SNPID)), SNPID := variant_id]
  }
  if ("rsID" %in% names(res)) {
    res[!is.na(rsID) & nzchar(as.character(rsID)) & grepl("^rs[0-9]+$", as.character(rsID), ignore.case = TRUE), SNPID := as.character(rsID)]
  }

  preferred <- c("SNPID", "rsID", "variant_id", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "P", "N")
  data.table::setcolorder(res, c(intersect(preferred, names(res)), setdiff(names(res), preferred)))
  res
}

easycoloc_harmonized_gwas_output_cols <- function() {
  c("SNPID", "variant_id", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "P", "N")
}

easycoloc_prepare_harmonized_gwas_output <- function(dt, sample_size_n = NULL) {
  res <- easycoloc_standardize_harmonized_gwas(dt, sample_size_n = sample_size_n)
  if ("rsID" %in% names(res)) {
    valid_rsid <- !is.na(res$rsID) &
      nzchar(as.character(res$rsID)) &
      grepl("^rs[0-9]+$", as.character(res$rsID), ignore.case = TRUE)
    res[valid_rsid, SNPID := as.character(rsID)]
  }
  output_cols <- easycoloc_harmonized_gwas_output_cols()
  missing_cols <- setdiff(output_cols, names(res))
  if (length(missing_cols) > 0) {
    stop(glue("Cannot write harmonized GWAS cache; missing required column(s): {paste(missing_cols, collapse = ', ')}"))
  }
  res <- res[, ..output_cols]
  res <- unique(res)
  data.table::setorder(res, CHR, POS, SNPID)
  res
}

easycoloc_normalize_build <- function(build) {
  build_chr <- as.character(build %||% "")
  build_chr <- tolower(gsub("^hg", "", build_chr))
  if (build_chr %in% c("19", "37", "grch37")) {
    return("19")
  }
  if (build_chr %in% c("38", "grch38")) {
    return("38")
  }
  build_chr
}

easycoloc_refseq_contig_map <- function(build = "19") {
  build <- easycoloc_normalize_build(build)
  chr <- c(as.character(1:22), "X", "Y", "MT")
  if (identical(build, "19")) {
    acc <- c(
      sprintf("NC_%06d.%d", 1:22, c(10, 11, 11, 11, 9, 11, 13, 10, 11, 10, 9, 11, 10, 8, 9, 9, 10, 9, 9, 10, 8, 10)),
      "NC_000023.10", "NC_000024.9", "NC_012920.1"
    )
  } else if (identical(build, "38")) {
    acc <- c(
      sprintf("NC_%06d.%d", 1:22, c(11, 12, 12, 12, 10, 12, 14, 11, 12, 11, 10, 12, 11, 9, 10, 10, 11, 10, 10, 11, 9, 11)),
      "NC_000023.11", "NC_000024.10", "NC_012920.1"
    )
  } else {
    acc <- chr
  }
  stats::setNames(acc, chr)
}

easycoloc_map_chr_for_vcf <- function(chr, vcf_file, build = "19") {
  chr_clean <- gsub("^chr", "", as.character(chr), ignore.case = TRUE)
  chr_clean[toupper(chr_clean) == "M"] <- "MT"
  if (is.null(vcf_file) || !nzchar(vcf_file) || !file.exists(vcf_file)) {
    return(chr_clean)
  }
  header <- tryCatch(system2(Sys.which("bcftools"), c("view", "-h", vcf_file), stdout = TRUE, stderr = FALSE), error = function(e) character())
  contig_lines <- grep("^##contig=<ID=", header, value = TRUE)
  contigs <- sub("^##contig=<ID=([^,>]+).*$", "\\1", contig_lines)
  if (length(contigs) == 0) {
    return(chr_clean)
  }
  if (any(chr_clean %in% contigs)) {
    return(chr_clean)
  }
  chr_prefixed <- paste0("chr", chr_clean)
  if (any(chr_prefixed %in% contigs)) {
    return(chr_prefixed)
  }
  refseq <- easycoloc_refseq_contig_map(build)
  mapped <- unname(refseq[chr_clean])
  mapped[is.na(mapped)] <- chr_clean[is.na(mapped)]
  mapped
}

easycoloc_is_queryable_vcf <- function(vcf_file) {
  if (is.null(vcf_file) || !nzchar(vcf_file) || !file.exists(vcf_file)) {
    return(FALSE)
  }
  if (grepl("\\.(vcf|bcf)(\\.gz)?$", vcf_file, ignore.case = TRUE)) {
    return(TRUE)
  }
  bcftools_bin <- Sys.which("bcftools")
  if (!nzchar(bcftools_bin)) {
    return(FALSE)
  }
  header <- tryCatch(system2(bcftools_bin, c("view", "-h", vcf_file), stdout = TRUE, stderr = FALSE), error = function(e) character())
  length(header) > 0 && any(grepl("^#CHROM", header))
}

easycoloc_complement_allele <- function(x) {
  x <- toupper(as.character(x))
  chartr("ACGT", "TGCA", x)
}

easycoloc_is_palindromic <- function(a1, a2) {
  a1 <- toupper(as.character(a1))
  a2 <- toupper(as.character(a2))
  paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
}

easycoloc_flip_effect_stats <- function(res, idx) {
  if (length(idx) == 0) {
    return(res)
  }
  old_ea <- copy(res$EA[idx])
  res[idx, EA := NEA]
  res[idx, NEA := old_ea]
  if ("BETA" %in% names(res)) res[idx, BETA := -BETA]
  if ("Z" %in% names(res)) res[idx, Z := -Z]
  if ("OR" %in% names(res)) res[idx, OR := 1 / OR]
  if ("EAF" %in% names(res)) res[idx[!is.na(res$EAF[idx])], EAF := 1 - EAF]
  res
}

easycoloc_validate_harmonized <- function(res_dt, verbose = TRUE, prefix = "[EasyColoc Harmonize]") {
  log_message <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) message(...)
  }
  n_invalid <- 0L
  if ("P" %in% names(res_dt)) {
    n_invalid_p <- sum(res_dt$P < 0 | res_dt$P > 1 | is.na(res_dt$P))
    if (n_invalid_p > 0) {
      warning(glue("{prefix} Found {n_invalid_p} invalid P-values (< 0, > 1, or NA)"))
      n_invalid <- n_invalid + n_invalid_p
    }
  }
  if ("SE" %in% names(res_dt)) {
    n_invalid_se <- sum(res_dt$SE <= 0 | is.na(res_dt$SE))
    if (n_invalid_se > 0) {
      warning(glue("{prefix} Found {n_invalid_se} invalid SE values (<= 0 or NA)"))
      n_invalid <- n_invalid + n_invalid_se
    }
  }
  if ("BETA" %in% names(res_dt) && "SE" %in% names(res_dt) && "P" %in% names(res_dt)) {
    z_from_beta <- res_dt$BETA / res_dt$SE
    z_from_p <- qnorm(res_dt$P / 2, lower.tail = FALSE) * sign(res_dt$BETA)
    n_inconsistent <- sum(abs(z_from_beta - z_from_p) > 1.0, na.rm = TRUE)
    if (n_inconsistent > nrow(res_dt) * 0.05) {
      warning(glue("{prefix} Found {n_inconsistent}/{nrow(res_dt)} variants with inconsistent BETA/SE/P"))
    }
  }
  if (n_invalid == 0) {
    log_message(glue("{prefix} Post-harmonization validation: PASS"), verbose = verbose)
  }
  invisible(n_invalid == 0)
}

easycoloc_fetch_ref_alleles <- function(dt, ref_fasta, chunk_size = 500000L, verbose = TRUE) {
  if (!requireNamespace("Rsamtools", quietly = TRUE) ||
    !requireNamespace("GenomicRanges", quietly = TRUE) ||
    !requireNamespace("Biostrings", quietly = TRUE)) {
    warning("Rsamtools/GenomicRanges/Biostrings not available; skipping FASTA allele check.")
    return(rep(NA_character_, nrow(dt)))
  }
  if (!file.exists(ref_fasta)) {
    warning(glue("Reference FASTA not found: {ref_fasta}; skipping FASTA allele check."))
    return(rep(NA_character_, nrow(dt)))
  }

  chr <- gsub("^chr", "", as.character(dt$CHR), ignore.case = TRUE)
  pos <- as.integer(dt$POS)
  ref <- rep(NA_character_, nrow(dt))
  ok <- !is.na(chr) & nzchar(chr) & !is.na(pos) & pos > 0
  idx_all <- which(ok)
  if (length(idx_all) == 0) {
    return(ref)
  }

  fai_file <- paste0(ref_fasta, ".fai")
  seqlevels_available <- character()
  if (file.exists(fai_file)) {
    fai <- tryCatch(data.table::fread(fai_file, header = FALSE, select = 1), error = function(e) NULL)
    if (!is.null(fai) && nrow(fai) > 0) seqlevels_available <- as.character(fai[[1]])
  }
  use_chr_prefix <- length(seqlevels_available) > 0 && any(grepl("^chr", seqlevels_available))
  seqnames <- if (use_chr_prefix) paste0("chr", chr) else chr

  fa <- Rsamtools::FaFile(ref_fasta)
  Rsamtools::open.FaFile(fa)
  on.exit(Rsamtools::close.FaFile(fa), add = TRUE)

  starts <- seq.int(1L, length(idx_all), by = chunk_size)
  for (start_i in starts) {
    idx <- idx_all[start_i:min(start_i + chunk_size - 1L, length(idx_all))]
    gr <- GenomicRanges::GRanges(seqnames = seqnames[idx], ranges = IRanges::IRanges(start = pos[idx], width = 1L))
    seqs <- tryCatch(as.character(Rsamtools::scanFa(fa, param = gr)), error = function(e) {
      warning(glue("FASTA lookup failed for chunk starting at row {idx[[1]]}: {conditionMessage(e)}"))
      rep(NA_character_, length(idx))
    })
    ref[idx] <- toupper(seqs)
  }
  ref
}

easycoloc_query_vcf_regions <- function(dt,
                                        vcf_file,
                                        fields,
                                        col_names,
                                        build = "19",
                                        verbose = TRUE,
                                        label = "VCF") {
  if (is.null(vcf_file) || !nzchar(vcf_file) || !file.exists(vcf_file)) {
    return(NULL)
  }
  if (!easycoloc_is_queryable_vcf(vcf_file)) {
    if (isTRUE(verbose)) message(glue("[EasyColoc Harmonize] {label} reference is not VCF/BCF; skipping: {vcf_file}"))
    return(NULL)
  }
  bcftools_bin <- Sys.which("bcftools")
  if (!nzchar(bcftools_bin)) {
    warning(glue("bcftools not found; skipping {label} query."))
    return(NULL)
  }
  if (!all(c("CHR", "POS") %in% names(dt))) {
    return(NULL)
  }
  ok <- !is.na(dt$CHR) & !is.na(dt$POS) & as.integer(dt$POS) > 0
  if (!any(ok)) {
    return(NULL)
  }

  tmp_regions <- tempfile(pattern = "easycoloc_vcf_regions_", fileext = ".bed")
  tmp_out <- tempfile(pattern = "easycoloc_vcf_query_", fileext = ".tsv")
  tmp_script <- tempfile(pattern = "easycoloc_vcf_query_", fileext = ".sh")
  on.exit(file.remove(c(tmp_regions, tmp_out, tmp_script)), add = TRUE)
  region_chr <- easycoloc_map_chr_for_vcf(dt$CHR[ok], vcf_file = vcf_file, build = build)
  region_dt <- data.table::data.table(
    CHROM = region_chr,
    POS = as.integer(dt$POS[ok])
  )
  region_dt <- unique(region_dt[!is.na(CHROM) & !is.na(POS)])
  merge_gap_bp <- suppressWarnings(as.integer(Sys.getenv("EASYCOLOC_VCF_REGION_MERGE_GAP_BP", unset = "100")))
  if (is.na(merge_gap_bp) || merge_gap_bp < 0L) merge_gap_bp <- 0L
  data.table::setorder(region_dt, CHROM, POS)
  region_dt[, GROUP := cumsum(c(TRUE, diff(POS) > merge_gap_bp + 1L)), by = CHROM]
  region_dt <- region_dt[, .(
    FROM = max(0L, min(POS) - 1L),
    TO = max(POS)
  ), by = .(CHROM, GROUP)]
  region_dt[, GROUP := NULL]
  region_dt <- unique(region_dt[!is.na(CHROM) & !is.na(FROM) & !is.na(TO)])
  data.table::setorder(region_dt, CHROM, FROM, TO)
  if (isTRUE(verbose)) {
    message(glue("[EasyColoc Harmonize] {label} query regions: {sum(ok)} variants -> {nrow(region_dt)} intervals (merge_gap_bp={merge_gap_bp})"))
  }
  data.table::fwrite(region_dt, tmp_regions, sep = "\t", col.names = FALSE, quote = FALSE)

  fmt <- paste0(paste(fields, collapse = "\\t"), "\\n")
  script_lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    paste(
      shQuote(bcftools_bin),
      "query",
      "-R", shQuote(tmp_regions),
      "-f", shQuote(fmt),
      shQuote(vcf_file),
      ">", shQuote(tmp_out)
    )
  )
  writeLines(script_lines, tmp_script)
  Sys.chmod(tmp_script, mode = "0700")
  exit_code <- system2("bash", tmp_script, stdout = FALSE, stderr = FALSE)
  if (!identical(as.integer(exit_code), 0L) || !file.exists(tmp_out) || file.info(tmp_out)$size == 0) {
    if (isTRUE(verbose)) message(glue("[EasyColoc Harmonize] {label} query returned no records."))
    return(NULL)
  }
  ref <- tryCatch(data.table::fread(tmp_out, header = FALSE, col.names = col_names), error = function(e) NULL)
  if (is.null(ref) || nrow(ref) == 0) {
    return(NULL)
  }
  if ("CHR" %in% names(ref)) {
    inverse_refseq <- easycoloc_refseq_contig_map(build)
    inverse_refseq <- stats::setNames(names(inverse_refseq), inverse_refseq)
    ref_chr <- as.character(ref$CHR)
    mapped_back <- unname(inverse_refseq[ref_chr])
    ref[, CHR := fifelse(!is.na(mapped_back), mapped_back, gsub("^chr", "", ref_chr, ignore.case = TRUE))]
  }
  if ("POS" %in% names(ref)) ref[, POS := as.integer(POS)]
  ref
}

easycoloc_query_ref_af <- function(dt, ref_vcf, ref_alt_freq = "AF", build = "19", verbose = TRUE) {
  if (is.null(ref_vcf) || !nzchar(ref_vcf) || !file.exists(ref_vcf)) {
    return(NULL)
  }
  ref <- easycoloc_query_vcf_regions(
    dt = dt,
    vcf_file = ref_vcf,
    fields = c("%CHROM", "%POS", "%REF", "%ALT", glue("%INFO/{ref_alt_freq}")),
    col_names = c("CHR", "POS", "REF", "ALT", "REF_AF"),
    build = build,
    verbose = verbose,
    label = "AF reference"
  )
  if (is.null(ref) || nrow(ref) == 0) {
    return(NULL)
  }
  ref[, REF := toupper(as.character(REF))]
  ref[, ALT := toupper(as.character(ALT))]
  ref[, REF_AF := suppressWarnings(as.numeric(sub(",.*$", "", as.character(REF_AF))))]
  ref <- ref[!duplicated(paste(CHR, POS, REF, ALT, sep = ":"))]
  ref
}

easycoloc_assign_rsid_fast <- function(dt, ref_dbsnp, build = "19", verbose = TRUE) {
  res <- as.data.table(copy(dt))
  if (is.null(ref_dbsnp) || !nzchar(ref_dbsnp) || !file.exists(ref_dbsnp)) {
    return(res)
  }
  if (!all(c("CHR", "POS", "EA", "NEA") %in% names(res))) {
    return(res)
  }
  ref <- easycoloc_query_vcf_regions(
    dt = res,
    vcf_file = ref_dbsnp,
    fields = c("%CHROM", "%POS", "%ID", "%REF", "%ALT"),
    col_names = c("CHR", "POS", "rsID", "REF", "ALT"),
    build = build,
    verbose = verbose,
    label = "dbSNP"
  )
  if (is.null(ref) || nrow(ref) == 0) {
    return(res)
  }
  ref[, REF := toupper(as.character(REF))]
  ref[, ALT := toupper(as.character(ALT))]
  ref[, ALT_LIST := strsplit(ALT, ",", fixed = TRUE)]
  ref <- ref[, .(ALT = unlist(ALT_LIST, use.names = FALSE)), by = .(CHR, POS, REF, rsID)]
  ref_forward <- ref[, .(CHR, POS, NEA = REF, EA = ALT, rsID)]
  ref_reverse <- ref[, .(CHR, POS, NEA = ALT, EA = REF, rsID)]
  ref_long <- unique(rbindlist(list(ref_forward, ref_reverse), use.names = TRUE))
  before_has <- "rsID" %in% names(res)
  merged <- merge(
    res[, .ROW_ID := .I][, .(.ROW_ID, CHR, POS, EA, NEA)],
    ref_long,
    by = c("CHR", "POS", "EA", "NEA"),
    all.x = FALSE,
    all.y = FALSE,
    allow.cartesian = TRUE
  )
  if (nrow(merged) == 0) {
    if (".ROW_ID" %in% names(res)) res[, .ROW_ID := NULL]
    return(res)
  }
  merged <- merged[!duplicated(.ROW_ID)]
  if (!before_has) res[, rsID := NA_character_]
  res[merged$.ROW_ID, rsID := merged$rsID]
  if (".ROW_ID" %in% names(res)) res[, .ROW_ID := NULL]
  if (isTRUE(verbose)) {
    message(glue("[EasyColoc Harmonize] dbSNP rsID fast assignment: matched={nrow(merged)}/{nrow(res)}"))
  }
  res
}

easycoloc_infer_palindromic_strand <- function(dt,
                                               ref_vcf,
                                               ref_alt_freq = "AF",
                                               build = "19",
                                               maf_threshold = 0.4,
                                               af_delta = 0.15,
                                               verbose = TRUE) {
  res <- as.data.table(copy(dt))
  if (!all(c("CHR", "POS", "EA", "NEA", "EAF") %in% names(res))) {
    return(res)
  }
  ref <- easycoloc_query_ref_af(res, ref_vcf = ref_vcf, ref_alt_freq = ref_alt_freq, build = build, verbose = verbose)
  if (is.null(ref) || nrow(ref) == 0) {
    return(res)
  }

  ref_forward <- ref[, .(CHR, POS, REF, ALT, REF_AF)]
  ref_reverse <- ref[, .(CHR, POS, REF = ALT, ALT = REF, REF_AF = 1 - REF_AF)]
  ref_long <- rbindlist(list(ref_forward, ref_reverse), use.names = TRUE)
  data.table::setnames(ref_long, c("REF", "ALT", "REF_AF"), c("NEA", "EA", "REF_EAF"))
  merged <- merge(
    res[, .ROW_ID := .I][, .(.ROW_ID, CHR, POS, EA, NEA, EAF)],
    ref_long,
    by = c("CHR", "POS", "EA", "NEA"),
    all.x = FALSE,
    all.y = FALSE,
    allow.cartesian = TRUE
  )
  if (nrow(merged) == 0) {
    return(res)
  }
  merged <- merged[!duplicated(.ROW_ID)]
  pal <- easycoloc_is_palindromic(merged$EA, merged$NEA)
  informative <- pal &
    !is.na(merged$EAF) &
    !is.na(merged$REF_EAF) &
    pmin(merged$EAF, 1 - merged$EAF) <= maf_threshold &
    pmin(merged$REF_EAF, 1 - merged$REF_EAF) <= maf_threshold
  flip_idx <- merged$.ROW_ID[informative & abs((1 - merged$EAF) - merged$REF_EAF) + af_delta < abs(merged$EAF - merged$REF_EAF)]
  if (length(flip_idx) > 0) {
    res <- easycoloc_flip_effect_stats(res, flip_idx)
    if ("ALLELE_FLIPPED" %in% names(res)) res[flip_idx, ALLELE_FLIPPED := TRUE]
  }
  res[, PALINDROMIC := easycoloc_is_palindromic(EA, NEA)]
  res[, AF_STRAND_INFERRED := FALSE]
  if (length(flip_idx) > 0) res[flip_idx, AF_STRAND_INFERRED := TRUE]
  if (".ROW_ID" %in% names(res)) res[, .ROW_ID := NULL]
  if (isTRUE(verbose)) {
    message(glue("[EasyColoc Harmonize] AF strand inference: matched={nrow(merged)}, palindromic={sum(pal)}, flipped={length(flip_idx)}"))
  }
  res
}

easycoloc_align_to_reference <- function(dt, ref_fasta, verbose = TRUE) {
  res <- as.data.table(copy(dt))
  if (!all(c("CHR", "POS", "EA", "NEA") %in% names(res))) {
    warning("Missing CHR/POS/EA/NEA; skipping FASTA allele alignment.")
    return(res)
  }
  ref_base <- easycoloc_fetch_ref_alleles(res, ref_fasta, verbose = verbose)
  ea <- toupper(as.character(res$EA))
  nea <- toupper(as.character(res$NEA))
  ref_ok <- ref_base %in% c("A", "C", "G", "T")
  single <- nchar(ea) == 1L & nchar(nea) == 1L
  can_check <- ref_ok & single

  flip <- can_check & ea == ref_base & nea != ref_base
  comp_nea <- easycoloc_complement_allele(nea)
  comp_ea <- easycoloc_complement_allele(ea)
  strand <- can_check & nea != ref_base & ea != ref_base & comp_nea == ref_base
  strand_flip <- can_check & nea != ref_base & ea != ref_base & comp_ea == ref_base

  if (any(strand, na.rm = TRUE)) {
    res[strand, `:=`(EA = comp_ea[strand], NEA = comp_nea[strand])]
  }
  if (any(strand_flip, na.rm = TRUE)) {
    res[strand_flip, `:=`(EA = comp_nea[strand_flip], NEA = comp_ea[strand_flip])]
  }
  flip_any <- flip | strand_flip
  if (any(flip_any, na.rm = TRUE)) {
    res <- easycoloc_flip_effect_stats(res, which(flip_any))
  }

  ref_match <- !can_check | toupper(as.character(res$NEA)) == ref_base
  res[, REF_ALLELE := ref_base]
  res[, ALLELE_FLIPPED := FALSE]
  res[flip_any, ALLELE_FLIPPED := TRUE]
  res[, STRAND_FLIPPED := FALSE]
  res[strand | strand_flip, STRAND_FLIPPED := TRUE]
  res[, REF_MATCH := ref_match]

  n_checked <- sum(can_check, na.rm = TRUE)
  n_flip <- sum(flip_any, na.rm = TRUE)
  n_strand <- sum(strand | strand_flip, na.rm = TRUE)
  n_mismatch <- sum(can_check & !ref_match, na.rm = TRUE)
  if (isTRUE(verbose)) {
    message(glue("[EasyColoc Harmonize] FASTA allele check: checked={n_checked}, flipped={n_flip}, strand_flipped={n_strand}, unresolved={n_mismatch}"))
  }
  res
}

easycoloc_liftover_table <- function(dt, liftover_chain, source_build, target_build, verbose = TRUE) {
  src <- easycoloc_normalize_build(source_build)
  tgt <- easycoloc_normalize_build(target_build)
  res <- as.data.table(copy(dt))
  if (identical(src, tgt) || !nzchar(src) || !nzchar(tgt)) {
    return(res)
  }
  if (is.null(liftover_chain) || !nzchar(liftover_chain) || !file.exists(liftover_chain)) {
    warning(glue("liftOver chain file not found: {liftover_chain}; coordinates remain on build {source_build}."))
    return(res)
  }
  lift_bin <- Sys.which("liftOver")
  if (!nzchar(lift_bin)) {
    warning("UCSC liftOver binary not found; coordinates remain unchanged.")
    return(res)
  }

  ok <- !is.na(res$CHR) & !is.na(res$POS) & as.integer(res$POS) > 0
  if (!any(ok)) {
    return(res)
  }
  tmp_prefix <- tempfile(pattern = "easycoloc_liftover_")
  bed_file <- paste0(tmp_prefix, ".bed")
  out_file <- paste0(tmp_prefix, ".lifted.bed")
  unmapped_file <- paste0(tmp_prefix, ".unmapped.bed")
  on.exit(file.remove(c(bed_file, out_file, unmapped_file)), add = TRUE)

  idx <- which(ok)
  bed <- data.table::data.table(
    chrom = paste0("chr", gsub("^chr", "", as.character(res$CHR[idx]), ignore.case = TRUE)),
    start = as.integer(res$POS[idx]) - 1L,
    end = as.integer(res$POS[idx]),
    row_id = idx
  )
  data.table::fwrite(bed, bed_file, sep = "\t", col.names = FALSE, quote = FALSE)
  exit_code <- system2(lift_bin, c(bed_file, liftover_chain, out_file, unmapped_file), stdout = FALSE, stderr = FALSE)
  if (!identical(as.integer(exit_code), 0L) || !file.exists(out_file)) {
    warning("liftOver failed; coordinates remain unchanged.")
    return(res)
  }

  lifted <- tryCatch(data.table::fread(out_file, header = FALSE), error = function(e) NULL)
  if (is.null(lifted) || nrow(lifted) == 0) {
    warning("liftOver produced no mapped variants; coordinates remain unchanged.")
    return(res)
  }
  data.table::setnames(lifted, c("chrom", "start", "end", "row_id"))
  lifted <- lifted[!duplicated(row_id)]
  row_id <- as.integer(lifted$row_id)
  res[row_id, CHR := gsub("^chr", "", as.character(lifted$chrom), ignore.case = TRUE)]
  res[row_id, POS := as.integer(lifted$end)]
  res[, LIFTOVER_MAPPED := FALSE]
  res[row_id, LIFTOVER_MAPPED := TRUE]
  if (isTRUE(verbose)) {
    message(glue("[EasyColoc Harmonize] liftOver {source_build}->{target_build}: mapped={length(row_id)}/{sum(ok)}"))
  }
  res
}

easycoloc_harmonize_table_full <- function(sumstats_dt,
                                           ref_fasta,
                                           ref_vcf = NULL,
                                           ref_dbsnp = NULL,
                                           ref_alt_freq = "AF",
                                           source_build = "19",
                                           target_build = "38",
                                           liftover_chain = NULL,
                                           sample_size_n = NULL,
                                           verbose = TRUE) {
  res_dt <- easycoloc_standardize_harmonized_gwas(sumstats_dt, sample_size_n = sample_size_n)
  if ("CHR" %in% names(res_dt)) res_dt[, CHR := gsub("^chr", "", as.character(CHR), ignore.case = TRUE)]
  if ("POS" %in% names(res_dt)) res_dt[, POS := as.integer(POS)]
  for (allele_col in intersect(c("EA", "NEA"), names(res_dt))) {
    res_dt[, (allele_col) := toupper(as.character(get(allele_col)))]
  }
  if ("variant_id" %notin% names(res_dt) && all(c("CHR", "POS", "NEA", "EA") %in% names(res_dt))) {
    res_dt <- add_variant_id(res_dt, chr_col = "CHR", pos_col = "POS", allele1_col = "NEA", allele2_col = "EA")
    if ("variant_id_flip" %in% names(res_dt)) res_dt[, variant_id_flip := NULL]
    if ("SNPID" %notin% names(res_dt)) res_dt[, SNPID := variant_id]
  }

  if (!is.null(ref_fasta) && nzchar(ref_fasta) && file.exists(ref_fasta)) {
    res_dt <- easycoloc_align_to_reference(res_dt, ref_fasta = ref_fasta, verbose = verbose)
  } else {
    warning(glue("Reference FASTA not found: {ref_fasta}; skipping allele alignment."))
  }
  res_dt <- easycoloc_infer_palindromic_strand(
    res_dt,
    ref_vcf = ref_vcf,
    ref_alt_freq = ref_alt_freq,
    build = source_build,
    verbose = verbose
  )
  assign_large_rsid <- tolower(Sys.getenv("EASYCOLOC_ASSIGN_RSID", unset = "auto")) %in% c("1", "true", "yes", "always")
  assign_auto_rsid <- tolower(Sys.getenv("EASYCOLOC_ASSIGN_RSID", unset = "auto")) %in% c("auto", "")
  should_assign_rsid <- !is.null(ref_dbsnp) &&
    nzchar(ref_dbsnp) &&
    file.exists(ref_dbsnp) &&
    (assign_large_rsid || (assign_auto_rsid && nrow(res_dt) <= 1000000L))
  if (should_assign_rsid) {
    res_dt <- easycoloc_assign_rsid_fast(
      res_dt,
      ref_dbsnp = ref_dbsnp,
      build = source_build,
      verbose = verbose
    )
  } else if (!is.null(ref_dbsnp) && nzchar(ref_dbsnp) && file.exists(ref_dbsnp) && isTRUE(verbose)) {
    message("[EasyColoc Harmonize] Skipping dbSNP rsID assignment for large input; set EASYCOLOC_ASSIGN_RSID=always to enable fast sweep.")
  }

  res_dt <- easycoloc_liftover_table(
    res_dt,
    liftover_chain = liftover_chain,
    source_build = source_build,
    target_build = target_build,
    verbose = verbose
  )

  if (("variant_id" %notin% names(res_dt) || all(is.na(res_dt$variant_id) | !nzchar(as.character(res_dt$variant_id)))) &&
      all(c("CHR", "POS", "NEA", "EA") %in% names(res_dt))) {
    res_dt[, variant_id := paste0("chr", gsub("^chr", "", CHR, ignore.case = TRUE), ":", POS, ":", NEA, ":", EA)]
  }
  res_dt <- easycoloc_standardize_harmonized_gwas(res_dt, sample_size_n = sample_size_n)
  easycoloc_validate_harmonized(res_dt, verbose = verbose)
  easycoloc_prepare_harmonized_gwas_output(res_dt, sample_size_n = sample_size_n)
}

easycoloc_harmony_part_dir <- function(final_output_file) {
  paste0(final_output_file, ".parts")
}

easycoloc_harmony_compress_output <- function(compress_output = NULL) {
  if (is.null(compress_output)) return(TRUE)
  isTRUE(compress_output)
}

easycoloc_harmonized_cache_path <- function(save_dir,
                                            dataset_id,
                                            source_build,
                                            target_build,
                                            compress_output = NULL) {
  ext <- if (easycoloc_harmony_compress_output(compress_output)) ".tsv.gz" else ".tsv"
  file.path(save_dir, glue("{dataset_id}_b{source_build}to{target_build}_harmonized{ext}"))
}

easycoloc_harmonized_cache_candidates <- function(save_dir,
                                                  dataset_id,
                                                  source_build,
                                                  target_build,
                                                  compress_output = NULL) {
  preferred <- easycoloc_harmonized_cache_path(
    save_dir = save_dir,
    dataset_id = dataset_id,
    source_build = source_build,
    target_build = target_build,
    compress_output = compress_output
  )
  legacy <- easycoloc_harmonized_cache_path(
    save_dir = save_dir,
    dataset_id = dataset_id,
    source_build = source_build,
    target_build = target_build,
    compress_output = FALSE
  )
  unique(c(preferred, legacy))
}

easycoloc_write_harmonized_cache <- function(dt, path) {
  compress <- if (grepl("\\.gz$", path, ignore.case = TRUE)) "gzip" else "none"
  data.table::fwrite(dt, path, sep = "\t", na = "NA", quote = FALSE, compress = compress)
}

easycoloc_harmony_cache_nonempty <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(FALSE)
  info <- file.info(path)
  !is.na(info$size) && info$size > 0
}

easycoloc_cleanup_harmony_parts <- function(final_output_file, verbose = TRUE) {
  part_dir <- easycoloc_harmony_part_dir(final_output_file)
  if (!dir.exists(part_dir)) return(invisible(FALSE))
  if (!easycoloc_harmony_cache_nonempty(final_output_file)) return(invisible(FALSE))
  unlink(part_dir, recursive = TRUE, force = TRUE)
  if (dir.exists(part_dir)) {
    warning(glue("Could not remove harmonization part directory: {part_dir}"), call. = FALSE)
    return(invisible(FALSE))
  }
  if (isTRUE(verbose)) message(glue("[EasyColoc Harmonize] Removed part directory after successful final cache write: {part_dir}"))
  invisible(TRUE)
}

easycoloc_migrate_legacy_harmony_cache <- function(cache_file,
                                                   preferred_file,
                                                   cached_dt = NULL,
                                                   verbose = TRUE) {
  if (is.null(cache_file) || is.null(preferred_file)) return(cache_file)
  if (!nzchar(cache_file) || !nzchar(preferred_file)) return(cache_file)
  if (identical(normalizePath(cache_file, mustWork = FALSE), normalizePath(preferred_file, mustWork = FALSE))) {
    return(cache_file)
  }
  if (!grepl("_harmonized\\.tsv$", cache_file) || !grepl("_harmonized\\.tsv\\.gz$", preferred_file)) {
    return(cache_file)
  }

  dt <- if (is.null(cached_dt)) {
    data.table::fread(cache_file, showProgress = FALSE)
  } else {
    cached_dt
  }
  dir.create(dirname(preferred_file), recursive = TRUE, showWarnings = FALSE)
  easycoloc_write_harmonized_cache(dt, preferred_file)
  if (!easycoloc_harmony_cache_nonempty(preferred_file)) {
    warning(glue("Could not migrate legacy harmonized cache to compressed output: {preferred_file}"), call. = FALSE)
    return(cache_file)
  }

  easycoloc_cleanup_harmony_parts(preferred_file, verbose = verbose)
  easycoloc_cleanup_harmony_parts(cache_file, verbose = verbose)
  removed_legacy <- file.remove(cache_file)
  if (!isTRUE(removed_legacy) && file.exists(cache_file)) {
    warning(glue("Could not remove legacy harmonized cache after compressed migration: {cache_file}"), call. = FALSE)
  }
  if (isTRUE(verbose)) message(glue("[EasyColoc Harmonize] Migrated legacy cache to compressed output: {preferred_file}"))
  preferred_file
}

easycoloc_file_signature <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return("missing")
  }
  info <- file.info(path)
  paste(
    normalizePath(path, mustWork = FALSE),
    ifelse(is.na(info$size), "NA", info$size),
    ifelse(is.na(info$mtime), "NA", format(info$mtime, "%Y-%m-%d %H:%M:%S %Z")),
    sep = "|"
  )
}

easycoloc_harmony_dependency_signature <- function(input_file = NULL,
                                                   ref_fasta = NULL,
                                                   ref_vcf = NULL,
                                                   ref_dbsnp = NULL,
                                                   liftover_chain = NULL,
                                                   source_build = "19",
                                                   target_build = "38",
                                                   sample_size_n = NULL) {
  pieces <- c(
    paste0("input=", easycoloc_file_signature(input_file)),
    paste0("ref_fasta=", easycoloc_file_signature(ref_fasta)),
    paste0("ref_vcf=", easycoloc_file_signature(ref_vcf)),
    paste0("ref_dbsnp=", easycoloc_file_signature(ref_dbsnp)),
    paste0("liftover_chain=", easycoloc_file_signature(liftover_chain)),
    paste0("source_build=", source_build),
    paste0("target_build=", target_build),
    paste0("sample_size_n=", if (is.null(sample_size_n)) "NA" else sample_size_n),
    paste0("assign_rsid=", Sys.getenv("EASYCOLOC_ASSIGN_RSID", unset = "auto")),
    paste0("schema=", paste(easycoloc_harmonized_gwas_output_cols(), collapse = ","))
  )
  paste(pieces, collapse = ";;")
}

easycoloc_done_value <- function(done_file, key) {
  if (!file.exists(done_file)) return(NA_character_)
  lines <- readLines(done_file, warn = FALSE)
  hit <- grep(paste0("^", key, "="), lines, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", key, "="), "", hit[[length(hit)]])
}

easycoloc_part_file_ok <- function(path, done_file = NULL, dependency_signature = NULL) {
  if (!file.exists(path)) return(FALSE)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) return(FALSE)
  header <- tryCatch(names(data.table::fread(path, nrows = 0, showProgress = FALSE)), error = function(e) character())
  if (!identical(header, easycoloc_harmonized_gwas_output_cols())) return(FALSE)
  if (!is.null(dependency_signature) && !is.null(done_file)) {
    done_sig <- easycoloc_done_value(done_file, "dependency_signature")
    return(identical(done_sig, dependency_signature))
  }
  TRUE
}

easycoloc_chr_sort_key <- function(chr) {
  chr <- gsub("^chr", "", as.character(chr), ignore.case = TRUE)
  fifelse(chr %in% as.character(1:22), sprintf("%02d", as.integer(chr)), fifelse(chr == "X", "23", fifelse(chr == "Y", "24", fifelse(chr %in% c("M", "MT"), "25", paste0("99_", chr)))))
}

easycoloc_harmonize_table_chunked <- function(sumstats_dt,
                                              ref_fasta,
                                              ref_vcf = NULL,
                                              ref_dbsnp = NULL,
                                              ref_alt_freq = "AF",
                                              source_build = "19",
                                              target_build = "38",
                                              liftover_chain = NULL,
                                              sample_size_n = NULL,
                                              final_output_file,
                                              chunk_by = "chr",
                                              chunk_min_rows = 1000000L,
                                              chunk_parallel_jobs = 1L,
                                              input_file = NULL,
                                              dependency_signature = NULL,
                                              verbose = TRUE) {
  if (is.null(final_output_file) || !nzchar(final_output_file)) {
    stop("final_output_file is required for chunked harmonization")
  }
  base_dt <- easycoloc_standardize_harmonized_gwas(sumstats_dt, sample_size_n = sample_size_n)
  if (!"CHR" %in% names(base_dt)) {
    stop("CHR column is required for chunked harmonization")
  }
  base_dt[, CHR := gsub("^chr", "", as.character(CHR), ignore.case = TRUE)]
  part_dir <- easycoloc_harmony_part_dir(final_output_file)
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(dependency_signature)) {
    dependency_signature <- easycoloc_harmony_dependency_signature(
      input_file = input_file,
      ref_fasta = ref_fasta,
      ref_vcf = ref_vcf,
      ref_dbsnp = ref_dbsnp,
      liftover_chain = liftover_chain,
      source_build = source_build,
      target_build = target_build,
      sample_size_n = sample_size_n
    )
  }

  chrs <- unique(base_dt$CHR)
  chrs <- chrs[!is.na(chrs) & nzchar(chrs)]
  chrs <- chrs[order(easycoloc_chr_sort_key(chrs))]
  if (length(chrs) == 0) {
    stop("No valid chromosomes found for chunked harmonization")
  }
  manifest_rows <- vector("list", length(chrs))
  part_files <- character(length(chrs))
  chunk_jobs <- list()
  log_message <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) message(...)
  }
  log_message(glue("[EasyColoc Harmonize] Chunked full harmonization: chunks={length(chrs)}, part_dir={part_dir}"), verbose = verbose)

  for (idx in seq_along(chrs)) {
    chr <- chrs[[idx]]
    safe_chr <- gsub("[^A-Za-z0-9_.-]+", "_", chr)
    part_file <- file.path(part_dir, glue("chr{safe_chr}.tsv"))
    done_file <- paste0(part_file, ".done")
    part_files[[idx]] <- part_file
    chr_n <- sum(base_dt$CHR == chr, na.rm = TRUE)
    if (file.exists(done_file) && easycoloc_part_file_ok(part_file, done_file = done_file, dependency_signature = dependency_signature)) {
      elapsed_sec <- suppressWarnings(as.numeric(easycoloc_done_value(done_file, "elapsed_sec")))
      output_rows <- suppressWarnings(as.integer(easycoloc_done_value(done_file, "output_rows")))
      log_message(glue("[EasyColoc Harmonize] Reusing completed chunk {idx}/{length(chrs)} chr{chr}: {part_file}"), verbose = verbose)
      manifest_rows[[idx]] <- data.table::data.table(chunk = idx, chr = chr, input_rows = chr_n, rows = output_rows, elapsed_sec = elapsed_sec, output = part_file, status = "cached")
      next
    }
    chunk_jobs[[length(chunk_jobs) + 1L]] <- list(
      idx = idx,
      chr = chr,
      chr_n = chr_n,
      part_file = part_file,
      done_file = done_file
    )
  }

  chunk_parallel_jobs <- suppressWarnings(as.integer(chunk_parallel_jobs))
  if (is.na(chunk_parallel_jobs) || chunk_parallel_jobs < 1L) chunk_parallel_jobs <- 1L
  chunk_parallel_jobs <- min(chunk_parallel_jobs, length(chunk_jobs))
  if (length(chunk_jobs) > 0L) {
    log_message(glue("[EasyColoc Harmonize] Processing {length(chunk_jobs)} chunk(s) with {chunk_parallel_jobs} parallel job(s)"), verbose = verbose)
  }

  run_chunk_job <- function(job) {
    idx <- job$idx
    chr <- job$chr
    chr_n <- job$chr_n
    part_file <- job$part_file
    done_file <- job$done_file
    chunk_start <- Sys.time()
    log_message(glue("[EasyColoc Harmonize] Chunk {idx}/{length(chrs)} chr{chr}: input_rows={chr_n} start={format(chunk_start, '%Y-%m-%d %H:%M:%S %Z')}"), verbose = verbose)
    chunk_dt <- base_dt[CHR == chr]
    output_dt <- easycoloc_harmonize_table_full(
      chunk_dt,
      ref_fasta = ref_fasta,
      ref_vcf = ref_vcf,
      ref_dbsnp = ref_dbsnp,
      ref_alt_freq = ref_alt_freq,
      source_build = source_build,
      target_build = target_build,
      liftover_chain = liftover_chain,
      sample_size_n = sample_size_n,
      verbose = verbose
    )
    tmp_part <- tempfile(pattern = paste0(basename(part_file), "."), tmpdir = dirname(part_file), fileext = ".tmp")
    data.table::fwrite(output_dt, tmp_part, sep = "\t", na = "NA", quote = FALSE)
    if (!file.rename(tmp_part, part_file)) {
      file.copy(tmp_part, part_file, overwrite = TRUE)
      file.remove(tmp_part)
    }
    chunk_end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(chunk_end, chunk_start, units = "secs"))
    tmp_done <- tempfile(pattern = paste0(basename(done_file), "."), tmpdir = dirname(done_file), fileext = ".tmp")
    writeLines(c(
      glue("chunk={idx}"),
      glue("chr={chr}"),
      glue("input_rows={chr_n}"),
      glue("output_rows={nrow(output_dt)}"),
      glue("started={format(chunk_start, '%Y-%m-%d %H:%M:%S %Z')}"),
      glue("completed={format(chunk_end, '%Y-%m-%d %H:%M:%S %Z')}"),
      glue("elapsed_sec={round(elapsed_sec, 3)}"),
      glue("dependency_signature={dependency_signature}")
    ), tmp_done)
    if (!file.rename(tmp_done, done_file)) {
      file.copy(tmp_done, done_file, overwrite = TRUE)
      file.remove(tmp_done)
    }
    log_message(glue("[EasyColoc Harmonize] Chunk {idx}/{length(chrs)} chr{chr} complete: output_rows={nrow(output_dt)}, elapsed_sec={round(elapsed_sec, 1)}"), verbose = verbose)
    data.table::data.table(chunk = idx, chr = chr, input_rows = chr_n, rows = nrow(output_dt), elapsed_sec = elapsed_sec, output = part_file, status = "completed")
  }

  if (length(chunk_jobs) > 0L) {
    completed_rows <- if (chunk_parallel_jobs > 1L && .Platform$OS.type == "unix") {
      parallel::mclapply(chunk_jobs, run_chunk_job, mc.cores = chunk_parallel_jobs, mc.set.seed = TRUE)
    } else {
      lapply(chunk_jobs, run_chunk_job)
    }
    for (row in completed_rows) {
      manifest_rows[[row$chunk[[1]]]] <- row
    }
  }

  part_done_files <- paste0(part_files, ".done")
  valid_part_files <- part_files[mapply(function(path, done) {
    easycoloc_part_file_ok(path, done_file = done, dependency_signature = dependency_signature)
  }, part_files, part_done_files)]
  if (length(valid_part_files) != length(part_files)) {
    missing_parts <- setdiff(part_files, valid_part_files)
    stop(glue("Chunked harmonization did not produce all part files: {paste(basename(missing_parts), collapse = ', ')}"))
  }
  manifest <- data.table::rbindlist(manifest_rows, fill = TRUE)
  manifest[, dependency_signature := dependency_signature]
  data.table::fwrite(manifest, file.path(part_dir, "manifest.tsv"), sep = "\t", quote = FALSE, na = "NA")
  log_message(glue("[EasyColoc Harmonize] Combining {length(valid_part_files)} chunk files"), verbose = verbose)
  output_dt <- data.table::rbindlist(lapply(valid_part_files, data.table::fread, showProgress = FALSE), use.names = TRUE)
  output_dt <- easycoloc_prepare_harmonized_gwas_output(output_dt, sample_size_n = sample_size_n)
  easycoloc_write_harmonized_cache(output_dt, final_output_file)
  easycoloc_cleanup_harmony_parts(final_output_file, verbose = verbose)
  log_message(glue("Saved result to: {final_output_file}"), verbose = verbose)
  output_dt
}

run_easycoloc_harmonization <- function(sumstats_dt,
                                        ref_fasta,
                                        ref_vcf = NULL,
                                        ref_dbsnp = NULL,
                                        ref_alt_freq = "AF",
                                        source_build = "19",
                                        target_build = "38",
                                        n_threads = 4,
                                        save_dir = NULL,
                                        dataset_id = NULL,
                                        input_file = NULL,
                                        verbose = TRUE,
                                        liftover_chain = NULL,
                                        sample_size_n = NULL,
                                        chunked = NULL,
                                        chunk_min_rows = 1000000L,
                                        chunk_parallel_jobs = 1L,
                                        compress_output = TRUE,
                                        reuse_cache = TRUE) {
  log_message <- function(..., verbose = TRUE) {
    if (isTRUE(verbose)) message(...)
  }
  log_message("--- Starting EasyColoc Native Harmonization & LiftOver ---", verbose = verbose)
  dependency_signature <- easycoloc_harmony_dependency_signature(
    input_file = input_file,
    ref_fasta = ref_fasta,
    ref_vcf = ref_vcf,
    ref_dbsnp = ref_dbsnp,
    liftover_chain = liftover_chain,
    source_build = source_build,
    target_build = target_build,
    sample_size_n = sample_size_n
  )

  source_file <- input_file
  use_cache <- !is.null(save_dir) && !is.null(dataset_id)
  final_output_file <- NULL
  if (use_cache) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    final_output_file <- easycoloc_harmonized_cache_path(
      save_dir = save_dir,
      dataset_id = dataset_id,
      source_build = source_build,
      target_build = target_build,
      compress_output = compress_output
    )
    cache_candidates <- easycoloc_harmonized_cache_candidates(
      save_dir = save_dir,
      dataset_id = dataset_id,
      source_build = source_build,
      target_build = target_build,
      compress_output = compress_output
    )
    existing_cache <- cache_candidates[file.exists(cache_candidates)]
    if (isTRUE(reuse_cache) && length(existing_cache) > 0L) {
      valid_cache_file <- NULL
      for (cache_file in existing_cache) {
        cache_info <- file.info(cache_file)
        cache_mtime <- cache_info$mtime
        cache_size_ok <- !is.na(cache_info$size) && cache_info$size > 0
        ref_mtime <- if (file.exists(ref_fasta)) file.info(ref_fasta)$mtime else NA
        input_mtime <- if (!is.null(source_file) && file.exists(source_file)) file.info(source_file)$mtime else NA
        dep_mtimes <- c(ref_mtime, input_mtime)
        dep_mtimes <- dep_mtimes[!is.na(dep_mtimes)]
        dep_mtime <- if (length(dep_mtimes) > 0) max(dep_mtimes) else NA
        required_cache_cols <- easycoloc_harmonized_gwas_output_cols()
        cache_header <- if (cache_size_ok) {
          tryCatch(names(data.table::fread(cache_file, nrows = 0, showProgress = FALSE)), error = function(e) character())
        } else {
          character()
        }
        cache_schema_ok <- identical(cache_header, required_cache_cols)
        cache_ok <- cache_size_ok && cache_schema_ok && !is.na(cache_mtime) && !is.na(dep_mtime) && cache_mtime >= dep_mtime
        if (cache_ok) {
          valid_cache_file <- cache_file
          break
        }
        if (!cache_schema_ok) {
          log_message(glue("Cached file lacks canonical harmony schema: {cache_file}"), verbose = verbose)
        } else {
          log_message(glue("Cached file is outdated or empty: {cache_file}"), verbose = verbose)
        }
      }
      if (!is.null(valid_cache_file)) {
        log_message(glue("Found cached harmonized file: {valid_cache_file}"), verbose = verbose)
        cached_dt <- data.table::fread(valid_cache_file)
        if (isTRUE(compress_output)) {
          valid_cache_file <- easycoloc_migrate_legacy_harmony_cache(
            cache_file = valid_cache_file,
            preferred_file = final_output_file,
            cached_dt = cached_dt,
            verbose = verbose
          )
        }
        return(easycoloc_standardize_harmonized_gwas(cached_dt, sample_size_n = sample_size_n))
      }
      log_message(glue("No valid cached harmonized file found; regenerating: {final_output_file}"), verbose = verbose)
    }
  }

  chunked_enabled <- if (is.null(chunked)) {
    use_cache && nrow(sumstats_dt) >= as.integer(chunk_min_rows)
  } else {
    isTRUE(chunked)
  }
  if (chunked_enabled && use_cache) {
    output_dt <- easycoloc_harmonize_table_chunked(
      sumstats_dt,
      ref_fasta = ref_fasta,
      ref_vcf = ref_vcf,
      ref_dbsnp = ref_dbsnp,
      ref_alt_freq = ref_alt_freq,
      source_build = source_build,
      target_build = target_build,
      liftover_chain = liftover_chain,
      sample_size_n = sample_size_n,
      final_output_file = final_output_file,
      chunk_min_rows = chunk_min_rows,
      chunk_parallel_jobs = chunk_parallel_jobs,
      input_file = input_file,
      dependency_signature = dependency_signature,
      verbose = verbose
    )
    log_message(glue("Processing complete. Input SNPs: {nrow(sumstats_dt)} -> Output SNPs: {nrow(output_dt)}"), verbose = verbose)
    return(output_dt)
  }

  output_dt <- easycoloc_harmonize_table_full(
    sumstats_dt,
    ref_fasta = ref_fasta,
    ref_vcf = ref_vcf,
    ref_dbsnp = ref_dbsnp,
    ref_alt_freq = ref_alt_freq,
    source_build = source_build,
    target_build = target_build,
    liftover_chain = liftover_chain,
    sample_size_n = sample_size_n,
    verbose = verbose
  )
  log_message(glue("Processing complete. Input SNPs: {nrow(sumstats_dt)} -> Output SNPs: {nrow(output_dt)}"), verbose = verbose)
  if (use_cache) {
    easycoloc_write_harmonized_cache(output_dt, final_output_file)
    log_message(glue("Saved result to: {final_output_file}"), verbose = verbose)
  }
  output_dt
}

# =============================================================================
# run_liftover_fallback: Iterative region contraction LiftOver fallback
# =============================================================================
# Optional liftOver helper using iterative region contraction. Gradually shrinks the region
# until successful mapping or region becomes too small.
#
# Arguments:
#   sumstats_dt: Input summary statistics data.table
#   chrom: Chromosome number
#   start_pos, end_pos: Region boundaries (will be modified in place)
#   liftOver_chain: Path to UCSC liftOver chain file
#   max_iterations: Maximum contraction iterations (default 20)
#   contraction_step: Base pairs to contract per iteration (default 5000)
#
# Returns:
#   list(success=BOOLEAN, positions=DATAFRAME, start=INT, end=INT)
# =============================================================================
run_liftover_fallback <- function(sumstats_dt,
                                  chrom,
                                  start_pos,
                                  end_pos,
                                  liftOver_chain,
                                  max_iterations = 20,
                                  contraction_step = 5000) {
  message(glue("[LiftOver Fallback] Starting iterative contraction: {start_pos}-{end_pos}"))

  collocStart <- as.integer(start_pos)
  collocStop <- as.integer(end_pos)

  hg38_positions <- NULL
  success <- FALSE

  for (iteration in 1:max_iterations) {
    message(glue("[LiftOver] Attempt {iteration}/{max_iterations}: {collocStart}-{collocStop}"))

    # Create BED file for liftOver
    bed_liftover <- data.frame(
      chrom = paste0("chr", chrom),
      start = c(collocStart - 1, collocStop - 1),
      end = c(collocStart, collocStop)
    )

    bed_file <- tempfile(pattern = "liftover_", fileext = ".bed")
    out_file <- paste0(bed_file, ".lifted")
    unmapped_file <- paste0(bed_file, ".unmapped")

    fwrite(bed_liftover, bed_file,
      sep = "\t",
      col.names = FALSE, quote = FALSE
    )

    # Run liftOver
    cmd <- paste("liftOver", bed_file, liftOver_chain, out_file, unmapped_file, "-bedPlus=3")
    exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    if (exit_code == 0 && file.exists(out_file)) {
      hg38_positions <- tryCatch(
        {
          fread(out_file, header = FALSE)
        },
        error = function(e) NULL
      )
    }

    # Check if mapping was successful (both boundaries mapped)
    if (!is.null(hg38_positions) && nrow(hg38_positions) >= 2) {
      if (!is.na(hg38_positions[2, 3])) {
        success <- TRUE
        message(glue("[LiftOver Fallback] SUCCESS at iteration {iteration}: {hg38_positions[1,3]}-{hg38_positions[2,3]}"))

        # Cleanup temp files
        file.remove(bed_file, out_file, unmapped_file)
        break
      }
    }

    # Contraction strategy: shrink region from both ends
    if (!is.null(hg38_positions) && nrow(hg38_positions) >= 1) {
      # Check which boundary failed
      if (is.na(hg38_positions[1, 3])) {
        message(glue("[LiftOver] Start position {collocStart} unmapped - contracting start"))
        collocStart <- collocStart + contraction_step
      }
      if (nrow(hg38_positions) >= 2 && is.na(hg38_positions[2, 3])) {
        message(glue("[LiftOver] End position {collocStop} unmapped - contracting end"))
        collocStop <- collocStop - contraction_step
      }
    } else {
      # Both failed - contract from both ends
      message(glue("[LiftOver] Both boundaries unmapped - contracting both ends by {contraction_step}bp"))
      collocStart <- collocStart + contraction_step
      collocStop <- collocStop - contraction_step
    }

    # Check if region is too small
    if (collocStart >= collocStop) {
      message("[LiftOver Fallback] Region collapsed - lifting failed")
      break
    }

    if (collocStop - collocStart < 10000) {
      message("[LiftOver Fallback] Region too small (< 10kb) - giving up")
      break
    }
  }

  if (!success) {
    message("[LiftOver Fallback] All iterations failed")
    return(list(
      success = FALSE,
      positions = NULL,
      start = as.integer(start_pos),
      end = as.integer(end_pos)
    ))
  }

  return(list(
    success = TRUE,
    positions = hg38_positions,
    start = hg38_positions[1, 3],
    end = hg38_positions[2, 3]
  ))
}

apply_liftover_positions <- function(sumstats_dt,
                                     chrom,
                                     start_pos,
                                     end_pos,
                                     hg38_positions) {
  dt <- as.data.table(copy(sumstats_dt))
  if (is.null(hg38_positions) || nrow(hg38_positions) < 2) {
    return(dt)
  }

  chr_col <- intersect(c("CHR", "chr", "CHROM"), names(dt))[1]
  pos_col <- intersect(c("POS", "pos", "BP", "POSITION"), names(dt))[1]
  if (is.na(chr_col) || is.na(pos_col)) {
    return(dt)
  }

  start_pos <- as.numeric(start_pos)
  end_pos <- as.numeric(end_pos)
  lifted_start <- as.numeric(hg38_positions[1, 3])
  lifted_end <- as.numeric(hg38_positions[2, 3])
  if (any(is.na(c(start_pos, end_pos, lifted_start, lifted_end)))) {
    return(dt)
  }

  source_span <- end_pos - start_pos
  lifted_span <- lifted_end - lifted_start
  if (source_span == 0) {
    scale_factor <- 1
  } else {
    scale_factor <- lifted_span / source_span
  }

  original_pos <- as.numeric(dt[[pos_col]])
  offset <- original_pos - start_pos
  dt[[pos_col]] <- as.integer(round(lifted_start + offset * scale_factor))
  dt[[chr_col]] <- gsub("^chr", "", as.character(hg38_positions[1, 1]), ignore.case = TRUE)
  dt
}

# =============================================================================
# apply_parse_geneID: Batch parse gene IDs and format for analysis
# =============================================================================
# Applies parse_geneID to a vector of gene IDs and returns formatted data.frame
# Useful for processing QTL phenotype IDs in batch
# =============================================================================
apply_parse_geneID <- function(gene_ids, org_db = "org.Hs.eg.db") {
  results <- lapply(gene_ids, function(gid) {
    parse_geneID(gid, org_db = org_db)
  })

  data.frame(
    original_id = gene_ids,
    geneSymbol = sapply(results, function(r) r$geneSymbol),
    is_apa = sapply(results, function(r) r$is_apa),
    is_ensembl = sapply(results, function(r) r$is_ensembl),
    transcript_id = sapply(results, function(r) r$transcript_id),
    ensembl_id = sapply(results, function(r) r$ensembl_id),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# prep_coloc_input_file: Merge GWAS and QTL summary statistics for colocalization
# =============================================================================
# Merges GWAS and QTL summary statistics by rsID (preferred) or position+allele.
# Handles strand flips, applies P-value floor for numerical stability.
#
# Arguments:
#   gwas_df, qtl_df: Data frames with GWAS/QTL summary statistics
#   gwas_cols, qtl_cols: Lists mapping standard names to actual column names
#   use_hash_table: If TRUE, use SNP hash tables for rsID->position conversion
#   min_snps: Minimum SNPs required for colocalization (default 30)
#   pvalue_floor: Floor value for P-values to prevent -Inf in -log10(P).
#                 Default 1e-300 follows GWAS catalog convention.
#                 Configurable via coloc_settings.pvalue_floor in config/global.yml.
#   verbose: If TRUE, print detailed progress messages
#
# Returns:
#   Merged data.table with columns for both GWAS and QTL (suffixes .gwas/.qtl)
# =============================================================================
prep_coloc_input_file <- function(gwas_df, qtl_df,
                                  gwas_cols = list(pval = "P", beta = "BETA", se = "SE", n = "N"),
                                  qtl_cols = list(pval = "pval_nominal", beta = "slope", se = "slope_se"),
                                  use_hash_table = TRUE,
                                  min_snps = 30,
                                  pvalue_floor = 1e-300,
                                  verbose = FALSE) {
  # localized message function to control verbosity
  msg <- function(...) {
    if (verbose) message(...)
  }

  msg("\n--- MERGING DATA ---")
  gwas_df <- as.data.table(gwas_df)
  qtl_df <- as.data.table(qtl_df)
  msg(glue("[MERGE] Input: GWAS={nrow(gwas_df)} SNPs, QTL={nrow(qtl_df)} SNPs"))
  # Standardize column names for merging
  # Ensure GWAS has 'rsid' (lowercase), QTL has 'variant_id' (lowercase)
  # First remove existing 'rsid' if 'rsID' exists, then rename
  if ("rsID" %in% names(gwas_df)) {
    if ("rsid" %in% names(gwas_df)) {
      gwas_df[, rsid := NULL]
    }
    setnames(gwas_df, "rsID", "rsid")
  }
  if ("variant_id" %in% names(qtl_df)) {
    # variant_id column already exists (lowercase)
  }
  # Check for rsID columns
  gwas_has_rsid <- any(tolower(names(gwas_df)) %in% c("snpid", "snp", "rsid"))
  qtl_has_rsid <- any(tolower(names(qtl_df)) %in% c("snpid", "variant_id", "snp", "rsid"))
  qtl_has_pos <- all(c("CHR", "POS") %in% names(qtl_df)) ||
    all(c("chr", "end") %in% names(qtl_df))
  gwas_has_pos <- all(c("CHR", "POS") %in% names(gwas_df))
  gwas_has_allele <- all(c("EA", "NEA") %in% names(gwas_df))
  qtl_has_allele <- all(c("ref", "alt") %in% names(qtl_df))
  if ("chr" %in% names(qtl_df) && "CHR" %notin% names(qtl_df)) {
    setnames(qtl_df, "chr", "CHR")
  }
  if ("end" %in% names(qtl_df) && "POS" %notin% names(qtl_df)) {
    setnames(qtl_df, "end", "POS")
  }
  if ("ref" %in% names(qtl_df) && "EA" %notin% names(qtl_df)) {
    if ("alt" %in% names(qtl_df)) {
      setnames(qtl_df, c("ref", "alt"), c("REF", "ALT"))
    }
  }
  if ("REF" %in% names(qtl_df) && "EA" %notin% names(qtl_df)) {
    qtl_df[, `:=`(EA = ALT, NEA = REF)]
  }
  merged_dt <- NULL
  merge_strategy <- "none"
  fill_snp_from_candidates <- function(dt, candidates, fallback = NULL) {
    if (is.null(dt) || nrow(dt) == 0) {
      return(dt)
    }
    dt[, snp := NA_character_]
    for (col in candidates) {
      if (col %in% names(dt)) {
        vals <- as.character(dt[[col]])
        fill_idx <- (is.na(dt$snp) | !nzchar(dt$snp)) & !is.na(vals) & nzchar(vals)
        if (any(fill_idx)) {
          dt[fill_idx, snp := vals[fill_idx]]
        }
      }
    }
    if (!is.null(fallback) && fallback %in% names(dt)) {
      vals <- as.character(dt[[fallback]])
      fill_idx <- (is.na(dt$snp) | !nzchar(dt$snp)) & !is.na(vals) & nzchar(vals)
      if (any(fill_idx)) {
        dt[fill_idx, snp := vals[fill_idx]]
      }
    }
    dt
  }

  msg("[MERGE] Attempting Strategy 1: rsID direct match...")
  # Use rsID format for matching (after standardization)
  gwas_id_col <- if ("rsid" %in% names(gwas_df)) "rsid" else if ("SNPID" %in% names(gwas_df)) "SNPID" else if ("SNP" %in% names(gwas_df)) "SNP" else "rsid"
  qtl_id_col <- if ("variant_id" %in% names(qtl_df)) "variant_id" else if ("SNPID" %in% names(qtl_df)) "SNPID" else if ("SNP" %in% names(qtl_df)) "SNP" else "snp"

  # Rename duplicate columns BEFORE merging
  dup_cols <- c("CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "N", "EAF", "AF")
  for (col in dup_cols) {
    if (col %in% names(gwas_df) && col %in% names(qtl_df)) {
      new_name_qtl <- paste0(col, ".qtl")
      new_name_gwas <- paste0(col, ".gwas")
      if (!(new_name_qtl %in% names(qtl_df))) {
        setnames(qtl_df, col, new_name_qtl)
      }
      if (!(new_name_gwas %in% names(gwas_df))) {
        setnames(gwas_df, col, new_name_gwas)
      }
    }
  }

  # Manual merge to avoid data.table segfault
  common_snps <- intersect(gwas_df[[gwas_id_col]], qtl_df[[qtl_id_col]])
  if (length(common_snps) > 0) {
    gwas_subset <- gwas_df[gwas_df[[gwas_id_col]] %in% common_snps, ]
    qtl_subset <- qtl_df[qtl_df[[qtl_id_col]] %in% common_snps, ]
    merged_dt <- merge(
      gwas_subset,
      qtl_subset,
      by.x = gwas_id_col,
      by.y = qtl_id_col,
      all = FALSE
    )
    merged_dt <- as.data.table(merged_dt)
    merge_strategy <- "rsid_direct"
    merged_dt <- fill_snp_from_candidates(
      merged_dt,
      c("rsid.gwas", "SNPID.gwas", "SNP.gwas", "rsid", "SNPID", "SNP"),
      fallback = gwas_id_col
    )
  } else {
    msg("[MERGE] ✗ Strategy 1 failed: No rsID overlap")
    merged_dt <- NULL
  }
  if ((is.null(merged_dt) || nrow(merged_dt) < min_snps) && use_hash_table && exists("convert_rsid_to_pos")) {
    msg("[MERGE] Attempting Strategy 2: hash_table conversion...")
    gwas_id_col <- if ("rsid" %in% names(gwas_df)) "rsid" else if ("SNPID" %in% names(gwas_df)) "SNPID" else if ("SNP" %in% names(gwas_df)) "SNP" else NULL
    if (!is.null(gwas_id_col) && gwas_has_pos) {
      gwas_chrpos <- convert_rsid_to_pos(
        gwas_df[[gwas_id_col]],
        gwas_df$CHR
      )
      gwas_df$gwas_chrpos_id <- gwas_chrpos
      n_converted <- sum(!is.na(gwas_chrpos))
      msg(glue("[MERGE]   Converted {n_converted}/{nrow(gwas_df)} GWAS rsIDs ({round(n_converted/nrow(gwas_df)*100, 1)}%)"))
      qtl_id_col <- if ("variant_id" %in% names(qtl_df)) "variant_id" else if ("snp" %in% names(qtl_df)) "snp" else NULL
      if (!is.null(qtl_id_col)) {
        # Use dplyr to avoid segfault
        merged_dt <- inner_join(
          gwas_df,
          qtl_df,
          by = c("gwas_chrpos_id" = qtl_id_col),
          suffix = c(".gwas", ".qtl"),
          relationship = "many-to-many"
        )
        merged_dt <- as.data.table(merged_dt)
        msg(glue("[MERGE] ✓ Strategy 2: {nrow(merged_dt)} SNPs matched"))
        merge_strategy <- "hash_conversion"
        merged_dt <- fill_snp_from_candidates(
          merged_dt,
          c("rsid.gwas", "SNPID.gwas", "SNP.gwas", "rsid", "SNPID", "SNP"),
          fallback = "gwas_chrpos_id"
        )
        snp_sample <- head(merged_dt$snp[!is.na(merged_dt$snp)], 3)
        if (length(snp_sample) > 0 && !all(grepl("^rs", snp_sample))) {
          msg(glue("[MERGE]   Warning: Non-rsID format detected: {paste(snp_sample, collapse=', ')}"))
        }
      } else {
        msg("[MERGE] ✗ Strategy 2 failed: No matches after hash conversion")
      }
    }
  }
  if ((is.null(merged_dt) || nrow(merged_dt) < min_snps) && gwas_has_pos && qtl_has_pos) {
    msg("[MERGE] Attempting Strategy 3: position + allele matching...")
    # Use renamed columns if they exist, otherwise use original names
    gwas_chr_col <- if ("CHR.gwas" %in% names(gwas_df)) "CHR.gwas" else "CHR"
    gwas_pos_col <- if ("POS.gwas" %in% names(gwas_df)) "POS.gwas" else "POS"
    gwas_ea_col <- if ("EA.gwas" %in% names(gwas_df)) "EA.gwas" else "EA"
    gwas_nea_col <- if ("NEA.gwas" %in% names(gwas_df)) "NEA.gwas" else "NEA"

    qtl_chr_col <- if ("CHR.qtl" %in% names(qtl_df)) "CHR.qtl" else "CHR"
    qtl_pos_col <- if ("POS.qtl" %in% names(qtl_df)) "POS.qtl" else "POS"
    qtl_ea_col <- if ("EA.qtl" %in% names(qtl_df)) "EA.qtl" else "EA"
    qtl_nea_col <- if ("NEA.qtl" %in% names(qtl_df)) "NEA.qtl" else "NEA"

    if (gwas_has_allele) {
      # GWAS SNPID format is chr:pos:NEA:EA, so we use NEA:EA order
      gwas_df <- add_variant_id(gwas_df, gwas_chr_col, gwas_pos_col, gwas_nea_col, gwas_ea_col)
    }
    # QTL variant_id format is chr:pos:ref:alt = chr:pos:NEA:EA
    # Since QTL's NEA=REF and EA=ALT, we use NEA:EA = REF:ALT
    if (all(c(qtl_nea_col, qtl_ea_col) %in% names(qtl_df))) {
      qtl_df <- add_variant_id(qtl_df, qtl_chr_col, qtl_pos_col, qtl_nea_col, qtl_ea_col)
    }
    if ("variant_id" %in% names(gwas_df) && "variant_id" %in% names(qtl_df)) {
      # Forward Match
      merged_forward <- inner_join(
        gwas_df, qtl_df,
        by = "variant_id",
        suffix = c(".gwas", ".qtl"),
        relationship = "many-to-many"
      )
      merged_forward <- as.data.table(merged_forward)
      merged_flip <- inner_join(
        gwas_df, qtl_df,
        by = c("variant_id_flip" = "variant_id"),
        suffix = c(".gwas", ".qtl"),
        relationship = "many-to-many"
      )
      merged_flip <- as.data.table(merged_flip)
      # Mark flipped SNPs and flip BETA values.
      if (nrow(merged_flip) > 0) {
        merged_flip$allele_flipped <- TRUE
        if ("ALLELE_FLIPPED" %in% names(merged_flip)) {
          already_flipped <- any(merged_flip$ALLELE_FLIPPED == TRUE, na.rm = TRUE)
          if (already_flipped) {
            warning("Detected ALLELE_FLIPPED=TRUE - variants were already flipped during harmonization. Skipping downstream flip to prevent double-flipping.")
          } else {
            # ALLELE_FLIPPED=FALSE or all NA, safe to flip downstream
            if ("BETA.qtl" %in% names(merged_flip)) {
              merged_flip$BETA.qtl <- -merged_flip$BETA.qtl
            }
            if ("BETA.gwas" %in% names(merged_flip)) {
              merged_flip$BETA.gwas <- -merged_flip$BETA.gwas
            }
            merged_flip$downstream_flipped <- TRUE
          }
        } else {
          if ("BETA.qtl" %in% names(merged_flip)) {
            merged_flip$BETA.qtl <- -merged_flip$BETA.qtl
          }
          if ("BETA.gwas" %in% names(merged_flip)) {
            merged_flip$BETA.gwas <- -merged_flip$BETA.gwas
          }
          merged_flip$downstream_flipped <- TRUE
        }
      }
      if (nrow(merged_forward) > 0) {
        merged_forward$allele_flipped <- FALSE
      }

      merged_dt <- rbindlist(list(merged_forward, merged_flip), fill = TRUE, use.names = TRUE)
      if (nrow(merged_dt) > 0) {
        msg(glue("[MERGE] ✓ Strategy 3: {nrow(merged_dt)} SNPs matched (Forward: {nrow(merged_forward)}, Flipped: {nrow(merged_flip)})"))
        merge_strategy <- "pos_allele"
        merged_dt <- fill_snp_from_candidates(
          merged_dt,
          c("rsid.gwas", "SNPID.gwas", "SNP.gwas", "SNPID", "SNP", "variant_id"),
          fallback = "variant_id"
        )
        snp_sample <- head(merged_dt$snp[!is.na(merged_dt$snp)], 3)
        if (length(snp_sample) > 0 && !all(grepl("^rs", snp_sample))) {
          msg(glue("[MERGE]   Warning: Non-rsID format detected: {paste(snp_sample, collapse=', ')}"))
        }
      } else {
        msg("[MERGE] ✗ Strategy 3 failed: No position + allele matches")
      }
    }
  }
  if (is.null(merged_dt) || nrow(merged_dt) == 0) {
    stop("*** MERGE FAILED: No overlapping SNPs found using any strategy! ***")
  }

  # =============================================================================
  # Final Summary
  # =============================================================================
  # P-value Floor Application: Numerical Stability
  # =============================================================================
  # Rationale: P-values <= 0 or NA cause -Inf in -log10(P) transformations
  # used in visualization. The floor value 1e-300 follows GWAS catalog
  # convention and stays safely above R's .Machine$double.xmin (~2.2e-308).
  # This prevents numerical instability without introducing meaningful bias
  # since such extreme P-values are effectively zero for practical interpretation.
  # Configurable via coloc_settings.pvalue_floor in config/global.yml.
  p_gwas <- paste0(gwas_cols$pval, ".gwas")
  p_qtl <- paste0(qtl_cols$pval, ".qtl")

  for (p_col in c(p_gwas, p_qtl)) {
    if (p_col %in% names(merged_dt)) {
      p_values <- merged_dt[[p_col]]
      needs_floor <- is.na(p_values) | p_values <= 0 | p_values < pvalue_floor
      n_fix <- sum(needs_floor, na.rm = TRUE)
      if (n_fix > 0) {
        merged_dt[needs_floor, (p_col) := pvalue_floor]
        msg(glue("  Replaced {n_fix} P values (NA/<=0/<floor) with {pvalue_floor} in {p_col}"))
      }
    }
  }

  # Fix column names that might have double suffixes
  # Find and rename columns like "P.gwas.gwas" to "P.gwas"
  double_suffix_cols <- grep("\\.gwas\\.gwas$|\\.qtl\\.qtl$", names(merged_dt), value = TRUE)
  if (length(double_suffix_cols) > 0) {
    msg(glue("  Fixing {length(double_suffix_cols)} double-suffix columns..."))
    for (col in double_suffix_cols) {
      new_col <- gsub("\\.gwas\\.gwas$", ".gwas", col)
      new_col <- gsub("\\.qtl\\.qtl$", ".qtl", new_col)
      setnames(merged_dt, col, new_col)
    }
  }

  # Find actual P column (GWAS P should be "P", QTL P should be "P.qtl")
  p_gwas_col <- if ("P" %in% names(merged_dt)) "P" else grep("^P\\.gwas$", names(merged_dt), value = TRUE)[1]
  if (is.na(p_gwas_col)) {
    p_gwas_col <- "P" # Default fallback
  }
  msg(glue("  Using {p_gwas_col} for deduplication"))

  # Deduplication using data.table
  if (anyDuplicated(merged_dt$snp)) {
    n_dup <- sum(duplicated(merged_dt$snp))
    msg(glue("  Removing {n_dup} duplicate SNPs (keeping lowest {p_gwas_col})..."))

    # Use data.table for deduplication
    setorderv(merged_dt, p_gwas_col)
    merged_dt <- merged_dt[!duplicated(merged_dt$snp)]
  }

  # Quality Check
  required_cols <- c(
    "snp",
    paste0(c("BETA", "SE", "P"), ".gwas"),
    paste0(c("BETA", "SE", "P"), ".qtl")
  )

  missing_cols <- setdiff(required_cols, names(merged_dt))
  if (length(missing_cols) > 0) {
    warning(glue("Missing columns: {paste(missing_cols, collapse=', ')}"))
  }

  msg(glue("Final dataset ready: {nrow(merged_dt)} SNPs\n"))
  return(merged_dt)
}

query_tabix_region <- function(tabix_file, chrom, start, end, col_names = NULL) {
  if (!file.exists(tabix_file)) {
    warning(glue("Tabix file not found: {tabix_file}"))
    return(NULL)
  }
  cache_key <- tabix_cache_key(tabix_file, chrom, start, end)
  if (exists(cache_key, envir = .tabix_cache, inherits = FALSE)) {
    cached <- get(cache_key, envir = .tabix_cache, inherits = FALSE)
    return(if (is.null(cached)) NULL else copy(cached))
  }
  chrom_clean <- gsub("^chr", "", as.character(chrom))

  run_query <- function(query_chrom) {
    cmd <- glue("tabix {tabix_file} {query_chrom}:{start}-{end}")
    dt <- tryCatch(
      {
        suppressWarnings(fread(cmd = cmd, header = FALSE, showProgress = FALSE))
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (!is.null(dt) && nrow(dt) > 0) {
      return(dt)
    }
    return(NULL)
  }
  dt <- run_query(paste0("chr", chrom_clean))
  if (is.null(dt)) dt <- run_query(chrom_clean)
  assign(cache_key, if (is.null(dt)) NULL else copy(dt), envir = .tabix_cache)
  return(dt)
}

get_qtl_candidate_phenotypes <- function(meta_row, cfg_qtl, chrom, start, end,
                                         prefilter_sig_pairs = TRUE,
                                         verbose = FALSE) {
  if (!isTRUE(prefilter_sig_pairs)) {
    return(NULL)
  }
  sig_col <- cfg_qtl$qtl_info$columns$sig_filename
  if (is.null(sig_col) || !sig_col %in% names(meta_row)) {
    return(NULL)
  }

  sig_file <- as.character(meta_row[[sig_col]])
  if (is.na(sig_file) || !nzchar(sig_file) || !file.exists(sig_file)) {
    return(NULL)
  }

  sig_raw <- tryCatch(
    query_tabix_region(sig_file, chrom, start, end),
    error = function(e) NULL
  )
  if (is.null(sig_raw) || nrow(sig_raw) == 0) {
    return(character(0))
  }

  if (ncol(sig_raw) == length(cfg_qtl$QTL_all_header)) {
    colnames(sig_raw) <- cfg_qtl$QTL_all_header
  }

  pheno_col <- cfg_qtl$QTL_cols$phenotype
  if (is.null(pheno_col) || !pheno_col %in% names(sig_raw)) {
    if (isTRUE(verbose)) {
      message("[QTL] sigPairs phenotype column unavailable; falling back to allPairs query")
    }
    return(NULL)
  }

  phenos <- unique(as.character(sig_raw[[pheno_col]]))
  phenos <- phenos[!is.na(phenos) & nzchar(phenos)]
  phenos
}

get_ld_matrix <- function(variants, bfile, plink_bin, keep_file = NULL) {
  if (length(variants) == 0) {
    return(NULL)
  }
  variants <- sort(unique(as.character(variants)))
  cached_res <- ld_cache_lookup(variants, bfile, keep_file = keep_file)
  if (!is.null(cached_res)) {
    return(cached_res)
  }
  cached_res <- ld_disk_cache_read(variants, bfile, keep_file = keep_file)
  if (!is.null(cached_res)) {
    return(cached_res)
  }
  fn <- tempfile()
  on.exit(unlink(paste0(fn, "*")), add = TRUE)
  tryCatch(
    {
      fwrite(list(variants), file = fn, col.names = FALSE)
      plink_keep_file <- easycoloc_prepare_plink_keep_file(
        keep_file,
        bfile = bfile,
        temp_dir = dirname(fn),
        label = "ld",
        verbose = FALSE
      )
      keep_cmd <- if (!is.null(plink_keep_file) && file.exists(plink_keep_file)) paste0("--keep ", shQuote(plink_keep_file)) else ""
      plink_args <- c(
        "--bfile", bfile,
        "--extract", fn,
        if (!is.null(plink_keep_file) && file.exists(plink_keep_file)) c("--keep", plink_keep_file) else character(),
        "--allow-extra-chr",
        "--r", "square",
        "--make-just-bim",
        "--out", fn,
        "--silent"
      )
      if (exists("run_ld_plink", mode = "function")) {
        status <- run_ld_plink(plink_bin, plink_args, label = glue("matrix LD ({length(variants)} variants)"))
      } else {
        status <- system2(
          plink_bin,
          args = plink_args,
          stdout = FALSE,
          stderr = FALSE,
          timeout = 120L
        )
      }
      if (!identical(as.integer(status), 0L)) {
        message(glue("[LD] Matrix LD extraction skipped after PLINK status={as.integer(status)}"))
        return(NULL)
      }
      if (!file.exists(paste0(fn, ".ld"))) {
        return(NULL)
      }
      ld_mat <- as.matrix(fread(paste0(fn, ".ld")))
      bim <- fread(paste0(fn, ".bim"), col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))
      colnames(ld_mat) <- bim$SNP
      rownames(ld_mat) <- bim$SNP
      ld_res <- list(R = ld_mat, snp_info = bim)
      ld_cache_store(variants, bfile, keep_file = keep_file, ld_res = ld_res)
      ld_disk_cache_write(variants, bfile, keep_file = keep_file, ld_res = ld_res)
      return(ld_res)
    },
    error = function(e) {
      return(NULL)
    }
  )
}

annotate_coloc_locus_capture <- function(results_dt) {
  if (is.null(results_dt) || nrow(results_dt) == 0) {
    return(copy(results_dt))
  }

  dt <- as.data.table(copy(results_dt))
  if (!all(c("GWAS_ID", "QTL_ID", "Phenotype", "Locus") %in% names(dt))) {
    return(dt)
  }
  dt[
    ,
    `:=`(
      phenotype_locus_count = uniqueN(Locus),
      phenotype_loci = paste(sort(unique(as.character(Locus))), collapse = ";")
    ),
    by = .(GWAS_ID, QTL_ID, Phenotype)
  ]
  dt
}

# ------------------------------------------------------------------------------
# merge_all_results: Aggregate colocalization results across all loci
# ------------------------------------------------------------------------------
merge_all_results <- function(output_dir,
                              pp4_threshold = 0.8,
                              merge_susie = TRUE,
                              save_summary = TRUE) {
  abf_dir <- file.path(output_dir, "abf")
  susie_dir <- file.path(output_dir, "susie")

  # Find all locus result files
  abf_files <- list.files(abf_dir, pattern = "_locus_results\\.csv$", full.names = TRUE)

  if (length(abf_files) == 0) {
    message("[SUM] No ABF results found to merge.")
    return(invisible(NULL))
  }

  message(glue("[SUM] Found {length(abf_files)} locus result file(s) to merge"))

  # Read and combine all ABF results
  all_results <- rbindlist(lapply(abf_files, function(f) {
    tryCatch(
      {
        fread(f)
      },
      error = function(e) {
        warning(glue("Failed to read {basename(f)}: {e$message}"))
        return(NULL)
      }
    )
  }), fill = TRUE)

  if (nrow(all_results) == 0) {
    message("[SUM] No valid results found after merging.")
    return(invisible(NULL))
  }

  all_results <- annotate_coloc_locus_capture(all_results)
  all_results <- all_results[order(-PP4)]

  # Filter by PP4 threshold
  sig_results <- all_results[PP4 >= pp4_threshold]

  # Save merged results
  merged_file <- file.path(output_dir, "all_colocalization_results.csv")
  fwrite(all_results, merged_file)
  message(glue("[SUM] Saved merged results: {basename(merged_file)} ({nrow(all_results)} rows)"))

  stale_dedup_files <- c(
    file.path(output_dir, "deduplicated_colocalization_results.csv"),
    list.files(
      output_dir,
      pattern = "^significant_unique_trait_phenotype_PP4_.*\\.csv$",
      full.names = TRUE
    )
  )
  stale_dedup_files <- stale_dedup_files[file.exists(stale_dedup_files)]
  if (length(stale_dedup_files) > 0) {
    unlink(stale_dedup_files)
    message(glue("[SUM] Removed stale deduplicated output file(s): {length(stale_dedup_files)}"))
  }

  stale_significant_files <- list.files(
    output_dir,
    pattern = "^significant_colocalizations_PP4_.*\\.csv$",
    full.names = TRUE
  )
  if (length(stale_significant_files) > 0) {
    unlink(stale_significant_files)
    message(glue("[SUM] Removed stale significant output file(s): {length(stale_significant_files)}"))
  }

  # Save significant results
  if (nrow(sig_results) > 0) {
    sig_file <- file.path(output_dir, glue("significant_colocalizations_PP4_{pp4_threshold}.csv"))
    fwrite(sig_results, sig_file)
    message(glue("[SUM] Saved significant results (PP4 >= {pp4_threshold}): {basename(sig_file)} ({nrow(sig_results)} rows)"))
  } else {
    message(glue("[SUM] No significant results found with PP4 >= {pp4_threshold}"))
  }

  # Merge SuSiE results if requested
  if (merge_susie && dir.exists(susie_dir)) {
    susie_files <- list.files(susie_dir, pattern = "_susie\\.csv$", full.names = TRUE)

    if (length(susie_files) > 0) {
      message(glue("[SUM] Found {length(susie_files)} SuSiE result file(s)"))

      susie_results <- rbindlist(lapply(susie_files, function(f) {
        tryCatch(
          {
            dt <- fread(f)
            dt$source_file <- basename(f)
            return(dt)
          },
          error = function(e) {
            warning(glue("Failed to read {basename(f)}: {e$message}"))
            return(NULL)
          }
        )
      }), fill = TRUE)

      if (nrow(susie_results) > 0) {
        susie_merged_file <- file.path(output_dir, "all_susie_results.csv")
        fwrite(susie_results, susie_merged_file)
        message(glue("[SUM] Saved merged SuSiE results: {basename(susie_merged_file)} ({nrow(susie_results)} rows)"))
      }
    }
  }

  # Generate summary statistics if requested
  if (save_summary) {
    summary_stats <- list(
      total_tests = nrow(all_results),
      significant_colocalizations = nrow(sig_results),
      unique_trait_phenotype_tests = uniqueN(all_results[, .(GWAS_ID, QTL_ID, Phenotype)]),
      unique_significant_trait_phenotype = if (nrow(sig_results) > 0) uniqueN(sig_results[, .(GWAS_ID, QTL_ID, Phenotype)]) else 0L,
      unique_gwas_traits = length(unique(all_results$GWAS_ID)),
      unique_qtl_datasets = length(unique(all_results$QTL_ID)),
      unique_loci = length(unique(all_results$Locus)),
      unique_phenotypes = length(unique(all_results$Phenotype)),
      mean_pp4 = mean(all_results$PP4, na.rm = TRUE),
      median_pp4 = median(all_results$PP4, na.rm = TRUE),
      max_pp4 = max(all_results$PP4, na.rm = TRUE),
      mean_n_snps = mean(all_results$n_snps, na.rm = TRUE)
    )

    # Convert to data.frame for pretty printing
    summary_df <- data.frame(
      Metric = names(summary_stats),
      Value = unlist(summary_stats)
    )

    summary_file <- file.path(output_dir, "summary_statistics.txt")

    # Write formatted summary
    sink(summary_file)
    cat("=============================================================\n")
    cat("EasyColoc Analysis Summary\n")
    cat("=============================================================\n")
    cat(glue("Analysis completed: {Sys.time()}\n\n"))
    cat("Summary Statistics:\n")
    cat("-------------------------------------------------------------\n")
    print(summary_df, row.names = FALSE)
    cat("\n")

    cat("Locus Capture Note:\n")
    cat("-------------------------------------------------------------\n")
    cat("`significant_colocalizations` counts all locus-level rows.\n")
    cat("`unique_significant_trait_phenotype` counts unique GWAS/QTL/Phenotype combinations without choosing a representative locus.\n")
    cat("`phenotype_locus_count` and `phenotype_loci` annotate all rows with how many distinct loci captured the same GWAS/QTL/Phenotype.\n\n")

    if (nrow(sig_results) > 0) {
      cat("\nTop 10 Locus-Level Hits by PP4:\n")
      cat("-------------------------------------------------------------\n")
      top_events <- head(sig_results[, .(GWAS_ID, QTL_ID, Locus, Phenotype, PP4, n_snps, phenotype_locus_count, phenotype_loci)], 10)
      print(top_events, row.names = FALSE)
    }

    cat("\n=============================================================\n")
    sink()

    message(glue("[SUM] Saved summary statistics: {basename(summary_file)}"))
  }

  message("[SUM] Merge complete!")
  return(invisible(all_results))
}
