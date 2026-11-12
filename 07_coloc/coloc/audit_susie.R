#!/usr/bin/env Rscript

# Refit EasyColoc SuSiE comparisons with stage-specific QTL LD and max_iter = 1000.
# GWAS LD is the ancestry-matched 1000G reference stored in each RDS.
# SCAF1 has no historical SuSiE CSV, so it is run as a supplemental comparison.

# --- setup ---
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(parallel)
  library(susieR)
})

# --- config ---
APA_ROOT <- Sys.getenv("APA_ROOT")
EASYCOLOC_ROOT <- Sys.getenv("EASYCOLOC_ROOT", file.path(APA_ROOT, "release/07_coloc/easycoloc"))
DATA_DIR <- file.path(APA_ROOT, "data/07_coloc")
RESULTS_DIR <- file.path(DATA_DIR, "coloc_run")
SUSIE_DIR <- file.path(RESULTS_DIR, "susie")
RDS_DIR <- file.path(RESULTS_DIR, "rds")
OUTPUT_PREFIX <- file.path(DATA_DIR, "susie_maxiter1000")
OLD_AUDIT_FILE <- file.path(RESULTS_DIR, "susie_convergence_audit.tsv")
QTL_N <- c(prenatal = 189, postnatal = 508, interaction = 697)
QTL_PFILE <- c(
  prenatal = file.path(APA_ROOT, "data/04_qtl/pre_input/genotype.prenatal"),
  postnatal = file.path(APA_ROOT, "data/04_qtl/pre_input/genotype.postnatal"),
  interaction = file.path(APA_ROOT, "data/05_gxs/genotype")
)
PLINK2_BIN <- Sys.getenv("PLINK2_BIN", "plink2")
HIGHLIGHT_GENES <- c("ARL3", "MTRFR", "OGFOD2", "SCAF1", "ATP2B1")
L_VALUE <- 10L
MAX_ITER <- 1000L
GLOBAL_SEED <- 2026L
N_CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

source(file.path(EASYCOLOC_ROOT, "src/utils_format.R"))

output_files <- c(
  comparison = paste0(OUTPUT_PREFIX, "_comparison.tsv"),
  diagnostics = paste0(OUTPUT_PREFIX, "_diagnostics.tsv"),
  coloc = paste0(OUTPUT_PREFIX, "_coloc.tsv"),
  summary = paste0(OUTPUT_PREFIX, "_old_new_summary.tsv"),
  highlighted = paste0(OUTPUT_PREFIX, "_highlighted_genes.tsv")
)

# --- helpers ---
empty_coloc_row <- function(comparison_file, rds_file, analysis_scope,
                            writeback_status) {
  data.table(
    comparison_file = comparison_file, rds_file = basename(rds_file),
    analysis_scope = analysis_scope, nsnps = NA_integer_,
    hit1 = NA_character_, hit2 = NA_character_,
    PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
    PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
    idx1 = NA_integer_, idx2 = NA_integer_,
    writeback_status = writeback_status
  )
}

diagnose_fit <- function(z, R, n, side, common_snps, fit_converged) {
  s_value <- tryCatch(
    estimate_s_rss(z = z, R = R, n = n),
    error = function(e) e
  )
  if (inherits(s_value, "error")) {
    return(data.table(
      side = side, s = NA_real_, n_kriging_outliers = NA_integer_,
      kriging_outlier_snps = NA_character_, diagnostic_status = "estimate_s_error",
      diagnostic_error = conditionMessage(s_value)
    ))
  }
  kriging <- tryCatch(
    kriging_rss(z = z, R = R, n = n, s = s_value),
    error = function(e) e
  )
  if (inherits(kriging, "error")) {
    return(data.table(
      side = side, s = s_value, n_kriging_outliers = NA_integer_,
      kriging_outlier_snps = NA_character_, diagnostic_status = "kriging_error",
      diagnostic_error = conditionMessage(kriging)
    ))
  }
  conditional <- as.data.table(kriging$conditional_dist)
  outlier_index <- which(conditional$logLR > 2 & abs(conditional$z) > 2)
  data.table(
    side = side,
    s = s_value,
    n_kriging_outliers = length(outlier_index),
    kriging_outlier_snps = if (length(outlier_index)) {
      paste(common_snps[outlier_index], collapse = ";")
    } else {
      ""
    },
    diagnostic_status = "complete",
    diagnostic_error = NA_character_
  )
}

fit_one <- function(task, task_index) {
  comparison_file <- task$comparison_file
  rds_file <- task$rds_file
  scope <- task$analysis_scope
  base_record <- data.table(
    comparison_file = comparison_file,
    rds_file = basename(rds_file),
    analysis_scope = scope,
    status = "pending"
  )
  fail <- function(fail_status, message = NA_character_) {
    base_record[, `:=`(
      status = fail_status, error_message = message,
      comparison_converged = FALSE, new_max_pp4 = NA_real_
    )]
    list(
      comparison = base_record,
      diagnostics = data.table(),
      coloc = empty_coloc_row(comparison_file, rds_file, scope, fail_status)
    )
  }

  if (!file.exists(rds_file)) return(fail("missing_rds"))
  object <- tryCatch(readRDS(rds_file), error = function(e) e)
  if (inherits(object, "error")) return(fail("rds_error", conditionMessage(object)))

  data <- as.data.table(object$merged_data)
  qtl_id <- as.character(object$qtl_id)
  qtl_n <- unname(QTL_N[qtl_id])
  if (is.na(qtl_n)) return(fail("unknown_qtl_sample_size", qtl_id))

  # historical SuSiE max PP4; NA when the CSV is absent (SCAF1) or has no PP.H4.abf
  old_path <- file.path(SUSIE_DIR, comparison_file)
  old_tab <- if (file.exists(old_path)) {
    tryCatch(fread(old_path, showProgress = FALSE), error = function(e) NULL)
  }
  old_pp4 <- if (!is.null(old_tab) && "PP.H4.abf" %in% names(old_tab) &&
                 nrow(old_tab) > 0L && any(is.finite(old_tab[["PP.H4.abf"]]))) {
    max(old_tab[["PP.H4.abf"]], na.rm = TRUE)
  } else {
    NA_real_
  }

  base_record[, `:=`(
    gene_symbol = as.character(object$gene_symbol),
    gwas_id = as.character(object$gwas_id),
    qtl_id = qtl_id,
    n_input_rows = nrow(data),
    n_qtl_samples = qtl_n,
    old_max_pp4 = old_pp4
  )]

  required <- c(
    "BETA.gwas", "SE.gwas", "BETA.qtl", "SE.qtl", "N", "snp", "variant_id"
  )
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    return(fail("missing_columns", paste(missing_columns, collapse = ",")))
  }
  data <- data[
    is.finite(BETA.gwas) & is.finite(SE.gwas) &
      is.finite(BETA.qtl) & is.finite(SE.qtl) & is.finite(N) &
      SE.gwas > 0 & SE.qtl > 0
  ]
  data <- data[!duplicated(snp) & !duplicated(variant_id)]
  gwas_ld <- tryCatch(
    get_ld_matrix(data$snp, object$plink_bfile, "plink"),
    error = function(e) e
  )
  if (inherits(gwas_ld, "error") || is.null(gwas_ld$R)) {
    message_text <- if (inherits(gwas_ld, "error")) conditionMessage(gwas_ld) else "NULL GWAS LD matrix"
    return(fail("ld_error", message_text))
  }
  common <- intersect(data$snp, colnames(gwas_ld$R))
  base_record[, n_common_snps := length(common)]
  if (length(common) < 5L) return(fail("too_few_common_snps"))

  # stage QTL LD: PLINK2 allele dosages from the QTL PGEN, coded to QTL EA/NEA
  qtl_pfile <- unname(QTL_PFILE[qtl_id])
  qtl_rows <- match(common, data$snp)
  qtl_variant_ids <- data$variant_id[qtl_rows]
  qtl_effect_alleles <- data$EA.qtl[qtl_rows]
  qtl_non_effect_alleles <- data$NEA.qtl[qtl_rows]
  ld_prefix <- tempfile("susie_qtl_ld_", tmpdir = tempdir())
  ld_extract <- paste0(ld_prefix, ".extract")
  ld_log <- paste0(ld_prefix, ".log.capture")
  qtl_ld <- tryCatch({
    variant_ids <- unique(as.character(qtl_variant_ids))
    variant_ids <- variant_ids[!is.na(variant_ids) & nzchar(variant_ids)]
    if (length(qtl_effect_alleles) != length(variant_ids) ||
        length(qtl_non_effect_alleles) != length(variant_ids)) {
      stop("Stage LD allele vectors do not match variant IDs")
    }
    fwrite(data.table(variant_id = variant_ids), ld_extract,
           col.names = FALSE, quote = FALSE)
    plink_status <- system2(
      PLINK2_BIN,
      args = c(
        "--pfile", qtl_pfile,
        "--extract", ld_extract,
        "--export", "A",
        "--allow-extra-chr",
        "--threads", "1",
        "--out", ld_prefix
      ),
      stdout = ld_log,
      stderr = ld_log
    )
    if (!identical(plink_status, 0L)) {
      log_tail <- if (file.exists(ld_log)) {
        paste(tail(readLines(ld_log, warn = FALSE), 20L), collapse = " | ")
      } else ""
      stop(sprintf("PLINK2 stage LD failed (exit=%s): %s", plink_status, log_tail))
    }
    dosage <- fread(paste0(ld_prefix, ".raw"), showProgress = FALSE)
    raw_names <- names(dosage)[-(1:6)]
    raw_ids <- sub("_[^_]+$", "", raw_names)
    raw_alleles <- toupper(sub("^.*_", "", raw_names))
    dosage_matrix <- as.matrix(dosage[, -(1:6), with = FALSE])
    storage.mode(dosage_matrix) <- "double"
    raw_index <- match(variant_ids, raw_ids)
    if (anyNA(raw_index) || anyDuplicated(raw_ids)) {
      stop("Stage dosage variant IDs do not match requested IDs for ", qtl_pfile)
    }
    dosage_matrix <- dosage_matrix[, raw_index, drop = FALSE]
    observed_alleles <- raw_alleles[raw_index]
    effect_alleles <- toupper(as.character(qtl_effect_alleles))
    non_effect_alleles <- toupper(as.character(qtl_non_effect_alleles))
    direct <- observed_alleles == effect_alleles
    flipped <- observed_alleles == non_effect_alleles
    if (any(!direct & !flipped)) {
      stop("Stage dosage allele coding does not match QTL EA/NEA for ", qtl_pfile)
    }
    if (any(flipped)) dosage_matrix[, flipped] <-
      2 - dosage_matrix[, flipped, drop = FALSE]
    R <- cor(dosage_matrix, use = "pairwise.complete.obs")
    dimnames(R) <- list(variant_ids, variant_ids)
    list(R = R, source = paste0("plink2_pgen_stage_effect_allele_dosage:", qtl_pfile))
  }, error = function(e) e, finally = unlink(c(
    ld_extract, ld_log, paste0(ld_prefix, ".log"),
    paste0(ld_prefix, ".unphased.vcor1"),
    paste0(ld_prefix, ".unphased.vcor1.vars"),
    paste0(ld_prefix, ".raw"), paste0(ld_prefix, ".export.log")
  )))
  if (inherits(qtl_ld, "error")) {
    return(fail("qtl_ld_error", conditionMessage(qtl_ld)))
  }
  qtl_keep <- qtl_variant_ids %in% rownames(qtl_ld$R)
  common <- common[qtl_keep]
  qtl_variant_ids <- qtl_variant_ids[qtl_keep]
  base_record[, `:=`(
    n_common_snps = length(common),
    gwas_ld_source = paste0("1000G_reference:", object$plink_bfile),
    qtl_ld_source = qtl_ld$source
  )]
  if (length(common) < 5L) return(fail("too_few_stage_common_snps"))

  data <- data[match(common, snp)]
  R_gwas <- gwas_ld$R[common, common, drop = FALSE]
  R_qtl <- qtl_ld$R[qtl_variant_ids, qtl_variant_ids, drop = FALSE]
  n_gwas <- as.numeric(data$N[1])

  set.seed(GLOBAL_SEED + 1L)
  fit_gwas <- tryCatch(
    susie_rss(
      bhat = data$BETA.gwas, shat = data$SE.gwas, R = R_gwas, n = n_gwas,
      L = L_VALUE, max_iter = MAX_ITER, verbose = FALSE
    ),
    error = function(e) e
  )
  set.seed(GLOBAL_SEED + 2L)
  fit_qtl <- tryCatch(
    susie_rss(
      bhat = data$BETA.qtl, shat = data$SE.qtl, R = R_qtl, n = qtl_n,
      L = L_VALUE, max_iter = MAX_ITER, verbose = FALSE
    ),
    error = function(e) e
  )

  gwas_error <- inherits(fit_gwas, "error")
  qtl_error <- inherits(fit_qtl, "error")
  gwas_converged <- !gwas_error && isTRUE(fit_gwas$converged)
  qtl_converged <- !qtl_error && isTRUE(fit_qtl$converged)
  gwas_cs_count <- if (gwas_error || is.null(fit_gwas$sets$cs)) {
    NA_integer_
  } else {
    length(fit_gwas$sets$cs)
  }
  qtl_cs_count <- if (qtl_error || is.null(fit_qtl$sets$cs)) {
    NA_integer_
  } else {
    length(fit_qtl$sets$cs)
  }
  base_record[, `:=`(
    gwas_converged = gwas_converged,
    qtl_converged = qtl_converged,
    comparison_converged = gwas_converged && qtl_converged,
    n_gwas_credible_sets = gwas_cs_count,
    n_qtl_credible_sets = qtl_cs_count,
    gwas_credible_set_status = if (gwas_error) "fit_error" else
      if (gwas_cs_count > 0L) "present" else "none",
    qtl_credible_set_status = if (qtl_error) "fit_error" else
      if (qtl_cs_count > 0L) "present" else "none",
    gwas_status = if (gwas_error) "error" else "complete",
    qtl_status = if (qtl_error) "error" else "complete",
    gwas_error = if (gwas_error) conditionMessage(fit_gwas) else NA_character_,
    qtl_error = if (qtl_error) conditionMessage(fit_qtl) else NA_character_
  )]

  diagnostics <- rbindlist(list(
    if (!gwas_error) diagnose_fit(
      data$BETA.gwas / data$SE.gwas, R_gwas, n_gwas, "GWAS", common, gwas_converged
    ) else data.table(
      side = "GWAS", s = NA_real_, n_kriging_outliers = NA_integer_,
      kriging_outlier_snps = NA_character_, diagnostic_status = "fit_error",
      diagnostic_error = conditionMessage(fit_gwas)
    ),
    if (!qtl_error) diagnose_fit(
      data$BETA.qtl / data$SE.qtl, R_qtl, qtl_n, "QTL", qtl_variant_ids, qtl_converged
    ) else data.table(
      side = "QTL", s = NA_real_, n_kriging_outliers = NA_integer_,
      kriging_outlier_snps = NA_character_, diagnostic_status = "fit_error",
      diagnostic_error = conditionMessage(fit_qtl)
    )
  ), fill = TRUE)
  if (nrow(diagnostics)) {
    diagnostics[, `:=`(
      comparison_file = comparison_file,
      rds_file = basename(rds_file),
      analysis_scope = scope,
      gene_symbol = as.character(object$gene_symbol),
      n_common_snps = length(common),
      outlier_definition = "kriging_rss conditional_dist: logLR > 2 and abs(z) > 2"
    )]
    setcolorder(diagnostics, c(
      "comparison_file", "rds_file", "analysis_scope", "gene_symbol", "side",
      "n_common_snps", "s", "n_kriging_outliers", "kriging_outlier_snps",
      "diagnostic_status", "diagnostic_error", "outlier_definition"
    ))
  }

  # coloc.bf_bf on the credible-set Bayes factors of the two SuSiE fits
  coloc_rows <- empty_coloc_row(comparison_file, rds_file, scope, "pending")
  if (!gwas_error && !qtl_error) {
    coloc_result <- tryCatch({
      cs1 <- fit_gwas$sets
      cs2 <- fit_qtl$sets
      if (is.null(cs1$cs) || is.null(cs2$cs) ||
          length(cs1$cs) == 0L || length(cs2$cs) == 0L) {
        list(summary = NULL, status = "empty_credible_set")
      } else {
        bf1 <- fit_gwas$lbf_variable[cs1$cs_index, , drop = FALSE]
        bf2 <- fit_qtl$lbf_variable[cs2$cs_index, , drop = FALSE]
        ret <- coloc::coloc.bf_bf(bf1, bf2)
        # coloc 5.2.3 needs a data.table summary before := ; coloc.bf_bf may return a data.frame
        if (!is.null(ret$summary)) {
          ret$summary <- as.data.table(ret$summary)
          ret$summary[, idx1 := cs1$cs_index[idx1]]
          ret$summary[, idx2 := cs2$cs_index[idx2]]
        }
        ret$status <- "complete"
        ret
      }
    }, error = function(e) e)
    if (inherits(coloc_result, "error")) {
      base_record[, `:=`(
        status = "coloc_error", error_message = conditionMessage(coloc_result),
        new_max_pp4 = NA_real_
      )]
    } else if (identical(coloc_result$status, "empty_credible_set")) {
      coloc_rows[, writeback_status := "empty_credible_set"]
      base_record[, `:=`(
        status = "empty_credible_set", error_message = NA_character_,
        new_max_pp4 = NA_real_, n_coloc_signal_pairs = 0L
      )]
    } else if (!is.null(coloc_result$summary) && nrow(coloc_result$summary)) {
      coloc_rows <- as.data.table(coloc_result$summary, keep.rownames = FALSE)
      coloc_rows[, `:=`(
        comparison_file = comparison_file,
        rds_file = basename(rds_file),
        analysis_scope = scope
      )]
      setcolorder(coloc_rows, c("comparison_file", "rds_file", "analysis_scope"))
      pp4 <- coloc_rows[["PP.H4.abf"]]
      base_record[, `:=`(
        status = "complete",
        error_message = NA_character_,
        new_max_pp4 = if (any(is.finite(pp4))) max(pp4, na.rm = TRUE) else NA_real_,
        n_coloc_signal_pairs = nrow(coloc_rows)
      )]
    } else {
      coloc_rows[, writeback_status := "empty_coloc_summary"]
      base_record[, `:=`(
        status = "empty_coloc_summary", error_message = NA_character_,
        new_max_pp4 = NA_real_, n_coloc_signal_pairs = 0L
      )]
    }
  } else {
    coloc_rows[, writeback_status := "fit_error"]
    base_record[, `:=`(
      status = "fit_error", error_message = NA_character_, new_max_pp4 = NA_real_
    )]
  }

  message(sprintf(
    "[SuSiE max_iter=1000] %d/%d %s GWAS=%s QTL=%s PP4=%s",
    task_index, length(tasks), comparison_file, gwas_converged, qtl_converged,
    format(base_record$new_max_pp4, digits = 4)
  ))
  list(comparison = base_record, diagnostics = diagnostics, coloc = coloc_rows)
}

# --- inputs ---
susie_files <- sort(list.files(SUSIE_DIR, pattern = "_susie[.]csv$", full.names = FALSE))
if (length(susie_files) != 125L) {
  stop(sprintf("Expected 125 existing SuSiE comparisons, found %d", length(susie_files)))
}
tasks <- lapply(susie_files, function(comparison_file) {
  list(
    comparison_file = comparison_file,
    rds_file = file.path(RDS_DIR, sub("_susie[.]csv$", "_coloc.rds", comparison_file)),
    analysis_scope = "primary_existing_susie_125"
  )
})
scaf1_rds <- list.files(RDS_DIR, pattern = "_SCAF1_coloc[.]rds$", full.names = TRUE)
if (length(scaf1_rds) != 1L) {
  stop(sprintf("Expected one supplemental SCAF1 RDS, found %d", length(scaf1_rds)))
}
tasks[[length(tasks) + 1L]] <- list(
  comparison_file = sub("_coloc[.]rds$", "_susie.csv", basename(scaf1_rds)),
  rds_file = scaf1_rds,
  analysis_scope = "supplemental_requested_SCAF1"
)

# --- fit ---
results <- if (.Platform$OS.type == "unix" && N_CORES > 1L) {
  mclapply(
    seq_along(tasks), function(i) fit_one(tasks[[i]], i),
    mc.cores = N_CORES, mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  lapply(seq_along(tasks), function(i) fit_one(tasks[[i]], i))
}

# --- collect ---
comparison <- rbindlist(lapply(results, `[[`, "comparison"), fill = TRUE)
diagnostics <- rbindlist(lapply(results, `[[`, "diagnostics"), fill = TRUE)
coloc_rows <- rbindlist(lapply(results, `[[`, "coloc"), fill = TRUE)

old_audit <- fread(OLD_AUDIT_FILE, showProgress = FALSE)
old_status <- old_audit[, .(
  comparison_file,
  old_gwas_converged = gwas_converged,
  old_qtl_converged = qtl_converged,
  old_comparison_converged = !comparison_unconverged
)]
comparison <- merge(comparison, old_status, by = "comparison_file", all.x = TRUE)
setorder(comparison, analysis_scope, comparison_file)
setorder(diagnostics, analysis_scope, comparison_file, side)
setorder(coloc_rows, analysis_scope, comparison_file, idx1, idx2)

# --- summary and audit tables ---
primary <- comparison[analysis_scope == "primary_existing_susie_125"]
summary <- data.table(
  metric = c(
    "n_comparisons", "n_converged_comparisons", "n_gwas_converged",
    "n_qtl_converged", "n_pp4_ge_0.7"
  ),
  definition = c(
    "Existing EasyColoc comparisons with a historical *_susie.csv",
    "Both GWAS and QTL susie_rss fits have $converged == TRUE",
    "GWAS susie_rss fit has $converged == TRUE",
    "QTL susie_rss fit has $converged == TRUE",
    "Comparison-level maximum coloc.susie PP.H4.abf >= 0.7"
  ),
  old_value = c(
    nrow(primary),
    sum(primary$old_comparison_converged %in% TRUE),
    sum(primary$old_gwas_converged %in% TRUE),
    sum(primary$old_qtl_converged %in% TRUE),
    sum(primary$old_max_pp4 >= 0.7, na.rm = TRUE)
  ),
  new_value = c(
    nrow(primary),
    sum(primary$comparison_converged %in% TRUE),
    sum(primary$gwas_converged %in% TRUE),
    sum(primary$qtl_converged %in% TRUE),
    sum(primary$new_max_pp4 >= 0.7, na.rm = TRUE)
  ),
  analysis_scope = "primary_existing_susie_125"
)

highlighted <- comparison[gene_symbol %in% HIGHLIGHT_GENES, .(
  gene_symbol, comparison_file, rds_file, analysis_scope, qtl_id, gwas_id,
  old_gwas_converged, old_qtl_converged, old_comparison_converged,
  gwas_converged, qtl_converged, comparison_converged,
  old_max_pp4, new_max_pp4,
  pp4_definition = "comparison-level maximum coloc.susie PP.H4.abf across signal pairs"
)]
diagnostic_s <- dcast(
  diagnostics[, .(comparison_file, side, s)],
  comparison_file ~ side, value.var = "s"
)
setnames(diagnostic_s, c("GWAS", "QTL"), c("gwas_s", "qtl_s"), skip_absent = TRUE)
highlighted <- merge(highlighted, diagnostic_s, by = "comparison_file", all.x = TRUE)
missing_highlights <- setdiff(HIGHLIGHT_GENES, highlighted$gene_symbol)
if (length(missing_highlights)) {
  stop(sprintf("Missing highlighted genes: %s", paste(missing_highlights, collapse = ", ")))
}
highlighted[, gene_order := match(gene_symbol, HIGHLIGHT_GENES)]
setorder(highlighted, gene_order, comparison_file)
highlighted[, gene_order := NULL]

fwrite(comparison, output_files[["comparison"]], sep = "\t", na = "NA")
fwrite(diagnostics, output_files[["diagnostics"]], sep = "\t", na = "NA")
fwrite(coloc_rows, output_files[["coloc"]], sep = "\t", na = "NA")
fwrite(summary, output_files[["summary"]], sep = "\t", na = "NA")
fwrite(highlighted, output_files[["highlighted"]], sep = "\t", na = "NA")

# --- writeback ---
# Replace the canonical per-comparison SuSiE tables and their aggregate in place.
# The RDS objects keep the ABF fields and receive the new SuSiE summary and PP4.
canonical_columns <- c(
  "nsnps", "hit1", "hit2", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
  "PP.H3.abf", "PP.H4.abf", "idx1", "idx2"
)
canonical_all <- list()
for (comparison_file in susie_files) {
  row_keep <- coloc_rows[["analysis_scope"]] == "primary_existing_susie_125" &
    coloc_rows[["comparison_file"]] == comparison_file
  rows <- coloc_rows[row_keep, canonical_columns, with = FALSE]
  fwrite(rows, file.path(SUSIE_DIR, comparison_file), sep = ",", na = "NA")
  if (nrow(rows)) {
    rows[, source_file := comparison_file]
    canonical_all[[length(canonical_all) + 1L]] <- rows
  }

  rds_path <- file.path(RDS_DIR, sub("_susie[.]csv$", "_coloc.rds", comparison_file))
  object <- readRDS(rds_path)
  object$susie_summary <- as.data.frame(rows)
  object$susie_best_pp4 <- if (nrow(rows)) max(rows$PP.H4.abf, na.rm = TRUE) else NA_real_
  saveRDS(object, rds_path, compress = TRUE)
}
canonical_all_dt <- if (length(canonical_all)) {
  rbindlist(canonical_all, fill = TRUE)
} else {
  data.table()
}
fwrite(canonical_all_dt, file.path(DATA_DIR, "all_susie_results.csv"),
       sep = ",", na = "NA")

print(summary)
