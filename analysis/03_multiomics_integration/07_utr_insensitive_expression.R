#!/usr/bin/env Rscript
# UTR-insensitive check: ComBat-corrected TPM and STAR junction-count expression.
# Then mediation reruns on coloc pairs (PP_H4 > 0.8) and cis APA-protein / APA-mRNA models.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(parallel)
  library(qvalue)
  library(rtracklayer)
  library(splines)
  library(sva)
})
setDTthreads(1L)

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(apa_root))
out_dir <- file.path(apa_root, "data/03_omics")
tpm_path <- file.path(apa_root, "raw/all.912.expr.tpm")
gtf_path <- file.path(apa_root, "ref/gencode.v40.annotation.gtf.gz")
meta_path <- file.path(apa_root, "data/04_qtl/CHB.all.meta")
current_expr_path <- file.path(apa_root, "data/01_quant/expr_after_combat.tsv")
protein_dir <- file.path(apa_root, "data/03_omics")

# cluster-only inputs (no copy under APA_ROOT)
junction_dir <- file.path(Sys.getenv("APA_HPC_HOME"), "project/work/ati/promoter/results/junction")
eqtl_dir <- file.path(Sys.getenv("APA_HPC_HOME"), "project/work/apa/11_eqtl")
expression_qtl_dir <- file.path(Sys.getenv("APA_HPC_HOME"), "project/work/eqtl")
apa_qtl_dir <- file.path(Sys.getenv("APA_HPC_HOME"), "project/work/apa/01_QTL")

n_boot <- 2000L
n_cores <- min(16L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16")))
set.seed(2026)

# --- helpers ---

strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

combat_expression <- function(x, meta) {
  batch <- factor(meta$batch)
  mod <- model.matrix(~ LifeStage, data = meta)
  within_var <- sapply(levels(batch), function(b) {
    apply(x[, batch == b, drop = FALSE], 1L, var)
  })
  keep <- !apply(within_var < 1e-8 | !is.finite(within_var), 1L, any)
  x <- x[keep, , drop = FALSE]
  corrected <- suppressWarnings(ComBat(
    dat = x, batch = batch, mod = mod, ref.batch = "batch9",
    par.prior = TRUE, prior.plots = FALSE
  ))
  stopifnot(all(is.finite(corrected)))
  corrected
}

write_matrix <- function(x, path, id_name = "gene_id") {
  fwrite(as.data.table(x, keep.rownames = id_name), path, sep = "\t", compress = "gzip")
}

read_covariates <- function(path, samples) {
  cov <- fread(path)
  cov_names <- cov[[1L]]
  x <- t(as.matrix(cov[, ..samples]))
  colnames(x) <- cov_names
  x
}

residualise_matrix <- function(y, covariates) {
  is_geno_pc <- grepl("^genoPC", colnames(covariates), ignore.case = TRUE)
  design <- cbind(Intercept = 1, covariates[, !is_geno_pc, drop = FALSE])
  t(qr.resid(qr(design), t(y)))
}

ols_fit <- function(x, y) {
  xtx <- crossprod(x)
  ch <- tryCatch(chol(xtx), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  inverse <- chol2inv(ch)
  coefficient <- as.vector(inverse %*% crossprod(x, y))
  residual <- y - x %*% coefficient
  df <- nrow(x) - ncol(x)
  if (df <= 1L) return(NULL)
  se <- sqrt(sum(residual^2) / df * diag(inverse))
  list(coef = coefficient, se = se, p = 2 * pt(-abs(coefficient / se), df))
}

ols_coef <- function(x, y) {
  fit <- tryCatch(.lm.fit(x, y), error = function(e) NULL)
  if (is.null(fit) || anyNA(fit$coefficients)) return(NULL)
  fit$coefficients
}

coefficient_test <- function(y, x, coefficient_index) {
  fit <- lm.fit(x, y)
  if (fit$rank < ncol(x)) return(c(t = NA_real_, p = NA_real_))
  residual_df <- nrow(x) - fit$rank
  covariance <- chol2inv(qr.R(fit$qr)) * sum(fit$residuals^2) / residual_df
  se <- sqrt(covariance[coefficient_index, coefficient_index])
  statistic <- fit$coefficients[coefficient_index] / se
  c(t = unname(statistic), p = unname(2 * pt(abs(statistic), residual_df, lower.tail = FALSE)))
}

run_mediation <- function(expression, measure) {
  stage_results <- lapply(names(stage_samples), function(stage_label) {
    stage_lower <- tolower(stage_label)
    pairs <- fread(file.path(eqtl_dir, "results", sprintf("coloc_%s_apa_eqtl.tsv", stage_lower)))[PP_H4 > 0.8]
    apa_file <- file.path(apa_qtl_dir, "pre_input", sprintf("phenotype.%s.GeneRange.bed.gz", stage_label))
    expression_cov_file <- file.path(expression_qtl_dir, "pre_input", sprintf("covariates_for_qtl.%s.txt", stage_label))
    apa_cov_file <- file.path(apa_qtl_dir, "pre_input", sprintf("covariates_for_qtl.%s.txt", stage_label))
    apa <- fread(apa_file)
    apa_ids <- apa[[4L]]
    apa_matrix <- as.matrix(apa[, -seq_len(4L)])
    rownames(apa_matrix) <- apa_ids
    samples <- Reduce(intersect, list(
      colnames(genotype), colnames(apa_matrix), colnames(expression),
      names(fread(expression_cov_file, nrows = 0L))[-1L],
      names(fread(apa_cov_file, nrows = 0L))[-1L]
    ))
    expression_covariates <- read_covariates(expression_cov_file, samples)
    apa_covariates <- read_covariates(apa_cov_file, samples)
    apa_residual <- residualise_matrix(apa_matrix[, samples, drop = FALSE], apa_covariates)
    expression_residual <- residualise_matrix(expression[, samples, drop = FALSE], expression_covariates)
    geno_pc <- expression_covariates[, grepl("^genoPC", colnames(expression_covariates), ignore.case = TRUE), drop = FALSE]
    gene_index <- match(strip_version(pairs$eqtl_pheno), strip_version(rownames(expression_residual)))
    pair_index <- match(pairs$apa_pheno, rownames(apa_residual))
    variant_index <- match(pairs$h4_snp, rownames(genotype))
    pairs[, `:=`(
      row = .I,
      expression_available = !is.na(gene_index),
      apa_available = !is.na(pair_index),
      genotype_available = !is.na(variant_index)
    )]
    test_index <- which(!is.na(gene_index) & !is.na(pair_index) & !is.na(variant_index))
    tested <- mclapply(test_index, function(i) {
      x <- genotype[variant_index[i], samples]
      m <- apa_residual[pair_index[i], ]
      y <- expression_residual[gene_index[i], ]
      covariates <- geno_pc
      ok <- is.finite(x) & is.finite(m) & is.finite(y) & apply(is.finite(covariates), 1L, all)
      x <- x[ok]
      m <- m[ok]
      y <- y[ok]
      covariates <- covariates[ok, , drop = FALSE]
      n <- length(x)
      alt_count <- sum(x)
      maf <- min(alt_count, 2 * n - alt_count) / (2 * n)
      n_carrier <- if (alt_count <= n) sum(x > 0.5) else sum(x < 1.5)
      if (n < 50L || sd(x) == 0 || sd(m) == 0 || sd(y) == 0 || n_carrier < 10L) return(NULL)

      d_x <- cbind(1, x, covariates)
      d_xm <- cbind(1, x, m, covariates)
      d_xy <- cbind(1, x, y, covariates)
      f_m <- ols_fit(d_x, m)
      f_ym <- ols_fit(d_xm, y)
      f_y <- ols_fit(d_x, y)
      f_mrev <- ols_fit(d_xy, m)
      if (any(vapply(list(f_m, f_ym, f_y, f_mrev), is.null, logical(1L)))) return(NULL)

      a <- f_m$coef[2L]
      b <- f_ym$coef[3L]
      c_total <- f_y$coef[2L]
      c_prime <- f_ym$coef[2L]
      b_rev <- f_mrev$coef[3L]
      d_prime <- f_mrev$coef[2L]
      ab_fwd <- a * b
      ab_rev <- c_total * b_rev
      se_ab_fwd <- sqrt(a^2 * f_ym$se[3L]^2 + b^2 * f_m$se[2L]^2 + f_m$se[2L]^2 * f_ym$se[3L]^2)
      se_ab_rev <- sqrt(c_total^2 * f_mrev$se[3L]^2 + b_rev^2 * f_y$se[2L]^2 + f_y$se[2L]^2 * f_mrev$se[3L]^2)
      boot <- matrix(NA_real_, nrow = n_boot, ncol = 3L)
      for (iteration in seq_len(n_boot)) {
        index <- sample.int(n, n, replace = TRUE)
        xb <- x[index]
        if (sd(xb) == 0) next
        cb <- covariates[index, , drop = FALSE]
        mb <- m[index]
        yb <- y[index]
        k_m <- ols_coef(cbind(1, xb, cb), mb)
        k_ym <- ols_coef(cbind(1, xb, mb, cb), yb)
        k_y <- ols_coef(cbind(1, xb, cb), yb)
        k_mr <- ols_coef(cbind(1, xb, yb, cb), mb)
        if (any(vapply(list(k_m, k_ym, k_y, k_mr), is.null, logical(1L)))) next
        boot[iteration, 1L] <- k_m[2L] * k_ym[3L]
        boot[iteration, 2L] <- k_y[2L] * k_mr[3L]
        boot[iteration, 3L] <- if (abs(k_y[2L]) > 1e-10) boot[iteration, 1L] / k_y[2L] else NA_real_
      }
      ci <- quantile(boot[, 1L], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      statistic <- data.table(
        n = n, maf = maf, a = a, p_a = f_m$p[2L], b = b, p_b = f_ym$p[3L],
        c_total = c_total, p_c = f_y$p[2L], c_prime = c_prime, p_c_prime = f_ym$p[2L],
        d_prime = d_prime, p_d_prime = f_mrev$p[2L],
        ab_fwd = ab_fwd, p_ab_fwd = 2 * pnorm(-abs(ab_fwd / se_ab_fwd)),
        ab_fwd_lo = ci[1L], ab_fwd_hi = ci[2L],
        ab_rev = ab_rev, p_ab_rev = 2 * pnorm(-abs(ab_rev / se_ab_rev)),
        prop_med = if (f_y$p[2L] < 0.05 && abs(c_total) > 1e-10) ab_fwd / c_total else NA_real_
      )
      cbind(data.table(row = i), statistic)
    }, mc.cores = n_cores, mc.set.seed = TRUE)
    tested <- rbindlist(tested, fill = TRUE)
    out <- merge(pairs, tested, by = "row", all.x = TRUE, sort = FALSE)
    setorder(out, row)
    out[, analysis_status := fcase(
      !expression_available, "expression unavailable after filtering",
      !apa_available, "APA phenotype unavailable",
      !genotype_available, "genotype unavailable",
      is.na(n), "model filter failed",
      default = "tested"
    )]
    out[, `:=`(stage = stage_lower, expression_measure = measure)]
    out
  })
  med <- rbindlist(stage_results, fill = TRUE)
  med[, fdr_ab_fwd := p.adjust(p_ab_fwd, "BH"), by = stage]
  med[, fdr_ab_rev := p.adjust(p_ab_rev, "BH"), by = stage]
  med[, fdr_c_prime := p.adjust(p_c_prime, "BH"), by = stage]
  med[, fdr_d_prime := p.adjust(p_d_prime, "BH"), by = stage]
  med[, mediation_call := fcase(
    analysis_status != "tested", "Not tested",
    fdr_ab_fwd < 0.05 & p_c_prime >= 0.05 & fdr_d_prime < 0.05, "APA -> EXP",
    fdr_ab_rev < 0.05 & p_d_prime >= 0.05 & fdr_c_prime < 0.05, "EXP -> APA",
    fdr_ab_fwd < 0.05 & fdr_c_prime < 0.05 & fdr_d_prime < 0.05, "Partial mediation",
    default = "No mediation"
  )]
  med[, apa_effect_on_expr := fcase(
    analysis_status != "tested", "Not tested",
    p_b >= 0.05, "not significant",
    b < 0, "longer 3'UTR -> lower expression",
    default = "longer 3'UTR -> higher expression"
  )]
  merge(med, alleles, by = "h4_snp", all.x = TRUE, sort = FALSE)
}

fit_cis <- function(expression, measure, include_protein) {
  expression <- expression[, samples, drop = FALSE]
  gene_index <- match(strip_version(universe$gencode_gene_base), strip_version(rownames(expression)))
  apa_index <- match(universe$site_id, rownames(apa))
  protein_index <- match(universe$protein_accession, rownames(protein))
  keep <- !is.na(gene_index) & !is.na(apa_index) & (!include_protein | !is.na(protein_index))
  u <- universe[keep]
  statistics <- t(vapply(which(keep), function(i) {
    a <- apa[apa_index[i], ]
    r <- expression[gene_index[i], ]
    rna_test <- coefficient_test(r, cbind(1, a, age_basis), 2L)
    if (include_protein) {
      protein_test <- coefficient_test(protein[protein_index[i], ], cbind(1, a, r, age_basis), 2L)
      c(t = protein_test["t"], p = protein_test["p"], t_rna = rna_test["t"], p_rna = rna_test["p"])
    } else {
      c(t = rna_test["t"], p = rna_test["p"])
    }
  }, numeric(if (include_protein) 4L else 2L)))
  statistics <- as.data.table(statistics)
  setnames(statistics, sub("\\..*$", "", names(statistics)))
  out <- cbind(u[, .(site_id, protein_accession, gene = gene_symbol, gencode_gene_base, cluster_id)], statistics)
  out[, `:=`(q = p.adjust(p, "BH"), expression_measure = measure)]
  if (include_protein) out[, q_rna := p.adjust(p_rna, "BH")]
  out
}

summarize_models <- function(x, direction) {
  rbindlist(lapply(c("overall", "C1", "C2"), function(stratum) {
    z <- if (stratum == "overall") x else x[cluster_id == as.integer(sub("C", "", stratum))]
    z[, .(
      n_pairs = .N,
      pi1 = if (any(is.finite(p))) 1 - qvalue::pi0est(p[is.finite(p)])$pi0 else NA_real_,
      n_q_lt_0_05 = sum(q < 0.05, na.rm = TRUE),
      positive_fraction_all = mean(t > 0, na.rm = TRUE),
      negative_fraction_all = mean(t < 0, na.rm = TRUE),
      positive_fraction_q_lt_0_05 = mean(t[q < 0.05] > 0, na.rm = TRUE),
      negative_fraction_q_lt_0_05 = mean(t[q < 0.05] < 0, na.rm = TRUE)
    ), by = expression_measure][, `:=`(stratum = stratum, model = direction)]
  }))
}

# --- retained samples and metadata ---

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
current_header <- names(fread(current_expr_path, nrows = 0L, check.names = FALSE))
current_samples <- current_header[-1L]
current_gene_ids <- fread(current_expr_path, select = 1L, header = TRUE)[[1L]]
meta <- fread(meta_path, data.table = FALSE)
rownames(meta) <- meta$sample
meta <- meta[current_samples, , drop = FALSE]
stopifnot(!anyNA(meta$batch))

# --- TPM expression ---

tpm <- fread(tpm_path, check.names = FALSE, showProgress = TRUE)
tpm_ids <- tpm[[1L]]
tpm[[1L]] <- NULL
tpm_samples <- names(tpm)
duplicate_20171125 <- which(tpm_samples == "XYA20171125")
if (length(duplicate_20171125) == 2L) tpm_samples[duplicate_20171125[2L]] <- "20171125C"
aliases <- c(
  XYA20180516B = "XYA20180516",
  XYA20210813B = "XYA20210813",
  XYA20230626B = "XYA20230626",
  XYA20231026B = "XYA20231026",
  XYA20240412A = "XYA20240412"
)
matched <- match(tpm_samples, names(aliases))
tpm_samples[!is.na(matched)] <- aliases[matched[!is.na(matched)]]
setnames(tpm, tpm_samples)
stopifnot(!anyDuplicated(names(tpm)), all(current_samples %chin% names(tpm)))
tpm_index <- match(strip_version(current_gene_ids), strip_version(tpm_ids))
stopifnot(!anyNA(tpm_index))
tpm_matrix <- as.matrix(tpm[tpm_index, ..current_samples])
rownames(tpm_matrix) <- current_gene_ids
expr_tpm <- combat_expression(log2(tpm_matrix + 1), meta)
write_matrix(expr_tpm, file.path(out_dir, "R1_expr_TPM.tsv.gz"))
rm(tpm, tpm_matrix)

# --- annotated non-terminal junction map ---

gtf <- rtracklayer::import(gtf_path)
exon <- gtf[gtf$type == "exon"]
exon <- data.table(
  chromosome = as.character(seqnames(exon)), start = start(exon), end = end(exon),
  strand = as.character(strand(exon)), gene_id = as.character(exon$gene_id),
  transcript_id = as.character(exon$transcript_id)
)
exon <- exon[chromosome %chin% paste0("chr", 1:22) & !is.na(gene_id) & !is.na(transcript_id)]
exon[, exon_rank := if (strand[1L] == "+") frank(start, ties.method = "first") else frank(-start, ties.method = "first"),
     by = transcript_id]
setorder(exon, transcript_id, exon_rank)
junction_map <- exon[, {
  if (.N >= 3L) {
    data.table(
      chromosome = chromosome[-.N],
      intron_start = pmin(end[-.N], end[-1L]) + 1L,
      intron_end = pmax(start[-.N], start[-1L]) - 1L,
      strand_code = if (strand[1L] == "+") 1L else 2L,
      gene_id = gene_id[1L],
      terminal = seq_len(.N - 1L) == .N - 1L
    )
  }
}, by = transcript_id]
junction_map <- junction_map[terminal == FALSE & intron_start <= intron_end]
junction_map[, key := paste(chromosome, intron_start, intron_end, strand_code, sep = ":")]
junction_map <- unique(junction_map[, .(key, gene_id)])
junction_map <- junction_map[strip_version(gene_id) %chin% strip_version(current_gene_ids)]
junction_genes <- unique(junction_map$gene_id)
rm(gtf, exon)

# --- STAR junction expression ---

junction_files <- file.path(junction_dir, paste0(current_samples, "_SJ.out.tab"))
stopifnot(all(file.exists(junction_files)))
junction_list <- mclapply(junction_files, function(path) {
  x <- fread(path, select = c(1L, 2L, 3L, 4L, 7L), col.names = c("chromosome", "intron_start", "intron_end", "strand_code", "unique_reads"))
  x[, key := paste(chromosome, intron_start, intron_end, strand_code, sep = ":")]
  x <- merge(x[key %chin% junction_map$key], junction_map, by = "key", allow.cartesian = TRUE)
  totals <- x[, .(unique_reads = sum(unique_reads)), by = gene_id]
  values <- numeric(length(junction_genes))
  values[match(totals$gene_id, junction_genes)] <- totals$unique_reads
  values
}, mc.cores = n_cores)
junction_counts <- do.call(cbind, junction_list)
dimnames(junction_counts) <- list(junction_genes, current_samples)
current_index <- match(strip_version(rownames(junction_counts)), strip_version(current_gene_ids))
keep <- !is.na(current_index) & rowSums(junction_counts >= 6) >= 0.20 * ncol(junction_counts)
junction_counts <- junction_counts[keep, , drop = FALSE]
rownames(junction_counts) <- current_gene_ids[current_index[keep]]
junction_counts <- junction_counts[!duplicated(rownames(junction_counts)), , drop = FALSE]
junction_dge <- calcNormFactors(DGEList(junction_counts), method = "TMM")
expr_junction <- combat_expression(cpm(junction_dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1), meta)
write_matrix(expr_junction, file.path(out_dir, "R1_expr_junction.tsv.gz"))
rm(junction_list, junction_counts, junction_dge)

# --- SJ sample coverage ---

stage_samples <- lapply(c(Prenatal = "Prenatal", Postnatal = "Postnatal"), function(stage) {
  names(fread(file.path(expression_qtl_dir, "pre_input", sprintf("phenotype.%s.bed.gz", stage)), nrows = 0L))[-seq_len(4L)]
})
qtl_samples <- unique(unlist(stage_samples, use.names = FALSE))
protein_manifest <- fread(file.path(protein_dir, "matched_sample_manifest.tsv"))
coverage <- rbindlist(list(
  data.table(cohort = "QTL", sample_id = qtl_samples),
  data.table(cohort = "protein", sample_id = protein_manifest$sample_id)
))
coverage[, sj_file := file.path(junction_dir, paste0(sample_id, "_SJ.out.tab"))]
coverage[, has_sj := file.exists(sj_file)]
fwrite(coverage, file.path(out_dir, "R1_SJ_sample_coverage.tsv"), sep = "\t")
fwrite(coverage[, .(n_samples = .N, n_with_sj = sum(has_sj)), by = cohort],
       file.path(out_dir, "R1_SJ_sample_coverage_summary.tsv"), sep = "\t")

# --- mediation reruns ---

traw <- fread(file.path(eqtl_dir, "coloc_target_geno.traw"))
genotype <- as.matrix(traw[, -(1:6)])
rownames(genotype) <- traw$SNP
colnames(genotype) <- sub("^[^_]*_", "", colnames(genotype))
flip <- traw$COUNTED != traw$ALT
genotype[flip, ] <- 2 - genotype[flip, ]
alleles <- data.table(h4_snp = traw$SNP, ref_allele = traw$COUNTED, alt_allele = traw$ALT)

med_tpm <- run_mediation(expr_tpm, "expr_TPM")
med_junction <- run_mediation(expr_junction, "junction")
fwrite(med_tpm, file.path(out_dir, "R1a_mediation_expr_TPM.tsv"), sep = "\t")
fwrite(med_junction, file.path(out_dir, "R1a_mediation_junction.tsv"), sep = "\t")
mediation_summary <- rbind(med_tpm, med_junction, fill = TRUE)[, .N,
  by = .(expression_measure, stage, mediation_call, apa_effect_on_expr)]
mediation_summary[, fraction := N / sum(N), by = .(expression_measure, stage)]
fwrite(mediation_summary, file.path(out_dir, "R1a_mediation_summary.tsv"), sep = "\t")

# --- 371-donor cis protein and APA-mRNA models ---

universe <- fread(file.path(protein_dir, "01_analysis_universe.tsv"))[primary_universe == TRUE]
apa <- readRDS(file.path(protein_dir, "apa_matched.rds"))
protein <- readRDS(file.path(protein_dir, "protein_matched.rds"))
current_rna <- readRDS(file.path(protein_dir, "rna_matched.rds"))
samples <- protein_manifest$sample_id
stopifnot(identical(colnames(apa), samples), identical(colnames(protein), samples), identical(colnames(current_rna), samples))
age_basis <- ns(log10(protein_manifest$age_years + 1), df = 2L)

protein_results <- rbindlist(list(
  fit_cis(expr_tpm, "expr_TPM", TRUE),
  fit_cis(expr_junction, "junction", TRUE)
), fill = TRUE)
rna_results <- rbindlist(list(
  fit_cis(current_rna, "current", FALSE),
  fit_cis(expr_tpm, "expr_TPM", FALSE),
  fit_cis(expr_junction, "junction", FALSE)
), fill = TRUE)
fwrite(protein_results[, .(site_id, protein_accession, gene, cluster_id, expression_measure, t, p, q)],
       file.path(out_dir, "R1b_cis_protein_models.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(rna_results[, .(site_id, protein_accession, gene, cluster_id, expression_measure, t, p, q)],
       file.path(out_dir, "R1c_cis_apa_mrna_models.tsv.gz"), sep = "\t", compress = "gzip")

model_summary <- rbind(
  summarize_models(protein_results, "protein ~ age spline + mRNA + PDUI"),
  summarize_models(rna_results, "mRNA ~ age spline + PDUI"),
  fill = TRUE
)
fwrite(model_summary, file.path(out_dir, "R1bc_cis_summary.tsv"), sep = "\t")

# --- unavailable exact null contract ---

null_contract <- data.table(
  requested_object = c("S1_pair_results_371", "same trans-gene null"),
  status = "missing",
  search_result = c(
    sprintf("No object or producer found; current traceable primary universe has %d pairs", nrow(universe)),
    "No pair-key table, producer, seed, or sampling rule found in workspace or Git history"
  ),
  consequence = "No substitute null was generated"
)
fwrite(null_contract, file.path(out_dir, "R1bc_missing_exact_contract.tsv"), sep = "\t")

cat(sprintf(
  "R1 complete: TPM genes=%d; junction genes=%d; mediation=%d/%d; cis universe=%d; QTL SJ=%d/%d; protein SJ=%d/%d\n",
  nrow(expr_tpm), nrow(expr_junction), nrow(med_tpm), nrow(med_junction), nrow(universe),
  coverage[cohort == "QTL", sum(has_sj)], coverage[cohort == "QTL", .N],
  coverage[cohort == "protein", sum(has_sj)], coverage[cohort == "protein", .N]
))
