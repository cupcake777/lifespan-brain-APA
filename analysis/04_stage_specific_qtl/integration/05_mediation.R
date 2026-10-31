# Mediation of significant coloc pairs (APA -> expression, expression -> APA) per life stage,
# with BH within stage and GO/KEGG enrichment of mediated genes. Writes data/04_qtl/04_mediation_*.
suppressMessages({
  library(data.table)
  library(parallel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# --- paths and settings ---
APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
D <- file.path(APA_ROOT, "data", "04_qtl")
PP_H4_MIN <- 0.8
N_BOOT <- 2000
FDR_CUT <- 0.05
CI_ALPHA <- 0.05
MIN_N <- 50
MIN_MINOR <- 10
SEED <- 2026
N_CORES <- max(1L, min(24L, suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))))
if (is.na(N_CORES)) N_CORES <- 8L
stage_key <- c(prenatal = "Prenatal", postnatal = "Postnatal")
call_levels <- c("APA -> EXP", "EXP -> APA", "Partial mediation", "No mediation")
set.seed(SEED)

# --- statistics (ols_fit, ols_coef, mediate_pair are also used by tests/qtl/test_mediation.R) ---
ols_fit <- function(X, y) {
  n <- nrow(X); p <- ncol(X)
  xtx <- crossprod(X)
  ch <- tryCatch(chol(xtx), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  xtx_inv <- chol2inv(ch)
  coef <- as.vector(xtx_inv %*% crossprod(X, y))
  resid <- y - X %*% coef
  df <- n - p
  if (df <= 1L) return(NULL)
  s2 <- sum(resid^2) / df
  se <- sqrt(s2 * diag(xtx_inv))
  list(coef = coef, se = se, df = df, p = 2 * pt(-abs(coef / se), df))
}

ols_coef <- function(X, y) {
  fit <- tryCatch(.lm.fit(X, y), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- fit$coefficients
  if (anyNA(cf)) return(NULL)
  cf
}

# residualise phenotypes on non-genotype covariates; genotype PCs returned separately
residualise_phenotypes <- function(bed_file, cov_file, geno_pc_pattern = "^genoPC") {
  bed <- fread(bed_file); cov <- fread(cov_file)
  cov_names <- cov[[1]]
  cov_samples <- setdiff(names(cov), names(cov)[1])
  bed_samples <- setdiff(names(bed), names(bed)[1:4])
  keep <- intersect(bed_samples, cov_samples)
  if (length(keep) < MIN_N) stop("Only ", length(keep), " shared samples for ", basename(bed_file))
  cov_mat <- t(as.matrix(cov[, keep, with = FALSE]))
  colnames(cov_mat) <- cov_names
  is_geno_pc <- grepl(geno_pc_pattern, colnames(cov_mat), ignore.case = TRUE)
  geno_pcs <- cov_mat[, is_geno_pc, drop = FALSE]
  other_cov <- cov_mat[, !is_geno_pc, drop = FALSE]
  pheno <- as.matrix(bed[, keep, with = FALSE])
  rownames(pheno) <- bed$pid
  design <- cbind(Intercept = 1, other_cov)
  resid <- t(qr.resid(qr(design), t(pheno)))
  dimnames(resid) <- list(rownames(pheno), keep)
  list(residuals = resid, geno_pcs = geno_pcs, samples = keep)
}

mediate_pair <- function(x, m, y, covar) {
  ok <- is.finite(x) & is.finite(m) & is.finite(y) & apply(is.finite(covar), 1, all)
  x <- x[ok]; m <- m[ok]; y <- y[ok]; covar <- covar[ok, , drop = FALSE]
  n <- length(x)
  alt_count <- sum(x)
  maf <- min(alt_count, 2 * n - alt_count) / (2 * n)
  n_carrier <- if (alt_count <= n) sum(x > 0.5) else sum(x < 1.5)

  if (n < MIN_N) return(list(status = "n_too_small", n = n))
  if (sd(x) == 0) return(list(status = "monomorphic", n = n))
  if (sd(m) == 0 || sd(y) == 0) return(list(status = "constant_phenotype", n = n))
  if (n_carrier < MIN_MINOR) return(list(status = "too_few_carriers", n = n))

  D_x <- cbind(1, x, covar)
  D_xm <- cbind(1, x, m, covar)
  D_xy <- cbind(1, x, y, covar)
  f_m <- ols_fit(D_x, m)
  f_ym <- ols_fit(D_xm, y)
  f_y <- ols_fit(D_x, y)
  f_mrev <- ols_fit(D_xy, m)
  if (any(vapply(list(f_m, f_ym, f_y, f_mrev), is.null, logical(1)))) return(list(status = "singular_design", n = n))

  a <- f_m$coef[2]; se_a <- f_m$se[2]; p_a <- f_m$p[2]
  cprime <- f_ym$coef[2]; se_cp <- f_ym$se[2]; p_cp <- f_ym$p[2]
  b <- f_ym$coef[3]; se_b <- f_ym$se[3]; p_b <- f_ym$p[3]
  ctot <- f_y$coef[2]; se_c <- f_y$se[2]; p_c <- f_y$p[2]
  dprime <- f_mrev$coef[2]; se_dp <- f_mrev$se[2]; p_dp <- f_mrev$p[2]
  b_rev <- f_mrev$coef[3]; se_br <- f_mrev$se[3]; p_br <- f_mrev$p[3]

  ab_fwd <- a * b
  ab_rev <- ctot * b_rev
  se_ab_fwd <- sqrt(a^2 * se_b^2 + b^2 * se_a^2 + se_a^2 * se_b^2)
  se_ab_rev <- sqrt(ctot^2 * se_br^2 + b_rev^2 * se_c^2 + se_c^2 * se_br^2)
  p_ab_fwd <- 2 * pnorm(-abs(ab_fwd / se_ab_fwd))
  p_ab_rev <- 2 * pnorm(-abs(ab_rev / se_ab_rev))

  # nonparametric bootstrap of the forward indirect effect
  boot <- matrix(NA_real_, nrow = N_BOOT, ncol = 3L)
  for (bi in seq_len(N_BOOT)) {
    idx <- sample.int(n, n, replace = TRUE)
    xb <- x[idx]
    if (sd(xb) == 0) next
    cb <- covar[idx, , drop = FALSE]; mb <- m[idx]; yb <- y[idx]
    k_m <- ols_coef(cbind(1, xb, cb), mb)
    k_ym <- ols_coef(cbind(1, xb, mb, cb), yb)
    k_y <- ols_coef(cbind(1, xb, cb), yb)
    k_mr <- ols_coef(cbind(1, xb, yb, cb), mb)
    if (is.null(k_m) || is.null(k_ym) || is.null(k_y) || is.null(k_mr)) next
    boot[bi, 1] <- k_m[2] * k_ym[3]
    boot[bi, 2] <- k_y[2] * k_mr[3]
    boot[bi, 3] <- if (abs(k_y[2]) > 1e-10) (k_m[2] * k_ym[3]) / k_y[2] else NA_real_
  }
  v <- boot[, 1]; v <- v[is.finite(v)]
  ci_fwd <- if (length(v) < 100L) c(NA_real_, NA_real_) else unname(quantile(v, c(0.025, 0.975)))
  prop_med <- if (p_c < 0.05 && abs(ctot) > 1e-10) ab_fwd / ctot else NA_real_

  list(status = "ok", n = n, maf = maf, n_carrier = n_carrier,
       a = a, se_a = se_a, p_a = p_a, b = b, se_b = se_b, p_b = p_b,
       c_total = ctot, se_c = se_c, p_c = p_c,
       c_prime = cprime, se_c_prime = se_cp, p_c_prime = p_cp,
       b_rev = b_rev, se_b_rev = se_br, p_b_rev = p_br,
       d_prime = dprime, se_d_prime = se_dp, p_d_prime = p_dp,
       ab_fwd = ab_fwd, se_ab_fwd = se_ab_fwd, p_ab_fwd = p_ab_fwd,
       ab_fwd_lo = ci_fwd[1], ab_fwd_hi = ci_fwd[2],
       ab_rev = ab_rev, se_ab_rev = se_ab_rev, p_ab_rev = p_ab_rev,
       prop_med = prop_med)
}

# --- genotype (coloc target SNPs, dosage on the alt allele) ---
traw <- fread(file.path(D, "coloc_target_geno.traw"))
geno_all <- as.matrix(traw[, -(1:6), with = FALSE])
rownames(geno_all) <- traw$SNP
colnames(geno_all) <- sub("^[^_]*_", "", colnames(geno_all))
flip <- traw$COUNTED != traw$ALT
if (any(flip)) geno_all[flip, ] <- 2 - geno_all[flip, ]
allele_dt <- data.table(h4_snp = traw$SNP, ref_allele = traw$COUNTED, alt_allele = traw$ALT)

# --- per stage: residualised phenotypes and mediation of coloc PP4 pairs ---
med_list <- lapply(names(stage_key), function(stage) {
  Stage <- stage_key[[stage]]
  pairs <- fread(file.path(D, sprintf("coloc_%s_apa_eqtl.tsv", stage)))[PP_H4 >= PP_H4_MIN]
  if (!nrow(pairs)) return(NULL)
  expr <- residualise_phenotypes(file.path(D, sprintf("eqtl.phenotype.%s.bed.gz", Stage)),
                                 file.path(D, sprintf("eqtl.covariates_for_qtl.%s.txt", Stage)))
  apa <- residualise_phenotypes(file.path(D, sprintf("phenotype.%s.GeneRange.bed.gz", Stage)),
                                file.path(D, sprintf("covariates_for_qtl.%s.txt", Stage)))
  samples <- Reduce(intersect, list(colnames(geno_all), colnames(expr$residuals), colnames(apa$residuals)))
  G <- geno_all[, samples, drop = FALSE]
  M <- apa$residuals[, samples, drop = FALSE]
  Y <- expr$residuals[, samples, drop = FALSE]
  C <- expr$geno_pcs[samples, , drop = FALSE]

  have <- pairs[, h4_snp %in% rownames(G) & apa_pheno %in% rownames(M) & eqtl_pheno %in% rownames(Y)]
  pairs[, resolvable := have]
  idx <- which(have)
  res <- mclapply(idx, function(i) {
    r <- mediate_pair(as.numeric(G[pairs$h4_snp[i], ]), as.numeric(M[pairs$apa_pheno[i], ]),
                      as.numeric(Y[pairs$eqtl_pheno[i], ]), C)
    c(list(row = i), r)
  }, mc.cores = N_CORES, mc.set.seed = TRUE)

  ok <- vapply(res, function(r) identical(r$status, "ok"), logical(1))
  if (!any(ok)) return(NULL)
  stat_dt <- rbindlist(lapply(res[ok], function(r) {
    r$status <- NULL
    as.data.table(r)
  }))
  out <- cbind(pairs[stat_dt$row], stat_dt[, -"row"])
  out[, stage := stage]
  out[]
})
med <- rbindlist(med_list, fill = TRUE)
if (!nrow(med)) stop("No pair could be tested")

# --- BH within stage, mediation calls and direction ---
med[, fdr_ab_fwd := p.adjust(p_ab_fwd, "BH"), by = stage]
med[, fdr_ab_rev := p.adjust(p_ab_rev, "BH"), by = stage]
med[, fdr_c_prime := p.adjust(p_c_prime, "BH"), by = stage]
med[, fdr_d_prime := p.adjust(p_d_prime, "BH"), by = stage]
med[, fwd_compatible := p_c_prime >= CI_ALPHA]
med[, rev_compatible := p_d_prime >= CI_ALPHA]
med[, fwd_residual := fdr_c_prime < FDR_CUT]
med[, rev_residual := fdr_d_prime < FDR_CUT]
med[, has_fwd_indirect := fdr_ab_fwd < FDR_CUT]
med[, has_rev_indirect := fdr_ab_rev < FDR_CUT]
med[, mediation_call := fcase(
  has_fwd_indirect & fwd_compatible & rev_residual, "APA -> EXP",
  has_rev_indirect & rev_compatible & fwd_residual, "EXP -> APA",
  has_fwd_indirect & fwd_residual & rev_residual, "Partial mediation",
  default = "No mediation")]
med[, mediation_call := factor(mediation_call, levels = call_levels)]
med[, apa_effect_on_expr := fcase(
  p_b >= 0.05, "not significant",
  b < 0, "longer 3'UTR -> lower expression",
  default = "longer 3'UTR -> higher expression")]
med <- merge(med, allele_dt, by = "h4_snp", all.x = TRUE, sort = FALSE)

core_cols <- c("stage", "gene_name", "ensg0", "h4_snp", "ref_allele", "alt_allele",
               "PP_H4", "mediation_call", "apa_effect_on_expr",
               "a", "p_a", "b", "p_b", "c_total", "p_c", "c_prime", "p_c_prime",
               "d_prime", "p_d_prime", "ab_fwd", "fdr_ab_fwd", "prop_med", "n", "maf")
core_cols <- intersect(core_cols, names(med))
setorder(med, stage, p_ab_fwd)
fwrite(med[, ..core_cols], file.path(D, "04_mediation_results.tsv"), sep = "\t", quote = FALSE)

summary_dt <- med[, .N, by = .(stage, mediation_call)]
summary_dt[, frac := N / sum(N), by = stage]
setorder(summary_dt, stage, mediation_call)
fwrite(summary_dt, file.path(D, "04_mediation_summary.tsv"), sep = "\t", quote = FALSE)

# --- GO / KEGG enrichment of APA-driven mediated genes against the shared QTL universe ---
universe_ensg <- unique(fread(file.path(D, "shared_qtl_universe_ensg.tsv"))$ensg)
target_ensg <- unique(med[mediation_call %in% c("APA -> EXP", "Partial mediation")]$ensg0)
target_ensg <- target_ensg[target_ensg %in% universe_ensg]
all_ensg <- unique(c(target_ensg, universe_ensg))
all_mapped <- suppressMessages(bitr(all_ensg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
setDT(all_mapped)
target_ids <- all_mapped[ENSEMBL %in% target_ensg, ENTREZID]
universe_ids <- all_mapped[ENSEMBL %in% universe_ensg, ENTREZID]

enrich_list <- list()
if (length(target_ids) >= 5 && length(universe_ids) > length(target_ids)) {
  for (ont in c("BP", "CC", "MF")) {
    ego <- enrichGO(gene = target_ids, universe = universe_ids, OrgDb = org.Hs.eg.db, ont = ont,
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    minGSSize = 10, maxGSSize = 500, readable = TRUE)
    dt <- as.data.table(ego)
    if (nrow(dt)) {
      dt[, ONTOLOGY := ont]
      enrich_list[[ont]] <- dt
    }
  }
  kegg_dt <- tryCatch({
    kk <- enrichKEGG(gene = target_ids, organism = "hsa", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10)
    dt <- as.data.table(kk)
    if (nrow(dt)) { dt[, ONTOLOGY := "KEGG"]; dt } else data.table()
  }, error = function(e) {
    message("KEGG enrichment skipped: ", e$message)
    data.table()
  })
  if (nrow(kegg_dt)) enrich_list[["KEGG"]] <- kegg_dt
}

if (length(enrich_list)) {
  enrich_dt <- rbindlist(enrich_list, fill = TRUE, use.names = TRUE)
  setorder(enrich_dt, p.adjust, na.last = TRUE)
  enrich_dt <- enrich_dt[, head(.SD, 1), by = Description]
  setorder(enrich_dt, p.adjust)
} else {
  enrich_dt <- data.table()
}
fwrite(enrich_dt, file.path(D, "04_mediation_enrichment.tsv"), sep = "\t", quote = FALSE)
