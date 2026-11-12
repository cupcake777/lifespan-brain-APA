#!/usr/bin/env Rscript
# ACSF3 mediation example (chr16:89160671 G>T): box statistics per genotype for expression and APA residuals.
# Phenotypes are residualised on all QTL covariates (incl. genoPCs), the tensorQTL convention.
# Summary only: box statistics per genotype; outlier values without sample IDs.
suppressMessages(library(data.table))

for (.v in c("APA_HPC_HOME", "APA_ROOT")) if (!nzchar(Sys.getenv(.v))) stop("Set ", .v, " (source config/paths.sh)")
H          <- Sys.getenv("APA_HPC_HOME")
eqtl_dir   <- file.path(H, "project/work/apa/11_eqtl")
apa_dir    <- file.path(H, "project/work/apa/01_QTL")
expr_dir   <- file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl", "eqtl")
D          <- file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl")
TARGET_SNP <- "chr16:89160671:G:T"
ENSG       <- "ENSG00000176715"           # ACSF3
APA_ID     <- "ENST00000317447.9|ACSF3"   # prefix of the APA pid in the GeneRange bed
OUT        <- file.path(D, "fig3_acsf3_boxstats.tsv")

# --- helpers ---
# residualised phenotype row for one pid
resid_row <- function(bed_file, cov_file, pid_pattern) {
  bed <- fread(cmd = paste("zcat", shQuote(bed_file)))
  row <- bed[grepl(pid_pattern, bed[[4]], fixed = TRUE)]
  if (nrow(row) != 1) stop(nrow(row), " rows match ", pid_pattern, " in ", basename(bed_file))
  cv  <- fread(cov_file)
  smp <- intersect(setdiff(names(bed), names(bed)[1:4]), setdiff(names(cv), names(cv)[1]))
  X   <- cbind(1, t(as.matrix(cv[, smp, with = FALSE])))
  y   <- as.numeric(row[, smp, with = FALSE])
  list(pid = row[[4]], r = setNames(as.numeric(qr.resid(qr(X), y)), smp), X = X, smp = smp)
}

# --- target SNP dosage (ALT = T copies) ---
traw <- fread(file.path(eqtl_dir, "coloc_target_geno.traw"))
g    <- traw[traw$SNP == TARGET_SNP]
if (nrow(g) != 1) stop("SNP not found in traw")
dos  <- as.numeric(g[, -(1:6), with = FALSE]); names(dos) <- sub("^[^_]*_", "", names(g)[-(1:6)])
if (g$COUNTED != g$ALT) dos <- 2 - dos            # ALT (T) copies, as in 04-1_mediate.R

# --- box statistics per stage, measure and genotype ---
out <- list()
for (Stage in c("Prenatal", "Postnatal")) {
  ph <- list(
    Expression = resid_row(file.path(expr_dir, "pre_input", sprintf("phenotype.%s.bed.gz", Stage)),
                           file.path(expr_dir, "pre_input", sprintf("covariates_for_qtl.%s.txt", Stage)), ENSG),
    APA        = resid_row(file.path(apa_dir, "pre_input", sprintf("phenotype.%s.GeneRange.bed.gz", Stage)),
                           file.path(apa_dir, "pre_input", sprintf("covariates_for_qtl.%s.txt", Stage)), APA_ID))
  shared <- Reduce(intersect, list(names(dos), ph$Expression$smp, ph$APA$smp))
  for (m in names(ph)) {
    r <- ph[[m]]$r[shared]; x <- dos[shared]
    rx <- qr.resid(qr(ph[[m]]$X[match(shared, ph[[m]]$smp), ]), x)   # tensorQTL-style slope check
    slope <- sum(rx * r) / sum(rx^2)
    for (k in 0:2) if (sum(x == k) > 0) {
      v <- r[x == k]
      q <- quantile(v, c(.25, .5, .75), names = FALSE); iqr <- q[3] - q[1]
      lo <- min(v[v >= q[1] - 1.5 * iqr]); hi <- max(v[v <= q[3] + 1.5 * iqr])
      bx <- data.table(n = length(v), mean = mean(v), q1 = q[1], median = q[2], q3 = q[3],
                       whisker_lo = lo, whisker_hi = hi,
                       outliers = paste(signif(sort(v[v < lo | v > hi]), 5), collapse = ","))
      out[[length(out) + 1]] <- cbind(data.table(stage = Stage, measure = m, pid = ph[[m]]$pid,
                                                 genotype = c("G/G", "G/T", "T/T")[k + 1], dosage_T = k,
                                                 slope_check = slope, n_stage = length(shared)), bx)
    }
  }
}
res <- rbindlist(out)
fwrite(res, OUT, sep = "\t")
