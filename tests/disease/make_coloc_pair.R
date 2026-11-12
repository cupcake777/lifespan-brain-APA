# Simulates two loci (AR(1) haplotype LD, rho = 0.9) for a GWAS and a QTL, then runs coloc.abf and coloc.susie.
# SIMSHARED: one causal SNP shared by GWAS and QTL (index 60), expect PP4 high.
# SIMDISTINCT: GWAS causal at index 40, QTL causal at index 100, expect PP3 high and PP4 low.
# Usage: Rscript make_coloc_pair.R <outdir>. Writes synthetic_truth.tsv for the test to check.

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
suppressMessages(library(coloc))

# --- settings ---

set.seed(20261019)
N_GWAS <- 1500
N_QTL  <- 1500
M      <- 120
RHO    <- 0.9
MAF    <- runif(M, 0.2, 0.45)
snp_ids <- sprintf("rsSIM%03d", seq_len(M))

# --- helpers ---

simulate_dosage <- function(n, maf, rho) {
  m <- length(maf)
  z <- matrix(0, 2 * n, m)
  e <- matrix(rnorm(2 * n * m), 2 * n, m)
  z[, 1] <- e[, 1]
  for (j in 2:m) z[, j] <- rho * z[, j - 1] + sqrt(1 - rho^2) * e[, j]
  h <- sweep(z, 2, qnorm(maf), FUN = "<") * 1L     # haplotype allele = latent below quantile
  h[1:n, ] + h[(n + 1):(2 * n), ]                  # dosage 0/1/2
}

univariate <- function(g, y) {
  n <- nrow(g)
  xc <- sweep(g, 2, colMeans(g))
  yc <- y - mean(y)
  sxx <- colSums(xc^2)
  b <- colSums(xc * yc) / sxx
  s2 <- (sum(yc^2) - b^2 * sxx) / (n - 2)
  list(beta = b, varbeta = s2 / sxx)
}

make_locus <- function(gene, causal_gwas, causal_qtl, b_gwas, b_qtl) {
  G <- simulate_dosage(N_GWAS, MAF, RHO)
  Q <- simulate_dosage(N_QTL, MAF, RHO)
  colnames(G) <- snp_ids; colnames(Q) <- snp_ids
  yg <- as.numeric(G[, causal_gwas] * b_gwas + rnorm(N_GWAS))
  yq <- as.numeric(Q[, causal_qtl] * b_qtl + rnorm(N_QTL))
  ug <- univariate(G, yg); uq <- univariate(Q, yq)
  LD <- cor(G); dimnames(LD) <- list(snp_ids, snp_ids)
  LDq <- cor(Q); dimnames(LDq) <- list(snp_ids, snp_ids)
  list(
    gene = gene,
    d_gwas = list(beta = ug$beta, varbeta = ug$varbeta, N = N_GWAS, type = "quant",
                  sdY = sd(yg), snp = snp_ids, MAF = MAF, LD = LD),
    d_qtl  = list(beta = uq$beta, varbeta = uq$varbeta, N = N_QTL, type = "quant",
                  sdY = sd(yq), snp = snp_ids, MAF = MAF, LD = LDq)
  )
}

# --- simulate and run coloc ---

loci <- list(
  make_locus("SIMSHARED",   causal_gwas = 60, causal_qtl = 60,  b_gwas = 0.25, b_qtl = 0.8),
  make_locus("SIMDISTINCT", causal_gwas = 40, causal_qtl = 100, b_gwas = 0.25, b_qtl = 0.8)
)

abf_rows <- list()
susie_rows <- list()
for (lc in loci) {
  ab <- coloc.abf(dataset1 = lc$d_gwas, dataset2 = lc$d_qtl, MAF = MAF)
  pp <- ab$summary
  abf_rows[[lc$gene]] <- data.frame(
    GWAS_ID = "SIM_TRAIT", QTL_ID = "SIM_QTL", Locus = lc$gene,
    Phenotype = paste0("ENSTSIM0001.1|", lc$gene, "|chr1:1000-2000|+"),
    PP0 = unname(pp["PP.H0.abf"]), PP1 = unname(pp["PP.H1.abf"]), PP2 = unname(pp["PP.H2.abf"]),
    PP3 = unname(pp["PP.H3.abf"]), PP4 = unname(pp["PP.H4.abf"]),
    n_snps = unname(pp["nsnps"]), identification_method = "coloc.abf", stringsAsFactors = FALSE)
  cs <- coloc.susie(dataset1 = lc$d_gwas, dataset2 = lc$d_qtl)
  susie_rows[[lc$gene]] <- data.frame(
    source_file = paste0("SIM_", lc$gene, "_susie"), nsnps = cs$summary$nsnps[1],
    PP.H4.abf = cs$summary$PP.H4.abf[1], stringsAsFactors = FALSE)
}

abf_df <- do.call(rbind, abf_rows)
susie_df <- do.call(rbind, susie_rows)
write.csv(abf_df, file.path(outdir, "synthetic_coloc_results.csv"), row.names = FALSE)
write.csv(susie_df, file.path(outdir, "synthetic_susie_results.csv"), row.names = FALSE)
truth <- data.frame(Locus = c("SIMSHARED", "SIMDISTINCT"),
                    expected = c("shared_causal_PP4_high", "distinct_causal_PP3_high_PP4_low"),
                    gwas_causal = c(snp_ids[60], snp_ids[40]), qtl_causal = c(snp_ids[60], snp_ids[100]))
write.table(truth, file.path(outdir, "synthetic_truth.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(abf_df[, c("Locus", "PP3", "PP4", "n_snps")])
print(susie_df)
