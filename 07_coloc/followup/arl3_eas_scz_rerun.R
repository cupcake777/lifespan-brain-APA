# ARL3 x EAS SCZ coloc.abf for apaQTL and eQTL, prenatal and postnatal; writes arl3_eas_scz_coloc_rerun.tsv.
# Inputs are built by arl3_prepare.py. GWAS: p-values, N = 14,006 cases + 16,755 controls, s = cases/N.
# apaQTL: p-values, N 189 prenatal / 508 postnatal. eQTL: beta and varbeta, N 182 / 497.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
suppressMessages(library(coloc))
d <- file.path(Sys.getenv("APA_ROOT"), "data/07_coloc")
apa <- read.delim(file.path(d, "in_apaqtl_eassc.tsv"))
eq  <- read.delim(file.path(d, "in_eqtl_regional_maf.tsv"))
# --- coloc per stage ---
N_GWAS <- 14006 + 16755; S_GWAS <- 14006 / N_GWAS
N_APA <- c(prenatal = 189, postnatal = 508); N_EQ <- c(prenatal = 182, postnatal = 497)
p1 <- 1e-4; p2 <- 1e-4; p12 <- 1e-5
out <- list()
for (st in c("prenatal", "postnatal")) {
  a <- subset(apa, stage == st & !is.na(af_apa) & af_apa > 0 & af_apa < 1)
  a <- a[!duplicated(a$snp), ]
  g_apa <- list(snp = a$snp, pvalues = 10^(-a$gwas_logp), N = N_GWAS, s = S_GWAS, type = "cc", MAF = a$af_apa)
  q_apa <- list(snp = a$snp, pvalues = 10^(-a$apa_logp), N = N_APA[[st]], type = "quant", MAF = a$af_apa)
  r1 <- coloc.abf(g_apa, q_apa, p1 = p1, p2 = p2, p12 = p12)
  out[[paste("apaQTL", st)]] <- r1$summary
  e <- subset(eq, stage == st & !is.na(maf_eq) & maf_eq > 0 & maf_eq < 1 & !is.na(se) & se > 0)
  e <- e[!duplicated(e$rsid), ]
  ai <- a[match(e$rsid, a$snp), ]   # GWAS p-values at the same SNPs
  keep <- !is.na(ai$snp)
  e <- e[keep, ]; ai <- ai[keep, ]
  g_eq <- list(snp = ai$snp, pvalues = 10^(-ai$gwas_logp), N = N_GWAS, s = S_GWAS, type = "cc", MAF = ai$af_apa)
  q_eq <- list(snp = e$rsid, beta = e$beta, varbeta = e$se^2, N = N_EQ[[st]], type = "quant", MAF = e$maf_eq)
  r2 <- coloc.abf(g_eq, q_eq, p1 = p1, p2 = p2, p12 = p12)
  out[[paste("eQTL", st)]] <- r2$summary
  cat(st, "apaQTL nsnps", r1$summary[["nsnps"]], "PP4", round(r1$summary[["PP.H4.abf"]], 4),
      "| eQTL nsnps", r2$summary[["nsnps"]], "PP4", round(r2$summary[["PP.H4.abf"]], 4), "\n")
}
# --- write ---
res <- do.call(rbind, lapply(names(out), function(k) data.frame(comparison = k, t(out[[k]]))))
write.table(res, file.path(d, "arl3_eas_scz_coloc_rerun.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
