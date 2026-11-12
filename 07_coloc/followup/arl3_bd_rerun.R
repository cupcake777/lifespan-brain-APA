# ARL3 x BD coloc.abf rerun for apaQTL and eQTL, prenatal and postnatal; prints PP4 and writes no files.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
suppressMessages(library(coloc))
bd <- read.delim(file.path(Sys.getenv("APA_ROOT"), "data/07_coloc/arl3_bd_coloc_input.tsv"))
# --- coloc per stage ---
N_BD <- 9643 + 89185; S_BD <- 9643 / N_BD
for (st in c("prenatal", "postnatal")) {
  b <- subset(bd, stage == st & !is.na(beta_bd) & !is.na(se_bd)); b <- b[!duplicated(b$rsid), ]
  a <- subset(b, !is.na(qtl_logp))
  r1 <- suppressWarnings(coloc.abf(
    list(snp = a$rsid, beta = a$beta_bd, varbeta = a$se_bd^2, N = N_BD, s = S_BD, type = "cc", MAF = a$af_apa),
    list(snp = a$rsid, pvalues = 10^(-a$qtl_logp), N = ifelse(st == "prenatal", 189, 508), type = "quant", MAF = a$af_apa),
    p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  e <- subset(b, !is.na(eq_beta) & !is.na(eq_se) & !is.na(af_exp))
  r2 <- suppressWarnings(coloc.abf(
    list(snp = e$rsid, beta = e$beta_bd, varbeta = e$se_bd^2, N = N_BD, s = S_BD, type = "cc", MAF = e$af_exp),
    list(snp = e$rsid, beta = e$eq_beta, varbeta = e$eq_se^2, N = ifelse(st == "prenatal", 182, 497), type = "quant", MAF = e$af_exp),
    p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  cat(st, "apaQTL PP4", round(r1$summary[["PP.H4.abf"]], 4), "eQTL PP4", round(r2$summary[["PP.H4.abf"]], 4), "\n")
}
