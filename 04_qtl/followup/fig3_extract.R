#!/usr/bin/env Rscript
# Fig 3 inputs that need cohort-level data (summary-level outputs only).
# A: eGene IDs per stage. B: covariate-adjusted expression by genotype for the PAS-variant candidates.
# C: group-mean 3'UTR coverage for ANKMY1 and ZSWIM7 from DaPars2 WIG files (only if WIG_DIR is set).
# Needs data.table and plink2 on PATH.
suppressMessages(library(data.table))

for (.v in c("APA_HPC_HOME", "APA_ROOT")) if (!nzchar(Sys.getenv(.v))) stop("Set ", .v, " (source config/paths.sh)")
H        <- Sys.getenv("APA_HPC_HOME")
EQ       <- file.path(H, "project/work/eqtl")
EQ_PRE   <- file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl", "eqtl", "pre_input")
APA_PI   <- file.path(H, "project/work/apa/01_QTL/pre_input")
D        <- file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl")
STAGES   <- c(prenatal = "Prenatal", postnatal = "Postnatal")

# --- helpers ---
dosage <- function(pfile, snps) {
  tmp <- tempfile(); writeLines(snps, paste0(tmp, ".snps"))
  st <- system(sprintf("plink2 --pfile %s --extract %s --export A --out %s >/dev/null 2>&1", shQuote(pfile), paste0(tmp, ".snps"), tmp))
  if (st != 0 || !file.exists(paste0(tmp, ".raw"))) stop("plink2 export failed for ", pfile)
  r <- fread(paste0(tmp, ".raw")); cols <- names(r)[-(1:6)]
  out <- lapply(cols, function(cn) {
    id <- sub("_[ACGT]+$", "", cn); counted <- sub(".*_", "", cn); alt <- tail(strsplit(id, ":")[[1]], 1)
    v <- as.numeric(r[[cn]]); if (counted != alt) v <- 2 - v            # ALT copies
    setNames(v, r$IID)
  }); names(out) <- sub("_[ACGT]+$", "", cols); out
}
pfile_of <- function(dir, st) {
  f <- Sys.glob(file.path(dir, sprintf("genotype.%s.pgen", st)))
  if (!length(f)) f <- Sys.glob(file.path(dir, sprintf("genotype*%s*.pgen", st)))
  if (!length(f)) stop("no genotype pgen in ", dir)
  sub("\\.pgen$", "", f[1])
}

# --- A. eGene sets ---
eg <- rbindlist(lapply(names(STAGES), function(st) {
  x <- fread(cmd = paste("zcat", shQuote(file.path(EQ, sprintf("%s.sig_QTL.txt.gz", st)))), select = "phenotype_id")
  x[, .(n_pairs = .N), by = .(ensg0 = sub("\\..*", "", phenotype_id))][, stage := st][]
}))
fwrite(eg, file.path(D, "fig3_egenes.tsv"), sep = "\t")

# --- B. expression by genotype for candidates ---
TG <- data.table(gene = c("CEP41", "ANKMY1", "HNRNPDL", "EARS2", "ZSWIM7", "SDHA", "UNC13A", "INPP1", "STAT6", "LIN54", "POP4", "COG3", "HID1", "ACSF3"),
                 ensg = c("ENSG00000106477", "ENSG00000144504", "ENSG00000152795", "ENSG00000103356", "ENSG00000214941", "ENSG00000073578", "ENSG00000130477",
                          "ENSG00000151689", "ENSG00000166888", "ENSG00000189308", "ENSG00000105171", "ENSG00000136152", "ENSG00000167861", "ENSG00000176715"),
                 snp  = c("chr7:130393789:C:T", "chr2:240521955:T:A", "chr4:82422592:T:C", "chr16:23525045:T:C", "chr17:15976596:C:T", "chr5:257077:A:G", "chr19:17626316:C:T",
                          "chr2:190362766:A:G", "chr12:57095865:C:T", "chr4:82927812:G:A", "chr19:29615752:G:A", "chr13:45493458:G:A", "chr17:74961620:T:C", "chr16:89160671:G:T"))
resB <- list()
for (st in names(STAGES)) {
  Stage <- STAGES[[st]]
  bed <- fread(cmd = paste("zcat", shQuote(file.path(EQ_PRE, sprintf("phenotype.%s.bed.gz", Stage)))))
  cv  <- fread(file.path(EQ_PRE, sprintf("covariates_for_qtl.%s.txt", Stage)))
  gd  <- tryCatch(dosage(pfile_of(EQ_PRE, st), TG$snp), error = function(e) dosage(pfile_of(APA_PI, st), TG$snp))
  smp0 <- intersect(setdiff(names(bed), names(bed)[1:4]), setdiff(names(cv), names(cv)[1]))
  for (i in seq_len(nrow(TG))) {
    row <- bed[startsWith(bed[[4]], TG$ensg[i])]; x0 <- gd[[TG$snp[i]]]
    if (nrow(row) != 1 || is.null(x0)) { message("B skip ", TG$gene[i], " ", st); next }
    smp <- intersect(smp0, names(x0)); X <- cbind(1, t(as.matrix(cv[, smp, with = FALSE])))
    y <- qr.resid(qr(X), as.numeric(row[, smp, with = FALSE])); x <- x0[smp]; rx <- qr.resid(qr(X), x)
    b <- sum(rx * y) / sum(rx^2); df <- length(y) - ncol(X) - 1; se <- sqrt(sum((y - b * rx)^2) / df / sum(rx^2))
    gl <- strsplit(TG$snp[i], ":")[[1]][3:4]; labs <- c(paste0(gl[1], "/", gl[1]), paste0(gl[1], "/", gl[2]), paste0(gl[2], "/", gl[2]))
    for (k in 0:2) if (sum(round(x) == k) > 0) {
      v <- y[round(x) == k]
      q <- quantile(v, c(.25, .5, .75), names = FALSE); iqr <- q[3] - q[1]
      lo <- min(v[v >= q[1] - 1.5 * iqr]); hi <- max(v[v <= q[3] + 1.5 * iqr])
      bx <- data.table(n = length(v), mean = mean(v), q1 = q[1], median = q[2], q3 = q[3], whisker_lo = lo, whisker_hi = hi,
                       outliers = paste(signif(sort(v[v < lo | v > hi]), 5), collapse = ","),
                       values = if (length(v) < 5) paste(signif(sort(v), 5), collapse = ",") else "")
      resB[[length(resB) + 1]] <- cbind(data.table(gene = TG$gene[i], snp = TG$snp[i], stage = Stage, genotype_genomic = labs[k + 1], dosage_alt = k,
                                                 slope = b, se = se, p = 2 * pt(-abs(b / se), df), n_stage = length(y)), bx)
    }
  }
}
resB <- rbindlist(resB); fwrite(resB, file.path(D, "fig3_expr_boxstats.tsv"), sep = "\t")

# --- C. group-mean 3'UTR coverage from DaPars2 WIGs (optional) ---
# minus-strand events; window = 3'UTR +/- flank. cUTR = proximal PAS .. UTR start (normaliser)
LOCI <- data.table(gene = c("ANKMY1", "ZSWIM7"), chr = c("chr2", "chr17"),
                   W0 = c(240521700L, 15976300L), W1 = c(240524600L, 15978400L),
                   utr_start = c(240524381L, 15978163L), ppas = c(240523884L, 15977740L), dpas = c(240521938L, 15976560L),
                   snp = c("chr2:240521955:T:A", "chr17:15976596:C:T"))
CAP <- 80L                                               # max samples per stage x genotype (random, seed 2026)
wd <- Sys.getenv("WIG_DIR")                              # optional: coverage only if WIG files are available
if (!nzchar(wd) || !dir.exists(wd)) { message("C skipped: WIG_DIR not set"); quit(save = "no", status = 0) }
wigs <- list.files(wd, pattern = "\\.wig(\\.gz)?$", full.names = TRUE); names(wigs) <- sub("\\.wig(\\.gz)?$", "", basename(wigs))
first <- readLines(if (grepl("gz$", wigs[1])) gzfile(wigs[1]) else wigs[1], n = 5); first <- first[!grepl("^(track|#)", first)][1]
LOCI[, chr_w := if (grepl("^chr", first)) chr else sub("^chr", "", chr)]

set.seed(2026); acc <- list()
for (st in names(STAGES)) {
  Stage <- STAGES[[st]]
  smp <- readLines(Sys.glob(file.path(APA_PI, sprintf("sample.%s.final", Stage)))[1])
  gd <- dosage(pfile_of(APA_PI, st), LOCI$snp); smp <- intersect(smp, names(wigs))
  picks <- lapply(seq_len(nrow(LOCI)), function(k) { x0 <- gd[[LOCI$snp[k]]][intersect(smp, names(gd[[LOCI$snp[k]]]))]
    lapply(0:2, function(g) { s <- names(x0)[round(x0) == g]; if (length(s) > CAP) sample(s, CAP) else s }) })
  todo <- unique(unlist(picks))
  V <- setNames(vector("list", length(todo)), todo)
  # one pass per sample over both regions; awk stops once every region has been passed
  for (ii in seq_along(todo)) {
    f <- wigs[[todo[ii]]]
    conds <- paste(sprintf("($1==\"%s\" && $3>%d && $2<%d){print %d\"\\t\"$2\"\\t\"$3\"\\t\"$4; h%d=1; next}",
                           LOCI$chr_w, LOCI$W0 - 1L, LOCI$W1, seq_len(nrow(LOCI)), seq_len(nrow(LOCI))), collapse = " ")
    done <- paste(sprintf("h%d", seq_len(nrow(LOCI))), collapse = " && ")
    cmd <- sprintf("%s %s | awk '%s %s {exit}'", if (grepl("gz$", f)) "zcat" else "cat", shQuote(f), conds, done)
    x <- tryCatch(fread(cmd = cmd, header = FALSE, col.names = c("k", "s", "e", "d")), error = function(e) data.table(k = integer(), s = integer(), e = integer(), d = numeric()))
    V[[ii]] <- lapply(seq_len(nrow(LOCI)), function(k) { L <- LOCI[k]; v <- numeric(L$W1 - L$W0 + 1L); xk <- x[x$k == k]
      if (nrow(xk)) for (j in seq_len(nrow(xk))) { a <- max(xk$s[j] + 1L, L$W0); b <- min(xk$e[j], L$W1); if (b >= a) v[(a:b) - L$W0 + 1L] <- xk$d[j] }
      v })
  }
  for (k in seq_len(nrow(LOCI))) { L <- LOCI[k]; pos <- L$W0:L$W1; cutr <- pos >= L$ppas & pos <= L$utr_start
    for (g in 0:2) { s_g <- picks[[k]][[g + 1]]; if (!length(s_g)) next
      M <- vapply(s_g, function(s) { v <- V[[s]][[k]]; m <- mean(v[cutr]); if (m > 0) v / m else rep(NA_real_, length(v)) }, numeric(length(pos)))
      M <- M[, colSums(is.na(M)) == 0, drop = FALSE]
      acc[[length(acc) + 1]] <- data.table(gene = L$gene, stage = Stage, dosage_alt = g, n = ncol(M), position = pos, rel_cov = rowMeans(M))
    } }
}
cv_out <- rbindlist(acc)[, .(rel_cov = mean(rel_cov)), by = .(gene, stage, dosage_alt, n, bin = (position %/% 10L) * 10L)]
fwrite(cv_out, file.path(D, "fig3_locus_cov.tsv"), sep = "\t")
