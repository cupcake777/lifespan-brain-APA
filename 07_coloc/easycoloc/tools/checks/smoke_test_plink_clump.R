#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

normalize_chr_for_plink <- function(x) {
  gsub("^chr", "", as.character(x), ignore.case = TRUE)
}

gwas_candidates <- c("./harmony/EAS_SCZ_b19to38_harmonized.tsv.gz", "./harmony/EAS_SCZ_b19to38_harmonized.tsv")
gwas_existing <- gwas_candidates[file.exists(gwas_candidates)]
gwas_file <- if (length(gwas_existing) > 0L) gwas_existing[[1L]] else gwas_candidates[[1L]]
bim_prefix <- Sys.getenv("EASYCOLOC_SMOKE_PLINK_PREFIX", unset = "refs/plink/hg38/1kg_hg38_filtered")
keep_file <- Sys.getenv("EASYCOLOC_SMOKE_KEEP_FILE", unset = "refs/plink/keep/EAS.sample")
bim_file <- paste0(bim_prefix, ".bim")
plink_prefix <- "temp/smoke_clump_EAS_SCZ"
input_file <- paste0(plink_prefix, ".qassoc")

if (!file.exists(gwas_file)) {
  cat("[SMOKE-SKIP] harmonized GWAS fixture missing:", paste(gwas_candidates, collapse = " or "), "\n")
  quit(save = "no", status = 0)
}

if (!file.exists(bim_file) || !file.exists(keep_file)) {
  cat("[SMOKE-SKIP] local PLINK references unavailable\n")
  cat("[SMOKE-SKIP] Set EASYCOLOC_SMOKE_PLINK_PREFIX and EASYCOLOC_SMOKE_KEEP_FILE to run this smoke test.\n")
  quit(save = "no", status = 0)
}

gwas <- fread(gwas_file)
gwas <- gwas[P <= 1e-7]

gwas[, CHR_CLUMP := normalize_chr_for_plink(CHR)]
gwas[, POS_CLUMP := as.numeric(POS)]
gwas[, EA_CLUMP := as.character(EA)]
gwas[, NEA_CLUMP := as.character(NEA)]

bim <- fread(
  bim_file,
  header = FALSE,
  select = c(1, 2, 4, 5, 6),
  col.names = c("CHR", "plink_snp_id", "POS", "V5", "V6")
)
bim[, CHR := normalize_chr_for_plink(CHR)]
bim[, POS := as.numeric(POS)]

merged <- merge(
  gwas,
  bim,
  by.x = c("CHR_CLUMP", "POS_CLUMP"),
  by.y = c("CHR", "POS"),
  all.x = FALSE,
  all.y = FALSE
)

matched <- merged[
  (EA_CLUMP == V5 & NEA_CLUMP == V6) |
  (EA_CLUMP == V6 & NEA_CLUMP == V5)
]

assert_true(nrow(merged) > 0, "CHR/POS merge returned zero rows")
assert_true(nrow(matched) > 0, "Dual-direction allele matching returned zero rows")

setorderv(matched, cols = "P", order = 1L, na.last = TRUE)
matched <- matched[!duplicated(plink_snp_id)]
clump_input <- matched[, .(SNP = plink_snp_id, P = P)]

assert_true(nrow(clump_input) > 0, "Translated PLINK clump input is empty")
fwrite(clump_input, input_file, sep = "\t", quote = FALSE)

plink_cmd <- paste(
  "plink",
  "--bfile", shQuote(bim_prefix),
  "--clump", shQuote(input_file),
  "--keep", shQuote(keep_file),
  "--clump-p1 1e-7",
  "--clump-p2 1e-7",
  "--clump-r2 0.2",
  "--clump-kb 1000",
  "--out", shQuote(plink_prefix)
)

exit_code <- system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
assert_true(exit_code == 0, "PLINK clump exited with non-zero status")
assert_true(file.exists(paste0(plink_prefix, ".clumped")), "PLINK did not write a .clumped file")

clumped <- fread(paste0(plink_prefix, ".clumped"))
assert_true(nrow(clumped) > 0, "PLINK .clumped file is empty")

cat("[SMOKE] PLINK clump smoke test passed\n")
cat("[SMOKE] merged_by_chr_pos:", nrow(merged), "\n")
cat("[SMOKE] dual_direction_matched:", nrow(matched), "\n")
cat("[SMOKE] clumps:", nrow(clumped), "\n")
