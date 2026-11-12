#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_format.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

sumstats_dt <- data.table(
  CHR = c(1, 1, 1),
  POS = c(100L, 120L, 150L),
  P = c(1e-8, 2e-7, 3e-6)
)

lifted_positions <- data.table(
  V1 = c("chr1", "chr1"),
  V2 = c(1000L, 1300L),
  V3 = c(1001L, 1301L)
)

updated_dt <- apply_liftover_positions(
  sumstats_dt = sumstats_dt,
  chrom = "1",
  start_pos = 100L,
  end_pos = 150L,
  hg38_positions = lifted_positions
)

assert_true(identical(updated_dt$CHR, c("1", "1", "1")), "CHR should be updated to lifted chromosome")
assert_true(identical(updated_dt$POS, c(1001L, 1121L, 1301L)), "POS should be remapped into the lifted interval")
assert_true(identical(sumstats_dt$POS, c(100L, 120L, 150L)), "input data should remain unchanged")

messy_harmonized <- data.table(
  SNPID = c("rs1", "rs1", "rs2"),
  rsID = c("rs1", "rs1", "rs2"),
  CHR = c("chr1", "chr1", "chr1"),
  POS = c(1001L, 1001L, 1121L),
  EA = c("a", "a", "g"),
  NEA = c("g", "g", "a"),
  EAF = c(0.2, 0.2, 0.7),
  BETA = c(0.1, 0.1, -0.2),
  SE = c(0.03, 0.03, 0.04),
  P = c(1e-8, 1e-8, 2e-7),
  N = c(1000, 1000, 1000),
  REF_ALLELE = c("G", "G", "A"),
  ALLELE_FLIPPED = c(FALSE, FALSE, TRUE),
  STRAND_FLIPPED = c(FALSE, FALSE, FALSE),
  REF_MATCH = c(TRUE, TRUE, TRUE),
  LIFTOVER_MAPPED = c(TRUE, TRUE, TRUE)
)

clean_harmonized <- easycoloc_prepare_harmonized_gwas_output(messy_harmonized)
expected_cols <- easycoloc_harmonized_gwas_output_cols()
assert_true(identical(names(clean_harmonized), expected_cols), "clean harmonized output schema should be exact and ordered")
assert_true(nrow(clean_harmonized) == 2L, "clean harmonized output should remove exact duplicate rows")
assert_true(!any(c("rsID", "REF_ALLELE", "ALLELE_FLIPPED", "STRAND_FLIPPED", "REF_MATCH", "LIFTOVER_MAPPED") %in% names(clean_harmonized)), "clean harmonized output should omit duplicate and process metadata columns")
assert_true(identical(clean_harmonized$SNPID, c("rs1", "rs2")), "clean harmonized output should keep rsid values in SNPID")
assert_true(identical(clean_harmonized$CHR, c("1", "1")), "clean harmonized output should normalize chromosome labels")
assert_true(identical(clean_harmonized$EA, c("A", "G")), "clean harmonized output should normalize alleles")

cat("[SMOKE] harmonization fallback smoke test passed\n")
