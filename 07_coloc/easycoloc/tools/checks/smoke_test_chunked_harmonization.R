#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_helpers.R")
source("src/utils_format.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

tmp_dir <- tempfile(pattern = "easycoloc_chunked_harmony_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(tmp_dir, "SMOKE_b38to38_harmonized.tsv.gz")
input_file <- file.path(tmp_dir, "input.tsv")

input_dt <- data.table(
  SNPID = c("rs1", "rs2", "rs3", "rs4"),
  CHR = c("1", "1", "2", "2"),
  POS = c(101L, 102L, 201L, 202L),
  EA = c("A", "C", "G", "T"),
  NEA = c("G", "T", "A", "C"),
  EAF = c(0.1, 0.2, 0.3, 0.4),
  BETA = c(0.1, -0.2, 0.3, -0.4),
  SE = c(0.05, 0.06, 0.07, 0.08),
  P = c(0.05, 0.01, 0.001, 0.2),
  N = c(100, 100, 100, 100)
)
data.table::fwrite(input_dt, input_file, sep = "\t", na = "NA", quote = FALSE)

res <- run_easycoloc_harmonization(
  input_dt,
  ref_fasta = file.path(tmp_dir, "missing.fa"),
  ref_vcf = NULL,
  ref_dbsnp = NULL,
  source_build = "38",
  target_build = "38",
  save_dir = tmp_dir,
  dataset_id = "SMOKE",
  input_file = input_file,
  sample_size_n = 100,
  chunked = TRUE,
  chunk_min_rows = 1L,
  chunk_parallel_jobs = 2L,
  verbose = FALSE
)

part_dir <- paste0(out_file, ".parts")
assert_true(file.exists(out_file), "final harmonized file missing")
assert_true(!dir.exists(part_dir), "part directory should be removed after successful final cache write")
assert_true(identical(names(res), easycoloc_harmonized_gwas_output_cols()), "output schema mismatch")
assert_true(nrow(res) == 4L, "output row count mismatch")

mtime_before <- file.info(out_file)$mtime
Sys.sleep(1)
res_cached <- run_easycoloc_harmonization(
  input_dt,
  ref_fasta = file.path(tmp_dir, "missing.fa"),
  ref_vcf = NULL,
  ref_dbsnp = NULL,
  source_build = "38",
  target_build = "38",
  save_dir = tmp_dir,
  dataset_id = "SMOKE",
  input_file = input_file,
  sample_size_n = 100,
  chunked = TRUE,
  chunk_min_rows = 1L,
  chunk_parallel_jobs = 2L,
  verbose = FALSE
)
mtime_after <- file.info(out_file)$mtime
assert_true(identical(mtime_before, mtime_after), "completed final cache should be reused")
assert_true(!dir.exists(part_dir), "part directory should remain absent when final cache is reused")
assert_true(nrow(res_cached) == 4L, "cached output row count mismatch")

unlink(out_file)
legacy_out_file <- file.path(tmp_dir, "SMOKE_b38to38_harmonized.tsv")
data.table::fwrite(res, legacy_out_file, sep = "\t", na = "NA", quote = FALSE)
Sys.sleep(1)
res_legacy_cached <- run_easycoloc_harmonization(
  input_dt,
  ref_fasta = file.path(tmp_dir, "missing.fa"),
  ref_vcf = NULL,
  ref_dbsnp = NULL,
  source_build = "38",
  target_build = "38",
  save_dir = tmp_dir,
  dataset_id = "SMOKE",
  input_file = input_file,
  sample_size_n = 100,
  chunked = TRUE,
  chunk_min_rows = 1L,
  chunk_parallel_jobs = 2L,
  verbose = FALSE
)
assert_true(nrow(res_legacy_cached) == 4L, "legacy tsv cache fallback row count mismatch")
assert_true(file.exists(out_file), "legacy cache fallback should migrate to preferred gz cache")
assert_true(!file.exists(legacy_out_file), "legacy cache should be removed after preferred gz migration")

cat("[SMOKE] chunked harmonization smoke test passed\n")
cat("[SMOKE] output dir:", tmp_dir, "\n")
