#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_format.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

variants_full <- c("rs1", "rs2", "rs3", "rs4")
ld_res_full <- list(
  R = matrix(
    c(
      1.0, 0.4, 0.3, 0.2,
      0.4, 1.0, 0.5, 0.1,
      0.3, 0.5, 1.0, 0.6,
      0.2, 0.1, 0.6, 1.0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(variants_full, variants_full)
  ),
  snp_info = data.table(
    CHR = c(1, 1, 1, 1),
    SNP = variants_full,
    CM = c(0, 0, 0, 0),
    POS = c(100, 200, 300, 400),
    A1 = c("A", "A", "A", "A"),
    A2 = c("G", "G", "G", "G")
  )
)

ld_cache_store(variants_full, "toy_bfile", keep_file = "toy_keep", ld_res = ld_res_full)

exact_hit <- ld_cache_lookup(c("rs1", "rs2", "rs3", "rs4"), "toy_bfile", keep_file = "toy_keep")
assert_true(!is.null(exact_hit), "exact LD cache hit missing")
assert_true(identical(dim(exact_hit$R), c(4L, 4L)), "exact cache matrix dimensions incorrect")

subset_hit <- ld_cache_lookup(c("rs2", "rs4"), "toy_bfile", keep_file = "toy_keep")
assert_true(!is.null(subset_hit), "subset LD cache hit missing")
assert_true(identical(colnames(subset_hit$R), c("rs2", "rs4")), "subset cache column names incorrect")
assert_true(isTRUE(all.equal(subset_hit$R["rs2", "rs4"], 0.1)), "subset cache value incorrect")

miss_hit <- ld_cache_lookup(c("rs9"), "toy_bfile", keep_file = "toy_keep")
assert_true(is.null(miss_hit), "unexpected LD cache miss result")

variants_requested <- c("rs10", "rs11", "rs12")
variants_returned <- c("rs10", "rs12")
ld_res_subset <- list(
  R = matrix(
    c(
      1.0, 0.7,
      0.7, 1.0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(variants_returned, variants_returned)
  ),
  snp_info = data.table(
    CHR = c(2, 2),
    SNP = variants_returned,
    CM = c(0, 0),
    POS = c(1000, 1200),
    A1 = c("A", "C"),
    A2 = c("G", "T")
  )
)

ld_cache_store(variants_requested, "subset_bfile", keep_file = NULL, ld_res = ld_res_subset)

dropped_hit <- ld_cache_lookup(c("rs11", "rs12"), "subset_bfile", keep_file = NULL)
assert_true(is.null(dropped_hit), "subset cache should miss when a requested SNP was dropped by PLINK")

all_dropped_hit <- ld_cache_lookup(c("rs11"), "subset_bfile", keep_file = NULL)
assert_true(is.null(all_dropped_hit), "all-dropped subset cache should return NULL")

actual_subset_hit <- ld_cache_lookup(c("rs10", "rs12"), "subset_bfile", keep_file = NULL)
assert_true(!is.null(actual_subset_hit), "subset cache should hit for actual PLINK-returned SNPs")
assert_true(identical(colnames(actual_subset_hit$R), c("rs10", "rs12")), "actual subset cache column names incorrect")

tmp_cache_dir <- tempfile("easycoloc_disk_ld_cache_")
dir.create(tmp_cache_dir, recursive = TRUE, showWarnings = FALSE)
old_disk_cache <- getOption("easycoloc.disk_ld_cache_dir")
on.exit(options(easycoloc.disk_ld_cache_dir = old_disk_cache), add = TRUE)
options(easycoloc.disk_ld_cache_dir = tmp_cache_dir)

tmp_ref <- tempfile("easycoloc_plink_ref_")
file.create(paste0(tmp_ref, ".bed"), paste0(tmp_ref, ".bim"), paste0(tmp_ref, ".fam"))
disk_ld_res <- list(
  R = ld_res_subset$R,
  snp_info = data.table(SNP = c("rs10", "rs12"))
)
ld_disk_cache_write(c("rs10", "rs12"), tmp_ref, keep_file = NULL, ld_res = disk_ld_res)
disk_hit <- ld_disk_cache_read(c("rs12", "rs10"), tmp_ref, keep_file = NULL)
assert_true(!is.null(disk_hit), "disk LD cache hit missing")
assert_true(identical(colnames(disk_hit$R), c("rs10", "rs12")), "disk LD cache matrix names incorrect")
writeLines("changed", paste0(tmp_ref, ".bim"))
rm(list = ls(envir = .ld_cache), envir = .ld_cache)
disk_miss_after_ref_change <- ld_disk_cache_read(c("rs12", "rs10"), tmp_ref, keep_file = NULL)
assert_true(is.null(disk_miss_after_ref_change), "disk LD cache should miss after reference signature changes")

cat("[SMOKE] LD cache smoke test passed\n")
