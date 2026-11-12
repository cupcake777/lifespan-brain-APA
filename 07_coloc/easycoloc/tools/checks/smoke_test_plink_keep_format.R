#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_format.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

tmp_dir <- tempfile(pattern = "easycoloc_plink_keep_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

bfile <- file.path(tmp_dir, "toy")
fam <- data.table(
  FID = c("0", "0", "F3"),
  IID = c("HG001", "HG002", "HG003"),
  father = 0,
  mother = 0,
  sex = 0,
  phenotype = -9
)
fwrite(fam, paste0(bfile, ".fam"), sep = "\t", col.names = FALSE)

one_col_keep <- file.path(tmp_dir, "one_col.sample")
writeLines(c("HG001", "HG003"), one_col_keep)

converted <- easycoloc_prepare_plink_keep_file(
  one_col_keep,
  bfile = bfile,
  temp_dir = tmp_dir,
  label = "smoke",
  verbose = FALSE
)
converted_dt <- fread(converted, header = FALSE, col.names = c("FID", "IID"))
assert_true(ncol(converted_dt) == 2L, "converted keep file must have two columns")
assert_true(identical(converted_dt$FID, c("0", "F3")), "converted keep file must preserve FID from .fam")
assert_true(identical(converted_dt$IID, c("HG001", "HG003")), "converted keep file IID mismatch")

two_col_keep <- file.path(tmp_dir, "two_col.sample")
fwrite(data.table(FID = c("0"), IID = c("HG002")), two_col_keep, sep = "\t", col.names = FALSE)
unchanged <- easycoloc_prepare_plink_keep_file(
  two_col_keep,
  bfile = bfile,
  temp_dir = tmp_dir,
  label = "smoke2",
  verbose = FALSE
)
assert_true(identical(unchanged, two_col_keep), "two-column keep file should be used unchanged")

cat("[SMOKE] PLINK keep format smoke test passed\n")
