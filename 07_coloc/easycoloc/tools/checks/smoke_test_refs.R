#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

source("tools/checks/smoke_project_fixture.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

fixture <- create_smoke_project_fixture(
  pattern = "easycoloc_refs_",
  gwas_build = "hg19",
  recom = "",
  qtl_paths = "relative"
)

cmd <- c(
  "tools/list_reference_requirements.R",
  "--global", fixture$global,
  "--gwas", fixture$gwas,
  "--qtl", fixture$qtl,
  "--include-qtl-files"
)
refs_out <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)

assert_true(any(grepl("\\[REF\\].*plink_hg38", refs_out)), "plink_hg38 requirement missing")
assert_true(any(grepl("\\[REF\\].*qtl_summary", refs_out)), "qtl_summary requirement missing")
assert_true(any(grepl("qtl_allpairs_example_eqtl", refs_out, fixed = TRUE)), "resolved QTL allPairs entry missing")
cat("[SMOKE] reference listing smoke test passed\n")
