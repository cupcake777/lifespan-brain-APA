#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

source("src/utils_config.R")
source("tools/checks/smoke_project_fixture.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

env_root <- tempfile(pattern = "easycoloc_env_root_")
Sys.setenv(EASYCOLOC_SMOKE_ROOT = env_root)
expanded_paths <- easycoloc_expand_env_vars(c(
  "$EASYCOLOC_SMOKE_ROOT/ref",
  "${EASYCOLOC_SMOKE_ROOT}/data",
  "plain/path"
))
assert_true(identical(expanded_paths[[1]], file.path(env_root, "ref")), "bare env var expansion failed")
assert_true(identical(expanded_paths[[2]], file.path(env_root, "data")), "braced env var expansion failed")
assert_true(identical(expanded_paths[[3]], "plain/path"), "plain path should not be rewritten")

fixture <- create_smoke_project_fixture(
  pattern = "easycoloc_doctor_",
  gwas_build = "hg38",
  recom = "ref/recom_map",
  include_results = TRUE,
  qtl_paths = "absolute"
)

gwas_cfg <- yaml::read_yaml(fixture$gwas)
second_gwas <- gwas_cfg$datasets[[1]]
second_gwas$id <- "SECOND"
second_gwas$name <- "Second"
gwas_cfg$datasets[[2]] <- second_gwas
yaml::write_yaml(gwas_cfg, fixture$gwas)

qtl_meta <- data.table::fread(file.path(fixture$root, "data", "qtl", "QTL_summary.csv"))
second_qtl <- data.table::copy(qtl_meta[1])
second_qtl[, Type := "second_eqtl"]
data.table::fwrite(rbind(qtl_meta, second_qtl), file.path(fixture$root, "data", "qtl", "QTL_summary.csv"))

cmd <- c(
  "tools/doctor_easycoloc.R",
  "--global", fixture$global,
  "--gwas", fixture$gwas,
  "--qtl", fixture$qtl
)
doctor_out <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)

assert_true(any(grepl("\\[DOCTOR\\] summary:", doctor_out)), "doctor summary line missing")
assert_true(!any(grepl("\\[FAIL\\]", doctor_out)), "doctor unexpectedly reported FAIL")

Sys.setenv(EASYCOLOC_GWAS_ID = "EXAMPLE", EASYCOLOC_QTL_ID = "example_eqtl")
Sys.setenv(
  EASYCOLOC_GLOBAL_CONFIG = fixture$global,
  EASYCOLOC_GWAS_CONFIG = fixture$gwas,
  EASYCOLOC_QTL_CONFIG = fixture$qtl
)
selected_cfg <- easycoloc_read_configs()
assert_true(length(selected_cfg$gwas$datasets) == 1L, "GWAS selector did not retain exactly one dataset")
assert_true(identical(selected_cfg$gwas$datasets[[1]]$id, "EXAMPLE"), "GWAS selector retained the wrong dataset")
selected_qtl <- easycoloc_filter_qtl_meta(qtl_meta, selected_cfg$qtl)
assert_true(nrow(selected_qtl) == 1L && identical(selected_qtl$Type[[1]], "example_eqtl"), "QTL selector retained the wrong dataset")

selected_doctor_out <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)
assert_true(!any(grepl("\\[FAIL\\]", selected_doctor_out)), "selected-pair doctor unexpectedly reported FAIL")
assert_true(any(grepl("selected_datasets", selected_doctor_out)), "selected QTL validation was not reported")

unknown_gwas_error <- tryCatch(
  {
    Sys.setenv(EASYCOLOC_GWAS_ID = "NOT_CONFIGURED")
    easycoloc_read_configs()
    ""
  },
  error = function(e) e$message
)
assert_true(grepl("Available IDs: EXAMPLE, SECOND", unknown_gwas_error, fixed = TRUE), "unknown GWAS ID did not fail with available IDs")
Sys.unsetenv(c(
  "EASYCOLOC_GWAS_ID", "EASYCOLOC_QTL_ID", "EASYCOLOC_GLOBAL_CONFIG",
  "EASYCOLOC_GWAS_CONFIG", "EASYCOLOC_QTL_CONFIG"
))
cat("[SMOKE] doctor smoke test passed\n")
