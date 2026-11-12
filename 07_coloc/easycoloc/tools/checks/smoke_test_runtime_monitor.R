#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(yaml)
})

source("src/utils_helpers.R")
source("src/utils_format.R")
source("src/utils_runtime.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

tmp_dir <- tempfile(pattern = "easycoloc_runtime_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "abf"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "rds"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "susie"), showWarnings = FALSE)

cfg_fp <- compute_runtime_config_fingerprint(c("config/global.yml", "run_coloc.R"))
initialize_runtime_tracker(tmp_dir, enabled = TRUE, config_fingerprint = cfg_fp)
register_active_run(log_file = file.path(tmp_dir, "run_easycoloc_full_20990101_000000.log"), run_label = "smoke")
write_runtime_heartbeat(stage = "smoke", message_text = "runtime tracker smoke test")

task_id <- build_task_id("SMOKE_GWAS", "chr1:1000:rs1", "prenatal", "GENE1")
runtime_update_task(
  task_id = task_id,
  status = "completed",
  gwas_id = "SMOKE_GWAS",
  locus_id = "chr1:1000:rs1",
  qtl_id = "prenatal",
  phenotype = "GENE1",
  identification_method = "PLINK_clump",
  result_file = file.path(tmp_dir, "abf", "SMOKE_GWAS_rs1_locus_results.csv"),
  pp4 = 0.91,
  n_snps = 42
)

fwrite(
  data.table(
    GWAS_ID = "SMOKE_GWAS",
    QTL_ID = "prenatal",
    Locus = "rs1",
    Phenotype = "GENE1",
    PP4 = 0.91,
    n_snps = 42,
    identification_method = "PLINK_clump"
  ),
  file.path(tmp_dir, "abf", "SMOKE_GWAS_rs1_locus_results.csv")
)

system2("Rscript", c("tools/monitor_easycoloc.R", tmp_dir), stdout = TRUE, stderr = TRUE)

snapshot_file <- file.path(tmp_dir, "runtime", "monitor_snapshot.json")
task_state_file <- file.path(tmp_dir, "runtime", "task_state.tsv")
assert_true(file.exists(snapshot_file), "monitor snapshot missing")
assert_true(file.exists(task_state_file), "task_state.tsv missing")

snapshot <- read_json(snapshot_file, simplifyVector = TRUE)
assert_true(snapshot$outputs$abf == 1, "ABF output count mismatch")
assert_true(snapshot$task_rows >= 1, "task row count mismatch")

cat("[SMOKE] runtime monitor smoke test passed\n")
cat("[SMOKE] output dir:", tmp_dir, "\n")
