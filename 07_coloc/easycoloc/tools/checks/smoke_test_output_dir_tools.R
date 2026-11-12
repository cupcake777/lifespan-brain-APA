#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

source("src/utils_runtime.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

tmp_dir <- tempfile(pattern = "easycoloc_output_tools_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "abf"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "rds"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "susie"), showWarnings = FALSE)

log_file <- file.path(tmp_dir, "run_easycoloc_full_20990101_000000.log")
writeLines(
  c(
    "Processing GWAS Dataset: SMOKE_GWAS",
    "[SUM] Merge complete!",
    "[REPORT] Interactive report saved: smoke_report.html"
  ),
  log_file
)

initialize_runtime_tracker(tmp_dir, enabled = TRUE, config_fingerprint = "smoke_fp")
register_active_run(log_file = log_file, run_label = "output_tools_smoke")
write_runtime_heartbeat(stage = "smoke", message_text = "output-dir tools smoke test")
runtime_update_task(
  task_id = build_task_id("SMOKE_GWAS", "chr1:1000:rs1", "prenatal", "GENE1"),
  status = "completed",
  gwas_id = "SMOKE_GWAS",
  locus_id = "chr1:1000:rs1",
  qtl_id = "prenatal",
  phenotype = "GENE1",
  result_file = file.path(tmp_dir, "abf", "SMOKE_GWAS_rs1_locus_results.csv"),
  pp4 = 0.87,
  n_snps = 31
)

fwrite(
  data.table(
    GWAS_ID = "SMOKE_GWAS",
    QTL_ID = "prenatal",
    Locus = "rs1",
    Phenotype = "GENE1",
    PP4 = 0.87,
    n_snps = 31
  ),
  file.path(tmp_dir, "abf", "SMOKE_GWAS_rs1_locus_results.csv")
)

Sys.setenv(
  EASYCOLOC_GLOBAL_CONFIG = file.path(tmp_dir, "missing_global.yml"),
  EASYCOLOC_GWAS_CONFIG = file.path(tmp_dir, "missing_gwas.yml"),
  EASYCOLOC_QTL_CONFIG = file.path(tmp_dir, "missing_qtl.yml")
)
on.exit(
  Sys.unsetenv(c("EASYCOLOC_GLOBAL_CONFIG", "EASYCOLOC_GWAS_CONFIG", "EASYCOLOC_QTL_CONFIG")),
  add = TRUE
)

snapshot_file <- tempfile(pattern = "easycoloc_snapshot_", fileext = ".json")
monitor_out <- system2(
  "Rscript",
  c("tools/monitor_easycoloc.R", "--snapshot-file", snapshot_file, tmp_dir),
  stdout = TRUE,
  stderr = TRUE
)
assert_true(any(grepl("\\[MONITOR\\] output_dir:", monitor_out)), "monitor output missing output_dir")
assert_true(file.exists(snapshot_file), "monitor snapshot override file missing")

status_out <- system2(
  "Rscript",
  c("tools/summarize_run_status.R", tmp_dir),
  stdout = TRUE,
  stderr = TRUE
)
assert_true(any(grepl("SMOKE_GWAS", status_out)), "status output missing inferred GWAS id")
assert_true(file.exists(file.path(tmp_dir, "run_status_summary.txt")), "status text summary missing")

manifest_out <- system2(
  "Rscript",
  c("tools/build_output_manifest.R", tmp_dir),
  stdout = TRUE,
  stderr = TRUE
)
assert_true(any(grepl("\\[MANIFEST\\] total files:", manifest_out)), "manifest output missing summary")
assert_true(file.exists(file.path(tmp_dir, "output_manifest.tsv")), "manifest file missing")

check_out <- suppressWarnings(system2(
  "Rscript",
  c("tools/check_run_completion.R", tmp_dir),
  stdout = TRUE,
  stderr = TRUE
))
check_status <- attr(check_out, "status")
if (is.null(check_status)) {
  check_status <- 0L
}
assert_true(check_status == 1L, "completion check should report incomplete status for smoke output")
assert_true(any(grepl("\\[CHECK\\] status: INCOMPLETE", check_out)), "completion check output missing status")

tmp_heartbeat_complete_dir <- tempfile(pattern = "easycoloc_output_tools_heartbeat_complete_")
dir.create(tmp_heartbeat_complete_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_heartbeat_complete_dir, "abf"), showWarnings = FALSE)

initialize_runtime_tracker(tmp_heartbeat_complete_dir, enabled = TRUE, config_fingerprint = "smoke_fp_heartbeat_complete")
register_active_run(log_file = NA_character_, run_label = "output_tools_heartbeat_complete")
write_runtime_heartbeat(stage = "pipeline_complete", message_text = "pipeline finished without managed log")
append_runtime_event(
  level = "INFO",
  stage = "summary_complete",
  message_text = "Merged ABF outputs"
)
append_runtime_event(
  level = "INFO",
  stage = "report_complete",
  message_text = "report.html"
)
cat(
  toJSON(
    list(
      timestamp = "2098-01-01 00:00:00",
      pid = 1L,
      level = "ERROR",
      stage = "gwas_failed",
      message = "stale failed event from a previous run",
      gwas_id = "STALE_GWAS",
      locus_id = NULL,
      qtl_id = NULL,
      phenotype = NULL,
      task_id = NULL,
      extras = list()
    ),
    auto_unbox = TRUE,
    null = "null"
  ),
  "\n",
  file = file.path(tmp_heartbeat_complete_dir, "runtime", "event_log.ndjson"),
  append = TRUE
)
fwrite(
  data.table(
    GWAS_ID = "HB_GWAS",
    QTL_ID = "prenatal",
    Locus = "rs1",
    Phenotype = "GENE1",
    PP4 = 0.91,
    n_snps = 42
  ),
  file.path(tmp_heartbeat_complete_dir, "abf", "HB_GWAS_rs1_locus_results.csv")
)

check_heartbeat_complete_out <- suppressWarnings(system2(
  "Rscript",
  c("tools/check_run_completion.R", tmp_heartbeat_complete_dir),
  stdout = TRUE,
  stderr = TRUE
))
check_heartbeat_complete_status <- attr(check_heartbeat_complete_out, "status")
if (is.null(check_heartbeat_complete_status)) {
  check_heartbeat_complete_status <- 0L
}
assert_true(check_heartbeat_complete_status == 0L, "completion check should pass from runtime heartbeat when run log is unavailable")
assert_true(any(grepl("\\[CHECK\\] status: COMPLETE", check_heartbeat_complete_out)), "heartbeat completion check output missing COMPLETE status")
assert_true(!any(grepl("rbindlist", check_heartbeat_complete_out, fixed = TRUE)), "completion check emitted NULL-field rbindlist warning")

tmp_complete_dir <- tempfile(pattern = "easycoloc_output_tools_complete_")
dir.create(tmp_complete_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_complete_dir, "abf"), showWarnings = FALSE)

complete_log_file <- file.path(tmp_complete_dir, "run_easycoloc_full_20990101_010101.log")
writeLines(
  c(
    "Processing GWAS Dataset: DONE_GWAS",
    "[SUM] Merge complete!",
    "[REPORT] Interactive report saved: done_report.html",
    "EasyColoc Analysis Complete!"
  ),
  complete_log_file
)

initialize_runtime_tracker(tmp_complete_dir, enabled = TRUE, config_fingerprint = "smoke_fp_complete")
register_active_run(log_file = complete_log_file, run_label = "output_tools_complete")
write_runtime_heartbeat(stage = "pipeline_complete", message_text = "pipeline finished")
append_runtime_event(
  level = "ERROR",
  stage = "gwas_failed",
  message_text = "subscript out of bounds",
  gwas_id = "DONE_GWAS"
)
runtime_update_task(
  task_id = build_task_id("DONE_GWAS", "chr2:2000:rs2", "postnatal", "GENE2"),
  status = "completed",
  gwas_id = "DONE_GWAS",
  locus_id = "chr2:2000:rs2",
  qtl_id = "postnatal",
  phenotype = "GENE2",
  result_file = file.path(tmp_complete_dir, "abf", "DONE_GWAS_rs2_locus_results.csv"),
  pp4 = 0.91,
  n_snps = 42
)

fwrite(
  data.table(
    GWAS_ID = "DONE_GWAS",
    QTL_ID = "postnatal",
    Locus = "rs2",
    Phenotype = "GENE2",
    PP4 = 0.91,
    n_snps = 42
  ),
  file.path(tmp_complete_dir, "abf", "DONE_GWAS_rs2_locus_results.csv")
)

status_complete_out <- system2(
  "Rscript",
  c("tools/summarize_run_status.R", tmp_complete_dir),
  stdout = TRUE,
  stderr = TRUE
)
assert_true(any(grepl("DONE_GWAS", status_complete_out)), "complete status output missing inferred GWAS id")
status_payload <- read_json(file.path(tmp_complete_dir, "run_status_summary.json"), simplifyVector = TRUE)
assert_true(identical(status_payload$heartbeat$stage, "pipeline_complete"), "status payload missing complete heartbeat")
assert_true(!isTRUE(status_payload$summary$current[[1]]), "completed run should not mark a GWAS as current")

check_complete_out <- suppressWarnings(system2(
  "Rscript",
  c("tools/check_run_completion.R", tmp_complete_dir),
  stdout = TRUE,
  stderr = TRUE
))
check_complete_status <- attr(check_complete_out, "status")
if (is.null(check_complete_status)) {
  check_complete_status <- 0L
}
assert_true(check_complete_status == 1L, "completion check should fail when runtime recorded gwas_failed")
assert_true(any(grepl("\\[CHECK\\] status: FAILED", check_complete_out)), "completion check should report FAILED when runtime recorded gwas_failed")

tmp_failed_dir <- tempfile(pattern = "easycoloc_output_tools_failed_")
dir.create(tmp_failed_dir, recursive = TRUE, showWarnings = FALSE)

failed_log_file <- file.path(tmp_failed_dir, "run_easycoloc_full_20990101_020202.log")
writeLines(
  c(
    "Processing GWAS Dataset: FAIL_GWAS",
    "EasyColoc pipeline finished with failures"
  ),
  failed_log_file
)

initialize_runtime_tracker(tmp_failed_dir, enabled = TRUE, config_fingerprint = "smoke_fp_failed")
register_active_run(log_file = failed_log_file, run_label = "output_tools_failed")
write_runtime_heartbeat(stage = "pipeline_failed", message_text = "pipeline finished with failures")
append_runtime_event(
  level = "ERROR",
  stage = "pipeline_failed",
  message_text = "EasyColoc pipeline finished with failures",
  gwas_id = "FAIL_GWAS"
)

status_failed_out <- system2(
  "Rscript",
  c("tools/summarize_run_status.R", tmp_failed_dir),
  stdout = TRUE,
  stderr = TRUE
)
assert_true(any(grepl("FAIL_GWAS", status_failed_out)), "failed status output missing inferred GWAS id")
failed_payload <- read_json(file.path(tmp_failed_dir, "run_status_summary.json"), simplifyVector = TRUE)
assert_true(identical(failed_payload$heartbeat$stage, "pipeline_failed"), "status payload missing failed heartbeat")
assert_true(!isTRUE(failed_payload$summary$current[[1]]), "failed run should not mark a GWAS as current")

cat("[SMOKE] output-dir tools smoke test passed\n")
cat("[SMOKE] output dir:", tmp_dir, "\n")
