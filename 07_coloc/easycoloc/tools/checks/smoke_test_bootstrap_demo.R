#!/usr/bin/env Rscript

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

tmp_dir <- tempfile(pattern = "easycoloc_demo_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
project_dir <- file.path(tmp_dir, "demo_project")

cmd <- c("easycoloc", "bootstrap-refs", "--demo", project_dir, "--force", "--run")
out <- system2("bash", cmd, stdout = TRUE, stderr = TRUE)
status <- attr(out, "status")
assert_true(is.null(status) || identical(status, 0L), paste(c("demo pipeline failed:", out), collapse = "\n"))

assert_true(file.exists(file.path(project_dir, "config", "global.yml")), "demo global config missing")
assert_true(file.exists(file.path(project_dir, "data", "gwas", "demo_gwas.tsv.gz")), "demo GWAS file missing")
assert_true(file.exists(file.path(project_dir, "data", "qtl", "all_pairs_demo.txt.gz")), "demo allPairs file missing")
assert_true(file.exists(file.path(project_dir, "data", "qtl", "all_pairs_demo.txt.gz.tbi")), "demo allPairs tabix missing")
assert_true(file.exists(file.path(project_dir, "refs", "plink", "hg38", "toy_chr22.bed")), "demo PLINK bed missing")
assert_true(file.exists(file.path(project_dir, "refs", "plink", "keep", "DEMO.sample")), "demo keep file missing")
assert_true(file.exists(file.path(project_dir, "refs", "annotation", "demo_chr22.gtf")), "demo gene annotation missing")
assert_true(file.exists(file.path(project_dir, "results", "coloc_report.html")), "demo HTML report missing")
assert_true(any(grepl("demo_project", out, fixed = TRUE)), "demo bootstrap output missing project path")
assert_true(any(grepl("[LD_plot] Using supplied LD table/matrix with 56 pairs", out, fixed = TRUE)), "demo coloc plot did not receive the complete LD matrix")
assert_true(any(grepl("PP.H4 = 0.9463", out, fixed = TRUE)), "demo coloc output missing expected strong PP.H4")

results_dir <- file.path(project_dir, "results")
all_results <- read.csv(file.path(results_dir, "all_colocalization_results.csv"))
sig_results <- read.csv(file.path(results_dir, "significant_colocalizations_PP4_0.7.csv"))

assert_true(nrow(all_results) == 1, "demo all results should contain one coloc test")
assert_true(nrow(sig_results) == 1, "demo significant results should contain one PP4 >= 0.7 hit")
assert_true(!file.exists(file.path(results_dir, "deduplicated_colocalization_results.csv")), "demo should not write deduplicated results")
assert_true(!file.exists(file.path(results_dir, "significant_unique_trait_phenotype_PP4_0.7.csv")), "demo should not write unique significant results")
assert_true(isTRUE(all_results$n_snps[1] == 8), "demo coloc test should include all 8 SNPs")
assert_true(isTRUE(all_results$PP4[1] > 0.7), "demo PP4 should exceed the configured significance threshold")
assert_true(isTRUE(abs(all_results$PP4[1] - sig_results$PP4[1]) < 1e-8), "demo significant PP4 should match all results")
assert_true(isTRUE(all_results$phenotype_locus_count[1] == 1), "demo should annotate one locus for the phenotype")

report_json <- readLines(file.path(results_dir, "coloc_report_data.json"), warn = FALSE)
assert_true(any(grepl('"total_tests": 1', report_json, fixed = TRUE)), "demo report summary should record one test")
assert_true(any(grepl('"significant_pp4_07": 1', report_json, fixed = TRUE)), "demo report summary should record one PP4 >= 0.7 hit")
assert_true(any(grepl('"mean_n_snps": 8', report_json, fixed = TRUE)), "demo report summary should record 8 mean SNPs")

event_log <- readLines(file.path(results_dir, "runtime", "event_log.ndjson"), warn = FALSE)
assert_true(!any(grepl('"level":"ERROR"', event_log, fixed = TRUE)), "demo runtime event log contains ERROR entries")
heartbeat <- readLines(file.path(results_dir, "runtime", "heartbeat.json"), warn = FALSE)
assert_true(any(grepl('"stage": "pipeline_complete"', heartbeat, fixed = TRUE)), "demo heartbeat should finish at pipeline_complete")
task_state <- read.delim(file.path(results_dir, "runtime", "task_state.tsv"), check.names = FALSE)
assert_true(nrow(task_state) == 1, "demo task state should contain one task")
assert_true(identical(task_state$status[1], "completed"), "demo task should be completed")
assert_true(isTRUE(task_state$pp4[1] > 0.7), "demo task PP4 should exceed the configured significance threshold")

plot_files <- list.files(file.path(results_dir, "plots"), pattern = "\\.pdf$", full.names = TRUE)
plot_png_files <- list.files(file.path(results_dir, "plots"), pattern = "\\.png$", full.names = TRUE)
rds_files <- list.files(file.path(results_dir, "rds"), pattern = "\\.rds$", full.names = TRUE)
layout_files <- list.files(file.path(results_dir, "plot_data"), pattern = "_locuscompare_layout_validation\\.json$", full.names = TRUE)
assert_true(length(plot_files) == 2, "demo should write coloc and locuscompare PDF plots")
assert_true(length(plot_png_files) == 2, "demo should write coloc and locuscompare PNG plots")
assert_true(length(rds_files) == 1, "demo should write one coloc RDS file")
assert_true(all(file.info(plot_files)$size > 5000), "demo PDF plot is unexpectedly small")
assert_true(all(file.info(plot_png_files)$size > 5000), "demo PNG plot is unexpectedly small")
assert_true(file.info(rds_files)$size > 500, "demo RDS file is unexpectedly small")
assert_true(length(layout_files) == 1, "demo should write one locuscompare layout-validation file")
layout_validation <- jsonlite::read_json(layout_files[[1]], simplifyVector = TRUE)
assert_true(identical(layout_validation$status, "PASS"), "demo locuscompare layout validation did not pass")
assert_true(isTRUE(layout_validation$gene_track$included), "demo locuscompare gene track is missing")
assert_true(isTRUE(layout_validation$gene_track$physical_column_aligned), "demo locuscompare gene track is misaligned")
assert_true(identical(layout_validation$ld$missing_ld_count, 0L), "demo locuscompare has missing LD values")
cat("[SMOKE] bootstrap demo smoke test passed\n")
