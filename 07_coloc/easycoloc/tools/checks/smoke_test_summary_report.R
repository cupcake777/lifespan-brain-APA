#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_format.R")
source("src/utils_report.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

tmp_dir <- tempfile(pattern = "easycoloc_summary_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tmp_dir, "abf"), showWarnings = FALSE)
dir.create(file.path(tmp_dir, "susie"), showWarnings = FALSE)

abf_dt1 <- data.table(
  GWAS_ID = "SMOKE_GWAS",
  QTL_ID = "postnatal",
  Locus = "rs10890255",
  Phenotype = "ENST00000654683.1|CCDC30|chr1:42482663-42484158|+",
  PP4 = 0.91,
  n_snps = 86,
  identification_method = "Distance_pruning"
)

abf_dt2 <- data.table(
  GWAS_ID = "SMOKE_GWAS",
  QTL_ID = "prenatal",
  Locus = "rs2819340",
  Phenotype = "ENST00000260746.6|ARL3|chr10:102673731-102676941|-",
  PP4 = 0.63,
  n_snps = 54,
  identification_method = "PLINK_clump"
)

abf_dt3 <- data.table(
  GWAS_ID = "SMOKE_GWAS",
  QTL_ID = "postnatal",
  Locus = "rs10890256",
  Phenotype = "ENST00000654683.1|CCDC30|chr1:42482663-42484158|+",
  PP4 = 0.88,
  n_snps = 79,
  identification_method = "PLINK_clump"
)

fwrite(abf_dt1, file.path(tmp_dir, "abf", "SMOKE_GWAS_rs10890255_locus_results.csv"))
fwrite(abf_dt2, file.path(tmp_dir, "abf", "SMOKE_GWAS_rs2819340_locus_results.csv"))
fwrite(abf_dt3, file.path(tmp_dir, "abf", "SMOKE_GWAS_rs10890256_locus_results.csv"))

susie_dt <- data.table(
  cs_id = 1L,
  variable = "rs10890255",
  variable_prob = 0.77,
  cs_log10bf = 12.4
)
fwrite(susie_dt, file.path(tmp_dir, "susie", "SMOKE_GWAS_postnatal_CCDC30_susie.csv"))

merged <- merge_all_results(
  output_dir = tmp_dir,
  pp4_threshold = 0.7,
  merge_susie = TRUE,
  save_summary = TRUE
)

assert_true(!is.null(merged), "merge_all_results() returned NULL")
assert_true(nrow(merged) == 3, "merge_all_results() row count mismatch")
assert_true(file.exists(file.path(tmp_dir, "all_colocalization_results.csv")), "merged results file missing")
assert_true(file.exists(file.path(tmp_dir, "significant_colocalizations_PP4_0.7.csv")), "significant results file missing")
assert_true(!file.exists(file.path(tmp_dir, "deduplicated_colocalization_results.csv")), "deduplicated results file should not be written")
assert_true(!file.exists(file.path(tmp_dir, "significant_unique_trait_phenotype_PP4_0.7.csv")), "deduplicated significant results file should not be written")
assert_true(file.exists(file.path(tmp_dir, "all_susie_results.csv")), "merged SuSiE file missing")
assert_true(file.exists(file.path(tmp_dir, "summary_statistics.txt")), "summary statistics file missing")

all_results <- fread(file.path(tmp_dir, "all_colocalization_results.csv"))
ccdc30_rows <- all_results[Phenotype == "ENST00000654683.1|CCDC30|chr1:42482663-42484158|+"]
assert_true(nrow(ccdc30_rows) == 2, "all results should retain both CCDC30 locus-level rows")
assert_true(all(ccdc30_rows$phenotype_locus_count == 2), "phenotype_locus_count mismatch for duplicated phenotype")
assert_true(all(ccdc30_rows$phenotype_loci == "rs10890255;rs10890256"), "phenotype_loci annotation mismatch")

report_ok <- generate_html_report(
  results_dir = tmp_dir,
  output_file = file.path(tmp_dir, "coloc_report.html"),
  project_name = "EasyColoc Smoke Test"
)

assert_true(isTRUE(report_ok), "generate_html_report() returned FALSE")
assert_true(file.exists(file.path(tmp_dir, "coloc_report.html")), "HTML report file missing")
assert_true(file.exists(file.path(tmp_dir, "coloc_report_data.json")), "report JSON data file missing")

cat("[SMOKE] summary/report smoke test passed\n")
cat("[SMOKE] output dir:", tmp_dir, "\n")
