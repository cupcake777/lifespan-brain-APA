#!/usr/bin/env Rscript

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

tmp_dir <- tempfile(pattern = "easycoloc_gtex_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
meta_dir <- file.path(tmp_dir, "meta")
eqtl_dir <- file.path(tmp_dir, "eqtl")
sqtl_dir <- file.path(tmp_dir, "sqtl")
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(eqtl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sqtl_dir, recursive = TRUE, showWarnings = FALSE)

sample_attr <- file.path(tmp_dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
writeLines(
  c(
    "SAMPID\tSMTSD",
    "GTEX-AAA\tWhole Blood",
    "GTEX-BBB\tWhole Blood",
    "GTEX-CCC\tMuscle - Skeletal"
  ),
  sample_attr
)

invisible(file.create(file.path(eqtl_dir, "Whole_Blood.allpairs.tab.gz")))
invisible(file.create(file.path(eqtl_dir, "Whole_Blood.v8.signif_variant_gene_pairs.tab.gz")))
invisible(file.create(file.path(sqtl_dir, "Muscle_Skeletal.v8.tab.gz")))
invisible(file.create(file.path(sqtl_dir, "Muscle_Skeletal.v8.signif_variant_gene_pairs.tab.gz")))

cmd <- c(
  "easycoloc",
  "bootstrap-refs",
  "--fetch-gtex-meta", meta_dir,
  "--gtex-url", paste0("file://", sample_attr),
  "--gtex-eqtl-dir", eqtl_dir,
  "--gtex-sqtl-dir", sqtl_dir
)
out <- system2("bash", cmd, stdout = TRUE, stderr = TRUE)

eqtl_summary <- file.path(meta_dir, "GTEx_v8_eQTL_summary.csv")
sqtl_summary <- file.path(meta_dir, "GTEx_v8_sQTL_summary.csv")
yaml_file <- file.path(meta_dir, "qtl_gtex_generated.yml")
assert_true(file.exists(file.path(meta_dir, basename(sample_attr))), "downloaded GTEx sample attributes missing")
assert_true(file.exists(eqtl_summary), "GTEx eQTL summary missing")
assert_true(file.exists(sqtl_summary), "GTEx sQTL summary missing")
assert_true(file.exists(yaml_file), "generated GTEx YAML missing")

eqtl_lines <- readLines(eqtl_summary)
sqtl_lines <- readLines(sqtl_summary)
yaml_lines <- readLines(yaml_file)
assert_true(any(grepl("Whole Blood", eqtl_lines, fixed = TRUE)), "Whole Blood row missing from eQTL summary")
assert_true(any(grepl("Muscle - Skeletal", sqtl_lines, fixed = TRUE)), "Muscle - Skeletal row missing from sQTL summary")
assert_true(any(grepl("GTEx_v8_eQTL_summary.csv", yaml_lines, fixed = TRUE)), "generated YAML does not point to eQTL summary")
assert_true(any(grepl("phenotype_id", yaml_lines, fixed = TRUE)), "generated YAML missing eQTL column mapping")
assert_true(any(grepl("build_gtex_summary", out, fixed = TRUE)), "GTEx bootstrap output missing summary step")
assert_true(any(grepl("write_qtl_yaml", out, fixed = TRUE)), "GTEx bootstrap output missing YAML step")
cat("[SMOKE] GTEx bootstrap smoke test passed\n")
