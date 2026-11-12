#!/usr/bin/env Rscript

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

tmp_dir <- tempfile(pattern = "easycoloc_1kg_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
ftp_root <- file.path(tmp_dir, "ftp")
release_dir <- file.path(ftp_root, "release", "20130502")
phase3_dir <- file.path(ftp_root, "phase3")
release38_dir <- file.path(ftp_root, "data_collections", "1000_genomes_project", "release", "20190312_biallelic_SNV_and_INDEL")
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(phase3_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(release38_dir, recursive = TRUE, showWarnings = FALSE)

panel_path <- file.path(phase3_dir, "20130502.integrated_call_samples_v3.20130502.ALL.panel")
panel_path_official <- file.path(release_dir, "integrated_call_samples_v3.20130502.ALL.panel")
writeLines(
  c(
    "sample\tpop\tsuper_pop\tgender",
    "NA1\tCHB\tEAS\tmale",
    "NA2\tJPT\tEAS\tfemale",
    "NA3\tCEU\tEUR\tmale",
    "NA4\tTSI\tEUR\tfemale"
  ),
  panel_path
)
invisible(file.copy(panel_path, panel_path_official, overwrite = TRUE))

vcf_path <- file.path(release_dir, "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
vcf_lines <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\tNA3\tNA4",
  "22\t200100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0",
  "22\t200200\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/1\t0/0\t1/1",
  "22\t200300\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/0"
)
plain_vcf <- sub("\\.gz$", "", vcf_path)
writeLines(vcf_lines, plain_vcf)
system2("bgzip", c("-f", plain_vcf), stdout = FALSE, stderr = FALSE)

vcf38_path <- file.path(release38_dir, "ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
plain_vcf38 <- sub("\\.gz$", "", vcf38_path)
writeLines(vcf_lines, plain_vcf38)
system2("bgzip", c("-f", plain_vcf38), stdout = FALSE, stderr = FALSE)

out_dir <- file.path(tmp_dir, "output_hg19")
cmd <- c(
  "easycoloc",
  "bootstrap-refs",
  "--setup-1kg", out_dir,
  "--build", "hg19",
  "--pop", "EAS",
  "--chromosomes", "22",
  "--base-url", paste0("file://", ftp_root)
)
out <- system2("bash", cmd, stdout = TRUE, stderr = TRUE)

assert_true(any(grepl("setup_1kg", out, fixed = TRUE)), "1KG setup output missing")
assert_true(file.exists(file.path(out_dir, "raw_vcf", basename(vcf_path))), "downloaded VCF missing")
assert_true(file.exists(file.path(out_dir, "plink", "1kg_phase3_hg19_chr22.bed")), "per-chromosome PLINK bed missing")
assert_true(file.exists(file.path(out_dir, "keep", "EAS.sample")), "population keep file missing")
assert_true(file.exists(file.path(out_dir, "plink", "1kg_phase3_hg19_EAS.bed")), "population PLINK panel missing")

out_dir38 <- file.path(tmp_dir, "output_hg38")
cmd38 <- c(
  "easycoloc",
  "bootstrap-refs",
  "--setup-1kg", out_dir38,
  "--build", "hg38",
  "--pop", "EAS",
  "--chromosomes", "22",
  "--base-url", paste0("file://", ftp_root)
)
out38 <- system2("bash", cmd38, stdout = TRUE, stderr = TRUE)
assert_true(any(grepl("hg38", out38, fixed = TRUE)), "hg38 setup output missing")
assert_true(file.exists(file.path(out_dir38, "raw_vcf", basename(vcf38_path))), "downloaded hg38 VCF missing")
assert_true(file.exists(file.path(out_dir38, "plink", "1kg_phase3_hg38_chr22.bed")), "hg38 per-chromosome PLINK bed missing")
assert_true(file.exists(file.path(out_dir38, "plink", "1kg_phase3_hg38_EAS.bed")), "hg38 population PLINK panel missing")

project_dir <- file.path(tmp_dir, "project")
invisible(system2("bash", c("easycoloc", "init", project_dir), stdout = TRUE, stderr = TRUE))
cmd_rewrite <- c(
  "easycoloc",
  "bootstrap-refs",
  "--setup-1kg", file.path(project_dir, "refs", "1kg_hg19"),
  "--build", "hg19",
  "--pop", "EAS",
  "--chromosomes", "22",
  "--base-url", paste0("file://", ftp_root),
  "--global", file.path(project_dir, "config", "global.yml"),
  "--gwas", file.path(project_dir, "config", "gwas.yml"),
  "--qtl", file.path(project_dir, "config", "qtl.yml"),
  "--rewrite-config"
)
invisible(system2("bash", cmd_rewrite, stdout = TRUE, stderr = TRUE))
cfg_lines <- readLines(file.path(project_dir, "config", "global.yml"))
assert_true(any(grepl("plink_hg19: refs/1kg_hg19/plink/1kg_phase3_hg19_EAS", cfg_lines, fixed = TRUE)), "plink_hg19 rewrite missing")
assert_true(any(grepl("plink_keep: refs/1kg_hg19/keep/EAS.sample", cfg_lines, fixed = TRUE)), "plink_keep rewrite missing")
cat("[SMOKE] 1KG setup smoke test passed\n")
