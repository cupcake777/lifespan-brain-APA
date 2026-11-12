#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

tmp_dir <- tempfile(pattern = "easycoloc_bootstrap_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
project_dir <- file.path(tmp_dir, "project")
shared_dir <- file.path(tmp_dir, "shared")
dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(shared_dir, "plink", "hg38"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(shared_dir, "plink", "keep"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(shared_dir, "af_ref"), recursive = TRUE, showWarnings = FALSE)

invisible(file.create(file.path(shared_dir, "plink", "hg38", paste0("1kg_hg38_filtered", c(".bed", ".bim", ".fam")))))
invisible(file.create(file.path(shared_dir, "plink", "keep", "EAS.sample")))
invisible(file.create(file.path(shared_dir, "af_ref", "1KG_hg19_EUR_AF.tsv.gz")))
invisible(file.create(file.path(shared_dir, "af_ref", "1KG_hg19_EAS_AF.tsv.gz")))
invisible(file.create(file.path(shared_dir, "af_ref", "1KG_hg38_EAS_AF.tsv.gz")))
invisible(file.create(file.path(shared_dir, c("GRCh37.fa.gz", "GRCh38.fa.gz", "dbsnp19.gz", "dbsnp38.gz", "genes.gtf.gz", "hg19ToHg38.over.chain.gz"))))
invisible(file.create(file.path(shared_dir, c("recomb_recombination_map_hapmap_format_hg38_chr_1.txt"))))

dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
init_res <- system2("bash", c("easycoloc", "init", project_dir), stdout = TRUE, stderr = TRUE)
assert_true(length(init_res) >= 1, "project init output missing")

global_path <- file.path(project_dir, "config", "global.yml")
gwas_path <- file.path(project_dir, "config", "gwas.yml")
qtl_path <- file.path(project_dir, "config", "qtl.yml")

gwas_cfg <- yaml::read_yaml(gwas_path)
gwas_cfg$datasets[[1]]$build <- "hg19"
gwas_cfg$datasets[[1]]$pop <- "EUR"
yaml::write_yaml(gwas_cfg, gwas_path)

source_cfg <- list(
  defaults = list(mode = "symlink"),
  resources = list(
    plink_hg38 = list(
      source = file.path(shared_dir, "plink", "hg38", "1kg_hg38_filtered"),
      target = "refs/plink/hg38/1kg_hg38_filtered"
    ),
    plink_keep = list(
      source = file.path(shared_dir, "plink", "keep", "EAS.sample"),
      target = "refs/plink/keep/EAS.sample"
    ),
    ref_genome_hg19 = list(
      source = file.path(shared_dir, "GRCh37.fa.gz"),
      target = "refs/fasta/GRCh37.fa.gz"
    ),
    ref_genome_hg38 = list(
      source = file.path(shared_dir, "GRCh38.fa.gz"),
      target = "refs/fasta/GRCh38.fa.gz"
    ),
    dbsnp_hg19 = list(
      source = file.path(shared_dir, "dbsnp19.gz"),
      target = "refs/dbsnp/dbsnp_hg19.gz"
    ),
    dbsnp_hg38 = list(
      source = file.path(shared_dir, "dbsnp38.gz"),
      target = "refs/dbsnp/dbsnp_hg38.gz"
    ),
    `1kg_af` = list(
      source = file.path(shared_dir, "af_ref"),
      target = "refs/af"
    ),
    gene_anno = list(
      source = file.path(shared_dir, "genes.gtf.gz"),
      target = "refs/gene/genes.gtf.gz"
    ),
    recom = list(
      source_prefix = file.path(shared_dir, "recomb"),
      target = "refs/recomb/recomb"
    ),
    liftover_chain = list(
      source = file.path(shared_dir, "hg19ToHg38.over.chain.gz"),
      target = "refs/liftover/hg19ToHg38.over.chain.gz"
    )
  )
)
source_yaml <- file.path(project_dir, "config", "reference_sources.test.yml")
yaml::write_yaml(source_cfg, source_yaml)

cmd <- c(
  "easycoloc",
  "bootstrap-refs",
  source_yaml,
  "--global", global_path,
  "--gwas", gwas_path,
  "--qtl", qtl_path,
  "--rewrite-config"
)
bootstrap_out <- system2("bash", cmd, stdout = TRUE, stderr = TRUE)
assert_true(any(grepl("\\[BOOTSTRAP\\].*plink_hg38", bootstrap_out)), "bootstrap did not materialize plink_hg38")

global_cfg <- yaml::read_yaml(global_path)
assert_true(identical(global_cfg$plink_hg38, "refs/plink/hg38/1kg_hg38_filtered"), "global.yml plink_hg38 not rewritten")
assert_true(identical(global_cfg$`1kg_af`, "refs/af"), "global.yml 1kg_af not rewritten")
assert_true(nzchar(Sys.readlink(file.path(project_dir, "refs", "plink", "hg38", "1kg_hg38_filtered.bed"))), "plink bed symlink missing")
assert_true(nzchar(Sys.readlink(file.path(project_dir, "refs", "af"))), "af directory should be a symlink")
assert_true(file.exists(file.path(project_dir, "reference_bootstrap_manifest.tsv")), "bootstrap manifest missing")
cat("[SMOKE] bootstrap reference smoke test passed\n")
