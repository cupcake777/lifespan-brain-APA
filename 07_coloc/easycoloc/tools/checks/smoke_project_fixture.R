create_smoke_project_fixture <- function(pattern,
                                         gwas_build = "hg38",
                                         recom = "ref/recom_map",
                                         include_results = FALSE,
                                         qtl_paths = c("relative", "absolute")) {
  qtl_paths <- match.arg(qtl_paths)

  tmp_dir <- tempfile(pattern = pattern)
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "config"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "data", "gwas"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "data", "qtl"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "ref"), recursive = TRUE, showWarnings = FALSE)
  if (include_results) {
    dir.create(file.path(tmp_dir, "results"), recursive = TRUE, showWarnings = FALSE)
  }

  plink_prefix <- file.path(tmp_dir, "ref", "toy_plink")
  invisible(file.create(paste0(plink_prefix, c(".bed", ".bim", ".fam"))))
  invisible(file.create(file.path(
    tmp_dir,
    "ref",
    c("keep.sample", "gene.gtf.gz", "ref19.fa.gz", "ref38.fa.gz", "dbsnp19.gz", "dbsnp38.gz", "hg19ToHg38.over.chain.gz")
  )))
  if (nzchar(recom)) {
    invisible(file.create(file.path(tmp_dir, recom)))
  }

  dir.create(file.path(tmp_dir, "ref", "af_ref"), recursive = TRUE, showWarnings = FALSE)
  invisible(file.create(file.path(tmp_dir, "ref", "af_ref", "1KG_hg19_EUR_AF.tsv.gz")))
  invisible(file.create(file.path(tmp_dir, "data", "gwas", "example.tsv.gz")))
  invisible(file.create(file.path(tmp_dir, "data", "qtl", c("all_pairs.txt.gz", "sig_pairs.txt.gz"))))

  all_pairs <- "all_pairs.txt.gz"
  sig_pairs <- "sig_pairs.txt.gz"
  if (identical(qtl_paths, "absolute")) {
    all_pairs <- file.path(tmp_dir, "data", "qtl", all_pairs)
    sig_pairs <- file.path(tmp_dir, "data", "qtl", sig_pairs)
  }
  writeLines(
    c(
      "Type,NumberRNASeqandGTSamples,allPairsTabixFilename,sigPairsTabixFilename",
      paste("example_eqtl,500", all_pairs, sig_pairs, sep = ",")
    ),
    file.path(tmp_dir, "data", "qtl", "QTL_summary.csv")
  )

  global_cfg <- list(
    project_root = "..",
    output_dir = "results",
    temp_dir = "temp",
    plink_hg38 = "ref/toy_plink",
    plink_keep = "ref/keep.sample",
    hash_table_dir = "",
    recom = recom,
    gene_anno = "ref/gene.gtf.gz",
    ref_genome_hg19 = "ref/ref19.fa.gz",
    ref_genome_hg38 = "ref/ref38.fa.gz",
    `1kg_af` = "ref/af_ref",
    dbsnp_hg19 = "ref/dbsnp19.gz",
    dbsnp_hg38 = "ref/dbsnp38.gz",
    harmonize_dir = "harmony",
    clump = list(p1 = 5e-8, p2 = 5e-8, r2 = 0.1, kb = 1000),
    coloc_settings = list(pp4_threshold = 0.7, susie_threshold = 0.5, min_snps = 10),
    harmonization_settings = list(liftover_chain = "ref/hg19ToHg38.over.chain.gz")
  )
  yaml::write_yaml(global_cfg, file.path(tmp_dir, "config", "global.yml"))

  gwas_cfg <- list(
    datasets = list(
      list(
        id = "EXAMPLE",
        name = "Example",
        file = "data/gwas/example.tsv.gz",
        build = gwas_build,
        pop = "EUR",
        type = "cc",
        sample_size_n = 10000,
        prop = 0.2,
        columns = list(
          snp = "SNP",
          chrom = "CHR",
          pos = "BP",
          a1 = "A1",
          a2 = "A2",
          pval = "P",
          beta = "BETA",
          se = "SE"
        )
      )
    )
  )
  yaml::write_yaml(gwas_cfg, file.path(tmp_dir, "config", "gwas.yml"))

  qtl_cfg <- list(
    qtl_info = list(
      build = "hg38",
      file = "data/qtl/QTL_summary.csv",
      columns = list(
        id = "Type",
        sample_size = "NumberRNASeqandGTSamples",
        all_filename = "allPairsTabixFilename",
        sig_filename = "sigPairsTabixFilename"
      )
    ),
    QTL_all_header = c("chr", "pos", "phenotype_id", "variant_id", "pval_nominal", "slope", "slope_se"),
    QTL_cols = list(
      phenotype = "phenotype_id",
      chrom = "chr",
      pos = "pos",
      pval = "pval_nominal",
      beta = "slope",
      se = "slope_se"
    )
  )
  yaml::write_yaml(qtl_cfg, file.path(tmp_dir, "config", "qtl.yml"))

  list(
    root = tmp_dir,
    global = file.path(tmp_dir, "config", "global.yml"),
    gwas = file.path(tmp_dir, "config", "gwas.yml"),
    qtl = file.path(tmp_dir, "config", "qtl.yml")
  )
}
