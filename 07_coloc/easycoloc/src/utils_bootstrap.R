suppressPackageStartupMessages({
  library(data.table)
  library(glue)
  library(yaml)
})

if (file.exists("src/utils_download.R")) {
  source("src/utils_download.R")
}

easycoloc_bootstrap_specs <- function(project_root) {
  data.table(
    resource = c(
      "plink_hg19", "plink_hg38", "plink_keep", "ref_genome_hg19", "ref_genome_hg38",
      "dbsnp_hg19", "dbsnp_hg38", "1kg_af", "gene_anno",
      "recom", "liftover_chain", "hash_table_dir"
    ),
    field = c(
      "plink_hg19", "plink_hg38", "plink_keep", "ref_genome_hg19", "ref_genome_hg38",
      "dbsnp_hg19", "dbsnp_hg38", "1kg_af", "gene_anno",
      "recom", "harmonization_settings.liftover_chain", "hash_table_dir"
    ),
    target_kind = c(
      "plink_prefix", "plink_prefix", "file", "file", "file",
      "file", "file", "directory", "file",
      "prefix_bundle", "file", "directory"
    ),
    target_rel = c(
      "refs/plink/hg19/1kg_hg19_filtered",
      "refs/plink/hg38/1kg_hg38_filtered",
      "refs/plink/keep/default.sample",
      "refs/fasta/GRCh37.fa.gz",
      "refs/fasta/GRCh38.fa.gz",
      "refs/dbsnp/dbsnp_hg19.gz",
      "refs/dbsnp/dbsnp_hg38.gz",
      "refs/af",
      "refs/gene/gene_annotation.gtf.gz",
      "refs/recomb/recomb",
      "refs/liftover/hg19ToHg38.over.chain.gz",
      "refs/hash_tables"
    ),
    project_root = project_root
  )
}

easycoloc_bootstrap_get_nested <- function(x, field) {
  parts <- strsplit(field, ".", fixed = TRUE)[[1]]
  value <- x
  for (part in parts) {
    if (is.null(value[[part]])) {
      return(NULL)
    }
    value <- value[[part]]
  }
  value
}

easycoloc_bootstrap_set_nested <- function(x, field, value) {
  parts <- strsplit(field, ".", fixed = TRUE)[[1]]
  if (length(parts) == 1) {
    x[[field]] <- value
    return(x)
  }
  current <- x[[parts[1]]]
  if (is.null(current)) {
    current <- list()
  }
  remaining <- paste(parts[-1], collapse = ".")
  current <- easycoloc_bootstrap_set_nested(current, remaining, value)
  x[[parts[1]]] <- current
  x
}

easycoloc_remove_path <- function(path) {
  if (!file.exists(path) && !dir.exists(path) && !nzchar(Sys.readlink(path))) {
    return(invisible(NULL))
  }
  info <- suppressWarnings(file.info(path))
  if (!is.na(info$isdir[1]) && isTRUE(info$isdir[1])) {
    unlink(path, recursive = TRUE, force = TRUE)
  } else {
    unlink(path, recursive = FALSE, force = TRUE)
  }
  invisible(NULL)
}

easycoloc_ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

easycoloc_materialize_file <- function(source, target, mode = "symlink", force = FALSE) {
  if (!file.exists(source)) {
    stop("Source file missing: ", source, call. = FALSE)
  }
  if (file.exists(target) || nzchar(Sys.readlink(target))) {
    if (!isTRUE(force)) {
      return("exists")
    }
    easycoloc_remove_path(target)
  }
  easycoloc_ensure_parent(target)
  ok <- if (identical(mode, "copy")) {
    file.copy(source, target, overwrite = FALSE)
  } else {
    file.symlink(source, target)
  }
  if (!isTRUE(ok)) {
    stop("Failed to materialize file: ", target, call. = FALSE)
  }
  "created"
}

easycoloc_materialize_directory <- function(source, target, mode = "symlink", force = FALSE) {
  if (!dir.exists(source)) {
    stop("Source directory missing: ", source, call. = FALSE)
  }
  if (file.exists(target) || dir.exists(target) || nzchar(Sys.readlink(target))) {
    if (!isTRUE(force)) {
      return("exists")
    }
    easycoloc_remove_path(target)
  }
  easycoloc_ensure_parent(target)
  if (identical(mode, "copy")) {
    stop("Directory copy is not supported; use mode=symlink for directories", call. = FALSE)
  }
  ok <- file.symlink(source, target)
  if (!isTRUE(ok)) {
    stop("Failed to create directory symlink: ", target, call. = FALSE)
  }
  "created"
}

easycoloc_materialize_plink_prefix <- function(source_prefix, target_prefix, mode = "symlink", force = FALSE) {
  suffixes <- c(".bed", ".bim", ".fam")
  statuses <- vapply(suffixes, function(suffix) {
    easycoloc_materialize_file(
      source = paste0(source_prefix, suffix),
      target = paste0(target_prefix, suffix),
      mode = mode,
      force = force
    )
  }, character(1))
  if (any(statuses == "created")) "created" else "exists"
}

easycoloc_materialize_prefix_bundle <- function(source_prefix, target_prefix, mode = "symlink", force = FALSE) {
  source_dir <- dirname(source_prefix)
  source_base <- basename(source_prefix)
  bundle_files <- list.files(source_dir, pattern = paste0("^", source_base), full.names = TRUE)
  if (length(bundle_files) == 0) {
    stop("No files found for source prefix: ", source_prefix, call. = FALSE)
  }
  target_dir <- dirname(target_prefix)
  target_base <- basename(target_prefix)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  statuses <- vapply(bundle_files, function(source_file) {
    suffix <- substring(basename(source_file), nchar(source_base) + 1L)
    target_file <- file.path(target_dir, paste0(target_base, suffix))
    easycoloc_materialize_file(source_file, target_file, mode = mode, force = force)
  }, character(1))
  if (any(statuses == "created")) "created" else "exists"
}

easycoloc_build_hash_tables <- function(bed_dir, target_dir, build = "hg38", force = FALSE) {
  if (!dir.exists(bed_dir)) {
    stop("dbSNP BED directory missing: ", bed_dir, call. = FALSE)
  }
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(force)) {
    unlink(file.path(target_dir, "chr_*_hash_table.rds"))
  }
  cmd <- c("tools/bed2rds.R", "-i", bed_dir, "-o", target_dir, "-b", build)
  exit_code <- system2("Rscript", cmd)
  if (!identical(exit_code, 0L)) {
    stop("Hash table build failed for: ", bed_dir, call. = FALSE)
  }
  "created"
}

easycoloc_bootstrap_apply <- function(cfg_bundle, source_cfg, rewrite_config = FALSE, force = FALSE) {
  project_root <- cfg_bundle$project_root
  specs <- easycoloc_bootstrap_specs(project_root)
  resources <- source_cfg$resources
  if (is.null(resources) || length(resources) == 0) {
    stop("No resources defined in source YAML", call. = FALSE)
  }

  default_mode <- if (!is.null(source_cfg$defaults$mode)) source_cfg$defaults$mode else "symlink"
  global_cfg <- yaml::read_yaml(cfg_bundle$paths$global)
  manifest <- data.table(
    resource = character(),
    field = character(),
    target = character(),
    target_rel = character(),
    action = character(),
    mode = character(),
    source = character()
  )

  for (i in seq_len(nrow(specs))) {
    spec <- specs[i]
    entry <- resources[[spec$resource]]
    if (is.null(entry)) {
      next
    }

    mode <- if (!is.null(entry$mode)) as.character(entry$mode) else default_mode
    target_rel <- if (!is.null(entry$target)) as.character(entry$target) else as.character(spec$target_rel)
    target <- normalizePath(file.path(project_root, target_rel), mustWork = FALSE)
    action <- "skipped"
    source_desc <- NA_character_

    if (!is.null(entry$build_from_bed_dir)) {
      source_desc <- as.character(entry$build_from_bed_dir)
      action <- easycoloc_build_hash_tables(
        bed_dir = source_desc,
        target_dir = target,
        build = if (!is.null(entry$build)) as.character(entry$build) else "hg38",
        force = force
      )
    } else if (identical(spec$target_kind, "directory")) {
      source_desc <- as.character(entry$source)
      action <- easycoloc_materialize_directory(source_desc, target, mode = mode, force = force)
    } else if (identical(spec$target_kind, "file")) {
      source_desc <- as.character(entry$source)
      action <- easycoloc_materialize_file(source_desc, target, mode = mode, force = force)
    } else if (identical(spec$target_kind, "plink_prefix")) {
      source_desc <- as.character(entry$source)
      action <- easycoloc_materialize_plink_prefix(source_desc, target, mode = mode, force = force)
    } else if (identical(spec$target_kind, "prefix_bundle")) {
      source_desc <- if (!is.null(entry$source_prefix)) as.character(entry$source_prefix) else as.character(entry$source)
      action <- easycoloc_materialize_prefix_bundle(source_desc, target, mode = mode, force = force)
    }

    manifest <- rbind(
      manifest,
      data.table(
        resource = as.character(spec$resource),
        field = as.character(spec$field),
        target = target,
        target_rel = target_rel,
        action = action,
        mode = mode,
        source = source_desc
      ),
      fill = TRUE
    )

    if (isTRUE(rewrite_config) && action %in% c("created", "exists")) {
      global_cfg <- easycoloc_bootstrap_set_nested(global_cfg, spec$field, target_rel)
    }
  }

  if (isTRUE(rewrite_config)) {
    yaml::write_yaml(global_cfg, cfg_bundle$paths$global)
  }

  manifest
}

easycoloc_create_project_skeleton <- function(dest_dir, force = FALSE) {
  if (dir.exists(dest_dir) && length(list.files(dest_dir, all.files = TRUE, no.. = TRUE)) > 0 && !isTRUE(force)) {
    stop("Destination exists and is not empty: ", dest_dir, call. = FALSE)
  }
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest_dir, "config"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest_dir, "data", "gwas"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest_dir, "data", "qtl"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest_dir, "results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest_dir, "temp"), recursive = TRUE, showWarnings = FALSE)
  normalizePath(dest_dir, mustWork = FALSE)
}

easycoloc_write_demo_ped_map <- function(prefix) {
  positions <- c(200100, 200200, 200300, 200400, 200500, 200600, 200700, 200800)
  map_dt <- data.table(
    CHR = rep(22, 8),
    SNP = paste0("chr22:", positions, ":A:G"),
    CM = 0,
    POS = positions
  )
  ped_rows <- list(
    c("D1", "D1", 0, 0, 1, 1, "A", "A", "A", "G", "G", "G", "A", "G", "A", "A", "A", "G", "G", "G", "A", "G"),
    c("D2", "D2", 0, 0, 2, 1, "A", "G", "A", "G", "G", "G", "A", "A", "A", "G", "A", "A", "G", "G", "A", "A"),
    c("D3", "D3", 0, 0, 1, 1, "G", "G", "A", "A", "A", "G", "A", "A", "A", "A", "A", "G", "A", "G", "A", "G"),
    c("D4", "D4", 0, 0, 2, 1, "A", "A", "A", "A", "A", "A", "A", "G", "A", "G", "A", "A", "A", "A", "A", "A"),
    c("D5", "D5", 0, 0, 1, 1, "A", "G", "G", "G", "G", "G", "A", "G", "A", "G", "A", "G", "G", "G", "A", "G"),
    c("D6", "D6", 0, 0, 2, 1, "A", "A", "A", "G", "A", "G", "A", "A", "A", "A", "A", "A", "A", "G", "A", "A")
  )
  ped_dt <- as.data.table(do.call(rbind, ped_rows))
  fwrite(map_dt, paste0(prefix, ".map"), sep = "\t", col.names = FALSE)
  fwrite(ped_dt, paste0(prefix, ".ped"), sep = "\t", col.names = FALSE)
}

easycoloc_create_demo_qtl <- function(project_dir) {
  qtl_dir <- file.path(project_dir, "data", "qtl")
  dir.create(qtl_dir, recursive = TRUE, showWarnings = FALSE)
  phenotype_id <- "ENSTDEMO001|GENEDEMO|chr22:200000-201000|+"
  positions <- c(200100, 200200, 200300, 200400, 200500, 200600, 200700, 200800)
  variant_ids <- paste0("chr22:", positions, ":A:G")
  qtl_beta <- c(0.90, 0.15, 0.08, 0.05, 0.03, 0.02, 0.01, 0.005)
  qtl_p <- c(5e-12, 2e-01, 3e-01, 4e-01, 5e-01, 6e-01, 7e-01, 8e-01)
  qtl_dt <- data.table(
    chr = rep("chr22", 8),
    pos = positions,
    ref = rep("A", 8),
    alt = rep("G", 8),
    phenotype_id = rep(phenotype_id, 8),
    variant_id = variant_ids,
    start_distance = c(100, 200, 300, 400, 500, 600, 700, 800),
    end_distance = c(900, 800, 700, 600, 500, 400, 300, 200),
    af = c(0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.24, 0.26),
    ma_samples = c(2, 2, 2, 3, 3, 3, 4, 4),
    ma_count = c(2, 2, 3, 3, 4, 4, 5, 5),
    pval_nominal = qtl_p,
    slope = qtl_beta,
    slope_se = abs(qtl_beta / qnorm(qtl_p / 2, lower.tail = FALSE))
  )
  all_pairs_tsv <- file.path(qtl_dir, "all_pairs_demo.txt")
  sig_pairs_tsv <- file.path(qtl_dir, "sig_pairs_demo.txt")
  fwrite(qtl_dt, all_pairs_tsv, sep = "\t", quote = FALSE, col.names = FALSE)
  fwrite(qtl_dt[pval_nominal <= 1e-4], sig_pairs_tsv, sep = "\t", quote = FALSE, col.names = FALSE)

  for (txt_file in c(all_pairs_tsv, sig_pairs_tsv)) {
    easycoloc_run_command("bgzip", c("-f", txt_file), fail_message = glue("bgzip failed for {txt_file}"))
    gz_file <- paste0(txt_file, ".gz")
    easycoloc_run_command("tabix", c("-f", "-s", "1", "-b", "2", "-e", "2", gz_file), fail_message = glue("tabix failed for {gz_file}"))
  }

  summary_dt <- data.table(
    Type = "demo_eqtl",
    NumberRNASeqandGTSamples = 1000,
    allPairsTabixFilename = "all_pairs_demo.txt.gz",
    sigPairsTabixFilename = "sig_pairs_demo.txt.gz"
  )
  fwrite(summary_dt, file.path(qtl_dir, "QTL_summary.demo.csv"))
}

easycoloc_write_demo_gwas <- function(project_dir) {
  gwas_dir <- file.path(project_dir, "data", "gwas")
  dir.create(gwas_dir, recursive = TRUE, showWarnings = FALSE)
  positions <- c(200100, 200200, 200300, 200400, 200500, 200600, 200700, 200800)
  gwas_beta <- c(0.90, 0.15, 0.08, 0.05, 0.03, 0.02, 0.01, 0.005)
  gwas_p <- c(5e-12, 2e-01, 3e-01, 4e-01, 5e-01, 6e-01, 7e-01, 8e-01)
  gwas_dt <- data.table(
    SNP = paste0("chr22:", positions, ":A:G"),
    CHR = rep(22, 8),
    BP = positions,
    A1 = rep("G", 8),
    A2 = rep("A", 8),
    P = gwas_p,
    BETA = gwas_beta,
    SE = abs(gwas_beta / qnorm(gwas_p / 2, lower.tail = FALSE)),
    EAF = c(0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.24, 0.26),
    N = 5000
  )
  fwrite(gwas_dt, file.path(gwas_dir, "demo_gwas.tsv.gz"), sep = "\t", quote = FALSE)
}

easycoloc_write_demo_gtf <- function(project_dir) {
  annotation_dir <- file.path(project_dir, "refs", "annotation")
  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
  gtf_path <- file.path(annotation_dir, "demo_chr22.gtf")
  gtf_lines <- c(
    'chr22\tEasyColoc\texon\t200050\t200180\t.\t+\t.\tgene_id "GENEDEMO1"; gene_name "GENEDEMO1"; gene_type "protein_coding";',
    'chr22\tEasyColoc\texon\t200300\t200460\t.\t+\t.\tgene_id "GENEDEMO1"; gene_name "GENEDEMO1"; gene_type "protein_coding";',
    'chr22\tEasyColoc\texon\t200560\t200680\t.\t-\t.\tgene_id "GENEDEMO2"; gene_name "GENEDEMO2"; gene_type "protein_coding";',
    'chr22\tEasyColoc\texon\t200730\t200900\t.\t-\t.\tgene_id "GENEDEMO2"; gene_name "GENEDEMO2"; gene_type "protein_coding";'
  )
  writeLines(gtf_lines, gtf_path)
  gtf_path
}

easycoloc_write_demo_configs <- function(project_dir) {
  global_cfg <- list(
    project_root = "..",
    output_dir = "results",
    temp_dir = "temp",
    sig_threshold = 0.0,
    use_parallel = FALSE,
    n_cores = 1,
    random_seed = 20240326,
    plink_hg38 = "refs/plink/hg38/toy_chr22",
    plink_keep = "refs/plink/keep/DEMO.sample",
    hash_table_dir = "",
    recom = "",
    gene_anno = "refs/annotation/demo_chr22.gtf",
    ref_genome_hg19 = "",
    ref_genome_hg38 = "",
    `1kg_af` = "refs/af",
    dbsnp_hg19 = "",
    dbsnp_hg38 = "",
    harmonize_dir = "harmony",
    clump = list(p1 = 1e-5, p2 = 1e-5, r2 = 0.2, kb = 500),
    coloc_settings = list(
      pp4_threshold = 0.7,
      susie_threshold = 0.0,
      min_snps = 5,
      min_snps_susie = 5,
      top_candidates = 50,
      maf_default = 0.1,
      maf_na_replacement = 0.05,
      maf_epsilon = 1e-6,
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = 1e-5,
      pvalue_floor = 1e-300
    ),
    analysis = list(
      flank_bp = 1000,
      plot_window_bp = 1000,
      prefilter_sig_pairs = TRUE
    ),
    runtime = list(
      enabled = TRUE,
      resume_completed_tasks = TRUE,
      skip_existing_locus_results = FALSE,
      write_task_events = TRUE,
      parallel_qtl_tasks = FALSE,
      precompute_locus_ld = FALSE,
      ld_plink_timeout = 120
    ),
    harmonization_settings = list(
      liftover_chain = "",
      chunked = TRUE,
      chunk_min_rows = 1000000,
      chunk_parallel_jobs = 1
    ),
    plot_settings = list(
      significance_threshold = 1e-5,
      plot_window_bp = 1000,
      plot_width = 6.0,
      plot_height = 4.5,
      title_phenotype_field = "gene",
      r2_breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
      r2_color_scale = "viridis",
      r2_colors = c("#000080", "#87CEEB", "#008000", "#FFA500", "#FF0000"),
      lead_snp_color = "#7F3C8D",
      credible_set_color = "#E67E22"
    )
  )
  yaml::write_yaml(global_cfg, file.path(project_dir, "config", "global.yml"))

  gwas_cfg <- list(
    datasets = list(
      list(
        id = "DEMO_GWAS",
        name = "Demo GWAS",
        file = "data/gwas/demo_gwas.tsv.gz",
        build = "hg38",
        pop = "DEMO",
        type = "cc",
        sample_size_n = 5000,
        prop = 0.40,
        columns = list(
          snp = "SNP",
          chrom = "CHR",
          pos = "BP",
          a1 = "A1",
          a2 = "A2",
          pval = "P",
          n = "N",
          af = "EAF",
          beta = "BETA",
          se = "SE"
        )
      )
    )
  )
  yaml::write_yaml(gwas_cfg, file.path(project_dir, "config", "gwas.yml"))

  qtl_cfg <- list(
    qtl_info = list(
      build = "hg38",
      file = "data/qtl/QTL_summary.demo.csv",
      columns = list(
        id = "Type",
        sample_size = "NumberRNASeqandGTSamples",
        all_filename = "allPairsTabixFilename",
        sig_filename = "sigPairsTabixFilename"
      )
    ),
    QTL_all_header = c(
      "chr", "pos", "ref", "alt", "phenotype_id", "variant_id",
      "start_distance", "end_distance", "af", "ma_samples", "ma_count",
      "pval_nominal", "slope", "slope_se"
    ),
    QTL_cols = list(
      phenotype = "phenotype_id",
      chrom = "chr",
      pos = "pos",
      pval = "pval_nominal",
      beta = "slope",
      se = "slope_se",
      af = "af"
    )
  )
  yaml::write_yaml(qtl_cfg, file.path(project_dir, "config", "qtl.yml"))
}

easycoloc_bootstrap_demo <- function(dest_dir, force = FALSE, run_pipeline = FALSE) {
  easycoloc_require_command(
    "plink",
    install_hint = "Install PLINK 1.9+ and ensure 'plink' is available in $PATH."
  )
  easycoloc_require_command(
    "bgzip",
    install_hint = "Install htslib/bgzip to create demo QTL tabix files."
  )
  easycoloc_require_command(
    "tabix",
    install_hint = "Install htslib/tabix to index demo QTL files."
  )

  project_dir <- easycoloc_create_project_skeleton(dest_dir, force = force)
  dir.create(file.path(project_dir, "refs", "plink", "hg38"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(project_dir, "refs", "plink", "keep"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(project_dir, "refs", "af"), recursive = TRUE, showWarnings = FALSE)

  demo_prefix <- file.path(project_dir, "temp", "toy_chr22_source")
  dir.create(dirname(demo_prefix), recursive = TRUE, showWarnings = FALSE)
  easycoloc_write_demo_ped_map(demo_prefix)

  out_prefix <- file.path(project_dir, "refs", "plink", "hg38", "toy_chr22")
  easycoloc_run_command(
    "plink",
    c("--file", demo_prefix, "--make-bed", "--out", out_prefix, "--silent"),
    fail_message = "PLINK failed while building the toy chr22 reference panel."
  )
  fam_dt <- fread(paste0(out_prefix, ".fam"), header = FALSE)
  keep_path <- file.path(project_dir, "refs", "plink", "keep", "DEMO.sample")
  fwrite(fam_dt[, .(V1, V2)], keep_path, sep = "\t", col.names = FALSE)

  easycoloc_write_demo_gwas(project_dir)
  easycoloc_create_demo_qtl(project_dir)
  easycoloc_write_demo_gtf(project_dir)
  easycoloc_write_demo_configs(project_dir)

  demo_readme <- c(
    "# EasyColoc Toy Demo",
    "",
    "This project is self-contained and should finish in under 2 minutes on a normal workstation.",
    "",
    "Run:",
    sprintf(
      "bash %s run --global %s --gwas %s --qtl %s",
      shQuote(file.path(normalizePath(getwd(), mustWork = FALSE), "easycoloc")),
      shQuote(file.path(project_dir, "config", "global.yml")),
      shQuote(file.path(project_dir, "config", "gwas.yml")),
      shQuote(file.path(project_dir, "config", "qtl.yml"))
    ),
    "",
    "Or from the EasyColoc repo root:",
    sprintf(
      "./easycoloc run --global %s --gwas %s --qtl %s",
      shQuote(file.path(project_dir, "config", "global.yml")),
      shQuote(file.path(project_dir, "config", "gwas.yml")),
      shQuote(file.path(project_dir, "config", "qtl.yml"))
    )
  )
  writeLines(demo_readme, file.path(project_dir, "DEMO.md"))

  if (isTRUE(run_pipeline)) {
    easycoloc_run_command(
      "bash",
      c(
        "easycoloc", "run",
        "--global", file.path(project_dir, "config", "global.yml"),
        "--gwas", file.path(project_dir, "config", "gwas.yml"),
        "--qtl", file.path(project_dir, "config", "qtl.yml")
      ),
      fail_message = "Toy demo pipeline execution failed."
    )
  }

  data.table(
    mode = "demo",
    project_dir = project_dir,
    global_config = file.path(project_dir, "config", "global.yml"),
    gwas_config = file.path(project_dir, "config", "gwas.yml"),
    qtl_config = file.path(project_dir, "config", "qtl.yml"),
    plink_prefix = out_prefix,
    keep_file = keep_path
  )
}

easycoloc_1kg_panel_url <- function(base_url) {
  paste0(base_url, "/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
}

easycoloc_1kg_vcf_url <- function(base_url, chromosome, build = "hg19") {
  if (identical(build, "hg38")) {
    return(
      paste0(
        base_url,
        "/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/",
        "ALL.chr", chromosome,
        ".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
      )
    )
  }
  paste0(
    base_url,
    "/release/20130502/ALL.chr", chromosome,
    ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  )
}

easycoloc_make_1kg_keep_file <- function(panel_file, fam_file, population, keep_file) {
  panel_dt <- fread(panel_file)
  fam_dt <- fread(fam_file, header = FALSE)
  setnames(fam_dt, c("FID", "IID", "PID", "MID", "SEX", "PHENO"))
  eligible <- panel_dt[panel_dt$pop == population | panel_dt$super_pop == population, .(sample)]
  if (nrow(eligible) == 0) {
    stop("No samples found for population: ", population, call. = FALSE)
  }
  keep_dt <- fam_dt[IID %in% eligible$sample, .(FID, IID)]
  if (nrow(keep_dt) == 0) {
    stop("Population samples do not overlap the PLINK .fam file for: ", population, call. = FALSE)
  }
  dir.create(dirname(keep_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(keep_dt, keep_file, sep = "\t", col.names = FALSE)
  keep_file
}

easycoloc_setup_1kg <- function(dest_dir,
                                pop = "ALL",
                                build = "hg19",
                                chromosomes = "1-22",
                                base_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp",
                                panel_url = NULL,
                                force = FALSE,
                                rewrite_config = FALSE,
                                cfg_bundle = NULL) {
  easycoloc_require_command(
    "plink",
    install_hint = "Install PLINK 1.9+ and ensure 'plink' is available in $PATH."
  )
  if (!build %in% c("hg19", "hg38")) {
    stop("Unsupported build: ", build, " (expected hg19 or hg38)", call. = FALSE)
  }
  downloader <- easycoloc_pick_downloader()
  message(glue("[1KG] Using downloader: {downloader}"))

  chroms <- easycoloc_parse_chromosomes(chromosomes)
  out_root <- normalizePath(dest_dir, mustWork = FALSE)
  raw_dir <- file.path(out_root, "raw_vcf")
  plink_dir <- file.path(out_root, "plink")
  keep_dir <- file.path(out_root, "keep")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plink_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(keep_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(panel_url)) {
    panel_url <- easycoloc_1kg_panel_url(base_url)
  }
  panel_file <- file.path(out_root, basename(panel_url))
  easycoloc_download_file(panel_url, panel_file, force = force)

  per_chr_prefixes <- character()
  for (chrom in chroms) {
    vcf_url <- easycoloc_1kg_vcf_url(base_url, chrom, build = build)
    vcf_file <- file.path(raw_dir, basename(vcf_url))
    easycoloc_download_file(vcf_url, vcf_file, force = force)
    chr_prefix <- file.path(plink_dir, paste0("1kg_phase3_", build, "_chr", chrom))
    if (isTRUE(force) || !all(file.exists(paste0(chr_prefix, c(".bed", ".bim", ".fam"))))) {
      easycoloc_run_command(
        "plink",
        c(
          "--vcf", vcf_file,
          "--double-id",
          "--set-missing-var-ids", "@:#:$1:$2",
          "--biallelic-only", "strict",
          "--snps-only", "just-acgt",
          "--make-bed",
          "--out", chr_prefix,
          "--silent"
        ),
        fail_message = glue("PLINK VCF conversion failed for chr{chrom}")
      )
    }
    per_chr_prefixes <- c(per_chr_prefixes, chr_prefix)
  }

  merged_prefix <- if (length(per_chr_prefixes) == 1) {
    per_chr_prefixes[[1]]
  } else {
    merge_list <- file.path(plink_dir, "merge_list.txt")
    writeLines(
      vapply(
        per_chr_prefixes[-1],
        function(prefix) paste0(prefix, ".bed ", prefix, ".bim ", prefix, ".fam"),
        character(1)
      ),
      con = merge_list
    )
    final_prefix <- file.path(plink_dir, paste0("1kg_phase3_", build, "_merged"))
    if (isTRUE(force) || !all(file.exists(paste0(final_prefix, c(".bed", ".bim", ".fam"))))) {
      easycoloc_run_command(
        "plink",
        c(
          "--bfile", per_chr_prefixes[[1]],
          "--merge-list", merge_list,
          "--make-bed",
          "--out", final_prefix,
          "--silent"
        ),
        fail_message = "PLINK merge failed for 1000 Genomes chromosomes."
      )
    }
    final_prefix
  }

  keep_file <- NA_character_
  pop_prefix <- merged_prefix
  if (!identical(pop, "ALL")) {
    keep_file <- file.path(keep_dir, paste0(pop, ".sample"))
    easycoloc_make_1kg_keep_file(
      panel_file = panel_file,
      fam_file = paste0(merged_prefix, ".fam"),
      population = pop,
      keep_file = keep_file
    )
    pop_prefix <- file.path(plink_dir, paste0("1kg_phase3_", build, "_", pop))
    if (isTRUE(force) || !all(file.exists(paste0(pop_prefix, c(".bed", ".bim", ".fam"))))) {
      easycoloc_run_command(
        "plink",
        c(
          "--bfile", merged_prefix,
          "--keep", keep_file,
          "--make-bed",
          "--out", pop_prefix,
          "--silent"
        ),
        fail_message = glue("PLINK keep/subset failed for population {pop}")
      )
    }
  }

  manifest <- data.table(
    mode = "setup_1kg",
    build = build,
    output_dir = out_root,
    panel_file = panel_file,
    merged_prefix = merged_prefix,
    population = pop,
    population_prefix = pop_prefix,
    keep_file = keep_file
  )

  if (isTRUE(rewrite_config) && !is.null(cfg_bundle)) {
    global_cfg <- yaml::read_yaml(cfg_bundle$paths$global)
    project_root <- cfg_bundle$project_root
    plink_field <- if (identical(build, "hg19")) "plink_hg19" else "plink_hg38"
    global_cfg[[plink_field]] <- sub(paste0("^", project_root, "/?"), "", pop_prefix)
    global_cfg$plink_keep <- if (!is.na(keep_file)) {
      sub(paste0("^", project_root, "/?"), "", keep_file)
    } else {
      global_cfg$plink_keep
    }
    yaml::write_yaml(global_cfg, cfg_bundle$paths$global)
  }

  manifest
}

easycoloc_slugify_tissue <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

easycoloc_rel_path <- function(path, start) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  start_norm <- normalizePath(start, winslash = "/", mustWork = FALSE)
  path_parts <- strsplit(sub("^/", "", path_norm), "/", fixed = TRUE)[[1]]
  start_parts <- strsplit(sub("^/", "", start_norm), "/", fixed = TRUE)[[1]]
  common <- 0L
  max_common <- min(length(path_parts), length(start_parts))
  while (common < max_common && identical(path_parts[common + 1L], start_parts[common + 1L])) {
    common <- common + 1L
  }
  up <- rep("..", length(start_parts) - common)
  down <- path_parts[(common + 1L):length(path_parts)]
  parts <- c(up, down)
  if (length(parts) == 0) "." else paste(parts, collapse = "/")
}

easycoloc_build_gtex_summary <- function(sample_attributes_file,
                                         qtl_dir,
                                         qtl_type = c("eqtl", "sqtl"),
                                         output_file) {
  qtl_type <- match.arg(qtl_type)
  if (!file.exists(sample_attributes_file)) {
    stop("GTEx sample attributes file not found: ", sample_attributes_file, call. = FALSE)
  }
  if (!dir.exists(qtl_dir)) {
    stop("GTEx QTL directory not found: ", qtl_dir, call. = FALSE)
  }

  sample_dt <- fread(sample_attributes_file)
  if (!"SMTSD" %in% names(sample_dt)) {
    stop("Sample attributes file does not contain SMTSD column", call. = FALSE)
  }
  tissue_counts <- sample_dt[
    ,
    .(NumberRNASeqSamples = .N, NumberRNASeqandGTSamples = .N),
    by = .(Tissue = SMTSD)
  ]
  tissue_counts[, tissue_slug := easycoloc_slugify_tissue(Tissue)]

  if (identical(qtl_type, "eqtl")) {
    all_files <- list.files(qtl_dir, pattern = "\\.allpairs\\.tab\\.gz$", full.names = TRUE)
    sig_files <- list.files(qtl_dir, pattern = "\\.v8\\.signif_variant_gene_pairs\\.tab\\.gz$", full.names = TRUE)
    extract_all_slug <- function(x) sub("\\.allpairs\\.tab\\.gz$", "", basename(x))
  } else {
    all_files <- list.files(qtl_dir, pattern = "\\.v8\\.tab\\.gz$", full.names = TRUE)
    all_files <- all_files[!grepl("\\.signif_variant_gene_pairs\\.tab\\.gz$", all_files)]
    sig_files <- list.files(qtl_dir, pattern = "\\.v8\\.signif_variant_gene_pairs\\.tab\\.gz$", full.names = TRUE)
    extract_all_slug <- function(x) sub("\\.v8\\.tab\\.gz$", "", basename(x))
  }
  extract_sig_slug <- function(x) sub("\\.v8\\.signif_variant_gene_pairs\\.tab\\.gz$", "", basename(x))

  file_dt <- data.table(
    tissue_slug = unique(c(vapply(all_files, extract_all_slug, character(1)), vapply(sig_files, extract_sig_slug, character(1))))
  )
  if (nrow(file_dt) == 0) {
    stop("No GTEx files found in directory: ", qtl_dir, call. = FALSE)
  }
  file_dt[, all_file := NA_character_]
  file_dt[, sig_file := NA_character_]
  if (length(all_files) > 0) {
    file_dt[match(vapply(all_files, extract_all_slug, character(1)), tissue_slug), all_file := all_files]
  }
  if (length(sig_files) > 0) {
    file_dt[match(vapply(sig_files, extract_sig_slug, character(1)), tissue_slug), sig_file := sig_files]
  }

  summary_dt <- merge(file_dt, tissue_counts, by = "tissue_slug", all.x = TRUE)
  summary_dt[, Tissue := ifelse(is.na(Tissue), gsub("_", " ", tissue_slug), Tissue)]
  summary_dir <- dirname(normalizePath(output_file, mustWork = FALSE))
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  summary_dt[, allPairsTabixFilename := ifelse(is.na(all_file), "NA", vapply(all_file, easycoloc_rel_path, character(1), start = summary_dir))]
  summary_dt[, sigPairsTabixFilename := ifelse(is.na(sig_file), "NA", vapply(sig_file, easycoloc_rel_path, character(1), start = summary_dir))]
  summary_dt[, Number_of_sigpheno := NA_integer_]
  setcolorder(
    summary_dt,
    c(
      "Tissue", "NumberRNASeqandGTSamples", "NumberRNASeqSamples", "Number_of_sigpheno",
      "allPairsTabixFilename", "sigPairsTabixFilename", "tissue_slug"
    )
  )
  fwrite(summary_dt[, !c("all_file", "sig_file")], output_file)
  output_file
}

easycoloc_write_gtex_qtl_yaml <- function(summary_csv,
                                          output_file,
                                          qtl_type = c("eqtl", "sqtl")) {
  qtl_type <- match.arg(qtl_type)
  summary_rel <- normalizePath(summary_csv, mustWork = FALSE)

  if (identical(qtl_type, "eqtl")) {
    payload <- list(
      qtl_info = list(
        build = "hg38",
        file = summary_rel,
        columns = list(
          id = "Tissue",
          sample_size = "NumberRNASeqandGTSamples",
          all_filename = "allPairsTabixFilename",
          sig_filename = "sigPairsTabixFilename"
        )
      ),
      QTL_all_header = c(
        "chr", "pos", "ref", "alt", "variant_id", "phenotype_id",
        "start_distance", "end_distance", "af", "ma_samples", "ma_count",
        "pval_nominal", "slope", "slope_se"
      ),
      QTL_cols = list(
        phenotype = "phenotype_id",
        chrom = "chr",
        pos = "pos",
        pval = "pval_nominal",
        beta = "slope",
        se = "slope_se"
      )
    )
  } else {
    payload <- list(
      qtl_info = list(
        build = "hg38",
        file = summary_rel,
        columns = list(
          id = "Tissue",
          sample_size = "NumberRNASeqandGTSamples",
          all_filename = "allPairsTabixFilename",
          sig_filename = "sigPairsTabixFilename"
        )
      ),
      QTL_all_header = c(
        "chrom_b38", "chromStart_b38", "chromEnd_b38", "eGeneID",
        "intron_chr", "intron_bp_first", "intron_bp_end", "intron_clu",
        "A1_sqtl", "A2_sqtl", "build", "tss_distance", "ma_samples",
        "ma_count", "maf", "pvalue_sQTL", "slope", "slope_se"
      ),
      QTL_cols = list(
        phenotype = "eGeneID",
        chrom = "chrom_b38",
        pos = "chromEnd_b38",
        pval = "pvalue_sQTL",
        beta = "slope",
        se = "slope_se"
      )
    )
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(payload, output_file)
  output_file
}

easycoloc_fetch_gtex_meta <- function(dest_dir,
                                      url = "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                                      force = FALSE,
                                      eqtl_dir = NULL,
                                      sqtl_dir = NULL,
                                      config_output_file = NULL,
                                      config_qtl_type = c("eqtl", "sqtl")) {
  config_qtl_type <- match.arg(config_qtl_type)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest_file <- file.path(dest_dir, basename(url))
  easycoloc_download_file(url, dest_file, force = force)

  outputs <- list()
  outputs[[1]] <- data.table(mode = "fetch_gtex_meta", file = dest_file, url = url)
  generated_eqtl_summary <- NULL
  generated_sqtl_summary <- NULL

  if (!is.null(eqtl_dir) && nzchar(eqtl_dir)) {
    eqtl_summary <- file.path(dest_dir, "GTEx_v8_eQTL_summary.csv")
    easycoloc_build_gtex_summary(dest_file, eqtl_dir, qtl_type = "eqtl", output_file = eqtl_summary)
    generated_eqtl_summary <- eqtl_summary
    outputs[[length(outputs) + 1L]] <- data.table(mode = "build_gtex_summary", qtl_type = "eqtl", file = eqtl_summary, source_dir = eqtl_dir)
  }
  if (!is.null(sqtl_dir) && nzchar(sqtl_dir)) {
    sqtl_summary <- file.path(dest_dir, "GTEx_v8_sQTL_summary.csv")
    easycoloc_build_gtex_summary(dest_file, sqtl_dir, qtl_type = "sqtl", output_file = sqtl_summary)
    generated_sqtl_summary <- sqtl_summary
    outputs[[length(outputs) + 1L]] <- data.table(mode = "build_gtex_summary", qtl_type = "sqtl", file = sqtl_summary, source_dir = sqtl_dir)
  }

  if (is.null(config_output_file) || !nzchar(config_output_file)) {
    config_output_file <- file.path(dest_dir, "qtl_gtex_generated.yml")
  }
  summary_for_yaml <- if (identical(config_qtl_type, "eqtl")) generated_eqtl_summary else generated_sqtl_summary
  if (is.null(summary_for_yaml)) {
    summary_for_yaml <- if (!is.null(generated_eqtl_summary)) generated_eqtl_summary else generated_sqtl_summary
  }
  yaml_type <- if (!is.null(generated_eqtl_summary) && identical(summary_for_yaml, generated_eqtl_summary)) "eqtl" else "sqtl"
  if (!is.null(summary_for_yaml)) {
    easycoloc_write_gtex_qtl_yaml(summary_for_yaml, config_output_file, qtl_type = yaml_type)
    outputs[[length(outputs) + 1L]] <- data.table(mode = "write_qtl_yaml", qtl_type = yaml_type, file = config_output_file, source_summary = summary_for_yaml)
  }

  rbindlist(outputs, fill = TRUE)
}
