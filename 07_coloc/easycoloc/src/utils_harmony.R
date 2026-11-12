read_valid_harmony_cache <- function(gwas_cfg, qtl_build) {
  if (is.null(cfg_global$harmonize_dir) || is.null(gwas_cfg$id)) {
    return(NULL)
  }
  source_build <- if (is.null(gwas_cfg$build)) "19" else gsub("hg", "", gwas_cfg$build)
  target_build <- gsub("hg", "", qtl_build)
  harmony_settings <- if (!is.null(cfg_global$harmonization_settings)) cfg_global$harmonization_settings else list()
  compress_output <- if (!is.null(harmony_settings$compress_output)) isTRUE(harmony_settings$compress_output) else TRUE
  cache_candidates <- easycoloc_harmonized_cache_candidates(
    save_dir = cfg_global$harmonize_dir,
    dataset_id = gwas_cfg$id,
    source_build = source_build,
    target_build = target_build,
    compress_output = compress_output
  )
  preferred_cache_file <- cache_candidates[[1L]]
  existing_cache <- cache_candidates[file.exists(cache_candidates)]
  if (length(existing_cache) == 0L) {
    return(NULL)
  }
  input_mtime <- if (!is.null(gwas_cfg$file) && file.exists(gwas_cfg$file)) file.info(gwas_cfg$file)$mtime else NA
  ref_fasta <- if (gwas_cfg$build == "hg19") cfg_global$ref_genome_hg19 else cfg_global$ref_genome_hg38
  ref_mtime <- if (!is.null(ref_fasta) && file.exists(ref_fasta)) file.info(ref_fasta)$mtime else NA
  dep_mtimes <- c(input_mtime, ref_mtime)
  dep_mtimes <- dep_mtimes[!is.na(dep_mtimes)]
  dep_mtime <- if (length(dep_mtimes) > 0) max(dep_mtimes) else NA
  valid_cache_file <- NULL
  for (cache_file in existing_cache) {
    cache_info <- file.info(cache_file)
    cache_size_ok <- !is.na(cache_info$size) && cache_info$size > 0
    cache_new_enough <- is.na(dep_mtime) || (!is.na(cache_info$mtime) && cache_info$mtime >= dep_mtime)
    cache_header <- if (cache_size_ok) {
      tryCatch(names(data.table::fread(cache_file, nrows = 0, showProgress = FALSE)), error = function(e) character())
    } else {
      character()
    }
    if (cache_size_ok && cache_new_enough && identical(cache_header, easycoloc_harmonized_gwas_output_cols())) {
      valid_cache_file <- cache_file
      break
    }
  }
  if (is.null(valid_cache_file)) {
    return(NULL)
  }

  message(glue("Found cached harmonized file: {valid_cache_file}"))
  cached_dt <- data.table::fread(valid_cache_file, showProgress = FALSE)
  if (isTRUE(compress_output)) {
    valid_cache_file <- easycoloc_migrate_legacy_harmony_cache(
      cache_file = valid_cache_file,
      preferred_file = preferred_cache_file,
      cached_dt = cached_dt,
      verbose = TRUE
    )
  }
  easycoloc_standardize_harmonized_gwas(cached_dt, sample_size_n = gwas_cfg$sample_size_n)
}

prepare_gwas_harmony <- function(gwas_cfg, qtl_build, n_threads = 1L, prefer_cache = TRUE) {
  n_threads <- max(1L, as.integer(n_threads))
  if (isTRUE(prefer_cache)) {
    cached <- read_valid_harmony_cache(gwas_cfg, qtl_build = qtl_build)
    if (!is.null(cached)) {
      return(cached)
    }
  }

  if (!file.exists(gwas_cfg$file)) {
    stop(glue("File missing: {gwas_cfg$file}"))
  }

  old_dt_threads <- data.table::getDTthreads()
  data.table::setDTthreads(n_threads)
  on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)

  gwas_select_cols <- unique(unname(unlist(gwas_cfg$columns, use.names = FALSE)))
  gwas_select_cols <- gwas_select_cols[!is.na(gwas_select_cols) & nzchar(gwas_select_cols)]
  gwas_raw <- data.table::fread(gwas_cfg$file, select = gwas_select_cols, showProgress = FALSE)
  gwas_std <- format_sumstats(
    gwas_raw,
    type = "gwas",
    col_map = as.list(gwas_cfg$columns),
    case_control = (gwas_cfg$type == "cc"),
    sample_size_n = gwas_cfg$sample_size_n
  )

  ref_fasta <- if (gwas_cfg$build == "hg19") cfg_global$ref_genome_hg19 else cfg_global$ref_genome_hg38
  pop <- if (is.null(gwas_cfg$pop)) "EUR" else gwas_cfg$pop
  build_tag <- tolower(gwas_cfg$build)
  af_root <- cfg_global[["1kg_af"]]
  if (is.null(af_root) || !nzchar(af_root)) {
    stop("1kg_af not configured in global config")
  }
  ref_vcf_1kg_candidates <- c(
    file.path(af_root, build_tag, glue("{pop}_AF.vcf.gz")),
    file.path(af_root, build_tag, glue("{pop}_AF.tsv.gz")),
    file.path(af_root, glue("1KG_{build_tag}_{pop}_AF.vcf.gz")),
    file.path(af_root, glue("1KG_{build_tag}_{pop}_AF.tsv.gz"))
  )
  existing_ref_vcf <- ref_vcf_1kg_candidates[file.exists(ref_vcf_1kg_candidates)]
  ref_vcf_1kg <- if (length(existing_ref_vcf) > 0) existing_ref_vcf[[1]] else ref_vcf_1kg_candidates[[1]]
  ref_dbsnp <- if (gwas_cfg$build == "hg19") cfg_global[["dbsnp_hg19"]] else cfg_global[["dbsnp_hg38"]]
  ref_alt_freq <- "AF"
  harmony_settings <- if (!is.null(cfg_global$harmonization_settings)) cfg_global$harmonization_settings else list()
  chunked <- if (!is.null(harmony_settings$chunked)) isTRUE(harmony_settings$chunked) else NULL
  compress_output <- if (!is.null(harmony_settings$compress_output)) isTRUE(harmony_settings$compress_output) else TRUE
  chunk_min_rows <- if (!is.null(harmony_settings$chunk_min_rows)) {
    as.integer(harmony_settings$chunk_min_rows)
  } else {
    1000000L
  }
  chunk_parallel_jobs <- if (!is.null(harmony_settings$chunk_parallel_jobs)) {
    as.integer(harmony_settings$chunk_parallel_jobs)
  } else {
    1L
  }

  gwas_harm <- run_easycoloc_harmonization(
    gwas_std,
    ref_fasta = ref_fasta,
    ref_vcf = ref_vcf_1kg,
    ref_dbsnp = ref_dbsnp,
    ref_alt_freq = ref_alt_freq,
    source_build = if (is.null(gwas_cfg$build)) "19" else gsub("hg", "", gwas_cfg$build),
    target_build = gsub("hg", "", qtl_build),
    n_threads = n_threads,
    save_dir = cfg_global$harmonize_dir,
    dataset_id = gwas_cfg$id,
    input_file = gwas_cfg$file,
    liftover_chain = harmony_settings$liftover_chain,
    sample_size_n = gwas_cfg$sample_size_n,
    chunked = chunked,
    chunk_min_rows = chunk_min_rows,
    chunk_parallel_jobs = chunk_parallel_jobs,
    compress_output = compress_output,
    reuse_cache = !isFALSE(prefer_cache)
  )

  easycoloc_standardize_harmonized_gwas(gwas_harm, sample_size_n = gwas_cfg$sample_size_n)
}

pre_harmonize_gwas_caches <- function(datasets, qtl_build, n_cores) {
  runtime_cfg <- if (!is.null(cfg_global$runtime)) cfg_global$runtime else list()
  enabled <- if (!is.null(runtime_cfg$pre_harmonize_gwas)) {
    isTRUE(runtime_cfg$pre_harmonize_gwas)
  } else {
    TRUE
  }
  datasets_existing <- Filter(function(gwas_cfg) !is.null(gwas_cfg$file) && file.exists(gwas_cfg$file), datasets)
  if (!enabled || length(datasets_existing) <= 1L) {
    return(invisible(TRUE))
  }

  max_jobs <- if (!is.null(runtime_cfg$harmony_parallel_jobs)) {
    as.integer(runtime_cfg$harmony_parallel_jobs)
  } else {
    min(2L, as.integer(n_cores))
  }
  max_jobs <- max(1L, min(max_jobs, length(datasets_existing)))
  if (max_jobs <= 1L) {
    return(invisible(TRUE))
  }

  per_job_threads <- max(1L, floor(as.integer(n_cores) / max_jobs))
  message(glue("[HARMONY] Pre-harmonizing GWAS caches with {max_jobs} parallel jobs ({per_job_threads} thread(s) per job)"))
  append_runtime_event(
    stage = "harmony_prefetch_start",
    message_text = "Pre-harmonizing GWAS caches",
    extras = list(total_gwas = length(datasets_existing), parallel_jobs = max_jobs, per_job_threads = per_job_threads)
  )
  write_runtime_heartbeat(
    stage = "harmony_prefetch_start",
    message_text = "Pre-harmonizing GWAS caches",
    counters = list(total_gwas = length(datasets_existing), parallel_jobs = max_jobs)
  )

  worker <- function(gwas_cfg) {
    tryCatch(
      {
        message(glue("[HARMONY] Start {gwas_cfg$id}"))
        Sys.setenv(
          OMP_NUM_THREADS = as.character(per_job_threads),
          MKL_NUM_THREADS = as.character(per_job_threads),
          OPENBLAS_NUM_THREADS = as.character(per_job_threads)
        )
        data.table::setDTthreads(per_job_threads)
        invisible(prepare_gwas_harmony(gwas_cfg, qtl_build = qtl_build, n_threads = per_job_threads))
        list(id = gwas_cfg$id, ok = TRUE, error = NA_character_)
      },
      error = function(e) {
        list(id = gwas_cfg$id, ok = FALSE, error = conditionMessage(e))
      }
    )
  }

  results <- if (.Platform$OS.type == "unix") {
    parallel::mclapply(datasets_existing, worker, mc.cores = max_jobs, mc.set.seed = TRUE)
  } else {
    lapply(datasets_existing, worker)
  }
  failed <- Filter(function(x) !isTRUE(x$ok), results)
  if (length(failed) > 0) {
    for (item in failed) {
      warning(glue("[HARMONY] {item$id} failed during pre-harmonization: {item$error}"))
      append_runtime_event(level = "ERROR", stage = "harmony_prefetch_failed", message_text = item$error, gwas_id = item$id)
    }
    stop(glue("[HARMONY] Pre-harmonization failed for {length(failed)} GWAS dataset(s)"))
  }
  append_runtime_event(stage = "harmony_prefetch_complete", message_text = "Pre-harmonized GWAS caches")
  invisible(TRUE)
}
