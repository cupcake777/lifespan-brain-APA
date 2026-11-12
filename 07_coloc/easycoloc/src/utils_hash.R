# ------------------------------------------------------------------------------
# src/utils_hash.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(glue)
})
.hash_cache <- new.env(parent = emptyenv())
.hash_dir <- NULL
initialize_hash_system <- function(hash_dir) {
    if (is.null(hash_dir) || !dir.exists(hash_dir)) {
        warning("Hash table directory not found. SNP ID conversion disabled.")
        return(FALSE)
    }
    .hash_dir <<- hash_dir
    hash_files <- list.files(hash_dir, pattern = "chr_.*_snp.*_hash_table\\.rds$", full.names = FALSE)
    if (length(hash_files) == 0) {
        warning("No hash table RDS files found.")
        return(FALSE)
    }
    .hash_cache$file_map <- list()
    for (f in hash_files) {
        chrom_match <- regmatches(f, regexpr("chr_(\\d+|X|Y|MT)_", f))
        if (length(chrom_match) == 0) {
            next
        }
        chrom <- gsub("chr_(\\d+|X|Y|MT)_.*", "\\1", chrom_match)
        .hash_cache$file_map[[chrom]] <- file.path(hash_dir, f)
    }
    message(glue("[HASH] Initialized system with {length(.hash_cache$file_map)} chromosomes available"))
    # message(glue("       Chromosomes: {paste(names(.hash_cache$file_map), collapse=', ')}"))
    message("       Note: Data will be loaded on-demand to save memory")

    return(TRUE)
}
load_chromosome_hash <- function(chrom) {
    chrom <- as.character(chrom)
    if (exists(chrom, envir = .hash_cache, inherits = FALSE)) {
        return(get(chrom, envir = .hash_cache))
    }

    if (is.null(.hash_cache$file_map) || !chrom %in% names(.hash_cache$file_map)) {
        return(NULL)
    }

    hash_file <- .hash_cache$file_map[[chrom]]

    message(glue("  [HASH] Loading chr{chrom} from {basename(hash_file)}..."))
    hash_data <- tryCatch({
        readRDS(hash_file)
    }, error = function(e) {
        warning(glue("Failed to load hash table for chr{chrom}: {e$message}"))
        return(NULL)
    })

    if (!is.null(hash_data)) {
        message(glue("  [HASH] chr{chrom} loaded: {format(length(hash_data), big.mark=',')} entries"))
        assign(chrom, hash_data, envir = .hash_cache)
    }

    return(hash_data)
}

convert_rsid_to_pos <- function(rsids, chromosomes) {
    result <- character(length(rsids))
    result[] <- NA_character_ # Default toNA

    if (is.null(.hash_cache$file_map)) {
        return(result)
    }
    unique_chroms <- unique(chromosomes[!is.na(chromosomes)])

    for (chrom in unique_chroms) {
        chrom_idx <- which(chromosomes == chrom)

        if (length(chrom_idx) == 0) next
        hash_data <- load_chromosome_hash(chrom)

        if (is.null(hash_data)) {
            next
        }
        chrom_rsids <- as.character(rsids[chrom_idx])
        matches <- chrom_rsids %in% names(hash_data)
        if (any(matches)) {
            result[chrom_idx[matches]] <- hash_data[chrom_rsids[matches]]
            if (is.list(result)) result <- unlist(result)
        }
    }

    return(result)
}