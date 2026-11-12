#!/usr/bin/env Rscript

# bed2rds.R - Convert dbSNP BED files directly to RDS hash tables
# Merged from dbsnp_hash_table_maker.py + hash2rds.R
# Author: EasyColoc

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

BED_DIR <- "."
OUT_DIR <- "."
BUILD <- "hg38"

usage <- function() {
  cat("Usage: Rscript bed2rds.R [OPTIONS]\n\n")
  cat("Convert dbSNP BED files directly to RDS hash tables.\n\n")
  cat("Options:\n")
  cat("  -i, --input DIR    Input directory containing BED files (default: current directory)\n")
  cat("  -o, --output DIR   Output directory for RDS files (default: same as input)\n")
  cat("  -b, --build BUILD  Genome build: hg38 or hg19 (default: hg38)\n")
  cat("  -h, --help         Show this help message\n")
  quit(status = 0)
}

args <- commandArgs(trailingOnly = TRUE)

if (any(c("-h", "--help") %in% args)) usage()

parse_arg <- function(flags) {
  for (flag in flags) {
    exact_idx <- which(args == flag)
    if (length(exact_idx) > 0) {
      idx <- exact_idx[1]
      if (idx < length(args) && !startsWith(args[idx + 1], "-")) {
        return(args[idx + 1])
      }
    }
    prefix <- paste0(flag, "=")
    value_idx <- grep(paste0("^", gsub("([\\-])", "\\\\\\1", prefix)), args)
    if (length(value_idx) > 0) {
      return(sub(paste0("^", gsub("([\\-])", "\\\\\\1", prefix)), "", args[value_idx[1]]))
    }
  }
  NULL
}

BED_DIR <- parse_arg(c("-i", "--input")) %||% BED_DIR
OUT_DIR <- parse_arg(c("-o", "--output")) %||% BED_DIR
BUILD <- parse_arg(c("-b", "--build")) %||% BUILD

cat("EasyColoc BED to RDS Converter\n")
cat("==============================\n")
cat(glue("Input directory: {BED_DIR}\n"))
cat(glue("Output directory: {OUT_DIR}\n"))
cat(glue("Genome build: {BUILD}\n\n"))

if (!dir.exists(BED_DIR)) {
  stop("Error: Input directory not found: ", BED_DIR)
}

dir.create(OUT_DIR, showWarnings = FALSE)

chr_list <- c(1:22, "X", "Y", "MT", "PAR", "Un", "Multi", "NotOn", "AltOnly")

cat("Processing chromosome files...\n\n")

for (chrom in chr_list) {
  bed_file <- file.path(BED_DIR, glue("bed_chr_{chrom}.bed.gz"))
  rds_file <- file.path(OUT_DIR, glue("chr_{chrom}_snp{BUILD}_hash_table.rds"))

  if (!file.exists(bed_file)) {
    cat(glue("[SKIP] chr{chrom}: BED file not found\n"))
    next
  }

  if (file.exists(rds_file)) {
    cat(glue("[SKIP] chr{chrom}: RDS already exists\n"))
    next
  }

  cat(glue("[PROCESS] chr{chrom}: "))

  t1 <- Sys.time()

  tryCatch({
    cmd <- glue("zcat {bed_file} | tail -n +2")
    con <- pipe(cmd, "r")
    dbsnp_data <- data.table(
      chrom = character(),
      chromStart = integer(),
      chromEnd = integer(),
      SNP = character(),
      score = character(),
      strand = character()
    )
    dbsnp_data <- fread(cmd, sep = "\t", header = FALSE, skip = 1,
                         colClasses = c("character", "integer", "integer",
                                       "character", "character", "character"))
    close(con)

    if (nrow(dbsnp_data) == 0) {
      cat("empty file\n")
      next
    }

    setnames(dbsnp_data, c("chrom", "chromStart", "chromEnd", "SNP", "score", "strand"))

    dbsnp_data[, chrom_chromEnd := paste0(chrom, "_", chromEnd)]

    rsToChrPosition <- dbsnp_data[, setNames(chrom_chromEnd, SNP)]

    rm(dbsnp_data)
    gc()

    saveRDS(rsToChrPosition, file = rds_file, compress = TRUE)

    t2 <- Sys.time()
    elapsed <- round(as.numeric(difftime(t2, t1, units = "secs")), 1)
    cat(glue("{formatC(length(rsToChrPosition), format='d')} SNPs saved in {elapsed}s\n"))

    rm(rsToChrPosition)
    gc()

  }, error = function(e) {
    cat(glue("ERROR: {e$message}\n"))
  })
}

cat("\nAll done!\n")
cat(glue("RDS files saved to: {OUT_DIR}\n"))
cat("You can now use these files with EasyColoc's hash table system.\n")
