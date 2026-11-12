#!/usr/bin/env Rscript
# =============================================================================
# rerun_plots.R - Regenerate plots from existing RDS files
# =============================================================================
# This script regenerates plots without running the full colocalization pipeline.
# It reads the RDS files that contain the saved data from previous runs.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glue)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
})

# Source utility functions
utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

# Read global config
cfg_global <- yaml::read_yaml("config/global.yml")

# Get plot settings from config
plot_sig_threshold <- if(!is.null(cfg_global$plot_settings$significance_threshold)) {
    as.numeric(cfg_global$plot_settings$significance_threshold)
} else {
    5.0e-8
}
plot_window_bp <- if(!is.null(cfg_global$plot_settings$plot_window_bp)) {
    as.integer(cfg_global$plot_settings$plot_window_bp)
} else {
    200000
}
plot_width <- if(!is.null(cfg_global$plot_settings$plot_width)) {
    as.numeric(cfg_global$plot_settings$plot_width)
} else {
    10
}
plot_height <- if(!is.null(cfg_global$plot_settings$plot_height)) {
    as.numeric(cfg_global$plot_settings$plot_height)
} else {
    8
}
title_phenotype_field <- if(!is.null(cfg_global$plot_settings$title_phenotype_field)) {
    cfg_global$plot_settings$title_phenotype_field
} else {
    "gene"
}
lead_snp_color <- if(!is.null(cfg_global$plot_settings$lead_snp_color)) {
    cfg_global$plot_settings$lead_snp_color
} else {
    "#7F3C8D"
}

# Directories
base_out_dir <- normalizePath(cfg_global$output_dir, mustWork = TRUE)
dir_plots <- file.path(base_out_dir, "plots")
dir_rds <- file.path(base_out_dir, "rds")

# List all RDS files
rds_files <- list.files(dir_rds, pattern = "\\.rds$", full.names = TRUE)

message(glue("Found {length(rds_files)} RDS files"))
message(glue("Regenerating plots with fix for is_credible column..."))

# Process each RDS file
success_count <- 0
error_count <- 0

for (rds_file in rds_files) {
    tryCatch({
        # Load data from RDS
        rds_data <- readRDS(rds_file)

        # Extract needed data
        merged_data <- rds_data$merged_data
        lead_snp <- rds_data$lead_snp
        gene_symbol <- rds_data$gene_symbol
        gwas_id <- rds_data$gwas_id
        qtl_id <- rds_data$qtl_id
        chrom <- rds_data$chrom
        plink_bfile <- rds_data$plink_bfile
        credible_set <- rds_data$credible_set

        # Assign to global environment for plot function
        assign("lead_SNP", lead_snp, envir = .GlobalEnv)
        assign("geneSymbol", gene_symbol, envir = .GlobalEnv)
        assign("plink_bfile", plink_bfile, envir = .GlobalEnv)
        assign("trait", gwas_id, envir = .GlobalEnv)

        # Load recombination data
        chrom_str <- paste0("chr", chrom)
        min_pos <- min(merged_data$POS.qtl, na.rm = TRUE)
        max_pos <- max(merged_data$POS.qtl, na.rm = TRUE)
        recomb_data <- load_recomb_map(chrom_str, min_pos, max_pos, cfg_global$recom)

        # Generate plot
        p <- plot_qtl_association(
            qtl_all_chrom = "CHR.qtl",
            qtl_all_pvalue = "P.qtl",
            leadSNP_DF = merged_data,
            ld_df = NULL,
            gtf_path = cfg_global$gene_anno,
            region_recomb = recomb_data,
            recomb_path = cfg_global$recom,
            show_lead_line = FALSE,
            qtl_type = qtl_id,
            phenotype_info = gene_symbol,
            significance_threshold = plot_sig_threshold,
            significance_label = NULL,
            title_phenotype_field = title_phenotype_field,
            plot_width = plot_width,
            plot_height = plot_height,
            credible_set = credible_set,
            plot_window_bp = plot_window_bp
        )

        # Save plot
        fname <- basename(gsub("\\.rds$", "", rds_file))
        fname <- gsub("_coloc$", "", fname)  # Remove trailing _coloc if present (avoid duplication)
        pname_pdf <- paste0(fname, "_coloc.pdf")
        pname_png <- paste0(fname, "_coloc.png")
        pdf_path <- file.path(dir_plots, pname_pdf)
        png_path <- file.path(dir_plots, pname_png)

        save_res <- save_plot_with_fallback(
            plot_obj = p,
            pdf_path = pdf_path,
            png_path = png_path,
            plot_width = plot_width,
            plot_height = plot_height
        )

        if (identical(save_res$format, "png")) {
          message(glue("[PLOT] Saved as PNG: {basename(save_res$path)}"))
        } else {
          message(glue("[PLOT] Saved as PDF: {basename(save_res$path)}"))
        }

        success_count <- success_count + 1

    }, error = function(e) {
        message(glue("[ERROR] Failed to process {basename(rds_file)}: {e$message}"))
        error_count <- error_count + 1
    })
}

message("==============================================================")
message(glue("Plot regeneration complete!"))
message(glue("Successful: {success_count}"))
message(glue("Errors: {error_count}"))
message("==============================================================")
