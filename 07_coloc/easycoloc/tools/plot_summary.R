#!/usr/bin/env Rscript
# ==============================================================================
# EasyColoc Summary Plot Tool
# 
# Generates compact multi-locus summary figures for presentations:
#   1. PP4 Summary Dot Plot — all loci ranked by PP4, colored by QTL type
#   2. Mini LocusZoom Facet Grid — locuscomparer-style clean per-locus plots
#
# Usage:
#   Rscript tools/plot_summary.R --top_n 6 --ncol 3
#   Rscript tools/plot_summary.R --genes "APOE,ARL3,SBNO1" --ncol 3
#   Rscript tools/plot_summary.R --summary_only --top_n 20
#   Rscript tools/plot_summary.R --help
#
# Requires: ggplot2, data.table, cowplot, ggrepel
# Pipeline must have been run with RDS output (results/rds/) for facet grid.
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(ggrepel)
})

# Source ld_extract() from pipeline utilities
# Detect script location via commandArgs (works with Rscript)
script_path <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        normalizePath(sub("^--file=", "", file_arg[1]))
    } else {
        normalizePath(".")
    }
}, error = function(e) ".")
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
utils_plot_path <- file.path(project_root, "src", "utils_plot.R")
utils_helpers_path <- file.path(project_root, "src", "utils_helpers.R")
if (file.exists(utils_plot_path)) {
    source(utils_plot_path, local = TRUE)
    message("[INIT] Loaded ld_extract() from ", utils_plot_path)
} else {
    message("[WARN] Could not find src/utils_plot.R at ", utils_plot_path)
    message("[WARN] LD computation unavailable — run from project root or set correct path")
    ld_extract <- NULL
}
if (file.exists(utils_helpers_path)) {
    source(utils_helpers_path, local = TRUE)
    message("[INIT] Loaded helpers from ", utils_helpers_path)
} else {
    stop("Could not find src/utils_helpers.R at ", utils_helpers_path, call. = FALSE)
}

# ==============================================================================
# CLI Argument Parsing (base R — no external dependency)
# ==============================================================================

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    opts <- list(
        results_dir = "results",
        output_dir  = "results/summary",
        top_n       = NULL,
        genes       = NULL,
        gwas_id     = NULL,
        ncol        = 3,
        width       = NULL,
        height      = NULL,
        dpi         = 300,
        summary_only = FALSE,
        facet_only   = FALSE,
        rds_only     = FALSE,
        pp4_threshold = 0.7,
        help        = FALSE
    )
    
    if (length(args) == 0) {
        opts$top_n <- 6
        return(opts)
    }
    
    i <- 1
    while (i <= length(args)) {
        arg <- args[i]
        if (arg == "--help" || arg == "-h") {
            opts$help <- TRUE
        } else if (arg == "--summary_only") {
            opts$summary_only <- TRUE
        } else if (arg == "--facet_only") {
            opts$facet_only <- TRUE
        } else if (arg == "--rds_only") {
            opts$rds_only <- TRUE
        } else if (arg %in% c("--results_dir", "--output_dir", "--top_n", "--genes",
                               "--gwas_id", "--ncol", "--width", "--height", "--dpi",
                               "--pp4_threshold")) {
            if (i + 1 > length(args)) stop(paste("Missing value for", arg))
            key <- sub("^--", "", arg)
            val <- args[i + 1]
            if (key %in% c("top_n", "ncol", "width", "height", "dpi")) {
                opts[[key]] <- as.numeric(val)
            } else if (key == "pp4_threshold") {
                opts[[key]] <- as.numeric(val)
            } else {
                opts[[key]] <- val
            }
            i <- i + 1
        } else {
            warning(paste("Unknown argument:", arg))
        }
        i <- i + 1
    }
    
    if (is.null(opts$top_n) && is.null(opts$genes)) {
        opts$top_n <- 6
    }
    
    return(opts)
}

print_usage <- function() {
    cat("
EasyColoc Summary Plot Tool
===========================

Usage:
  Rscript tools/plot_summary.R [OPTIONS]

Options:
  --results_dir DIR     Results directory (default: results)
  --output_dir DIR      Output directory for summary plots (default: results/summary)
  --top_n N             Select top N loci by PP4 (default: 6)
  --genes \"G1,G2,...\"   Select specific genes (comma-separated, matches gene symbols)
  --gwas_id ID          Filter by GWAS ID (e.g., 'EAS_SCZ')
  --ncol N              Number of columns in facet grid (default: 3)
  --width W             Output width in inches (auto-calculated if omitted)
  --height H            Output height in inches (auto-calculated if omitted)
  --dpi N               DPI for raster output (default: 300)
  --pp4_threshold T     PP4 threshold line on dot plot (default: 0.7)
  --summary_only        Only generate PP4 summary dot plot (no facet grid)
  --facet_only          Only generate facet grid (no summary dot plot)
  --rds_only            Scan RDS files directly (skip CSV, use during pipeline run)
  --help, -h            Show this help message

Examples:
  # Top 6 loci, 3 columns
  Rscript tools/plot_summary.R --top_n 6 --ncol 3

  # Specific genes
  Rscript tools/plot_summary.R --genes \"APOE,ARL3,SBNO1\"

  # Filter by GWAS, top 9 in 3x3 grid
  Rscript tools/plot_summary.R --gwas_id EAS_SCZ --top_n 9 --ncol 3

  # Only the PP4 summary dot plot, top 20
  Rscript tools/plot_summary.R --summary_only --top_n 20

  # Use RDS files directly (during pipeline run, no CSV needed)
  Rscript tools/plot_summary.R --rds_only --top_n 9 --ncol 3
")
}

# ==============================================================================
# Data Loading & Filtering
# ==============================================================================

load_results <- function(results_dir) {
    abf_dir <- file.path(results_dir, "abf")
    
    merged_file <- file.path(results_dir, "all_colocalization_results.csv")
    if (file.exists(merged_file)) {
        message("[INFO] Loading merged results: ", merged_file)
        dt <- fread(merged_file)
    } else {
        files <- Sys.glob(file.path(abf_dir, "*_results.csv"))
        if (length(files) == 0) stop("No results found in ", abf_dir)
        message("[INFO] Loading ", length(files), " result files from ", abf_dir)
        dt <- rbindlist(lapply(files, fread), fill = TRUE)
    }
    
    dt <- dt[!is.na(PP4)]
    dt[, gene_symbol := sub("^[^|]*\\|([^|]+)\\|.*$", "\\1", Gene)]
    dt <- dt[order(-PP4)]
    
    message("[INFO] Loaded ", nrow(dt), " results (", 
            sum(dt$PP4 >= 0.7), " with PP4 >= 0.7)")
    return(dt)
}

load_results_from_rds <- function(rds_dir) {
    rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(rds_files) == 0) stop("No RDS files found in ", rds_dir)
    message("[RDS_ONLY] Scanning ", length(rds_files), " RDS files...")
    
    # SuSiE CSV directory (sibling of rds/)
    susie_dir <- file.path(dirname(rds_dir), "susie")
    has_susie_dir <- dir.exists(susie_dir)
    if (has_susie_dir) {
        message("[RDS_ONLY] Found susie/ directory — will match SuSiE results")
    }
    
    rows <- lapply(rds_files, function(f) {
        d <- tryCatch(readRDS(f), error = function(e) NULL)
        if (is.null(d)) return(NULL)
        
        # Try to get SuSiE PP4:
        # 1. From RDS directly (new pipeline runs include it)
        susie_pp4 <- NA_real_
        if (!is.null(d$susie_best_pp4) && !is.na(d$susie_best_pp4)) {
            susie_pp4 <- d$susie_best_pp4
        } else if (has_susie_dir) {
            # 2. Fallback: match susie/ CSV by filename pattern
            rds_base <- gsub("_coloc\\.rds$", "", basename(f))
            susie_csv <- file.path(susie_dir, paste0(rds_base, "_susie.csv"))
            if (file.exists(susie_csv)) {
                susie_dt <- tryCatch(fread(susie_csv), error = function(e) NULL)
                if (!is.null(susie_dt) && "PP.H4.abf" %in% names(susie_dt) && nrow(susie_dt) > 0) {
                    susie_pp4 <- max(susie_dt$PP.H4.abf, na.rm = TRUE)
                }
            }
        }
        
        data.table(
            GWAS_ID    = d$gwas_id,
            QTL_ID     = d$qtl_id,
            Locus      = d$lead_snp,
            Gene       = d$gene_symbol,
            PP4        = d$pp4,
            susie_pp4  = susie_pp4,
            n_snps     = nrow(d$merged_data),
            rds_path   = f
        )
    })
    
    dt <- rbindlist(rows[!sapply(rows, is.null)], fill = TRUE)
    dt <- dt[!is.na(PP4)]
    dt[, gene_symbol := sub("^[^|]*\\|([^|]+)\\|.*$", "\\1", Gene)]
    dt <- dt[order(-PP4)]
    
    n_susie <- sum(!is.na(dt$susie_pp4))
    message("[RDS_ONLY] Found ", nrow(dt), " loci (", 
            sum(dt$PP4 >= 0.7), " with PP4 >= 0.7, ",
            n_susie, " with SuSiE results)")
    return(dt)
}


filter_results <- function(dt, opts) {
    filtered <- copy(dt)
    
    if (!is.null(opts$gwas_id)) {
        filtered <- filtered[GWAS_ID == opts$gwas_id]
        if (nrow(filtered) == 0) {
            available <- unique(dt$GWAS_ID)
            stop("No results for GWAS_ID '", opts$gwas_id, 
                 "'. Available: ", paste(available, collapse = ", "))
        }
        message("[FILTER] GWAS_ID = '", opts$gwas_id, "': ", nrow(filtered), " rows")
    }
    
    if (!is.null(opts$genes)) {
        gene_list <- trimws(unlist(strsplit(opts$genes, ",")))
        filtered <- filtered[gene_symbol %in% gene_list]
        if (nrow(filtered) == 0) {
            candidates <- unique(dt$gene_symbol)
            close_matches <- candidates[sapply(gene_list, function(g) {
                grepl(g, candidates, ignore.case = TRUE)
            }, simplify = FALSE) |> unlist() |> which()]
            stop("No results for genes: ", paste(gene_list, collapse = ", "),
                 if (length(close_matches) > 0) paste0("\n  Did you mean: ", paste(head(close_matches, 5), collapse = ", ")) else "")
        }
        message("[FILTER] Genes: ", paste(gene_list, collapse = ", "), ": ", nrow(filtered), " rows")
    }
    
    # Deduplicate: keep best PP4 per gene_symbol (across all GWAS/QTL)
    filtered <- filtered[, .SD[which.max(PP4)], by = gene_symbol]
    filtered <- filtered[order(-PP4)]
    
    if (!is.null(opts$top_n)) {
        n <- min(opts$top_n, nrow(filtered))
        filtered <- filtered[1:n]
        message("[FILTER] Top ", n, " loci selected")
    }
    
    return(filtered)
}

# ==============================================================================
# Filename Mapping: Results Row → RDS Path
# ==============================================================================

find_plot_rds <- function(row, rds_dir) {
    # Reconstruct filename as pipeline does: {GWAS_ID}_{QTL_ID}_{Gene}_coloc.rds
    fname <- sanitize_filename(paste0(row$GWAS_ID, "_", row$QTL_ID, "_", row$Gene, "_coloc.rds"))
    rds_path <- file.path(rds_dir, fname)
    
    if (file.exists(rds_path)) return(rds_path)
    
    # Fallback: search by gene symbol
    pattern <- paste0(".*", row$gene_symbol, ".*_coloc\\.rds$")
    candidates <- list.files(rds_dir, pattern = pattern, full.names = TRUE)
    
    if (length(candidates) > 0) {
        # Prefer exact GWAS_ID match
        gwas_match <- grep(row$GWAS_ID, candidates, value = TRUE)
        if (length(gwas_match) > 0) return(gwas_match[1])
        return(candidates[1])
    }
    
    return(NULL)
}

# ==============================================================================
# PP4 Summary Dot Plot
# ==============================================================================

plot_pp4_summary <- function(dt, pp4_threshold = 0.7) {
    plot_dt <- copy(dt)
    
    # Create display label: Gene (GWAS)
    plot_dt[, display_label := paste0(gene_symbol, " (", GWAS_ID, ")")]
    
    # Order by PP4
    plot_dt[, display_label := factor(display_label, 
                                       levels = rev(plot_dt$display_label))]
    
    # Color palette for QTL types
    qtl_types <- unique(plot_dt$QTL_ID)
    if (length(qtl_types) <= 3) {
        qtl_colors <- c("#E67E22", "#3498DB", "#2ECC71")[seq_along(qtl_types)]
    } else {
        qtl_colors <- scales::hue_pal()(length(qtl_types))
    }
    names(qtl_colors) <- qtl_types
    
    p <- ggplot(plot_dt, aes(x = PP4, y = display_label, color = QTL_ID)) +
        geom_vline(xintercept = pp4_threshold, linetype = "dashed", 
                   color = "grey40", linewidth = 0.5) +
        geom_segment(aes(x = 0, xend = PP4, yend = display_label),
                     linewidth = 0.4, color = "grey70") +
        geom_point(size = 3.5, alpha = 0.9) +
        scale_color_manual(values = qtl_colors, name = "QTL Type") +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                           expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = "Colocalization Summary",
            subtitle = paste0("PP.H4 (posterior probability of shared causal variant) | ",
                              "dashed line = ", pp4_threshold),
            x = "PP.H4",
            y = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 9),
            legend.position = "bottom",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            plot.margin = margin(10, 15, 10, 10)
        )
    
    return(p)
}

# ==============================================================================
# ABF vs SuSiE Comparison Dot Plot
#
# For loci that have both ABF and SuSiE PP4 values, creates a paired dot plot:
#   - Each locus shown as a row
#   - ABF PP4 (circle) and SuSiE PP4 (triangle) connected by a line
#   - Color indicates agreement (both high, both low, or discordant)
# ==============================================================================

plot_abf_vs_susie <- function(dt, pp4_threshold = 0.7) {
    # Filter to loci that have SuSiE results
    plot_dt <- copy(dt[!is.na(susie_pp4)])
    
    if (nrow(plot_dt) == 0) {
        message("[INFO] No loci with SuSiE results — skipping comparison plot")
        return(NULL)
    }
    
    # Create display label
    plot_dt[, display_label := paste0(gene_symbol, " (", GWAS_ID, ")")]
    
    # Classify agreement
    plot_dt[, agreement := fifelse(
        PP4 >= pp4_threshold & susie_pp4 >= pp4_threshold, "Both confirm",
        fifelse(PP4 < pp4_threshold & susie_pp4 < pp4_threshold, "Both reject",
                "Discordant")
    )]
    
    # Order by ABF PP4
    plot_dt[, display_label := factor(display_label, 
                                       levels = rev(plot_dt[order(PP4)]$display_label))]
    
    # Melt to long format for plotting both points
    long_dt <- melt(plot_dt, id.vars = c("display_label", "agreement"),
                    measure.vars = c("PP4", "susie_pp4"),
                    variable.name = "method", value.name = "pp4_value")
    long_dt[, method := fifelse(method == "PP4", "ABF", "SuSiE")]
    
    # Agreement colors
    agree_colors <- c("Both confirm" = "#2ECC71", "Both reject" = "#E74C3C", 
                       "Discordant" = "#F39C12")
    
    p <- ggplot() +
        # Threshold line
        geom_vline(xintercept = pp4_threshold, linetype = "dashed", 
                   color = "grey40", linewidth = 0.5) +
        # Connecting segment between ABF and SuSiE
        geom_segment(
            data = plot_dt,
            aes(x = PP4, xend = susie_pp4, 
                y = display_label, yend = display_label,
                color = agreement),
            linewidth = 0.6, alpha = 0.6
        ) +
        # ABF points (circle)
        geom_point(
            data = long_dt[method == "ABF"],
            aes(x = pp4_value, y = display_label, color = agreement),
            shape = 16, size = 3, alpha = 0.9
        ) +
        # SuSiE points (triangle)
        geom_point(
            data = long_dt[method == "SuSiE"],
            aes(x = pp4_value, y = display_label, color = agreement),
            shape = 17, size = 3, alpha = 0.9
        ) +
        scale_color_manual(values = agree_colors, name = "Agreement") +
        scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                           expand = expansion(mult = c(0, 0))) +
        # Manual legend for method shapes
        annotate("point", x = 0.03, y = 1.2, shape = 16, size = 2.5, color = "grey30") +
        annotate("text", x = 0.07, y = 1.2, label = "ABF", size = 2.8, hjust = 0, color = "grey30") +
        annotate("point", x = 0.15, y = 1.2, shape = 17, size = 2.5, color = "grey30") +
        annotate("text", x = 0.19, y = 1.2, label = "SuSiE", size = 2.8, hjust = 0, color = "grey30") +
        labs(
            title = "ABF vs SuSiE Colocalization",
            subtitle = paste0("Circle = ABF PP.H4, Triangle = SuSiE best PP.H4 | ",
                              "dashed line = ", pp4_threshold),
            x = "PP.H4",
            y = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = 9),
            legend.position = "bottom",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            plot.margin = margin(10, 15, 10, 10)
        )
    
    return(p)
}

# ==============================================================================
# Mini LocusZoom Plot (locuscomparer style)
#
# Creates a clean, compact LocusZoom-style plot from saved RDS data.
# Style reference: locuscomparer::make_locuszoom()
#   - theme_classic(), no gridlines, no recomb ribbon, no gene track
#   - Lead SNP: purple diamond (shape=23), labeled via ggrepel
#   - LD-colored points: r² bins mapped to blue→red palette
#   - X-axis in Mb, Y-axis: -log10(P)
# ==============================================================================

mini_locuszoom <- function(rds_data, plink_bin = "plink", susie_pp4 = NA_real_) {
    
    merged <- rds_data$merged_data
    lead_snp <- rds_data$lead_snp
    gene_sym <- sub("^[^|]*\\|([^|]+)\\|.*$", "\\1", rds_data$gene_symbol)
    chrom <- rds_data$chrom
    bfile <- rds_data$plink_bfile
    pp4 <- rds_data$pp4
    credible_set <- rds_data$credible_set
    
    # --- Build plot data frame ---
    # Resolve column names flexibly
    col_snp <- intersect(c("rsid", "snp", "SNP", "SNPID.gwas"), names(merged))[1]
    col_pos <- intersect(c("POS.qtl", "POS.gwas", "POS", "pos", "BP"), names(merged))[1]
    col_p   <- intersect(c("P.qtl", "P.gwas", "P", "pval"), names(merged))[1]
    
    if (is.na(col_snp) || is.na(col_pos) || is.na(col_p)) {
        message("[WARN] Missing columns in RDS for ", gene_sym)
        return(ggplot() + theme_void() + 
                   annotate("text", x = 0.5, y = 0.5, label = gene_sym, size = 3))
    }
    
    plot_df <- data.frame(
        rsid = as.character(merged[[col_snp]]),
        position = as.numeric(merged[[col_pos]]),
        p_value = as.numeric(merged[[col_p]]),
        stringsAsFactors = FALSE
    )
    
    # Clean
    plot_df <- plot_df[!is.na(plot_df$rsid) & plot_df$rsid != "" &
                           !is.na(plot_df$position) &
                           !is.na(plot_df$p_value) & plot_df$p_value > 0, ]
    
    if (nrow(plot_df) == 0) {
        return(ggplot() + theme_void() + 
                   annotate("text", x = 0.5, y = 0.5, label = gene_sym, size = 3))
    }
    
    plot_df$logp <- -log10(plot_df$p_value)
    
    # --- Determine lead SNP ---
    if (!lead_snp %in% plot_df$rsid) {
        lead_snp <- plot_df$rsid[which.min(plot_df$p_value)]
    }
    
    # --- Compute LD ---
    plot_df$r2 <- NA_real_
    
    if (!is.null(bfile) && file.exists(paste0(bfile, ".bed")) && !is.null(ld_extract)) {
        snps_for_ld <- unique(c(lead_snp, plot_df$rsid))
        ld_result <- tryCatch(
            ld_extract(snps_for_ld, bfile, plink_bin),
            error = function(e) NULL
        )
        
        if (!is.null(ld_result)) {
            # LD to lead SNP
            ld_to_lead <- ld_result[ld_result$SNP_A == lead_snp, ]
            
            if (nrow(ld_to_lead) == 0) {
                # Try as SNP_B
                ld_to_lead <- ld_result[ld_result$SNP_B == lead_snp, ]
                if (nrow(ld_to_lead) > 0) {
                    ld_to_lead <- data.frame(
                        SNP_A = ld_to_lead$SNP_B,
                        SNP_B = ld_to_lead$SNP_A,
                        R = ld_to_lead$R,
                        stringsAsFactors = FALSE
                    )
                }
            }
            
            if (nrow(ld_to_lead) > 0) {
                ld_to_lead$r2_calc <- abs(ld_to_lead$R)^2
                plot_df <- merge(
                    plot_df,
                    ld_to_lead[, c("SNP_B", "r2_calc")],
                    by.x = "rsid", by.y = "SNP_B",
                    all.x = TRUE
                )
                plot_df$r2 <- ifelse(!is.na(plot_df$r2_calc), plot_df$r2_calc, plot_df$r2)
                plot_df$r2_calc <- NULL
            }
        }
    }
    
    # Lead SNP r² = 1.0, missing = 0.05
    plot_df$r2[plot_df$rsid == lead_snp] <- 1.0
    plot_df$r2[is.na(plot_df$r2)] <- 0.05
    
    # --- Assign colors (locuscomparer palette) ---
    # r² bins: [0, 0.2) blue4, [0.2, 0.4) skyblue, [0.4, 0.6) darkgreen,
    #          [0.6, 0.8) orange, [0.8, 1.0] red
    ld_colors <- c("blue4" = "blue4", "skyblue" = "skyblue", 
                    "darkgreen" = "darkgreen", "orange" = "orange", "red" = "red")
    
    plot_df$ld_bin <- cut(
        plot_df$r2,
        breaks = c(-0.01, 0.2, 0.4, 0.6, 0.8, 1.01),
        labels = c("blue4", "skyblue", "darkgreen", "orange", "red"),
        include.lowest = TRUE
    )
    plot_df$ld_bin <- as.character(plot_df$ld_bin)
    
    # --- Shapes & sizes ---
    plot_df$is_lead <- plot_df$rsid == lead_snp
    plot_df$point_shape <- ifelse(plot_df$is_lead, 23, 21)  # diamond vs circle
    plot_df$point_size  <- ifelse(plot_df$is_lead, 3, 1.5)
    
    # Lead SNP color override: purple
    plot_df$fill_color <- plot_df$ld_bin
    plot_df$fill_color[plot_df$is_lead] <- "purple"
    
    # --- Label: only lead SNP ---
    plot_df$label <- ""
    plot_df$label[plot_df$is_lead] <- lead_snp
    
    # --- Sort: lead SNP on top ---
    plot_df <- plot_df[order(plot_df$is_lead), ]
    
    # --- Build title (include SuSiE PP4 if available) ---
    if (!is.na(susie_pp4)) {
        title_text <- paste0(gene_sym, " (ABF=", sprintf("%.2f", pp4), 
                             ", SuSiE=", sprintf("%.2f", susie_pp4), ")")
    } else {
        title_text <- paste0(gene_sym, " (PP4=", sprintf("%.2f", pp4), ")")
    }
    
    # --- Plot (locuscomparer style) ---
    p <- ggplot(plot_df, aes(x = position, y = logp)) +
        # Background points
        geom_point(
            data = plot_df[!plot_df$is_lead, ],
            aes(fill = fill_color),
            shape = 21, size = 1.5, alpha = 0.8, stroke = 0,
            color = "grey30"
        ) +
        # Lead SNP diamond
        geom_point(
            data = plot_df[plot_df$is_lead, ],
            aes(fill = fill_color),
            shape = 23, size = 3, alpha = 1, stroke = 0.3,
            color = "black"
        ) +
        # Lead SNP label
        geom_text_repel(
            data = plot_df[plot_df$label != "", ],
            aes(label = label),
            size = 2.2, fontface = "italic",
            max.overlaps = 20,
            segment.size = 0.3, segment.linetype = "dotted",
            segment.color = "grey40",
            min.segment.length = 0,
            box.padding = 0.3, point.padding = 0.2,
            nudge_y = 0.5
        ) +
        scale_fill_identity() +
        scale_x_continuous(
            labels = function(x) sprintf("%.1f", x / 1e6),
            expand = expansion(mult = 0.02)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
        labs(
            title = title_text,
            x = paste0("chr", chrom, " (Mb)"),
            y = bquote(-log[10] * "(P)")
        ) +
        theme_classic(base_size = 9) +
        theme(
            plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6),
            plot.margin = unit(c(0.3, 0.5, 0.3, 0.3), "lines")
        )
    
    return(p)
}

# ==============================================================================
# Mini LocusZoom Facet Grid
# ==============================================================================

build_facet_grid <- function(selected_dt, rds_dir, ncol_grid = 3, plink_bin = "plink") {
    if (nrow(selected_dt) == 0) {
        message("[WARN] No loci to plot in facet grid")
        return(NULL)
    }
    
    panel_list <- list()
    
    for (i in seq_len(nrow(selected_dt))) {
        row <- selected_dt[i]
        rds_path <- find_plot_rds(row, rds_dir)
        
        if (is.null(rds_path)) {
            message("[WARN] RDS not found for ", row$gene_symbol, 
                    " (", row$GWAS_ID, "/", row$QTL_ID, ")")
            # Create placeholder panel
            panel_list[[length(panel_list) + 1]] <- ggplot() + 
                theme_void() +
                annotate("text", x = 0.5, y = 0.5, 
                         label = paste0(row$gene_symbol, "\n(no data)"),
                         size = 3, color = "grey60")
            next
        }
        
        message("[PLOT] ", i, "/", nrow(selected_dt), " — ", row$gene_symbol, 
                " (PP4=", sprintf("%.2f", row$PP4), ")")
        
        rds_data <- tryCatch(readRDS(rds_path), error = function(e) NULL)
        if (is.null(rds_data)) {
            message("[WARN] Failed to read RDS: ", basename(rds_path))
            panel_list[[length(panel_list) + 1]] <- ggplot() + theme_void()
            next
        }
        
        panel <- tryCatch(
            mini_locuszoom(rds_data, plink_bin = plink_bin,
                           susie_pp4 = if ("susie_pp4" %in% names(row)) row$susie_pp4 else NA_real_),
            error = function(e) {
                message("[WARN] mini_locuszoom failed for ", row$gene_symbol, ": ", e$message)
                ggplot() + theme_void() + 
                    annotate("text", x = 0.5, y = 0.5, 
                             label = paste0(row$gene_symbol, "\n(plot error)"),
                             size = 3, color = "grey60")
            }
        )
        
        panel_list[[length(panel_list) + 1]] <- panel
    }
    
    if (length(panel_list) == 0) {
        message("[WARN] No plots generated — facet grid skipped")
        return(NULL)
    }
    
    n_panels <- length(panel_list)
    nrow_grid <- ceiling(n_panels / ncol_grid)
    
    # Pad with empty panels if needed
    while (length(panel_list) < nrow_grid * ncol_grid) {
        panel_list[[length(panel_list) + 1]] <- ggplot() + theme_void()
    }
    
    grid_plot <- plot_grid(
        plotlist = panel_list,
        ncol = ncol_grid,
        nrow = nrow_grid,
        align = "hv"
    )
    
    return(grid_plot)
}

# ==============================================================================
# Main
# ==============================================================================

main <- function() {
    opts <- parse_args()
    
    if (opts$help) {
        print_usage()
        return(invisible(NULL))
    }
    
    # Print configuration
    message("==============================================================")
    message("EasyColoc Summary Plot Tool")
    message("==============================================================")
    message("[CONFIG] Results dir:  ", opts$results_dir)
    message("[CONFIG] Output dir:   ", opts$output_dir)
    if (!is.null(opts$top_n))  message("[CONFIG] Top N:        ", opts$top_n)
    if (!is.null(opts$genes))  message("[CONFIG] Genes:        ", opts$genes)
    if (!is.null(opts$gwas_id)) message("[CONFIG] GWAS ID:     ", opts$gwas_id)
    message("[CONFIG] Columns:      ", opts$ncol)
    message("[CONFIG] PP4 threshold:", opts$pp4_threshold)
    message("[CONFIG] Summary only: ", opts$summary_only)
    message("[CONFIG] Facet only:   ", opts$facet_only)
    message("[CONFIG] RDS only:     ", opts$rds_only)
    message("")
    
    # Load and filter results
    if (opts$rds_only) {
        rds_dir <- file.path(opts$results_dir, "rds")
        if (!dir.exists(rds_dir)) stop("RDS directory not found: ", rds_dir)
        dt <- load_results_from_rds(rds_dir)
    } else {
        dt <- load_results(opts$results_dir)
    }
    selected <- filter_results(dt, opts)
    if (nrow(selected) == 0) {
        stop("No loci matched your filter criteria")
    }
    
    message("\n[SELECTED] ", nrow(selected), " loci:")
    for (i in seq_len(nrow(selected))) {
        r <- selected[i]
        message(sprintf("  %d. %s (PP4=%.3f, %s/%s, Locus=%s)", 
                        i, r$gene_symbol, r$PP4, r$GWAS_ID, r$QTL_ID, r$Locus))
    }
    message("")
    
    # Create output directory
    dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 1. PP4 Summary Dot Plot
    if (!opts$facet_only) {
        message("[PLOT] Generating PP4 summary dot plot...")
        
        # For summary plot, use more loci if available
        summary_dt <- if (!is.null(opts$genes) || opts$rds_only) {
            selected
        } else {
            n_summary <- max(nrow(selected), 15)
            dt_for_summary <- filter_results(dt, list(
                gwas_id = opts$gwas_id, genes = opts$genes, 
                top_n = n_summary
            ))
            dt_for_summary
        }
        
        p_summary <- plot_pp4_summary(summary_dt, pp4_threshold = opts$pp4_threshold)
        
        # Calculate dimensions
        n_loci <- nrow(summary_dt)
        summary_h <- max(3, n_loci * 0.35 + 1.5)
        summary_w <- 7
        
        summary_pdf <- file.path(opts$output_dir, "pp4_summary.pdf")
        summary_png <- file.path(opts$output_dir, "pp4_summary.png")
        
        ggsave(summary_pdf, p_summary, width = summary_w, height = summary_h,
               device = cairo_pdf)
        ggsave(summary_png, p_summary, width = summary_w, height = summary_h,
               dpi = opts$dpi)
        
        message("[DONE] PP4 summary: ", summary_pdf)
        message("[DONE] PP4 summary: ", summary_png)
        
        # 1b. ABF vs SuSiE Comparison Plot (if SuSiE results available)
        if ("susie_pp4" %in% names(summary_dt) && any(!is.na(summary_dt$susie_pp4))) {
            message("[PLOT] Generating ABF vs SuSiE comparison plot...")
            p_compare <- plot_abf_vs_susie(summary_dt, pp4_threshold = opts$pp4_threshold)
            if (!is.null(p_compare)) {
                n_susie_loci <- sum(!is.na(summary_dt$susie_pp4))
                compare_h <- max(3, n_susie_loci * 0.45 + 1.5)
                compare_w <- 7.5
                
                compare_pdf <- file.path(opts$output_dir, "abf_vs_susie.pdf")
                compare_png <- file.path(opts$output_dir, "abf_vs_susie.png")
                
                ggsave(compare_pdf, p_compare, width = compare_w, height = compare_h,
                       device = cairo_pdf)
                ggsave(compare_png, p_compare, width = compare_w, height = compare_h,
                       dpi = opts$dpi)
                
                message("[DONE] ABF vs SuSiE: ", compare_pdf)
                message("[DONE] ABF vs SuSiE: ", compare_png)
            }
        }
    }
    
    # 2. Mini LocusZoom Facet Grid (from RDS data)
    if (!opts$summary_only) {
        message("\n[PLOT] Generating mini LocusZoom facet grid...")
        
        rds_dir <- file.path(opts$results_dir, "rds")
        if (!dir.exists(rds_dir)) {
            message("[WARN] RDS directory not found: ", rds_dir)
            message("[HINT] Re-run the pipeline to generate RDS files (run_coloc.R now saves them automatically)")
            message("[HINT] RDS files are required for the locuscomparer-style mini plots")
        } else {
            n_rds <- length(list.files(rds_dir, pattern = "\\.rds$"))
            message("[INFO] Found ", n_rds, " RDS files in ", rds_dir)
            
            grid_plot <- build_facet_grid(
                selected, rds_dir,
                ncol_grid = opts$ncol,
                plink_bin = "plink"
            )
            
            if (!is.null(grid_plot)) {
                n_panels <- nrow(selected)
                nrow_grid <- ceiling(n_panels / opts$ncol)
                
                # Auto-calculate dimensions: compact mini plots
                # Each panel ~3.0 inches wide, ~2.5 inches tall
                facet_w <- if (!is.null(opts$width)) opts$width else opts$ncol * 3.0
                facet_h <- if (!is.null(opts$height)) opts$height else nrow_grid * 2.5
                
                facet_pdf <- file.path(opts$output_dir, "facet_locuszoom.pdf")
                facet_png <- file.path(opts$output_dir, "facet_locuszoom.png")
                
                ggsave(facet_pdf, grid_plot, width = facet_w, height = facet_h,
                       device = cairo_pdf)
                ggsave(facet_png, grid_plot, width = facet_w, height = facet_h,
                       dpi = opts$dpi)
                
                message("[DONE] Facet grid: ", facet_pdf)
                message("[DONE] Facet grid: ", facet_png)
            }
        }
    }
    
    message("\n==============================================================")
    message("Summary plots saved to: ", opts$output_dir)
    message("==============================================================")
}

main()
