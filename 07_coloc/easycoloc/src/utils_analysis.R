# ------------------------------------------------------------------------------
# src/utils_analysis.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(susieR)
})

get_coloc_results <- function(colocInputFile, gwas_type = "cc", gwas_prop = 0.5, gwas_N = NULL, qtl_N = NULL, use_susie = TRUE, susie_threshold = 0.75, plink_bfile = NULL, plink_bin = "plink", p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, maf_default = 0.1, maf_na_replacement = 0.05, maf_epsilon = 1e-6, keep_file = NULL, ld_res = NULL) {
    # ==========================================================================
    # get_clean_maf: Extract or impute minor allele frequencies
    # ==========================================================================
    # Priority: (1) EAF/MAF/af column if present
    #           (2) Default maf_default (configurable, default 0.1)
    # Missing values replaced with maf_na_replacement (default 0.05)
    # Values clamped to [maf_epsilon, 1-maf_epsilon] for numerical stability
    get_clean_maf <- function(df, suffix) {
        candidates <- c(paste0("EAF", suffix), paste0("MAF", suffix), paste0("af", suffix))
        if (suffix == ".gwas") candidates <- c(candidates, "EAF", "MAF", "af")
        raw_freq <- NULL
        for (col in candidates) if (col %in% names(df)) { raw_freq <- df[[col]]; break }

        # Case 1: No frequency column found - use default
        if (is.null(raw_freq)) {
            warning(sprintf("[MAF] No allele frequency column found; using default MAF=%.2f for all %d SNPs",
                           maf_default, nrow(df)))
            return(rep(maf_default, nrow(df)))
        }

        freq <- as.numeric(raw_freq); maf <- pmin(freq, 1 - freq)

        # Case 2: Missing values in frequency column
        n_na_maf <- sum(is.na(maf))
        if (n_na_maf > 0) {
            message(sprintf("[MAF] %d SNPs (%.1f%%) had missing MAF - replaced with %.3f",
                           n_na_maf, n_na_maf/nrow(df)*100, maf_na_replacement))
            maf[is.na(maf)] <- maf_na_replacement
        }

        # Clamp to valid range to prevent numerical issues
        maf[maf < maf_epsilon] <- maf_epsilon
        maf[maf > (1 - maf_epsilon)] <- 1 - maf_epsilon
        return(maf)
    }

    N_gwas_val <- if("N.gwas" %in% names(colocInputFile)) colocInputFile$N.gwas else if ("N" %in% names(colocInputFile)) colocInputFile$N else gwas_N
    if (is.null(N_gwas_val)) stop("GWAS N missing.")
    N_qtl_val <- if("N.qtl" %in% names(colocInputFile)) colocInputFile$N.qtl else qtl_N
    if (is.null(N_qtl_val)) stop("QTL N missing.")
    
    d1 <- list(snp = as.character(colocInputFile$snp), beta = as.numeric(colocInputFile$BETA.gwas), varbeta = as.numeric(colocInputFile$SE.gwas)^2, pvalues = as.numeric(colocInputFile$P.gwas), type = gwas_type, N = as.numeric(N_gwas_val), MAF = get_clean_maf(colocInputFile, ".gwas"))
    if (gwas_type == "cc") { if (is.null(gwas_prop)) stop("Case prop required."); d1$s <- as.numeric(gwas_prop) }
    d2 <- list(snp = as.character(colocInputFile$snp), beta = as.numeric(colocInputFile$BETA.qtl), varbeta = as.numeric(colocInputFile$SE.qtl)^2, pvalues = as.numeric(colocInputFile$P.qtl), type = "quant", N = as.numeric(N_qtl_val), MAF = get_clean_maf(colocInputFile, ".qtl"))
    
    message("Running Coloc ABF...")
    res_abf <- coloc.abf(d1, d2, p1 = p1, p2 = p2, p12 = p12)
    pp4 <- as.numeric(res_abf$summary["PP.H4.abf"])
    
    if (!use_susie || is.na(pp4) || pp4 < susie_threshold) {
        return(res_abf)
    }
    
    message(sprintf("  -> Strong signal detected (PP.H4 = %.4f). Running SuSiE...", pp4))
    
    # Check if SNPs are in rsID format (required for PLINK)
    snp_sample <- head(d1$snp, 5)
    is_rsid_format <- all(grepl("^rs", snp_sample), na.rm = TRUE)
    if (!is_rsid_format) {
        message("SNPs not in rsID format - skipping SuSiE (PLINK requires rsIDs)")
        return(res_abf)
    }
    
    if (is.null(ld_res) && (!exists("get_ld_matrix") || is.null(plink_bfile))) { 
        warning("LD tools missing. Skipping SuSiE."); 
        return(res_abf) 
    }

    if (is.null(ld_res)) {
        ld_res <- get_ld_matrix(d1$snp, plink_bfile, plink_bin, keep_file = keep_file)
    }
    
    if (is.null(ld_res) || is.null(ld_res$R)) { 
        warning(glue("LD calculation failed (NULL matrix). Check if SNP IDs (e.g. {d1$snp[1]}) match PLINK reference."))
        return(res_abf) 
    }
    
    common_snps <- intersect(d1$snp, colnames(ld_res$R))
    if (length(common_snps) < 5) { 
        warning(glue("Too few overlapping SNPs ({length(common_snps)}) between data and LD ref. Skipping SuSiE."))
        return(res_abf) 
    }
    
    filter_data <- function(d, snps) {
        idx <- match(snps, d$snp)
        d$snp <- d$snp[idx]; d$beta <- d$beta[idx]; d$varbeta <- d$varbeta[idx]
        if(!is.null(d$MAF)) d$MAF <- d$MAF[idx]
        return(d)
    }
    
    d1_susie <- filter_data(d1, common_snps); d2_susie <- filter_data(d2, common_snps)
    R_susie <- ld_res$R[common_snps, common_snps]

    # ==========================================================================
    # SuSiE-RSS Fine-Mapping: Set seeds for reproducibility
    # ==========================================================================
    # susie_rss() uses iterative variational inference which involves
    # stochastic initialization. We use offset seeds from the global seed
    # to ensure: (1) reproducibility across runs, (2) different seeds for
    # GWAS vs QTL to avoid artificial correlation.
    # Offsets: +1 for GWAS, +2 for QTL
    if (exists("global_seed", mode = "numeric")) {
        set.seed(global_seed + 1)
        message("[SuSiE] Seed set for GWAS fine-mapping (global_seed + 1)")
    }
    susie_1 <- try(susie_rss(bhat = d1_susie$beta, shat = sqrt(d1_susie$varbeta), R = R_susie, n = d1_susie$N[1], verbose = FALSE), silent = TRUE)

    if (exists("global_seed", mode = "numeric")) {
        set.seed(global_seed + 2)
        message("[SuSiE] Seed set for QTL fine-mapping (global_seed + 2)")
    }
    susie_2 <- try(susie_rss(bhat = d2_susie$beta, shat = sqrt(d2_susie$varbeta), R = R_susie, n = d2_susie$N[1], verbose = FALSE), silent = TRUE)
    
    if (inherits(susie_1, "try-error") || inherits(susie_2, "try-error")) { 
        warning("SuSiE failed to converge."); 
        return(res_abf) 
    }
    
    res_susie <- coloc.susie(susie_1, susie_2)
    message("  SuSiE Coloc complete.")

    # ==========================================================================
    # Extract Credible Set from SuSiE Results
    # ==========================================================================
    tryCatch({
        credible_set <- extract_credible_set(
            susie_obj = susie_1,
            susie_res = susie_2,
            credible_level = 0.95,
            snp_ids = common_snps
        )

        if (!is.null(credible_set) && nrow(credible_set) > 0) {
            res_abf$credible_set <- credible_set
            message(sprintf("  [Credible Set] %d SNPs in 95%% credible set",
                           sum(credible_set$is_credible)))
        }
    }, error = function(e) {
        warning(sprintf("[Credible Set] Extraction failed: %s", e$message))
    })

    res_abf$susie_result <- res_susie
    return(res_abf)
}

# =============================================================================
# extract_credible_set: Extract credible set SNPs from SuSiE results
# =============================================================================
# Extracts 95% credible set SNPs with posterior inclusion probabilities (PIP)
# and cluster assignments from SuSiE fine-mapping results.
#
# Arguments:
#   susie_obj: SuSiE model object (from susie_rss or susie)
#   susie_res: SuSiE RSS result with lbf_all
#   credible_level: Credible level (default 0.95 for 95% credible set)
#
# Returns:
#   data.frame with columns:
#     - snp: SNP identifier
#     - cluster: SuSiE effect cluster (1, 2, 3...)
#     - PIP: Posterior inclusion probability
#     - SNP_PPs: SNP-level posterior probabilities
#     - is_credible: Whether in 95% credible set
#     - lead_in_cluster: Whether this is the lead SNP in its cluster
# =============================================================================
extract_credible_set <- function(susie_obj, susie_res,
                                  credible_level = 0.95,
                                  snp_ids = NULL) {

    if (is.null(susie_obj) || is.null(susie_res)) {
        return(NULL)
    }

    # Extract PIPs (Posterior Inclusion Probabilities)
    pip_values <- susie_obj$pip

    # Extract cluster assignments
    cluster_assign <- susie_obj$sets$cs_index

    # Get number of effects (clusters)
    L <- if(!is.null(susie_obj$alpha)) {
        dim(susie_obj$alpha)[2]
    } else if(!is.null(susie_res$alpha)) {
        dim(susie_res$alpha)[2]
    } else {
        length(unique(cluster_assign[!is.na(cluster_assign)]))
    }

    # Calculate SNP-level posterior probabilities
    if (!is.null(susie_res$lbf_all)) {
        # Convert log Bayes factors to posterior probabilities
        lbf_all <- susie_res$lbf_all
        max_lbf <- max(lbf_all, na.rm = TRUE)
        pp_raw <- exp(lbf_all - max_lbf)  # Numerical stability
        pp_total <- sum(pp_raw, na.rm = TRUE)
        snp_pps <- pp_raw / pp_total
    } else if (!is.null(susie_obj$alpha)) {
        # Use alpha matrix directly
        alpha_mat <- susie_obj$alpha
        snp_pps <- rowSums(alpha_mat, na.rm = TRUE)
    } else {
        snp_pps <- pip_values
    }

    # Build credible set data frame
    n_snps <- length(pip_values)

    if (is.null(snp_ids)) {
        snp_ids <- paste0("SNP_", seq_len(n_snps))
    }

    credible_df <- data.frame(
        snp = snp_ids,
        cluster = cluster_assign,
        PIP = pip_values,
        SNP_PPs = snp_pps,
        is_credible = FALSE,
        lead_in_cluster = FALSE,
        stringsAsFactors = FALSE
    )

    # Identify credible set SNPs per cluster
    for (l in 1:L) {
        cluster_snps <- which(cluster_assign == l)

        if (length(cluster_snps) == 0) next

        # Sort by PIP descending
        sorted_idx <- order(credible_df$PIP[cluster_snps], decreasing = TRUE)
        cluster_snps_sorted <- cluster_snps[sorted_idx]

        # Accumulate PPs until reaching credible level
        cumsum_pp <- 0
        credible_snps <- c()

        for (idx in cluster_snps_sorted) {
            cumsum_pp <- cumsum_pp + credible_df$SNP_PPs[idx]
            credible_snps <- c(credible_snps, idx)

            if (cumsum_pp >= credible_level) break
        }

        # Mark credible set SNPs
        credible_df$is_credible[credible_snps] <- TRUE

        # Mark lead SNP (highest PIP in cluster)
        if (length(cluster_snps_sorted) > 0) {
            credible_df$lead_in_cluster[cluster_snps_sorted[1]] <- TRUE
        }
    }

    # Add summary metadata
    attr(credible_df, "n_clusters") <- L
    attr(credible_df, "n_credible_snps") <- sum(credible_df$is_credible)
    attr(credible_df, "credible_level") <- credible_level

    message(sprintf("[Credible Set] Extracted %d SNPs in %d clusters (%d in 95%% CS)",
                    n_snps, L, sum(credible_df$is_credible)))

    return(credible_df)
}

# =============================================================================
# annotate_credible_set: Annotate credible set SNPs with functional information
# =============================================================================
# Adds functional annotations to credible set SNPs including:
# - Gene overlap (from GTF)
# - Consequence type (missense, synonymous, etc.)
# - Regulome annotations (if available)
#
# Arguments:
#   credible_df: Data frame from extract_credible_set
#   gtf_path: Path to GTF annotation file
#   dbsnp_path: Path to dbSNP VCF for RSID lookup (optional)
#
# Returns:
#   Augmented data.frame with annotation columns
# =============================================================================
annotate_credible_set <- function(credible_df,
                                   gtf_path = NULL,
                                   snp_info = NULL) {

    if (is.null(credible_df) || nrow(credible_df) == 0) {
        return(NULL)
    }

    # Initialize annotation columns
    credible_df$gene_symbol <- NA_character_
    credible_df$consequence <- "intergenic"
    credible_df$distance_to_TSS <- NA_integer_

    # If SNP info available (from PLINK BIM), add position-based annotations
    if (!is.null(snp_info) && is.data.frame(snp_info)) {
        # Merge with SNP info
        if (all(c("SNP", "POS", "CHR") %in% names(snp_info))) {
            credible_df <- merge(
                credible_df,
                snp_info[, c("SNP", "POS", "CHR")],
                by.x = "snp",
                by.y = "SNP",
                all.x = TRUE
            )
        }
    }

    # Add GTF-based gene annotations if available
    if (!is.null(gtf_path) && file.exists(gtf_path)) {
        tryCatch({
            library(rtracklayer)
            library(GenomicFeatures)
            library(GenomicRanges)

            # Create GRanges for credible SNPs
            if ("POS" %in% names(credible_df) && "CHR" %in% names(credible_df)) {
                snp_gr <- GRanges(
                    seqnames = paste0("chr", credible_df$CHR),
                    ranges = IRanges(credible_df$POS, credible_df$POS)
                )

                # Import GTF
                gtf <- import(gtf_path)

                # Find overlapping genes
                overlaps <- findOverlaps(snp_gr, gtf[gtf$type == "gene"])

                if (length(overlaps) > 0) {
                    for (i in seq_along(overlaps)) {
                        snp_idx <- queryHits(overlaps)[i]
                        gene_idx <- subjectHits(overlaps)[i]

                        gene_sym <- gtf$gene_name[gene_idx]
                        if (is.na(gene_sym)) {
                            gene_sym <- gtf$gene_id[gene_idx]
                        }

                        credible_df$gene_symbol[snp_idx] <- as.character(gene_sym)

                        # Calculate distance to TSS
                        gene_range <- gtf[gene_idx]
                        tss_pos <- if(strand(gene_range) == "+") {
                            start(gene_range)
                        } else {
                            end(gene_range)
                        }

                        credible_df$distance_to_TSS[snp_idx] <- as.integer(
                            credible_df$POS[snp_idx] - tss_pos
                        )
                    }
                }

                # Classify consequence
                exon_overlaps <- findOverlaps(snp_gr, gtf[gtf$type == "exon"])
                if (length(exon_overlaps) > 0) {
                    credible_df$consequence[queryHits(exon_overlaps)] <- "exonic"
                }

                utr_overlaps <- findOverlaps(snp_gr, gtf[gtf$type == "exon" & grepl("UTR", gtf$gene_type)])
                if (length(utr_overlaps) > 0) {
                    credible_df$consequence[queryHits(utr_overlaps)] <- "UTR"
                }
            }

            message("[Annotation] Added functional annotations from GTF")
        }, error = function(e) {
            warning(sprintf("[Annotation] GTF annotation failed: %s", e$message))
        })
    }

    # Summary statistics
    n_credible <- sum(credible_df$is_credible, na.rm = TRUE)
    n_annotated <- sum(!is.na(credible_df$gene_symbol), na.rm = TRUE)

    message(sprintf("[Annotation] %d credible SNPs, %d annotated with genes",
                    n_credible, n_annotated))

    return(credible_df)
}

# =============================================================================
# summarize_susie_results: Generate summary table for SuSiE fine-mapping
# =============================================================================
# Creates a publication-ready summary table of SuSiE results across all loci
#
# Arguments:
#   all_results: List of coloc results with SuSiE
#   output_file: Path to save summary CSV (optional)
#
# Returns:
#   data.frame with summary statistics per locus
# =============================================================================
summarize_susie_results <- function(all_results, output_file = NULL) {

    summary_list <- lapply(names(all_results), function(locus_id) {
        res <- all_results[[locus_id]]

        if (is.null(res$susie_result)) {
            return(NULL)
        }

        susie_res <- res$susie_result
        abf_summary <- res$summary

        # Extract SuSiE summary statistics
        n_clusters <- if(!is.null(susie_res$summary)) {
            susie_res$summary["sets"]
        } else {
            NA_integer_
        }

        pp4_conditional <- abf_summary["PP.H4.abf"] /
                           (abf_summary["PP.H3.abf"] + abf_summary["PP.H4.abf"])

        data.frame(
            Locus = locus_id,
            PP4 = abf_summary["PP.H4.abf"],
            PP3 = abf_summary["PP.H3.abf"],
            PP4_conditional = pp4_conditional,
            n_credible_clusters = n_clusters,
            has_coloc_signal = abf_summary["PP.H4.abf"] >= 0.5,
            stringsAsFactors = FALSE
        )
    })

    summary_df <- do.call(rbind, summary_list)

    if (!is.null(output_file) && !is.null(summary_df)) {
        fwrite(summary_df, output_file)
        message(sprintf("[Summary] Saved SuSiE summary: %s", output_file))
    }

    return(summary_df)
}
