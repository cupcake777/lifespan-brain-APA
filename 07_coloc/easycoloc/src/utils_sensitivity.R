# =============================================================================
# src/utils_sensitivity.R
# =============================================================================
# Prior Sensitivity Analysis for Colocalization
# =============================================================================
# Tests how colocalization results (PP4) change under different prior
# probability settings. This is critical for demonstrating result robustness.
#
# Reference: Giambartolomei et al. (2014) PLoS Genet
# =============================================================================

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# =============================================================================
# run_sensitivity_analysis: Test PP4 stability across prior settings
# =============================================================================
# Arguments:
#   colocInputFile: Merged GWAS+QTL data
#   gwas_type: "cc" or "quant"
#   gwas_prop: Case proportion (for cc)
#   gwas_N, qtl_N: Sample sizes
#   prior_grid: List of (p1, p2, p12) tuples to test
#   use_susie: Whether to also run SuSiE for each setting
#
# Returns:
#   data.frame with PP4 values for each prior setting
# =============================================================================
run_sensitivity_analysis <- function(colocInputFile,
                                      gwas_type = "cc",
                                      gwas_prop = 0.5,
                                      gwas_N = NULL,
                                      qtl_N = NULL,
                                      prior_grid = NULL,
                                      use_susie = FALSE,
                                      plink_bfile = NULL) {

  # Default prior grid (follows coloc manual recommendations)
  if (is.null(prior_grid)) {
    prior_grid <- list(
      list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),   # Default
      list(p1 = 1e-3, p2 = 1e-4, p12 = 1e-5),   # Higher P(GWAS)
      list(p1 = 1e-4, p2 = 1e-3, p12 = 1e-5),   # Higher P(QTL)
      list(p1 = 1e-3, p2 = 1e-3, p12 = 1e-5),   # Both higher
      list(p1 = 1e-5, p2 = 1e-5, p12 = 1e-6),   # More conservative
      list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-6),   # Lower shared prior
      list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-4)    # Higher shared prior
    )
  }

  message("[Sensitivity] Testing ", length(prior_grid), " prior settings...")

  # Prepare datasets
  N_gwas_val <- if("N.gwas" %in% names(colocInputFile)) colocInputFile$N.gwas else gwas_N
  N_qtl_val <- if("N.qtl" %in% names(colocInputFile)) colocInputFile$N.qtl else qtl_N

  d1 <- list(
    snp = as.character(colocInputFile$snp),
    beta = as.numeric(colocInputFile$BETA.gwas),
    varbeta = as.numeric(colocInputFile$SE.gwas)^2,
    pvalues = as.numeric(colocInputFile$P.gwas),
    type = gwas_type,
    N = as.numeric(N_gwas_val)
  )
  if (gwas_type == "cc") {
    d1$s <- as.numeric(gwas_prop)
  }

  d2 <- list(
    snp = as.character(colocInputFile$snp),
    beta = as.numeric(colocInputFile$BETA.qtl),
    varbeta = as.numeric(colocInputFile$SE.qtl)^2,
    pvalues = as.numeric(colocInputFile$P.qtl),
    type = "quant",
    N = as.numeric(N_qtl_val)
  )

  # Add MAF if available
  if ("MAF.gwas" %in% names(colocInputFile)) {
    d1$MAF <- as.numeric(colocInputFile$MAF.gwas)
    d2$MAF <- as.numeric(colocInputFile$MAF.qtl)
  }

  # Run coloc for each prior setting
  results <- lapply(seq_along(prior_grid), function(i) {
    prior <- prior_grid[[i]]

    res <- tryCatch({
      coloc.abf(d1, d2,
                p1 = prior$p1,
                p2 = prior$p2,
                p12 = prior$p12)
    }, error = function(e) {
      message(sprintf("[Sensitivity] Failed at setting %d: %s", i, e$message))
      return(NULL)
    })

    if (is.null(res)) {
      return(data.frame(
        setting = i,
        p1 = prior$p1,
        p2 = prior$p2,
        p12 = prior$p12,
        PP.H0 = NA,
        PP.H1 = NA,
        PP.H2 = NA,
        PP.H3 = NA,
        PP.H4 = NA,
        logABF = NA
      ))
    }

    data.frame(
      setting = i,
      p1 = prior$p1,
      p2 = prior$p2,
      p12 = prior$p12,
      PP.H0 = res$summary["PP.H0.abf"],
      PP.H1 = res$summary["PP.H1.abf"],
      PP.H2 = res$summary["PP.H2.abf"],
      PP.H3 = res$summary["PP.H3.abf"],
      PP.H4 = res$summary["PP.H4.abf"],
      logABF = res$summary["logABF"]
    )
  })

  results_df <- bind_rows(results)

  # Add interpretive labels
  results_df <- results_df %>%
    mutate(
      interpretation = case_when(
        PP.H4 >= 0.8 ~ "strong_coloc",
        PP.H4 >= 0.5 ~ "suggestive_coloc",
        PP.H3 + PP.H4 >= 0.5 ~ "shared_signal",
        PP.H1 >= 0.5 ~ "gwas_only",
        PP.H2 >= 0.5 ~ "qtl_only",
        TRUE ~ "no_signal"
      )
    )

  message(sprintf("[Sensitivity] PP4 range: %.3f - %.3f",
                  min(results_df$PP.H4, na.rm = TRUE),
                  max(results_df$PP.H4, na.rm = TRUE)))

  # Calculate sensitivity metrics
  pp4_values <- results_df$PP.H4[!is.na(results_df$PP.H4)]
  if (length(pp4_values) > 1) {
    results_df$pp4_range <- max(pp4_values) - min(pp4_values)
    results_df$pp4_cv <- sd(pp4_values) / mean(pp4_values)

    message(sprintf("[Sensitivity] PP4 range: %.3f, CV: %.2f",
                    results_df$pp4_range[1],
                    results_df$pp4_cv[1]))
  }

  return(results_df)
}

# =============================================================================
# plot_sensitivity_results: Visualize prior sensitivity
# =============================================================================
# Creates a bar plot showing PP4 values across prior settings
#
# Arguments:
#   sensitivity_df: Results from run_sensitivity_analysis
#   output_file: Path to save plot (PDF or PNG)
#
# Returns:
#   ggplot object (also saved to file if output_file specified)
# =============================================================================
plot_sensitivity_results <- function(sensitivity_df,
                                      output_file = NULL,
                                      title = "Prior Sensitivity Analysis") {

  # Reshape for plotting
  plot_df <- sensitivity_df %>%
    select(setting, p1, p2, p12, PP.H4) %>%
    mutate(
      setting_label = sprintf("S%d\n(p1=%.0e, p2=%.0e, p12=%.0e)",
                               setting, p1, p2, p12),
      pp4_category = case_when(
        PP.H4 >= 0.8 ~ "Strong (≥0.8)",
        PP.H4 >= 0.5 ~ "Suggestive (0.5-0.8)",
        TRUE ~ "Weak (<0.5)"
      )
    )

  # Create bar plot
  p <- ggplot(plot_df, aes(x = setting_label, y = PP.H4, fill = pp4_category)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", PP.H4)),
              vjust = -0.5, size = 3, angle = 45) +
    scale_fill_manual(values = c("Strong (≥0.8)" = "#228B22",
                                  "Suggestive (0.5-0.8)" = "#FFA500",
                                  "Weak (<0.5)" = "#BDBDBD")) +
    labs(
      title = title,
      subtitle = "PP.H4 (colocalization probability) across different prior settings",
      x = "Prior Setting",
      y = "PP.H4 (Colocalization Probability)",
      fill = "Evidence Level"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 5, dpi = 300)
    message(sprintf("[Plot] Saved sensitivity plot: %s", output_file))
  }

  return(p)
}

# =============================================================================
# sensitivity_heatmap: 2D sensitivity across p1 and p2
# =============================================================================
# Creates a heatmap showing PP4 across a grid of p1 and p2 values
# (keeping p12/p1 ratio constant)
# =============================================================================
sensitivity_heatmap <- function(colocInputFile,
                                 gwas_type = "cc",
                                 gwas_prop = 0.5,
                                 gwas_N = NULL,
                                 qtl_N = NULL,
                                 p1_values = c(1e-5, 1e-4, 1e-3, 1e-2),
                                 p2_values = c(1e-5, 1e-4, 1e-3, 1e-2),
                                 p12_ratio = 0.1) {

  message("[Heatmap] Testing ", length(p1_values) * length(p2_values), " combinations...")

  N_gwas_val <- if("N.gwas" %in% names(colocInputFile)) colocInputFile$N.gwas else gwas_N
  N_qtl_val <- if("N.qtl" %in% names(colocInputFile)) colocInputFile$N.qtl else qtl_N

  d1 <- list(
    snp = as.character(colocInputFile$snp),
    beta = as.numeric(colocInputFile$BETA.gwas),
    varbeta = as.numeric(colocInputFile$SE.gwas)^2,
    pvalues = as.numeric(colocInputFile$P.gwas),
    type = gwas_type,
    N = as.numeric(N_gwas_val)
  )
  if (gwas_type == "cc") d1$s <- as.numeric(gwas_prop)

  d2 <- list(
    snp = as.character(colocInputFile$snp),
    beta = as.numeric(colocInputFile$BETA.qtl),
    varbeta = as.numeric(colocInputFile$SE.qtl)^2,
    pvalues = as.numeric(colocInputFile$P.qtl),
    type = "quant",
    N = as.numeric(N_qtl_val)
  )

  if ("MAF.gwas" %in% names(colocInputFile)) {
    d1$MAF <- as.numeric(colocInputFile$MAF.gwas)
    d2$MAF <- as.numeric(colocInputFile$MAF.qtl)
  }

  # Full grid search
  results <- expand.grid(p1 = p1_values, p2 = p2_values)
  results$PP4 <- NA_real_
  results$PP3 <- NA_real_

  for (i in 1:nrow(results)) {
    p1 <- results$p1[i]
    p2 <- results$p2[i]
    p12 <- p1 * p12_ratio  # Keep p12 proportional to p1

    res <- tryCatch({
      coloc.abf(d1, d2, p1 = p1, p2 = p2, p12 = p12)
    }, error = function(e) NULL)

    if (!is.null(res)) {
      results$PP4[i] <- res$summary["PP.H4.abf"]
      results$PP3[i] <- res$summary["PP.H3.abf"]
    }
  }

  # Create heatmap
  plot_df <- results %>%
    mutate(
      p1_label = format(p1, scientific = TRUE),
      p2_label = format(p2, scientific = TRUE),
      pp4_category = case_when(
        PP4 >= 0.8 ~ "Strong",
        PP4 >= 0.5 ~ "Suggestive",
        TRUE ~ "Weak"
      )
    )

  p <- ggplot(plot_df, aes(x = p2_label, y = p1_label, fill = PP4)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", PP4)), size = 4) +
    scale_fill_gradientn(
      colors = c("#BDBDBD", "#90EE90", "#228B22", "#006400"),
      values = scales::rescale(c(0, 0.5, 0.8, 1)),
      name = "PP4"
    ) +
    labs(
      title = "Prior Sensitivity Heatmap",
      subtitle = sprintf("PP4 across p1 × p2 grid (p12/p1 = %.2f)", p12_ratio),
      x = "p2 (QTL prior)",
      y = "p1 (GWAS prior)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  return(list(plot = p, data = plot_df))
}

# =============================================================================
# write_sensitivity_report: Generate text summary
# =============================================================================
write_sensitivity_report <- function(sensitivity_df, output_file) {

  sink(output_file)

  cat("=============================================================\n")
  cat("Prior Sensitivity Analysis Report\n")
  cat("=============================================================\n")
  cat(sprintf("Generated: %s\n\n", Sys.time()))

  cat("SETTINGS TESTED\n")
  cat("-------------------------------------------------------------\n")

  for (i in 1:nrow(sensitivity_df)) {
    row <- sensitivity_df[i, ]
    cat(sprintf("\nSetting %d:\n", row$setting))
    cat(sprintf("  p1 = %.2e, p2 = %.2e, p12 = %.2e\n", row$p1, row$p2, row$p12))
    cat(sprintf("  PP.H4 = %.4f (%s)\n", row$PP.H4, row$interpretation))
  }

  # Summary statistics
  pp4_vals <- sensitivity_df$PP.H4[!is.na(sensitivity_df$PP.H4)]

  cat("\n\nSUMMARY STATISTICS\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  PP4 range: %.4f - %.4f (Δ = %.4f)\n",
              min(pp4_vals), max(pp4_vals), max(pp4_vals) - min(pp4_vals)))
  cat(sprintf("  PP4 mean: %.4f\n", mean(pp4_vals)))
  cat(sprintf("  PP4 median: %.4f\n", median(pp4_vals)))
  cat(sprintf("  PP4 CV: %.2f\n", sd(pp4_vals) / mean(pp4_vals)))

  # Robustness assessment
  cat("\n\nROBUSTNESS ASSESSMENT\n")
  cat("-------------------------------------------------------------\n")

  n_strong <- sum(sensitivity_df$PP.H4 >= 0.8, na.rm = TRUE)
  n_suggestive <- sum(sensitivity_df$PP.H4 >= 0.5 & sensitivity_df$PP.H4 < 0.8, na.rm = TRUE)
  n_weak <- sum(sensitivity_df$PP.H4 < 0.5, na.rm = TRUE)
  n_total <- n_strong + n_suggestive + n_weak

  cat(sprintf("  Strong evidence (PP4 ≥ 0.8): %d/%d settings (%.0f%%)\n",
              n_strong, n_total, n_strong/n_total*100))
  cat(sprintf("  Suggestive evidence (0.5 ≤ PP4 < 0.8): %d/%d (%.0f%%)\n",
              n_suggestive, n_total, n_suggestive/n_total*100))
  cat(sprintf("  Weak evidence (PP4 < 0.5): %d/%d (%.0f%%)\n",
              n_weak, n_total, n_weak/n_total*100))

  # Conclusion
  cat("\n\nCONCLUSION\n")
  cat("-------------------------------------------------------------\n")

  if (n_strong / n_total >= 0.8) {
    cat("  ✓ Results are ROBUST across prior settings\n")
    cat("  Strong colocalization evidence in most scenarios\n")
  } else if (n_strong + n_suggestive >= n_total * 0.5) {
    cat("  ~ Results show MODERATE robustness\n")
    cat("  Evidence varies with prior choice\n")
  } else {
    cat("  ✗ Results are SENSITIVE to prior choice\n")
    cat("  Interpret with caution\n")
  }

  sink()

  message(sprintf("[Report] Saved sensitivity report: %s", output_file))
}
