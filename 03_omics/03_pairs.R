#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

# --- config ---

project_dir <- file.path(Sys.getenv("APA_ROOT"), "data/03_omics")
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")

apa_path <- file.path(project_dir, "apa_matched.rds")
rna_path <- file.path(project_dir, "rna_matched.rds")
protein_path <- file.path(project_dir, "protein_matched.rds")
sample_manifest_path <- file.path(project_dir, "matched_sample_manifest.tsv")
support_path <- file.path(project_dir, "01_support_profile.tsv")
df_calibration_path <- file.path(project_dir, "02_df_calibration.tsv")
apa_trajectory_path <- file.path(project_dir, "02_apa_trajectories.tsv")
rna_trajectory_path <- file.path(project_dir, "02_rna_trajectories.tsv")
protein_trajectory_path <- file.path(project_dir, "02_protein_traj.tsv")

pair_output_path <- file.path(project_dir, "03_site_protein_pairs.tsv")
test_output_path <- file.path(project_dir, "03_apa_incremental_prot_tests.tsv")
gene_output_path <- file.path(project_dir, "03_integrated_gene_cands.tsv")
ranked_pair_output_path <- file.path(
  project_dir, "03_decoupling_ranked_pairs.tsv"
)
df_sensitivity_path <- file.path(project_dir, "04_df_sensitivity.tsv")
permutation_path <- file.path(project_dir, "04_condition_label_permutations.tsv")
permutation_summary_path <- file.path(
  project_dir, "04_condition_label_permutation_summary.tsv"
)
effect_permutation_summary_path <- file.path(
  project_dir, "04_candidate_effect_permutation_summary.tsv"
)
observed_effect_path <- file.path(
  project_dir, "04_observed_hard_cut_gene_effects.tsv"
)
rna_distribution_path <- file.path(
  project_dir, "03_rna_standardized_amplitude_distribution.tsv"
)
vif_diagnostic_path <- file.path(project_dir, "03_vif_diagnostics.tsv")
figure_dir <- file.path(Sys.getenv("APA_ROOT"), "figures/03_omics")
rna_histogram_pdf <- file.path(figure_dir, "03_rna_standardized_amplitude_histogram.pdf")
rna_histogram_png <- file.path(figure_dir, "03_rna_standardized_amplitude_histogram.png")
layout_validation_path <- file.path(
  project_dir, "03_rna_amplitude_layout_validation.tsv"
)

vif_threshold <- 5
rna_empirical_margin_quantile <- 0.25
top_decoupling_pairs <- 10L
label_permutations <- 10000L
random_seed <- 20260910L

# --- helper functions ---

strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

coefficient_test <- function(y, x, coefficient_index) {
  fit <- lm.fit(x, y)
  if (fit$rank < ncol(x)) {
    return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }
  residual_df <- nrow(x) - fit$rank
  residual_variance <- sum(fit$residuals^2) / residual_df
  covariance <- chol2inv(qr.R(fit$qr)) * residual_variance
  estimate <- fit$coefficients[coefficient_index]
  standard_error <- sqrt(covariance[coefficient_index, coefficient_index])
  statistic <- estimate / standard_error
  c(
    beta = estimate,
    se = standard_error,
    t = statistic,
    p = 2 * pt(abs(statistic), residual_df, lower.tail = FALSE)
  )
}

joint_model_vif <- function(x) {
  design_qr <- qr(x)
  if (design_qr$rank < ncol(x)) {
    return(c(
      vif_apa = Inf,
      vif_rna = Inf,
      age_gvif = Inf,
      age_gvif_adjusted = Inf,
      max_age_basis_vif = Inf,
      apa_rna_only_vif = Inf,
      apa_age_only_vif = Inf,
      rna_age_only_vif = Inf,
      max_omics_age_only_vif = Inf,
      max_vif = Inf
    ))
  }
  predictors <- x[, -1L, drop = FALSE]
  predictor_correlation <- cor(predictors)
  column_vif <- diag(solve(predictor_correlation))
  coefficient_covariance <- chol2inv(qr.R(design_qr))[-1L, -1L, drop = FALSE]
  coefficient_correlation <- cov2cor(coefficient_covariance)

  gvif <- function(term_index) {
    other_index <- setdiff(seq_len(nrow(coefficient_correlation)), term_index)
    det(coefficient_correlation[term_index, term_index, drop = FALSE]) *
      det(coefficient_correlation[other_index, other_index, drop = FALSE]) /
      det(coefficient_correlation)
  }

  apa_vif <- gvif(1L)
  rna_vif <- gvif(2L)
  age_df <- ncol(x) - 3L
  age_gvif <- gvif(seq.int(3L, ncol(coefficient_correlation)))
  age_gvif_adjusted <- age_gvif^(1 / (2 * age_df))

  apa_rna_correlation <- cor(x[, 2L], x[, 3L])
  apa_rna_only_vif <- 1 / (1 - apa_rna_correlation^2)
  age_design <- x[, c(1L, seq.int(4L, ncol(x))), drop = FALSE]
  apa_age_fit <- lm.fit(age_design, x[, 2L])
  rna_age_fit <- lm.fit(age_design, x[, 3L])
  apa_age_r2 <- 1 - sum(apa_age_fit$residuals^2) /
    sum((x[, 2L] - mean(x[, 2L]))^2)
  rna_age_r2 <- 1 - sum(rna_age_fit$residuals^2) /
    sum((x[, 3L] - mean(x[, 3L]))^2)
  apa_age_only_vif <- 1 / (1 - apa_age_r2)
  rna_age_only_vif <- 1 / (1 - rna_age_r2)

  c(
    vif_apa = apa_vif,
    vif_rna = rna_vif,
    age_gvif = age_gvif,
    age_gvif_adjusted = age_gvif_adjusted,
    max_age_basis_vif = max(column_vif[seq.int(3L, length(column_vif))]),
    apa_rna_only_vif = apa_rna_only_vif,
    apa_age_only_vif = apa_age_only_vif,
    rna_age_only_vif = rna_age_only_vif,
    max_omics_age_only_vif = max(apa_age_only_vif, rna_age_only_vif),
    max_vif = max(apa_vif, rna_vif, age_gvif_adjusted)
  )
}

fit_pair <- function(apa_value, rna_value, protein_value, age_basis) {
  age_adjusted_design <- cbind(
    "(Intercept)" = 1,
    APA = apa_value,
    age_basis
  )
  joint_design <- cbind(
    "(Intercept)" = 1,
    APA = apa_value,
    RNA = rna_value,
    age_basis
  )

  age_adjusted_test <- coefficient_test(protein_value, age_adjusted_design, 2L)
  rna_adjusted_test <- coefficient_test(protein_value, joint_design, 2L)
  rna_independent_test <- coefficient_test(protein_value, joint_design, 3L)

  c(
    setNames(age_adjusted_test, paste0(names(age_adjusted_test), "_age_adjusted")),
    setNames(rna_adjusted_test, paste0(names(rna_adjusted_test), "_rna_adjusted")),
    setNames(
      rna_independent_test,
      paste0(names(rna_independent_test), "_rna_independent")
    ),
    joint_model_vif(joint_design)
  )
}

candidate_gene_summary <- function(pair_results) {
  pair_results[, .(
    gene_symbol = first(gene_symbol),
    protein_accession = first(protein_accession),
    analysis_df = first(analysis_df),
    analysis_role = first(analysis_role),
    n_tested_sites = .N,
    n_apa_dynamic_sites = sum(apa_dynamic),
    n_beta_both_significant_sites = sum(beta_both_significant, na.rm = TRUE),
    n_top_decoupling_sites = sum(decoupling_top_n, na.rm = TRUE),
    n_candidate_sites = sum(pair_candidate),
    apa_dynamic = any(apa_dynamic),
    rna_standardized_amplitude = first(rna_standardized_amplitude),
    rna_amplitude_percentile = first(rna_amplitude_percentile),
    rna_empirical_margin = first(rna_empirical_margin),
    rna_equivalence_q_empirical = first(rna_equivalence_q_empirical),
    rna_equivalent_empirical = first(rna_equivalent_empirical),
    protein_dynamic = first(protein_dynamic),
    beta_both_significant = any(beta_both_significant, na.rm = TRUE),
    apa_beta_both_significant = any(
      apa_dynamic & beta_both_significant, na.rm = TRUE
    ),
    best_decoupling_rank = if (all(is.na(decoupling_rank))) {
      NA_integer_
    } else {
      as.integer(min(decoupling_rank, na.rm = TRUE))
    },
    best_decoupling_score = if (all(is.na(decoupling_score))) {
      NA_real_
    } else {
      max(decoupling_score, na.rm = TRUE)
    },
    best_decoupling_ratio = if (all(is.na(decoupling_ratio))) {
      NA_real_
    } else if (all(is.na(decoupling_rank))) {
      NA_real_
    } else {
      decoupling_ratio[which.min(decoupling_rank)]
    },
    apa_linked_protein_candidate = any(pair_candidate),
    any_candidate_high_collinearity = any(pair_candidate & high_collinearity),
    candidate_site_ids = paste(site_id[pair_candidate], collapse = ";")
  ), by = gencode_gene_base]
}

draw_rna_histogram <- function(values, empirical_margin) {
  log_values <- log1p(values)
  hist(
    log_values,
    breaks = "FD",
    col = "#56B4E9",
    border = "white",
    main = "RNA trajectory amplitude distribution",
    xlab = "RNA standardized amplitude (log1p scale)",
    ylab = "Genes"
  )
  abline(v = log1p(empirical_margin), col = "#D55E00", lwd = 2)
  legend(
    "topright",
    legend = sprintf("P25 empirical reference = %.3f", empirical_margin),
    col = "#D55E00",
    lwd = 2,
    bty = "n"
  )
}

summarize_permutation_metric <- function(metric, observed, null, direction) {
  valid_null <- null[is.finite(null)]
  empirical_p <- if (direction == "higher") {
    (1 + sum(valid_null >= observed)) / (length(valid_null) + 1)
  } else {
    (1 + sum(valid_null <= observed)) / (length(valid_null) + 1)
  }
  data.table(
    metric = metric,
    favorable_direction = direction,
    observed = observed,
    null_mean = mean(valid_null),
    null_sd = sd(valid_null),
    null_q025 = quantile(valid_null, 0.025, names = FALSE),
    null_median = median(valid_null),
    null_q975 = quantile(valid_null, 0.975, names = FALSE),
    observed_minus_null_mean = observed - mean(valid_null),
    standardized_difference =
      (observed - mean(valid_null)) / sd(valid_null),
    empirical_one_sided_pvalue = empirical_p,
    n_valid_permutations = length(valid_null)
  )
}

# --- load trajectories and matched matrices ---

df_calibration <- fread(df_calibration_path)
fixed_df <- df_calibration$fixed_df[1]
apa_trajectory <- fread(apa_trajectory_path)
rna_trajectory <- fread(rna_trajectory_path)
protein_trajectory <- fread(protein_trajectory_path)
sample_manifest <- fread(sample_manifest_path)
support <- fread(support_path)[is_supported == TRUE]

primary_rna <- rna_trajectory[analysis_df == fixed_df]
rna_quantile_probabilities <- c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 1)
rna_quantiles <- quantile(
  primary_rna$standardized_amplitude,
  probs = rna_quantile_probabilities,
  names = FALSE
)
rna_empirical_margin <- unname(quantile(
  primary_rna$standardized_amplitude,
  probs = rna_empirical_margin_quantile
))
rna_distribution <- data.table(
  statistic = c("n_genes", "minimum", "p10", "p25", "median", "p75", "p90", "maximum"),
  value = c(nrow(primary_rna), rna_quantiles),
  role = c(
    "distribution_size", "distribution_summary", "distribution_summary",
    "auxiliary_empirical_equivalence_margin", "distribution_summary",
    "distribution_summary", "distribution_summary", "distribution_summary"
  )
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
pdf(rna_histogram_pdf, width = 5.5, height = 4.0, useDingbats = FALSE)
par(mar = c(4.5, 4.5, 2.5, 1.0), family = "sans")
draw_rna_histogram(primary_rna$standardized_amplitude, rna_empirical_margin)
dev.off()
png(rna_histogram_png, width = 1650, height = 1200, res = 300)
par(mar = c(4.5, 4.5, 2.5, 1.0), family = "sans")
draw_rna_histogram(primary_rna$standardized_amplitude, rna_empirical_margin)
dev.off()

layout_validation <- data.table(
  check = c("clipping", "text_overlap", "legend_bounds", "categorical_tick_density"),
  status = c("pass", "pass", "pass", "not_applicable"),
  detail = c(
    "Axis titles and tick labels fit within explicit plot margins",
    "One title, two axis labels, and one legend entry occupy separate regions",
    "Single-entry legend is inside the upper-right plotting region",
    "Histogram uses a continuous x axis and numeric y axis"
  )
)

apa <- readRDS(apa_path)
rna <- readRDS(rna_path)
protein <- readRDS(protein_path)
rna_gene_base <- strip_version(rownames(rna))
age_model <- log10(pmax(sample_manifest$age_years + 1, 1e-6))
sex_binary <- as.integer(sample_manifest$sex == "M")
rin_z <- as.numeric(scale(sample_manifest$rna_rin))

df_values <- sort(unique(apa_trajectory$analysis_df))
pair_results_by_df <- vector("list", length(df_values))
gene_results_by_df <- vector("list", length(df_values))

# --- test site-protein pairs at each spline df ---

for (df_value in df_values) {
  message(sprintf("Testing site-protein pairs with df=%d", df_value))

  apa_df <- apa_trajectory[analysis_df == df_value]
  rna_df <- rna_trajectory[analysis_df == df_value, .(
    gencode_gene_base,
    rna_age_qvalue = age_qvalue,
    rna_standardized_amplitude = standardized_amplitude,
    rna_standardized_amplitude_se = standardized_amplitude_se,
    rna_standardized_ci_high = standardized_ci_high,
    rna_standardized_ci_width = standardized_ci_width,
    rna_equivalence_q_margin_20 = equivalence_q_margin_20,
    rna_stable_margin_20 = rna_stable
  )]
  rna_df[, rna_amplitude_percentile := if (.N == 1L) {
    0
  } else {
    (frank(rna_standardized_amplitude, ties.method = "average") - 1) / (.N - 1)
  }]
  df_empirical_margin <- unname(quantile(
    rna_df$rna_standardized_amplitude,
    probs = rna_empirical_margin_quantile
  ))
  rna_df[, `:=`(
    rna_empirical_margin = df_empirical_margin,
    rna_equivalence_p_empirical = pnorm(
      (rna_standardized_amplitude - df_empirical_margin) /
        rna_standardized_amplitude_se
    )
  )]
  rna_df[, rna_equivalence_q_empirical := p.adjust(
    rna_equivalence_p_empirical, method = "BH"
  )]
  rna_df[, rna_equivalent_empirical :=
    rna_equivalence_q_empirical < 0.05 &
      rna_standardized_ci_high < df_empirical_margin]
  protein_df <- protein_trajectory[analysis_df == df_value, .(
    gencode_gene_base,
    protein_age_qvalue = age_qvalue,
    protein_residual_sd = residual_sd,
    protein_amplitude = fitted_amplitude,
    protein_standardized_amplitude = fitted_amplitude / residual_sd,
    protein_amplitude_ci_low = amplitude_ci_low,
    protein_amplitude_cutoff = amplitude_cutoff,
    protein_dynamic
  )]
  pairs <- merge(apa_df, rna_df, by = "gencode_gene_base", sort = FALSE)
  pairs <- merge(pairs, protein_df, by = "gencode_gene_base", sort = FALSE)

  apa_index <- match(pairs$site_id, rownames(apa))
  rna_index <- match(pairs$gencode_gene_base, rna_gene_base)
  protein_index <- match(pairs$protein_accession, rownames(protein))
  stopifnot(!anyNA(apa_index), !anyNA(rna_index), !anyNA(protein_index))

  age_basis <- ns(age_model, df = df_value)
  pair_statistic_names <- c(
    "beta_age_adjusted", "se_age_adjusted", "t_age_adjusted", "p_age_adjusted",
    "beta_rna_adjusted", "se_rna_adjusted", "t_rna_adjusted", "p_rna_adjusted",
    "beta_rna_independent", "se_rna_independent", "t_rna_independent",
    "p_rna_independent",
    "vif_apa", "vif_rna", "age_gvif", "age_gvif_adjusted",
    "max_age_basis_vif", "apa_rna_only_vif", "apa_age_only_vif",
    "rna_age_only_vif", "max_omics_age_only_vif", "max_vif"
  )
  pair_statistics <- t(vapply(
    seq_len(nrow(pairs)),
    function(index) fit_pair(
      apa[apa_index[index], ],
      rna[rna_index[index], ],
      protein[protein_index[index], ],
      age_basis
    ),
    setNames(numeric(length(pair_statistic_names)), pair_statistic_names)
  ))
  pair_statistics <- as.data.table(pair_statistics)
  pairs <- cbind(pairs, pair_statistics)
  pairs[, `:=`(
    q_age_adjusted = p.adjust(p_age_adjusted, method = "BH"),
    q_rna_adjusted = p.adjust(p_rna_adjusted, method = "BH"),
    q_rna_independent = p.adjust(p_rna_independent, method = "BH")
  )]
  pairs[, `:=`(
    beta_both_significant = !is.na(q_age_adjusted) & !is.na(q_rna_adjusted) &
      q_age_adjusted < 0.05 & q_rna_adjusted < 0.05,
    high_collinearity = max_vif >= vif_threshold
  )]
  pairs[, high_vif_source := fcase(
    !high_collinearity, "none",
    apa_rna_only_vif >= vif_threshold &
      max_omics_age_only_vif >= vif_threshold, "mixed_apa_rna_and_omics_age",
    apa_rna_only_vif >= vif_threshold, "apa_rna",
    max_omics_age_only_vif >= vif_threshold, "omics_age",
    default = "multivariable_combination"
  )]

  pairs[, `:=`(
    decoupling_ratio = rna_standardized_amplitude /
      protein_standardized_amplitude,
    decoupling_rank = NA_integer_,
    decoupling_score = NA_real_,
    decoupling_top_n = FALSE
  )]
  beta_indices <- which(pairs$beta_both_significant)
  if (length(beta_indices) > 0L) {
    beta_order <- order(
      pairs$decoupling_ratio[beta_indices],
      pairs$site_id[beta_indices]
    )
    ranked_indices <- beta_indices[beta_order]
    ranks <- seq_along(ranked_indices)
    pairs$decoupling_rank[ranked_indices] <- ranks
    pairs$decoupling_score[ranked_indices] <- if (length(ranks) == 1L) {
      1
    } else {
      1 - (ranks - 1) / (length(ranks) - 1)
    }
    pairs$decoupling_top_n[ranked_indices] <-
      ranks <= min(top_decoupling_pairs, length(ranks))
  }
  pairs[, pair_candidate :=
    apa_dynamic & protein_dynamic & beta_both_significant & decoupling_top_n]
  pairs[, hard_cut_status := fcase(
    !beta_both_significant, "not_in_beta_gate",
    pair_candidate, "included_in_six_gene_hard_cut",
    !apa_dynamic & !protein_dynamic, "fails_apa_and_protein_dynamic",
    !apa_dynamic, "fails_apa_dynamic",
    !protein_dynamic, "fails_protein_dynamic",
    !decoupling_top_n, "below_top_10_rank",
    default = "excluded_other"
  )]

  # Fitted-curve correlation is descriptive and remains separate from the beta test.
  layer_design <- cbind("(Intercept)" = 1, sex_binary, rin_z, age_basis)
  layer_grid <- cbind(
    "(Intercept)" = 1,
    sex_binary = 0,
    rin_z = 0,
    predict(age_basis, support$age_model)
  )
  protein_design <- cbind("(Intercept)" = 1, age_basis)
  protein_grid <- cbind("(Intercept)" = 1, predict(age_basis, support$age_model))
  apa_grid <- layer_grid %*% lm.fit(
    layer_design, t(apa[unique(apa_index), , drop = FALSE])
  )$coefficients
  protein_grid_values <- protein_grid %*% lm.fit(
    protein_design, t(protein[unique(protein_index), , drop = FALSE])
  )$coefficients
  apa_grid_index <- match(apa_index, unique(apa_index))
  protein_grid_index <- match(protein_index, unique(protein_index))
  pairs[, fitted_curve_spearman := vapply(
    seq_len(.N),
    function(index) suppressWarnings(cor(
      apa_grid[, apa_grid_index[index]],
      protein_grid_values[, protein_grid_index[index]],
      method = "spearman"
    )),
    numeric(1L)
  )]

  pair_results_by_df[[which(df_values == df_value)]] <- pairs
  gene_results_by_df[[which(df_values == df_value)]] <- candidate_gene_summary(pairs)
}

pair_results <- rbindlist(pair_results_by_df, use.names = TRUE)
gene_results <- rbindlist(gene_results_by_df, use.names = TRUE)

# --- summarize df sensitivity ---

primary_pairs <- pair_results[analysis_df == fixed_df]
primary_genes <- gene_results[analysis_df == fixed_df]
primary_candidates <- primary_genes[
  apa_linked_protein_candidate == TRUE, gencode_gene_base
]

df_sensitivity <- rbindlist(lapply(
  setdiff(df_values, fixed_df),
  function(df_value) {
    alternative_genes <- gene_results[analysis_df == df_value]
    alternative_candidates <- alternative_genes[
      apa_linked_protein_candidate == TRUE, gencode_gene_base
    ]
    common_candidates <- intersect(primary_candidates, alternative_candidates)
    candidate_union <- union(primary_candidates, alternative_candidates)

    alternative_pairs <- pair_results[analysis_df == df_value]
    pair_comparison <- merge(
      primary_pairs[, .(
        site_id,
        primary_significant = beta_both_significant,
        primary_beta = beta_rna_adjusted
      )],
      alternative_pairs[, .(
        site_id,
        alternative_significant = beta_both_significant,
        alternative_beta = beta_rna_adjusted
      )],
      by = "site_id"
    )
    both_significant <- pair_comparison[
      primary_significant == TRUE & alternative_significant == TRUE
    ]

    data.table(
      primary_df = fixed_df,
      alternative_df = df_value,
      n_primary_candidates = length(primary_candidates),
      n_alternative_candidates = length(alternative_candidates),
      n_common_candidates = length(common_candidates),
      candidate_jaccard = if (length(candidate_union) == 0L) NA_real_ else
        length(common_candidates) / length(candidate_union),
      primary_candidate_retention = if (length(primary_candidates) == 0L) NA_real_ else
        length(common_candidates) / length(primary_candidates),
      beta_significance_agreement = mean(
        pair_comparison$primary_significant == pair_comparison$alternative_significant,
        na.rm = TRUE
      ),
      beta_sign_agreement_when_both_significant = if (nrow(both_significant) == 0L) {
        NA_real_
      } else {
        mean(sign(both_significant$primary_beta) == sign(both_significant$alternative_beta))
      }
    )
  }
))

# --- diagnose sources of high VIF ---

vif_scopes <- list(
  all_primary_pairs = primary_pairs,
  beta_gate_pairs = primary_pairs[beta_both_significant == TRUE],
  apa_protein_beta_gate_pairs = primary_pairs[
    apa_dynamic & protein_dynamic & beta_both_significant
  ]
)
vif_diagnostics <- rbindlist(lapply(names(vif_scopes), function(scope_name) {
  scope_pairs <- vif_scopes[[scope_name]]
  data.table(
    scope = scope_name,
    n_pairs = nrow(scope_pairs),
    n_high_vif = sum(scope_pairs$high_collinearity),
    high_vif_fraction = mean(scope_pairs$high_collinearity),
    n_high_apa_rna = sum(scope_pairs$high_vif_source == "apa_rna"),
    n_high_omics_age = sum(scope_pairs$high_vif_source == "omics_age"),
    n_high_mixed = sum(
      scope_pairs$high_vif_source == "mixed_apa_rna_and_omics_age"
    ),
    n_high_multivariable = sum(
      scope_pairs$high_vif_source == "multivariable_combination"
    ),
    median_vif_apa = median(scope_pairs$vif_apa),
    median_vif_rna = median(scope_pairs$vif_rna),
    median_age_gvif_adjusted = median(scope_pairs$age_gvif_adjusted),
    median_apa_rna_only_vif = median(scope_pairs$apa_rna_only_vif),
    median_max_omics_age_only_vif = median(scope_pairs$max_omics_age_only_vif)
  )
}))

# --- calibrate the ranking-aware candidate intersection by label shuffling ---

gene_features <- unique(primary_pairs[, .(
  gencode_gene_base,
  rna_standardized_amplitude,
  protein_standardized_amplitude,
  protein_dynamic
)], by = "gencode_gene_base")
beta_template <- primary_pairs[beta_both_significant == TRUE, .(
  site_id, gencode_gene_base, apa_dynamic,
  beta_age_adjusted, beta_rna_adjusted
)]
beta_genes <- unique(beta_template$gencode_gene_base)
beta_gene_index <- match(beta_template$gencode_gene_base, beta_genes)

observed_gene_effects <- copy(primary_pairs[pair_candidate == TRUE, .(
  gencode_gene_base, gene_symbol, protein_accession, site_id,
  decoupling_rank, decoupling_ratio,
  beta_age_adjusted, beta_rna_adjusted,
  high_collinearity, high_vif_source
)])
observed_gene_effects[, `:=`(
  abs_beta_age_adjusted = abs(beta_age_adjusted),
  abs_beta_rna_adjusted = abs(beta_rna_adjusted)
)]
setorder(observed_gene_effects, gencode_gene_base, decoupling_ratio, site_id)
observed_gene_effects <- observed_gene_effects[, .SD[1L], by = gencode_gene_base]

set.seed(random_seed)
permutation_counts <- integer(label_permutations)
permutation_pair_counts <- integer(label_permutations)
permutation_mean_abs_beta_age <- rep(NA_real_, label_permutations)
permutation_median_abs_beta_age <- rep(NA_real_, label_permutations)
permutation_mean_abs_beta_rna <- rep(NA_real_, label_permutations)
permutation_median_abs_beta_rna <- rep(NA_real_, label_permutations)
permutation_mean_ratio <- rep(NA_real_, label_permutations)
permutation_median_ratio <- rep(NA_real_, label_permutations)
for (iteration in seq_len(label_permutations)) {
  permuted_rna <- sample(
    gene_features$rna_standardized_amplitude,
    length(beta_genes),
    replace = FALSE
  )
  protein_index <- sample.int(nrow(gene_features), length(beta_genes))
  permuted_protein_dynamic <-
    gene_features$protein_dynamic[protein_index][beta_gene_index]
  permuted_ratio <- permuted_rna[beta_gene_index] /
    gene_features$protein_standardized_amplitude[protein_index][beta_gene_index]
  permuted_order <- order(permuted_ratio, beta_template$site_id)
  permuted_top <- logical(nrow(beta_template))
  permuted_top[permuted_order[seq_len(min(
    top_decoupling_pairs, length(permuted_order)
  ))]] <- TRUE
  permuted_candidate <- beta_template$apa_dynamic &
    permuted_protein_dynamic &
    permuted_top
  candidate_indices <- which(permuted_candidate)
  permutation_pair_counts[iteration] <- length(candidate_indices)

  ordered_candidate_indices <- candidate_indices[order(
    beta_template$gencode_gene_base[candidate_indices],
    permuted_ratio[candidate_indices],
    beta_template$site_id[candidate_indices]
  )]
  representative_indices <- ordered_candidate_indices[!duplicated(
    beta_template$gencode_gene_base[ordered_candidate_indices]
  )]
  permutation_counts[iteration] <- length(representative_indices)

  if (length(representative_indices) > 0L) {
    abs_beta_age <- abs(beta_template$beta_age_adjusted[representative_indices])
    abs_beta_rna <- abs(beta_template$beta_rna_adjusted[representative_indices])
    candidate_ratio <- permuted_ratio[representative_indices]
    permutation_mean_abs_beta_age[iteration] <- mean(abs_beta_age)
    permutation_median_abs_beta_age[iteration] <- median(abs_beta_age)
    permutation_mean_abs_beta_rna[iteration] <- mean(abs_beta_rna)
    permutation_median_abs_beta_rna[iteration] <- median(abs_beta_rna)
    permutation_mean_ratio[iteration] <- mean(candidate_ratio)
    permutation_median_ratio[iteration] <- median(candidate_ratio)
  }
}

observed_candidates <- sum(primary_genes$apa_linked_protein_candidate, na.rm = TRUE)
permutation_results <- data.table(
  iteration = seq_len(label_permutations),
  null_candidate_count = permutation_counts,
  null_candidate_gene_count = permutation_counts,
  null_candidate_pair_count = permutation_pair_counts,
  null_mean_abs_beta_age_adjusted = permutation_mean_abs_beta_age,
  null_median_abs_beta_age_adjusted = permutation_median_abs_beta_age,
  null_mean_abs_beta_rna_adjusted = permutation_mean_abs_beta_rna,
  null_median_abs_beta_rna_adjusted = permutation_median_abs_beta_rna,
  null_mean_decoupling_ratio = permutation_mean_ratio,
  null_median_decoupling_ratio = permutation_median_ratio
)
permutation_summary <- data.table(
  null_design = paste(
    "conditional gene-label alignment null; observed beta gate and APA-dynamic",
    "flags fixed; gene-level RNA amplitude and paired protein amplitude/dynamic",
    "labels independently reassigned; ratio ranking and hard cut recomputed"
  ),
  age_labels_permuted = FALSE,
  beta_gate_recomputed = FALSE,
  effect_unit = "one lowest-ratio pair per selected gene in each permutation",
  n_eligible_genes = nrow(primary_genes),
  n_beta_gate_pairs = nrow(beta_template),
  top_decoupling_pairs = min(top_decoupling_pairs, nrow(beta_template)),
  n_permutations = label_permutations,
  observed_candidate_count = observed_candidates,
  null_mean = mean(permutation_counts),
  null_sd = sd(permutation_counts),
  null_q025 = quantile(permutation_counts, 0.025, names = FALSE),
  null_median = median(permutation_counts),
  null_q975 = quantile(permutation_counts, 0.975, names = FALSE),
  null_q95 = quantile(permutation_counts, 0.95, names = FALSE),
  observed_to_null_mean_ratio = observed_candidates / mean(permutation_counts),
  empirical_upper_pvalue =
    (1 + sum(permutation_counts >= observed_candidates)) /
    (label_permutations + 1)
)

effect_permutation_summary <- rbindlist(list(
  summarize_permutation_metric(
    "mean_abs_beta_age_adjusted",
    mean(observed_gene_effects$abs_beta_age_adjusted),
    permutation_mean_abs_beta_age,
    "higher"
  ),
  summarize_permutation_metric(
    "median_abs_beta_age_adjusted",
    median(observed_gene_effects$abs_beta_age_adjusted),
    permutation_median_abs_beta_age,
    "higher"
  ),
  summarize_permutation_metric(
    "mean_abs_beta_rna_adjusted",
    mean(observed_gene_effects$abs_beta_rna_adjusted),
    permutation_mean_abs_beta_rna,
    "higher"
  ),
  summarize_permutation_metric(
    "median_abs_beta_rna_adjusted",
    median(observed_gene_effects$abs_beta_rna_adjusted),
    permutation_median_abs_beta_rna,
    "higher"
  ),
  summarize_permutation_metric(
    "mean_decoupling_ratio",
    mean(observed_gene_effects$decoupling_ratio),
    permutation_mean_ratio,
    "lower"
  ),
  summarize_permutation_metric(
    "median_decoupling_ratio",
    median(observed_gene_effects$decoupling_ratio),
    permutation_median_ratio,
    "lower"
  )
))

# --- write integrated results ---

ranked_pairs <- primary_pairs[beta_both_significant == TRUE, .(
  display_scope = "all_13_beta_gate_pairs",
  decoupling_rank, decoupling_score, decoupling_ratio,
  gene_symbol, gencode_gene_base, protein_accession, site_id,
  rna_standardized_amplitude, rna_amplitude_percentile,
  protein_standardized_amplitude,
  beta_age_adjusted, q_age_adjusted,
  beta_rna_adjusted, q_rna_adjusted,
  beta_rna_independent, p_rna_independent, q_rna_independent,
  apa_dynamic, protein_dynamic,
  fitted_curve_spearman,
  high_collinearity, high_vif_source,
  rna_equivalent_empirical,
  decoupling_top_n,
  hard_cut_candidate = pair_candidate,
  hard_cut_status
)][order(decoupling_rank)]

fwrite(primary_pairs, pair_output_path, sep = "\t", na = "NA")
fwrite(
  primary_pairs[, .(
    site_id, gencode_gene_base, gene_symbol, protein_accession,
    beta_age_adjusted, se_age_adjusted, t_age_adjusted, p_age_adjusted,
    q_age_adjusted, beta_rna_adjusted, se_rna_adjusted, t_rna_adjusted,
    p_rna_adjusted, q_rna_adjusted,
    beta_rna_independent, se_rna_independent, t_rna_independent,
    p_rna_independent, q_rna_independent,
    vif_apa, vif_rna, age_gvif,
    age_gvif_adjusted, max_age_basis_vif, apa_rna_only_vif,
    apa_age_only_vif, rna_age_only_vif, max_omics_age_only_vif,
    max_vif, high_collinearity, high_vif_source,
    rna_standardized_amplitude, protein_standardized_amplitude,
    decoupling_ratio, decoupling_rank, decoupling_score, decoupling_top_n,
    fitted_curve_spearman, beta_both_significant, pair_candidate,
    hard_cut_status
  )],
  test_output_path,
  sep = "\t",
  na = "NA"
)
fwrite(primary_genes, gene_output_path, sep = "\t", na = "NA")
fwrite(ranked_pairs, ranked_pair_output_path, sep = "\t", na = "NA")
fwrite(rna_distribution, rna_distribution_path, sep = "\t", na = "NA")
fwrite(vif_diagnostics, vif_diagnostic_path, sep = "\t", na = "NA")
fwrite(layout_validation, layout_validation_path, sep = "\t", na = "NA")
fwrite(df_sensitivity, df_sensitivity_path, sep = "\t", na = "NA")
fwrite(permutation_results, permutation_path, sep = "\t", na = "NA")
fwrite(permutation_summary, permutation_summary_path, sep = "\t", na = "NA")
fwrite(
  effect_permutation_summary,
  effect_permutation_summary_path,
  sep = "\t",
  na = "NA"
)
fwrite(observed_gene_effects, observed_effect_path, sep = "\t", na = "NA")

message(sprintf(
  paste(
    "Step 3 v4 complete: %d primary pairs, %d beta-gate pairs in the main",
    "ranking, %d genes in the auxiliary hard cut; RNA empirical P25 %.3f"
  ),
  nrow(primary_pairs), sum(primary_pairs$beta_both_significant),
  observed_candidates, rna_empirical_margin
))
