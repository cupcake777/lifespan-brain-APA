# Sensitivity of the ranked APA-RNA-protein candidates to the trajectory model: fixed spline vs adaptive GAM.
# Compares trajectory shapes, pair coefficients and k-check gates for the ranked pairs.

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(splines)
})
setDTthreads(1L)

# --- config ---

if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
mod <- file.path(Sys.getenv("APA_ROOT"), "data/03_omics")

fixed_df <- fread(file.path(mod, "02_df_calibration.tsv"))$fixed_df[1]
k_ladder <- c(20L, 30L, 40L, 60L, 80L, 100L, 120L)
k_index_p_threshold <- 0.05
edf_margin <- 1
edf_review_threshold <- 4.5
shape_correlation_threshold <- 0.90
random_seed <- 20260910L

# --- helpers ---

build_formula <- function(linear_terms, smooth_term = NULL) {
  rhs <- c(linear_terms, smooth_term)
  as.formula(if (length(rhs) == 0L) {
    "y ~ 1"
  } else {
    paste("y ~", paste(rhs, collapse = " + "))
  })
}

fit_adaptive_gam <- function(dat, linear_terms, seed_offset) {
  set.seed(random_seed + seed_offset)
  eligible_k <- k_ladder[k_ladder < length(unique(dat$age_model))]
  for (k in eligible_k) {
    fit <- gam(
      build_formula(
        linear_terms,
        sprintf("s(age_model, k = %d, bs = 'tp')", k)
      ),
      data = dat,
      method = "REML"
    )
    smooth_table <- summary(fit)$s.table
    k_result <- k.check(fit, subsample = 5000, n.rep = 400)
    edf <- as.numeric(smooth_table[1L, "edf"])
    k_index <- as.numeric(k_result[1L, "k-index"])
    k_pvalue <- as.numeric(k_result[1L, "p-value"])
    k_pass <- is.finite(k_pvalue) && k_pvalue > k_index_p_threshold &&
      edf < k - edf_margin
    if (k_pass || k == tail(eligible_k, 1L)) {
      return(list(
        fit = fit,
        final_k = k,
        edf = edf,
        smooth_pvalue = as.numeric(smooth_table[1L, "p-value"]),
        k_index = k_index,
        k_pvalue = k_pvalue,
        k_check_pass = k_pass
      ))
    }
  }
}

extract_coefficient <- function(fit, term) {
  coefficient_table <- if (inherits(fit, "gam")) {
    summary(fit)$p.table
  } else {
    summary(fit)$coefficients
  }
  data.table(
    beta = as.numeric(coefficient_table[term, "Estimate"]),
    se = as.numeric(coefficient_table[term, "Std. Error"]),
    statistic = as.numeric(coefficient_table[term, "t value"]),
    pvalue = as.numeric(coefficient_table[term, "Pr(>|t|)"])
  )
}

fit_trajectory <- function(
    y, layer, gene_symbol, feature_id, rank_label,
    dat, grid_dat, linear_terms, seed_offset) {
  model_dat <- copy(dat)
  model_dat[, y := as.numeric(y)]

  fixed_null <- lm(build_formula(linear_terms), data = model_dat)
  fixed_fit <- lm(
    build_formula(
      linear_terms,
      sprintf("ns(age_model, df = %d)", fixed_df)
    ),
    data = model_dat
  )
  fixed_age_pvalue <- anova(fixed_null, fixed_fit)$`Pr(>F)`[2L]

  gam_fit <- fit_adaptive_gam(model_dat, linear_terms, seed_offset)

  fixed_prediction <- predict(fixed_fit, newdata = grid_dat, se.fit = TRUE)
  gam_prediction <- predict(
    gam_fit$fit,
    newdata = grid_dat,
    se.fit = TRUE,
    unconditional = TRUE
  )

  fixed_values <- as.numeric(fixed_prediction$fit)
  gam_values <- as.numeric(gam_prediction$fit)
  fixed_z <- as.numeric(scale(fixed_values))
  gam_z <- as.numeric(scale(gam_values))

  comparison <- data.table(
    gene_symbol = gene_symbol,
    layer = layer,
    feature_id = feature_id,
    decoupling_ranks = rank_label,
    n_samples = nrow(model_dat),
    fixed_age_df = fixed_df,
    fixed_age_pvalue = fixed_age_pvalue,
    gam_final_k = gam_fit$final_k,
    gam_edf = gam_fit$edf,
    gam_edf_minus_fixed_df = gam_fit$edf - fixed_df,
    gam_age_pvalue = gam_fit$smooth_pvalue,
    gam_k_index = gam_fit$k_index,
    gam_k_pvalue = gam_fit$k_pvalue,
    gam_k_check_pass = gam_fit$k_check_pass,
    shape_pearson = cor(fixed_values, gam_values, method = "pearson"),
    shape_spearman = cor(fixed_values, gam_values, method = "spearman"),
    standardized_curve_rmse = sqrt(mean((fixed_z - gam_z)^2)),
    endpoint_direction_concordant =
      sign(tail(fixed_values, 1L) - fixed_values[1L]) ==
      sign(tail(gam_values, 1L) - gam_values[1L]),
    peak_age_difference_years = abs(
      grid_dat$age_years[which.max(fixed_values)] -
        grid_dat$age_years[which.max(gam_values)]
    ),
    trough_age_difference_years = abs(
      grid_dat$age_years[which.min(fixed_values)] -
        grid_dat$age_years[which.min(gam_values)]
    )
  )
  comparison[, `:=`(
    gam_edf_materially_above_fixed = gam_edf >= edf_review_threshold,
    shape_review_flag = shape_pearson < shape_correlation_threshold |
      !endpoint_direction_concordant,
    localized_nonlinearity_review = gam_edf >= edf_review_threshold |
      shape_pearson < shape_correlation_threshold |
      !endpoint_direction_concordant
  )]

  z_value <- qnorm(0.975)
  curves <- rbindlist(list(
    data.table(
      gene_symbol = gene_symbol,
      layer = layer,
      feature_id = feature_id,
      decoupling_ranks = rank_label,
      grid_dat[, .(grid_index, age_model, age_years)],
      model = sprintf("fixed_ns_df%d", fixed_df),
      fitted = fixed_values,
      se = as.numeric(fixed_prediction$se.fit)
    ),
    data.table(
      gene_symbol = gene_symbol,
      layer = layer,
      feature_id = feature_id,
      decoupling_ranks = rank_label,
      grid_dat[, .(grid_index, age_model, age_years)],
      model = "adaptive_gam",
      fitted = gam_values,
      se = as.numeric(gam_prediction$se.fit)
    )
  ))
  curves[, `:=`(
    lower_95 = fitted - z_value * se,
    upper_95 = fitted + z_value * se
  )]
  curves[, fitted_z := as.numeric(scale(fitted)), by = model]

  list(comparison = comparison, curves = curves)
}

# --- load frozen candidates and matched data ---

ranked_pairs <- fread(file.path(mod, "03_decoupling_ranked_pairs.tsv"))[order(decoupling_rank)]
apa <- readRDS(file.path(mod, "apa_matched.rds"))
rna <- readRDS(file.path(mod, "rna_matched.rds"))
protein <- readRDS(file.path(mod, "protein_matched.rds"))
sample_manifest <- fread(file.path(mod, "matched_sample_manifest.tsv"))
support <- fread(file.path(mod, "01_support_profile.tsv"))[is_supported == TRUE]

rna_gene_base <- sub("\\.[0-9]+$", "", rownames(rna))
age_model <- log10(pmax(sample_manifest$age_years + 1, 1e-6))
model_dat <- data.table(
  age_model = age_model,
  sex_binary = as.integer(sample_manifest$sex == "M"),
  rin_z = as.numeric(scale(sample_manifest$rna_rin))
)
grid_dat <- data.table(
  grid_index = support$grid_index,
  age_model = support$age_model,
  age_years = support$grid_age,
  sex_binary = 0,
  rin_z = 0
)

gene_rank_labels <- ranked_pairs[, .(
  decoupling_ranks = paste(decoupling_rank, collapse = ";")
), by = .(gene_symbol, gencode_gene_base, protein_accession)]

# --- compare fixed and adaptive trajectories ---

trajectory_results <- list()
curve_results <- list()
result_index <- 0L

for (index in seq_len(nrow(ranked_pairs))) {
  pair <- ranked_pairs[index]
  result_index <- result_index + 1L
  fitted <- fit_trajectory(
    apa[match(pair$site_id, rownames(apa)), ],
    "APA", pair$gene_symbol, pair$site_id,
    as.character(pair$decoupling_rank),
    model_dat, grid_dat, c("sex_binary", "rin_z"), result_index
  )
  trajectory_results[[result_index]] <- fitted$comparison
  curve_results[[result_index]] <- fitted$curves
}

for (index in seq_len(nrow(gene_rank_labels))) {
  gene <- gene_rank_labels[index]
  result_index <- result_index + 1L
  fitted <- fit_trajectory(
    rna[match(gene$gencode_gene_base, rna_gene_base), ],
    "RNA", gene$gene_symbol, gene$gencode_gene_base,
    gene$decoupling_ranks,
    model_dat, grid_dat, c("sex_binary", "rin_z"), result_index
  )
  trajectory_results[[result_index]] <- fitted$comparison
  curve_results[[result_index]] <- fitted$curves

  result_index <- result_index + 1L
  fitted <- fit_trajectory(
    protein[match(gene$protein_accession, rownames(protein)), ],
    "protein", gene$gene_symbol, gene$protein_accession,
    gene$decoupling_ranks,
    model_dat, grid_dat, character(), result_index
  )
  trajectory_results[[result_index]] <- fitted$comparison
  curve_results[[result_index]] <- fitted$curves
}

trajectory_comparison <- rbindlist(trajectory_results)
trajectory_comparison[, `:=`(
  fixed_age_q_scoped = p.adjust(fixed_age_pvalue, method = "BH"),
  gam_age_q_scoped = p.adjust(gam_age_pvalue, method = "BH")
), by = layer]
trajectory_comparison[, `:=`(
  fixed_age_significant_scoped = fixed_age_q_scoped < 0.05,
  gam_age_significant_scoped = gam_age_q_scoped < 0.05
)]
trajectory_comparison[, age_significance_concordant :=
  fixed_age_significant_scoped == gam_age_significant_scoped]
setorder(trajectory_comparison, gene_symbol, layer, feature_id)

curve_data <- rbindlist(curve_results)
setorder(curve_data, gene_symbol, layer, feature_id, model, grid_index)

# --- compare fixed and GAM-adjusted pair coefficients ---

coefficient_results <- vector("list", nrow(ranked_pairs) * 3L)
coefficient_index <- 0L

for (index in seq_len(nrow(ranked_pairs))) {
  pair <- ranked_pairs[index]
  pair_dat <- copy(model_dat)
  pair_dat[, `:=`(
    APA = as.numeric(apa[match(pair$site_id, rownames(apa)), ]),
    RNA = as.numeric(rna[match(pair$gencode_gene_base, rna_gene_base), ]),
    y = as.numeric(protein[match(pair$protein_accession, rownames(protein)), ])
  )]

  fixed_m0 <- lm(
    build_formula("APA", sprintf("ns(age_model, df = %d)", fixed_df)),
    data = pair_dat
  )
  fixed_m1 <- lm(
    build_formula(c("APA", "RNA"), sprintf("ns(age_model, df = %d)", fixed_df)),
    data = pair_dat
  )
  gam_m0 <- fit_adaptive_gam(pair_dat, "APA", 1000L + index)
  gam_m1 <- fit_adaptive_gam(pair_dat, c("APA", "RNA"), 2000L + index)

  test_definitions <- list(
    list(model = "M0", term = "APA", fixed = fixed_m0, gam = gam_m0,
         global_q = pair$q_age_adjusted),
    list(model = "M1", term = "APA", fixed = fixed_m1, gam = gam_m1,
         global_q = pair$q_rna_adjusted),
    list(model = "M1", term = "RNA", fixed = fixed_m1, gam = gam_m1,
         global_q = pair$q_rna_independent)
  )

  for (definition in test_definitions) {
    coefficient_index <- coefficient_index + 1L
    fixed_coefficient <- extract_coefficient(definition$fixed, definition$term)
    gam_coefficient <- extract_coefficient(
      definition$gam$fit, definition$term
    )
    coefficient_results[[coefficient_index]] <- data.table(
      decoupling_rank = pair$decoupling_rank,
      gene_symbol = pair$gene_symbol,
      gencode_gene_base = pair$gencode_gene_base,
      protein_accession = pair$protein_accession,
      site_id = pair$site_id,
      hard_cut_candidate = pair$hard_cut_candidate,
      high_vif_source = pair$high_vif_source,
      model = definition$model,
      term = definition$term,
      fixed_beta = fixed_coefficient$beta,
      fixed_se = fixed_coefficient$se,
      fixed_t = fixed_coefficient$statistic,
      fixed_pvalue = fixed_coefficient$pvalue,
      fixed_global_qvalue = definition$global_q,
      fixed_global_significant = definition$global_q < 0.05,
      gam_beta = gam_coefficient$beta,
      gam_se = gam_coefficient$se,
      gam_t = gam_coefficient$statistic,
      gam_pvalue = gam_coefficient$pvalue,
      gam_age_edf = definition$gam$edf,
      gam_final_k = definition$gam$final_k,
      gam_k_index = definition$gam$k_index,
      gam_k_pvalue = definition$gam$k_pvalue,
      gam_k_check_pass = definition$gam$k_check_pass
    )
  }
}

coefficient_comparison <- rbindlist(coefficient_results)
coefficient_comparison[, `:=`(
  fixed_qvalue_scoped = p.adjust(fixed_pvalue, method = "BH"),
  gam_qvalue_scoped = p.adjust(gam_pvalue, method = "BH")
), by = .(model, term)]
coefficient_comparison[, `:=`(
  beta_direction_concordant = sign(fixed_beta) == sign(gam_beta),
  fixed_raw_significant = fixed_pvalue < 0.05,
  gam_raw_significant = gam_pvalue < 0.05,
  fixed_scoped_significant = fixed_qvalue_scoped < 0.05,
  gam_scoped_significant = gam_qvalue_scoped < 0.05
)]
coefficient_comparison[, `:=`(
  raw_significance_concordant =
    fixed_raw_significant == gam_raw_significant,
  scoped_significance_concordant =
    fixed_scoped_significant == gam_scoped_significant
)]
setorder(coefficient_comparison, decoupling_rank, model, term)

# --- summarize sensitivity results ---

sensitivity_summary <- rbindlist(list(
  data.table(
    section = "analysis_scope",
    metric = c(
      "ranked_pairs", "unique_genes", "trajectory_fits",
      "curve_rows", "matched_samples"
    ),
    value = as.character(c(
      nrow(ranked_pairs), uniqueN(ranked_pairs$gencode_gene_base),
      nrow(trajectory_comparison), nrow(curve_data), nrow(model_dat)
    )),
    detail = c(
      sprintf("Current %d-pair beta-gate ranking", nrow(ranked_pairs)),
      "Genes represented in the frozen ranking",
      sprintf("%d APA sites plus one RNA and one protein trajectory per gene",
              nrow(ranked_pairs)),
      "Fixed and adaptive predictions over supported age grid",
      "Matched donors used for every model"
    )
  ),
  data.table(
    section = "trajectory_robustness",
    metric = c(
      "edf_at_least_4.5", "shape_review_flags",
      "age_significance_disagreements", "failed_k_checks"
    ),
    value = as.character(c(
      sum(trajectory_comparison$gam_edf_materially_above_fixed),
      sum(trajectory_comparison$shape_review_flag),
      sum(!trajectory_comparison$age_significance_concordant),
      sum(!trajectory_comparison$gam_k_check_pass)
    )),
    detail = c(
      sprintf("Adaptive GAM EDF materially above fixed df=%d", fixed_df),
      "Pearson shape correlation below 0.90 or endpoint direction mismatch",
      "Scoped BH age-effect call differs between fixed and GAM trajectory",
      "Final selected k did not pass the prespecified k.check gate"
    )
  ),
  coefficient_comparison[, .(
    section = "pair_coefficient_robustness",
    metric = paste0(
      tolower(term), "_", tolower(model), "_scoped_fdr_disagreements"
    ),
    value = as.character(sum(!scoped_significance_concordant)),
    detail = paste0(
      "Fixed df=", fixed_df, " versus adaptive GAM; BH within ",
      nrow(ranked_pairs), " pairs for ", term,
      " in ", model
    )
  ), by = .(model, term)][, c("model", "term") := NULL],
  data.table(
    section = "pair_coefficient_robustness",
    metric = "failed_pair_model_k_checks",
    value = as.character(uniqueN(
      coefficient_comparison[gam_k_check_pass == FALSE, .(site_id, model)]
    )),
    detail = "Unique M0/M1 models; RNA and APA rows from the same M1 fit count once"
  )
), use.names = TRUE)

# --- write GAM sensitivity outputs ---

fwrite(trajectory_comparison, file.path(mod, "04_gam_trajectory_comparison_v2.tsv"), sep = "\t", na = "NA")
fwrite(coefficient_comparison, file.path(mod, "04_gam_pair_coefficient_comparison_v2.tsv"), sep = "\t", na = "NA")
fwrite(curve_data, file.path(mod, "04_gam_curve_data_v2.tsv"), sep = "\t", na = "NA")
fwrite(sensitivity_summary, file.path(mod, "04_gam_sensitivity_summary_v2.tsv"), sep = "\t", na = "NA")
