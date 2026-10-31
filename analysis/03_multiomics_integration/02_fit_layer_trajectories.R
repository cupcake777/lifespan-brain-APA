# Fixed-df natural-spline trajectories for the matched APA, RNA and protein layers.
# Degrees of freedom come from the upstream dynamic APA events; RNA and protein use the same design.

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

# --- config ---

if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
mod <- file.path(Sys.getenv("APA_ROOT"), "data/03_omics")

rna_equivalence_margin <- 0.20
rna_margin_sensitivity <- c(0.10, 0.30)
protein_absolute_margin <- 0.25
protein_age_permutations <- 1000L
random_seed <- 20260910L

# --- helpers ---

fit_shared_trajectory <- function(y, x_full, x_null, x_grid, confidence) {
  full_fit <- lm.fit(x_full, y)
  null_fit <- lm.fit(x_null, y)

  full_rss <- colSums(full_fit$residuals^2)
  null_rss <- colSums(null_fit$residuals^2)
  residual_df <- nrow(x_full) - full_fit$rank
  age_df <- full_fit$rank - null_fit$rank
  residual_sd <- sqrt(full_rss / residual_df)

  f_stat <- pmax((null_rss - full_rss) / age_df, 0) /
    (full_rss / residual_df)
  age_pvalue <- pf(f_stat, age_df, residual_df, lower.tail = FALSE)

  fitted_grid <- x_grid %*% full_fit$coefficients
  max_index <- max.col(t(fitted_grid), ties.method = "first")
  min_index <- max.col(t(-fitted_grid), ties.method = "first")
  feature_index <- seq_len(ncol(y))
  amplitude <- fitted_grid[cbind(max_index, feature_index)] -
    fitted_grid[cbind(min_index, feature_index)]
  endpoint_delta <- fitted_grid[nrow(fitted_grid), ] - fitted_grid[1L, ]

  contrast <- x_grid[max_index, , drop = FALSE] -
    x_grid[min_index, , drop = FALSE]
  xtx_inverse <- chol2inv(qr.R(qr(x_full)))
  contrast_variance <- rowSums((contrast %*% xtx_inverse) * contrast)
  amplitude_se <- residual_sd * sqrt(pmax(contrast_variance, 0))
  z_value <- qnorm(1 - (1 - confidence) / 2)

  data.table(
    age_pvalue = age_pvalue,
    residual_sd = residual_sd,
    endpoint_delta = endpoint_delta,
    fitted_amplitude = amplitude,
    amplitude_se = amplitude_se,
    amplitude_ci_low = pmax(amplitude - z_value * amplitude_se, 0),
    amplitude_ci_high = amplitude + z_value * amplitude_se
  )
}

# --- fixed spline df from upstream dynamic genes ---

upstream_gam_path <- file.path(Sys.getenv("APA_ROOT"), "data/02_traj/02_sites.tsv")
upstream_gam <- fread(
  upstream_gam_path,
  select = c("gene_id", "smooth_edf", "directional_pass")
)[directional_pass == TRUE & is.finite(smooth_edf)]

gene_edf <- upstream_gam[, .(gene_median_edf = median(smooth_edf)), by = gene_id]
fixed_df <- as.integer(round(median(gene_edf$gene_median_edf)))
df_values <- sort(unique(c(fixed_df - 1L, fixed_df, fixed_df + 1L)))

df_calibration <- data.table(
  source_file = upstream_gam_path,
  dynamic_definition = "directional_pass == TRUE",
  gene_summary = "median event smooth_edf within gene",
  n_dynamic_events = nrow(upstream_gam),
  n_dynamic_genes = nrow(gene_edf),
  gene_edf_min = min(gene_edf$gene_median_edf),
  gene_edf_q25 = quantile(gene_edf$gene_median_edf, 0.25),
  gene_edf_median = median(gene_edf$gene_median_edf),
  gene_edf_q75 = quantile(gene_edf$gene_median_edf, 0.75),
  gene_edf_max = max(gene_edf$gene_median_edf),
  fixed_df = fixed_df,
  sensitivity_df_low = fixed_df - 1L,
  sensitivity_df_high = fixed_df + 1L
)

# --- load matched data and primary analysis universe ---

apa <- readRDS(file.path(mod, "apa_matched.rds"))
rna <- readRDS(file.path(mod, "rna_matched.rds"))
protein <- readRDS(file.path(mod, "protein_matched.rds"))
sample_manifest <- fread(file.path(mod, "matched_sample_manifest.tsv"))
support <- fread(file.path(mod, "01_support_profile.tsv"))[is_supported == TRUE]
universe <- fread(file.path(mod, "01_analysis_universe.tsv"))[primary_universe == TRUE]

site_map <- unique(universe[, .(
  site_id, gencode_gene_base, gene_symbol, protein_accession,
  cluster_id, discovery_delta_logit
)], by = "site_id")
gene_map <- unique(universe[, .(
  gencode_gene_base, gene_symbol, protein_accession
)], by = "gencode_gene_base")

rna_gene_base <- sub("\\.[0-9]+$", "", rownames(rna))
apa_index <- match(site_map$site_id, rownames(apa))
rna_index <- match(gene_map$gencode_gene_base, rna_gene_base)
protein_index <- match(gene_map$protein_accession, rownames(protein))
stopifnot(
  !anyNA(apa_index), !anyNA(rna_index), !anyNA(protein_index),
  identical(colnames(apa), sample_manifest$sample_id),
  identical(colnames(rna), sample_manifest$sample_id),
  identical(colnames(protein), sample_manifest$sample_id)
)

apa_y <- t(apa[apa_index, , drop = FALSE])
rna_y <- t(rna[rna_index, , drop = FALSE])
protein_y <- t(protein[protein_index, , drop = FALSE])

age_model <- log10(pmax(sample_manifest$age_years + 1, 1e-6))
sex_binary <- as.integer(sample_manifest$sex == "M")
rin_z <- as.numeric(scale(sample_manifest$rna_rin))
supported_age <- support$age_model

# --- fit primary and df sensitivity trajectories ---

apa_results <- vector("list", length(df_values))
rna_results <- vector("list", length(df_values))
protein_results <- vector("list", length(df_values))

for (df_value in df_values) {
  age_basis <- ns(age_model, df = df_value)
  grid_basis <- predict(age_basis, supported_age)
  covariate_design <- cbind(
    "(Intercept)" = 1,
    sex_binary = sex_binary,
    rin_z = rin_z
  )
  covariate_grid <- cbind(
    "(Intercept)" = rep(1, length(supported_age)),
    sex_binary = rep(0, length(supported_age)),
    rin_z = rep(0, length(supported_age))
  )
  layer_design <- cbind(covariate_design, age_basis)
  layer_grid <- cbind(covariate_grid, grid_basis)
  protein_null <- matrix(1, nrow = length(age_model), ncol = 1L,
                         dimnames = list(NULL, "(Intercept)"))
  protein_design <- cbind("(Intercept)" = 1, age_basis)
  protein_grid <- cbind("(Intercept)" = 1, grid_basis)

  apa_fit <- fit_shared_trajectory(
    apa_y, layer_design, covariate_design, layer_grid, confidence = 0.95
  )
  apa_fit[, age_qvalue := p.adjust(age_pvalue, method = "BH")]
  apa_fit[, direction_concordant :=
    sign(endpoint_delta) == sign(site_map$discovery_delta_logit)]
  apa_fit[, `:=`(
    apa_dynamic = age_qvalue < 0.05 & direction_concordant,
    analysis_df = df_value,
    analysis_role = ifelse(df_value == fixed_df, "primary", "df_sensitivity")
  )]
  apa_results[[which(df_values == df_value)]] <- cbind(site_map, apa_fit)

  rna_fit <- fit_shared_trajectory(
    rna_y, layer_design, covariate_design, layer_grid, confidence = 0.90
  )
  rna_fit[, age_qvalue := p.adjust(age_pvalue, method = "BH")]
  rna_fit[, `:=`(
    standardized_amplitude = fitted_amplitude / residual_sd,
    standardized_amplitude_se = amplitude_se / residual_sd,
    standardized_ci_low = amplitude_ci_low / residual_sd,
    standardized_ci_high = amplitude_ci_high / residual_sd
  )]
  rna_fit[, standardized_ci_width := standardized_ci_high - standardized_ci_low]
  all_margins <- c(rna_equivalence_margin, rna_margin_sensitivity)
  for (margin in all_margins) {
    margin_label <- sprintf("%02d", round(100 * margin))
    p_column <- paste0("equivalence_p_margin_", margin_label)
    q_column <- paste0("equivalence_q_margin_", margin_label)
    stable_column <- paste0("rna_stable_margin_", margin_label)
    rna_fit[, (p_column) := pnorm(
      (standardized_amplitude - margin) / standardized_amplitude_se
    )]
    rna_fit[, (q_column) := p.adjust(get(p_column), method = "BH")]
    rna_fit[, (stable_column) :=
      get(q_column) < 0.05 & standardized_ci_high < margin]
  }
  rna_fit[, `:=`(
    rna_equivalence_margin = rna_equivalence_margin,
    rna_stable = get(sprintf(
      "rna_stable_margin_%02d", round(100 * rna_equivalence_margin)
    )),
    analysis_df = df_value,
    analysis_role = ifelse(df_value == fixed_df, "primary", "df_sensitivity")
  )]
  rna_results[[which(df_values == df_value)]] <- cbind(gene_map, rna_fit)

  protein_fit <- fit_shared_trajectory(
    protein_y, protein_design, protein_null, protein_grid, confidence = 0.95
  )

  # permutation null for the protein amplitude, seeded per df
  set.seed(random_seed + df_value)
  n_samples <- nrow(protein_y)
  n_features <- ncol(protein_y)
  grid_design <- cbind("(Intercept)" = 1, grid_basis)
  null_amplitude <- matrix(NA_real_, nrow = n_features, ncol = protein_age_permutations)
  for (iteration in seq_len(protein_age_permutations)) {
    permuted_design <- cbind(
      "(Intercept)" = 1,
      age_basis[sample.int(n_samples), , drop = FALSE]
    )
    coefficients <- qr.coef(qr(permuted_design), protein_y)
    fitted_grid <- grid_design %*% coefficients
    null_amplitude[, iteration] <-
      apply(fitted_grid, 2L, max) - apply(fitted_grid, 2L, min)
  }
  permutation_cutoff <- apply(null_amplitude, 1L, quantile, probs = 0.95, names = FALSE)

  protein_fit[, age_qvalue := p.adjust(age_pvalue, method = "BH")]
  protein_fit[, `:=`(
    permutation_amplitude_cutoff = permutation_cutoff,
    amplitude_cutoff = pmax(protein_absolute_margin, permutation_cutoff),
    analysis_df = df_value,
    analysis_role = ifelse(df_value == fixed_df, "primary", "df_sensitivity")
  )]
  protein_fit[, protein_dynamic :=
    age_qvalue < 0.05 & amplitude_ci_low > amplitude_cutoff]
  protein_results[[which(df_values == df_value)]] <- cbind(gene_map, protein_fit)
}

# --- write trajectory results ---

fwrite(df_calibration, file.path(mod, "02_df_calibration.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(apa_results), file.path(mod, "02_apa_trajectories.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(rna_results), file.path(mod, "02_rna_trajectories.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(protein_results), file.path(mod, "02_protein_trajectories.tsv"), sep = "\t", na = "NA")
