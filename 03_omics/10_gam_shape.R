#!/usr/bin/env Rscript
# Stage-stratified donor bootstrap of GAM-vs-linear trajectory shape (Spearman rho).

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(splines)
})
setDTthreads(1L)

# --- config ---
apa_root <- Sys.getenv("APA_ROOT")
omics_dir <- file.path(apa_root, "data", "03_omics")
apa_path <- file.path(omics_dir, "apa_matched.rds")
rna_path <- file.path(omics_dir, "rna_matched.rds")
protein_path <- file.path(omics_dir, "protein_matched.rds")
manifest_path <- file.path(omics_dir, "matched_sample_manifest.tsv")
support_path <- file.path(omics_dir, "01_support_profile.tsv")
comparison_path <- file.path(omics_dir, "04_gam_traj_comp.tsv")
df_calibration_path <- file.path(omics_dir, "02_df_calibration.tsv")
output_path <- file.path(omics_dir, "06_gam_shape_spearman_boot.tsv")

n_bootstrap <- 1000L
random_seed <- 20260923L
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
fixed_df <- fread(df_calibration_path)$fixed_df[1]

# --- helpers ---
build_formula <- function(covariates, smooth, k = NULL) {
  age_term <- if (smooth) {
    sprintf("s(age_model, k = %d, bs = 'tp')", k)
  } else {
    sprintf("ns(age_model, df = %d)", fixed_df)
  }
  as.formula(paste("y ~", paste(c(covariates, age_term), collapse = " + ")))
}

# --- loading data and trajectory definitions ---
apa <- readRDS(apa_path)
rna <- readRDS(rna_path)
protein <- readRDS(protein_path)
manifest <- fread(manifest_path)
support <- fread(support_path)[is_supported == TRUE]
comparison <- fread(comparison_path)

stopifnot(identical(colnames(apa), colnames(rna)),
          identical(colnames(apa), colnames(protein)))
manifest <- manifest[match(colnames(apa), sample_id)]
stopifnot(!anyNA(manifest$sample_id), nrow(comparison) > 0L,
          !anyDuplicated(comparison[, .(layer, feature_id)]),
          all(comparison$gam_final_k %in% c(20L, 30L, 40L, 60L, 80L, 100L, 120L)))

rna_gene_base <- sub("\\.[0-9]+$", "", rownames(rna))
feature_values <- lapply(seq_len(nrow(comparison)), function(i) {
  row <- comparison[i]
  if (row$layer == "APA") {
    as.numeric(apa[match(row$feature_id, rownames(apa)), ])
  } else if (row$layer == "RNA") {
    as.numeric(rna[match(row$feature_id, rna_gene_base), ])
  } else {
    as.numeric(protein[match(row$feature_id, rownames(protein)), ])
  }
})
stopifnot(all(vapply(feature_values, function(x) all(is.finite(x)), logical(1))))

model_data <- data.table(
  age_model = log10(manifest$age_years + 1),
  sex_binary = as.integer(manifest$sex == "M"),
  rin_z = as.numeric(scale(manifest$rna_rin))
)
grid_data <- data.table(
  age_model = support$age_model,
  sex_binary = 0,
  rin_z = 0
)
strata <- split(seq_len(nrow(manifest)), manifest$life_stage)

# --- stage-stratified donor bootstrap ---
bootstrap_one <- function(draw) {
  set.seed(random_seed + draw)
  donor_index <- unlist(lapply(strata, function(i) sample(i, length(i), TRUE)),
                        use.names = FALSE)
  dat <- model_data[donor_index]
  k_used <- pmin(comparison$gam_final_k, length(unique(dat$age_model)) - 1L)

  rho <- vapply(seq_len(nrow(comparison)), function(i) {
    dat[, y := feature_values[[i]][donor_index]]
    covariates <- if (comparison$layer[i] %in% c("APA", "RNA")) {
      c("sex_binary", "rin_z")
    } else {
      character()
    }
    tryCatch({
      fixed_fit <- lm(build_formula(covariates, FALSE), data = dat)
      gam_fit <- gam(
        build_formula(covariates, TRUE, k_used[i]),
        data = dat, method = "REML"
      )
      fixed_curve <- predict(fixed_fit, newdata = grid_data)
      gam_curve <- predict(gam_fit, newdata = grid_data)
      cor(fixed_curve, gam_curve, method = "spearman")
    }, error = function(e) NA_real_)
  }, numeric(1))
  list(rho = rho, k_used = k_used)
}

draws <- parallel::mclapply(
  seq_len(n_bootstrap), bootstrap_one,
  mc.cores = n_cores, mc.preschedule = TRUE
)
draw_matrix <- do.call(rbind, lapply(draws, `[[`, "rho"))
k_matrix <- do.call(rbind, lapply(draws, `[[`, "k_used"))

# --- summarize bootstrap uncertainty ---
bootstrap_summary <- rbindlist(lapply(seq_len(nrow(comparison)), function(i) {
  rho <- draw_matrix[, i]
  valid <- is.finite(rho)
  data.table(
    gene_symbol = comparison$gene_symbol[i],
    layer = comparison$layer[i],
    feature_id = comparison$feature_id[i],
    decoupling_ranks = comparison$decoupling_ranks[i],
    point_shape_spearman = comparison$shape_spearman[i],
    bootstrap_ci_lower_95 = quantile(rho[valid], 0.025, names = FALSE),
    bootstrap_ci_upper_95 = quantile(rho[valid], 0.975, names = FALSE),
    bootstrap_prob_positive = mean(rho[valid] > 0),
    bootstrap_prob_negative = mean(rho[valid] < 0),
    bootstrap_prob_gt_0_9 = mean(rho[valid] > 0.9),
    n_reduced_k_draws = sum(k_matrix[, i] < comparison$gam_final_k[i]),
    minimum_bootstrap_gam_k = min(k_matrix[, i]),
    n_successful_draws = sum(valid),
    n_requested_draws = n_bootstrap
  )
}))
setorder(bootstrap_summary, gene_symbol, layer, feature_id)
fwrite(bootstrap_summary, output_path, sep = "\t", na = "NA")
