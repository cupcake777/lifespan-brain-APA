# Integrate APA, RNA and protein for the sites in the analysis universe.
# Fits donor-level trajectories, decoupling pairs, robustness cases, donor bootstrap and repeated cross-validation.

suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(matrixStats)
})
setDTthreads(1L)

# --- config ---

if (!nzchar(Sys.getenv("APA_ROOT")) || !nzchar(Sys.getenv("APA_PROTEOMICS_HOME"))) {
  stop("Set APA_ROOT and APA_PROTEOMICS_HOME (source config/paths.sh)")
}
mod <- file.path(Sys.getenv("APA_ROOT"), "data/03_omics")
fixed_df <- fread(file.path(mod, "02_df_calibration.tsv"))$fixed_df[1]
bootstrap_iterations <- 1000L
cv_folds <- 5L
cv_repeats <- 20L
protein_age_permutations <- 1000L
bootstrap_workers <- 4L
random_seed <- 20260912L
top_n_values <- c(5L, 10L, 20L)

# --- helpers ---

fit_pairs <- function(samples, pair_index, protein_values = P, extra = NULL) {
  basis <- ns(age_model[samples], df = fixed_df)
  base <- cbind(1, basis, if (!is.null(extra)) extra[samples, , drop = FALSE])
  q <- qr(base)
  pair_genes <- unique(universe$gene_index[pair_index])
  gi <- match(universe$gene_index[pair_index], pair_genes)
  a <- qr.resid(q, A[samples, pair_index, drop = FALSE])
  r <- qr.resid(q, R[samples, pair_genes, drop = FALSE])
  p <- qr.resid(q, protein_values[samples, pair_genes, drop = FALSE])
  aa <- colSums(a*a)
  rr <- colSums(r*r)[gi]
  pp <- colSums(p*p)[gi]
  ar <- colSums(a*r[, gi, drop = FALSE])
  ap <- colSums(a*p[, gi, drop = FALSE])
  rp <- colSums(r*p)[gi]
  determinant <- aa*rr-ar^2
  valid_m0 <- aa > 1e-14*colSums(A[samples, pair_index, drop = FALSE]^2)
  valid_m1 <- valid_m0 & (rr-ar^2/aa) >
    1e-14*colSums(R[samples, pair_genes, drop = FALSE]^2)[gi]
  beta_m0 <- ap/aa
  beta_m1 <- (ap*rr-rp*ar)/determinant
  beta_rna <- (rp*aa-ap*ar)/determinant
  rss0 <- pmax(pp-ap^2/aa, 0)
  rss1 <- pmax(pp-beta_m1*ap-beta_rna*rp, 0)
  df0 <- length(samples)-q$rank-1L
  df1 <- length(samples)-q$rank-2L
  se0 <- sqrt(rss0/df0/aa)
  se1 <- sqrt(rss1/df1*rr/determinant)
  se_rna <- sqrt(rss1/df1*aa/determinant)
  out <- data.table(
    beta_age_adjusted = beta_m0, beta_rna_adjusted = beta_m1,
    beta_rna_independent = beta_rna,
    p_age_adjusted = 2*pt(abs(beta_m0/se0), df0, lower.tail = FALSE),
    p_rna_adjusted = 2*pt(abs(beta_m1/se1), df1, lower.tail = FALSE),
    p_rna_independent = 2*pt(abs(beta_rna/se_rna), df1, lower.tail = FALSE),
    se_rna_adjusted = se1)
  out[which(!valid_m0), c("beta_age_adjusted", "p_age_adjusted") := NA_real_]
  out[which(!valid_m1), c("beta_rna_adjusted", "beta_rna_independent",
    "p_rna_adjusted", "p_rna_independent", "se_rna_adjusted") := NA_real_]
  out[, `:=`(q_age_adjusted = p.adjust(p_age_adjusted, "BH"),
             q_rna_adjusted = p.adjust(p_rna_adjusted, "BH"),
             q_rna_independent = p.adjust(p_rna_independent, "BH"))]
  out[, beta_both_significant := !is.na(q_age_adjusted) & !is.na(q_rna_adjusted) &
        q_age_adjusted < .05 & q_rna_adjusted < .05]
  out
}

trajectory_fit <- function(y, x_null, basis, grid_basis, confidence = .95) {
  null_qr <- qr(x_null)
  x_null <- x_null[, null_qr$pivot[seq_len(null_qr$rank)], drop = FALSE]
  x <- cbind(x_null, basis)
  x_grid <- cbind(matrix(0, nrow(grid_basis), ncol(x_null)), grid_basis)
  x_grid[, 1] <- 1
  full <- lm.fit(x, y)
  null <- lm.fit(x_null, y)
  stopifnot(full$rank == ncol(x))
  rss <- colSums(full$residuals^2)
  rss_null <- colSums(null$residuals^2)
  rdf <- nrow(x)-full$rank
  residual_sd <- sqrt(rss/rdf)
  fitted <- x_grid %*% full$coefficients
  maxima <- max.col(t(fitted), ties.method = "first")
  minima <- max.col(t(-fitted), ties.method = "first")
  amplitude <- fitted[cbind(maxima, seq_len(ncol(y)))]-
    fitted[cbind(minima, seq_len(ncol(y)))]
  contrast <- x_grid[maxima, , drop = FALSE]-x_grid[minima, , drop = FALSE]
  inverse <- chol2inv(qr.R(full$qr))
  amplitude_se <- residual_sd*sqrt(pmax(rowSums((contrast %*% inverse)*contrast), 0))
  age_p <- pf(pmax((rss_null-rss)/fixed_df, 0)/(rss/rdf), fixed_df, rdf,
              lower.tail = FALSE)
  data.table(age_q = p.adjust(age_p, "BH"), residual_sd = residual_sd,
    amplitude = amplitude, standardized_amplitude = amplitude/residual_sd,
    amplitude_low = pmax(amplitude-qnorm((1+confidence)/2)*amplitude_se, 0),
    endpoint_delta = fitted[nrow(fitted), ]-fitted[1, ])
}

fit_workflow <- function(samples, pair_index, protein_values = P, extra = NULL,
                         batch = NULL, cutoff = NULL, seed = random_seed) {
  basis <- ns(age_model[samples], df = fixed_df)
  unique_age <- age_model[unique(samples)]
  age_range <- diff(range(age_model[samples]))
  supported <- vapply(support$age_model, function(z)
    sum(abs(unique_age-z) <= .08*age_range) >= 5 &&
      sum(abs(unique_age-z) <= .04*age_range) >= 3,
    logical(1))
  grid_basis <- predict(basis, support$age_model[supported])
  layer_null <- cbind(1, sex_binary[samples], rin_z[samples],
    if (!is.null(batch)) batch[samples, , drop = FALSE])
  protein_null <- cbind(rep(1, length(samples)),
    if (!is.null(extra)) extra[samples, , drop = FALSE])
  gene_index <- unique(universe$gene_index[pair_index])
  pair_gene <- match(universe$gene_index[pair_index], gene_index)
  at <- trajectory_fit(A[samples, pair_index, drop = FALSE], layer_null, basis, grid_basis)
  rt <- trajectory_fit(R[samples, gene_index, drop = FALSE], layer_null, basis, grid_basis, .90)
  pt <- trajectory_fit(protein_values[samples, gene_index, drop = FALSE], protein_null, basis, grid_basis)
  if (is.null(cutoff)) {
    # permutation cutoff for the protein amplitude, seeded per call
    y <- protein_values[samples, gene_index, drop = FALSE]
    set.seed(seed)
    null_qr <- qr(protein_null)
    x_null <- protein_null[, null_qr$pivot[seq_len(null_qr$rank)], drop = FALSE]
    x <- cbind(x_null, basis)
    grid <- cbind(matrix(0, nrow(grid_basis), ncol(x_null)), grid_basis)
    grid[, 1] <- 1
    op <- qr.coef(qr(x), diag(nrow(y)))
    null <- lm.fit(x_null, y)
    amplitudes <- matrix(NA_real_, ncol(y), protein_age_permutations)
    for (iteration in seq_len(protein_age_permutations)) {
      permuted <- null$fitted.values+null$residuals[sample.int(nrow(y)), , drop = FALSE]
      fitted <- grid %*% (op %*% permuted)
      amplitudes[, iteration] <- colMaxs(fitted)-colMins(fitted)
    }
    cutoff <- pmax(.25, rowQuantiles(amplitudes, probs = .95, drop = TRUE))
  } else cutoff <- cutoff[gene_index]
  out <- cbind(universe[pair_index, .(site_id, gene_symbol, gencode_gene_base,
                                    protein_accession, cluster_id)],
               fit_pairs(samples, pair_index, protein_values, extra))
  out[, `:=`(
    apa_dynamic = at$age_q < .05 & sign(at$endpoint_delta) == sign(universe$discovery_delta_logit[pair_index]),
    protein_dynamic = (pt$age_q < .05 & pt$amplitude_low > cutoff)[pair_gene],
    rna_standardized_amplitude = rt$standardized_amplitude[pair_gene],
    protein_standardized_amplitude = pt$standardized_amplitude[pair_gene],
    rna_percentile = ((frank(rt$standardized_amplitude, ties.method = "average")-1)/(nrow(rt)-1))[pair_gene],
    decoupling_rank = NA_integer_)]
  out[, decoupling_ratio := rna_standardized_amplitude/protein_standardized_amplitude]
  selected <- which(out$beta_both_significant)
  ordered <- selected[order(out$decoupling_ratio[selected], out$site_id[selected])]
  out$decoupling_rank[ordered] <- seq_along(ordered)
  out[, hard_cut_candidate := beta_both_significant & apa_dynamic & protein_dynamic &
        !is.na(decoupling_rank) & decoupling_rank <= 10L]
  out
}

# --- load manifest, fixed universe and matched matrices ---

manifest <- fread(file.path(mod, "matched_sample_manifest.tsv"))
universe <- fread(file.path(mod, "01_analysis_universe.tsv"))[
  sensitivity_universe == TRUE]
setorder(universe, site_id)
genes <- unique(universe[, .(gencode_gene_base, gene_symbol, protein_accession)],
                by = "gencode_gene_base")
universe[, gene_index := match(gencode_gene_base, genes$gencode_gene_base)]
primary_index <- which(universe$primary_universe)
ranked <- fread(file.path(mod, "03_decoupling_ranked_pairs.tsv"))
setorder(ranked, decoupling_rank)
primary_ranked_index <- match(ranked$site_id, universe$site_id[primary_index])
apa <- readRDS(file.path(mod, "apa_matched.rds"))
rna <- readRDS(file.path(mod, "rna_matched.rds"))
protein <- readRDS(file.path(mod, "protein_matched.rds"))
stopifnot(identical(colnames(apa), manifest$sample_id),
          identical(colnames(rna), manifest$sample_id),
          identical(colnames(protein), manifest$sample_id))
A <- t(apa[universe$site_id, , drop = FALSE])
R <- t(rna[match(genes$gencode_gene_base, sub("\\.[0-9]+$", "", rownames(rna))), , drop = FALSE])
P <- t(protein[genes$protein_accession, , drop = FALSE])
rm(apa, rna, protein)
protein_objects <- new.env()
load(file.path(mod, "protein_all_modifiedresults_covadj_adj_final.RData"),
     envir = protein_objects)
P_imputed <- t(protein_objects$data_rmOut_imp$log2_value[
  genes$protein_accession, manifest$protein_sample_id, drop = FALSE])

# --- source metadata and pre-imputation observation mask ---

source_dir <- Sys.getenv("APA_PROTEOMICS_HOME")
if (!nzchar(source_dir)) stop("Set APA_PROTEOMICS_HOME to the authorized proteomics source directory")
source_files <- c("CHB_meta_1959_250827.RData", "CHB_meta_1960_251108.RData")
provenance <- rbindlist(lapply(source_files, function(filename) {
  metadata <- new.env()
  load(file.path(mod, filename), envir = metadata)
  x <- as.data.table(metadata$meta_all)
  x[, match_id := gsub("-", "_", ID)]
  selected <- x[match(manifest$protein_sample_id, match_id)]
  stopifnot(!anyNA(selected$ID), !anyDuplicated(selected$match_id))
  data.table(sample_id = manifest$sample_id, protein_sample_id = manifest$protein_sample_id,
    source_file = paste0(source_dir, "/", filename), source_id = selected$ID,
    bank = selected$bank, source_age = selected$age, ms_analysis_name = selected$MS_anaysisName,
    ms_sample_name = selected$MS_SampleName, ptm_project_id = selected$PTMprojectID,
    pfc_in_sample_name = grepl("PFC", selected$ID, ignore.case = TRUE) |
      grepl("PFC", selected$MS_anaysisName, ignore.case = TRUE) |
      grepl("PFC", selected$MS_SampleName, ignore.case = TRUE),
    brain_region = "PFC", brain_region_source = "investigator_confirmation_2026-09-12",
    tissue_metadata_status = "no_explicit_tissue_or_brain_region_field")
}))
fwrite(provenance, file.path(mod, "protein_source_provenance.tsv"), sep = "\t")
pre <- new.env()
load(file.path(mod, "protein_all_modified_results_rmOut.RData"), envir = pre)
feature <- fread(file.path(mod, "protein_feature_annotation.tsv"))
stopifnot(all(feature$Protein.accession %in% rownames(pre$data_rmOut$intensity)),
          all(manifest$protein_sample_id %in% colnames(pre$data_rmOut$intensity)))
pre_intensity <- as.matrix(pre$data_rmOut$intensity[feature$Protein.accession,
                                                 manifest$protein_sample_id, drop = FALSE])
observed_mask <- !is.na(pre_intensity)
colnames(observed_mask) <- manifest$sample_id
saveRDS(observed_mask, file.path(mod, "protein_observed_mask.rds"))
imputation <- feature[, .(protein_accession = Protein.accession, gene_symbol = Gene.name)]
imputation[, `:=`(n_donors = nrow(manifest), n_observed = rowSums(observed_mask),
  imputed_fraction = rowMeans(!observed_mask),
  prenatal_imputed_fraction = rowMeans(!observed_mask[, manifest$age_years < 0, drop = FALSE]),
  postnatal_imputed_fraction = rowMeans(!observed_mask[, manifest$age_years >= 0, drop = FALSE]))]
fwrite(imputation, file.path(mod, "04_protein_imputation.tsv"), sep = "\t")
manifest_names <- c(source_files, "protein_all_modified_results_rmOut.RData",
                    "protein_all_modifiedresults_covadj_adj_final.RData")
manifest_paths <- file.path(mod, manifest_names)
source_manifest <- data.table(artifact = c("source_metadata_20250827", "source_metadata_20251108",
    "source_preimputation_protein", "source_final_protein"), path = manifest_paths,
  source = paste0(source_dir, "/", manifest_names),
  n_features = c(18L, 19L, nrow(pre$data_rmOut$intensity), nrow(feature)),
  n_samples = c(1960L, 1960L, ncol(pre$data_rmOut$intensity), ncol(protein_objects$quant_final)),
  notes = c(rep("PFC confirmed by investigator on 2026-09-12; source metadata has no explicit region field", 2),
    "Observed intensities exactly match the final object's imputed intensity on retained cells",
    "MD5 matched the genetics source: e6723f96913e3ad004e631bec2a1a2a0"),
  md5 = unname(tools::md5sum(manifest_paths)))
input_manifest <- fread(file.path(mod, "input_manifest.tsv"))
fwrite(rbindlist(list(input_manifest[!artifact %chin% source_manifest$artifact], source_manifest)),
       file.path(mod, "input_manifest.tsv"), sep = "\t")
rm(pre, protein_objects, pre_intensity)

# --- covariates and fixed references ---

age_model <- log10(manifest$age_years + 1)
sex_binary <- as.integer(manifest$sex == "M")
rin_z <- as.numeric(scale(manifest$rna_rin))
observed_covariates <- model.matrix(~ sex_binary + rin_z +
  factor(manifest$protein_digest_batch) + factor(manifest$protein_disease_class))[, -1, drop = FALSE]
batch_covariates <- model.matrix(~ factor(manifest$rna_batch))[, -1, drop = FALSE]
support <- fread(file.path(mod, "01_support_profile.tsv"))
protein_reference <- fread(file.path(mod, "02_protein_trajectories.tsv"))[
  analysis_df == fixed_df]
fixed_cutoff <- protein_reference$amplitude_cutoff[
  match(genes$gencode_gene_base, protein_reference$gencode_gene_base)]
all_samples <- seq_len(nrow(manifest))

# --- primary coefficient family, checked against the archived run ---

primary <- fit_workflow(all_samples, primary_index, cutoff = fixed_cutoff)
cached <- fread(file.path(mod, "03_apa_incremental_protein_tests.tsv"))
comparison <- merge(primary, cached, by = "site_id", suffixes = c("", "_cached"))
expected_primary_pairs <- length(primary_index)
stopifnot(nrow(primary) == expected_primary_pairs,
  nrow(cached) == expected_primary_pairs,
  nrow(comparison) == expected_primary_pairs,
  !anyDuplicated(primary$site_id), !anyDuplicated(cached$site_id),
  setequal(primary$site_id, universe$site_id[primary_index]),
  setequal(cached$site_id, universe$site_id[primary_index]),
  identical(is.na(comparison$beta_rna_adjusted), is.na(comparison$beta_rna_adjusted_cached)),
  max(abs(comparison$beta_rna_adjusted-comparison$beta_rna_adjusted_cached), na.rm = TRUE) < 1e-8,
  max(abs(comparison$q_rna_adjusted-comparison$q_rna_adjusted_cached), na.rm = TRUE) < 1e-8,
  setequal(primary[beta_both_significant == TRUE, site_id], ranked$site_id),
  setequal(primary[hard_cut_candidate == TRUE, site_id], ranked[hard_cut_candidate == TRUE, site_id]))
fwrite(primary, file.path(mod, "04_robustness_primary.tsv.gz"), sep = "\t")

# --- fixed-universe and protein-processing sensitivities ---

cases <- list(
  exclude_PTB752 = list(samples = which(manifest$canonical_id != "PTB752"), index = primary_index),
  imputed_observed_covariates = list(samples = all_samples, index = primary_index,
                                    protein = P_imputed, extra = observed_covariates),
  all_unambiguous_proteins = list(samples = all_samples, index = seq_len(nrow(universe))),
  prenatal_continuous_age = list(samples = which(manifest$age_years < 0), index = primary_index),
  postnatal_continuous_age = list(samples = which(manifest$age_years >= 0), index = primary_index),
  post_combat_batch = list(samples = all_samples, index = primary_index,
                          extra = batch_covariates, batch = batch_covariates))
case_summary <- list()
for (name in names(cases)) {
  spec <- cases[[name]]
  elapsed <- proc.time()[3]
  out <- fit_workflow(spec$samples, spec$index,
    protein_values = if (is.null(spec$protein)) P else spec$protein,
    extra = spec$extra, batch = spec$batch, seed = random_seed)
  out[, sensitivity := name]
  original <- out[match(ranked$site_id, site_id)]
  selected_genes <- out[hard_cut_candidate == TRUE, unique(gencode_gene_base)]
  primary_genes <- ranked[hard_cut_candidate == TRUE, unique(gencode_gene_base)]
  case_summary[[name]] <- data.table(sensitivity = name, n_donors = length(spec$samples),
    n_pairs = nrow(out), n_genes = uniqueN(out$gencode_gene_base),
    beta_gate_pairs = sum(out$beta_both_significant), hard_cut_pairs = sum(out$hard_cut_candidate),
    hard_cut_genes = length(selected_genes), retained_reference_pairs = sum(original$beta_both_significant),
    reference_pairs_same_sign = sum(sign(original$beta_rna_adjusted) == sign(ranked$beta_rna_adjusted)),
    retained_reference_hard_cut_genes = length(intersect(selected_genes, primary_genes)),
    seconds = proc.time()[3]-elapsed)
  fwrite(out, file.path(mod, paste0("04_robustness_", name, ".tsv.gz")), sep = "\t")
}

# --- observed-only protein sensitivity, BH across the original pair family ---

elapsed <- proc.time()[3]
observed_tests <- rbindlist(lapply(split(primary_index, universe$gene_index[primary_index]), function(index) {
  accession <- universe$protein_accession[index[1]]
  samples <- which(observed_mask[accession, ])
  stopifnot(length(samples) > fixed_df+3L, length(unique(age_model[samples])) > fixed_df)
  fitted <- fit_pairs(samples, index)
  fitted[, `:=`(site_id = universe$site_id[index], n_observed = length(samples))]
  fitted
}))
observed_tests[, `:=`(q_age_adjusted = p.adjust(p_age_adjusted, "BH"),
  q_rna_adjusted = p.adjust(p_rna_adjusted, "BH"), q_rna_independent = p.adjust(p_rna_independent, "BH"))]
observed_tests[, beta_both_significant := !is.na(q_age_adjusted) & !is.na(q_rna_adjusted) &
  q_age_adjusted < .05 & q_rna_adjusted < .05]
observed_tests <- merge(primary[, .(site_id, gene_symbol, gencode_gene_base, protein_accession)],
                        observed_tests, by = "site_id")
observed_tests[, sensitivity := "observed_only_processed_protein"]
fwrite(observed_tests, file.path(mod, "04_robustness_observed_only.tsv.gz"), sep = "\t")
original <- observed_tests[match(ranked$site_id, site_id)]
case_summary[["observed_only_processed_protein"]] <- data.table(
  sensitivity = "observed_only_processed_protein", n_donors = NA_integer_,
  n_pairs = nrow(observed_tests), n_genes = uniqueN(observed_tests$gencode_gene_base),
  beta_gate_pairs = sum(observed_tests$beta_both_significant),
  hard_cut_pairs = NA_integer_, hard_cut_genes = NA_integer_,
  retained_reference_pairs = sum(original$beta_both_significant),
  reference_pairs_same_sign = sum(sign(original$beta_rna_adjusted) == sign(ranked$beta_rna_adjusted)),
  retained_reference_hard_cut_genes = NA_integer_, seconds = proc.time()[3]-elapsed,
  min_pair_n = min(observed_tests$n_observed), max_pair_n = max(observed_tests$n_observed))
fwrite(rbindlist(case_summary, fill = TRUE), file.path(mod, "04_robustness_summary.tsv"), sep = "\t")

# --- continuous ranking and heterogeneous multi-site genes ---

rank_sensitivity <- rbindlist(lapply(c("decoupling_ratio", "rna_percentile"), function(metric) {
  ordered <- primary[beta_both_significant == TRUE][order(get(metric), site_id)]
  rbindlist(lapply(top_n_values, function(n_top) {
    chosen <- head(ordered, n_top)[apa_dynamic == TRUE & protein_dynamic == TRUE]
    data.table(rank_metric = metric, top_n = n_top, candidate_pairs = nrow(chosen),
      candidate_genes = uniqueN(chosen$gencode_gene_base),
      genes = paste(sort(unique(chosen$gene_symbol)), collapse = ";"),
      reference_hard_cut_gene_retention = sum(unique(ranked[hard_cut_candidate == TRUE, gencode_gene_base]) %in%
                                              chosen$gencode_gene_base))
  }))
}))
fwrite(rank_sensitivity, file.path(mod, "04_ranking_sensitivity.tsv"), sep = "\t")
gene_types <- universe[primary_universe == TRUE,
  .(n_sites = .N, n_directions = uniqueN(cluster_id)), by = gencode_gene_base]
gene_types[, site_pattern := fifelse(n_sites == 1L, "single_site",
                              fifelse(n_directions == 1L, "concordant_multiple_sites", "opposite_direction_sites"))]
gene_hits <- primary[, .(beta_hit = any(beta_both_significant), candidate = any(hard_cut_candidate)),
                      by = gencode_gene_base]
gene_types <- merge(gene_types, gene_hits, by = "gencode_gene_base")
fwrite(gene_types[, .(n_genes = .N, n_beta_genes = sum(beta_hit), n_candidate_genes = sum(candidate)),
                  by = site_pattern], file.path(mod, "04_multisite_sensitivity.tsv"), sep = "\t")

# --- donor bootstrap with the full coefficient family and ranking refitted ---

bootstrap <- parallel::mclapply(seq_len(bootstrap_iterations), function(iteration) {
  set.seed(random_seed+iteration)
  samples <- sample.int(nrow(manifest), replace = TRUE)
  fitted <- fit_workflow(samples, primary_index, cutoff = fixed_cutoff)
  original <- fitted[primary_ranked_index]
  list(selected = which(fitted$beta_both_significant), hard = which(fitted$hard_cut_candidate),
       beta = original$beta_rna_adjusted, rank = original$decoupling_rank,
       summary = c(n_beta = sum(fitted$beta_both_significant), n_hard = sum(fitted$hard_cut_candidate),
                   n_genes = uniqueN(fitted[hard_cut_candidate == TRUE, gencode_gene_base])))
}, mc.cores = bootstrap_workers, mc.set.seed = FALSE)
stopifnot(all(vapply(bootstrap, is.list, logical(1))))
selection_frequency <- tabulate(unlist(lapply(bootstrap, `[[`, "selected")), nrow(primary))/bootstrap_iterations
hard_frequency <- tabulate(unlist(lapply(bootstrap, `[[`, "hard")), nrow(primary))/bootstrap_iterations
bootstrap_all <- primary[, .(site_id, gene_symbol, gencode_gene_base, beta_both_significant, hard_cut_candidate)]
bootstrap_all[, `:=`(beta_selection_frequency = selection_frequency,
                     hard_cut_frequency = hard_frequency, iterations = bootstrap_iterations)]
fwrite(bootstrap_all, file.path(mod, "04_donor_bootstrap_pairs.tsv"), sep = "\t")
beta_draws <- do.call(rbind, lapply(bootstrap, `[[`, "beta"))
rank_draws <- do.call(rbind, lapply(bootstrap, `[[`, "rank"))
bootstrap_ranked <- bootstrap_all[primary_ranked_index]
bootstrap_ranked[, `:=`(
  beta_median = colMedians(beta_draws),
  beta_low95 = colQuantiles(beta_draws, probs = .025, drop = TRUE),
  beta_high95 = colQuantiles(beta_draws, probs = .975, drop = TRUE),
  same_sign_fraction = colMeans(sign(beta_draws) == matrix(sign(ranked$beta_rna_adjusted),
                           nrow(beta_draws), nrow(ranked), byrow = TRUE)),
  median_rank_when_selected = apply(rank_draws, 2, median, na.rm = TRUE))]
fwrite(bootstrap_ranked, file.path(mod, "04_donor_bootstrap_ranked.tsv"), sep = "\t")
fwrite(as.data.table(do.call(rbind, lapply(bootstrap, `[[`, "summary"))),
       file.path(mod, "04_donor_bootstrap_counts.tsv"), sep = "\t")
rm(bootstrap)

# --- repeated donor cross-validation with training-fold coefficient selection ---

cv_rows <- vector("list", cv_repeats)
cv_selection_rows <- vector("list", cv_repeats)
for (repeat_id in seq_len(cv_repeats)) {
  set.seed(random_seed+10000L+repeat_id)
  fold <- integer(nrow(manifest))
  for (stage_value in unique(manifest$life_stage)) {
    ids <- which(manifest$life_stage == stage_value)
    fold[ids] <- sample(rep(seq_len(cv_folds), length.out = length(ids)))
  }
  sse_base <- sse_full <- sse_selected <- denominator <- numeric(nrow(primary))
  selected_folds <- numeric(nrow(primary))
  gi <- universe$gene_index[primary_index]
  for (fold_id in seq_len(cv_folds)) {
    training <- which(fold != fold_id)
    testing <- which(fold == fold_id)
    training_tests <- fit_pairs(training, primary_index)
    selected <- training_tests$beta_both_significant
    selected_folds <- selected_folds+selected
    basis <- ns(age_model[training], df = fixed_df)
    base_train <- cbind(1, basis)
    base_test <- cbind(1, predict(basis, age_model[testing]))
    base_qr <- qr(base_train)
    a_fit <- qr.coef(base_qr, A[training, primary_index, drop = FALSE])
    r_fit <- qr.coef(base_qr, R[training, , drop = FALSE])
    p_fit <- qr.coef(base_qr, P[training, , drop = FALSE])
    a <- A[training, primary_index, drop = FALSE]-base_train %*% a_fit
    r <- R[training, , drop = FALSE]-base_train %*% r_fit
    p <- P[training, , drop = FALSE]-base_train %*% p_fit
    aa <- colSums(a*a); rr <- colSums(r*r)[gi]
    ar <- colSums(a*r[, gi, drop = FALSE]); ap <- colSums(a*p[, gi, drop = FALSE])
    rp <- colSums(r*p)[gi]
    beta_a <- (ap*rr-rp*ar)/(aa*rr-ar^2)
    beta_r <- (rp*aa-ap*ar)/(aa*rr-ar^2)
    invalid <- !is.finite(training_tests$beta_rna_adjusted)
    beta_a[invalid] <- 0
    beta_r[invalid] <- (rp/rr)[invalid]
    a_test <- A[testing, primary_index, drop = FALSE]-base_test %*% a_fit
    r_test <- R[testing, , drop = FALSE]-base_test %*% r_fit
    p_base <- (base_test %*% p_fit)[, gi, drop = FALSE]
    prediction_base <- p_base+sweep(r_test[, gi, drop = FALSE], 2, rp/rr, `*`)
    prediction_full <- p_base+sweep(a_test, 2, beta_a, `*`)+
      sweep(r_test[, gi, drop = FALSE], 2, beta_r, `*`)
    observed <- P[testing, gi, drop = FALSE]
    error_base <- colSums((observed-prediction_base)^2)
    error_full <- colSums((observed-prediction_full)^2)
    sse_base <- sse_base+error_base
    sse_full <- sse_full+error_full
    sse_selected <- sse_selected+ifelse(selected, error_full, error_base)
    denominator <- denominator+colSums(sweep(observed, 2, colMeans(P[training, gi, drop = FALSE]), `-`)^2)
  }
  cv_rows[[repeat_id]] <- data.table(repeat_id = repeat_id, site_id = primary$site_id,
    gene_symbol = primary$gene_symbol, r2_rna_age = 1-sse_base/denominator,
    r2_apa_rna_age = 1-sse_full/denominator,
    delta_r2_fixed_candidate = (sse_base-sse_full)/denominator,
    delta_r2_training_selected = (sse_base-sse_selected)/denominator,
    training_selection_fraction = selected_folds/cv_folds)
  cv_selection_rows[[repeat_id]] <- data.table(repeat_id = repeat_id,
    n_predefined_pairs = nrow(primary), mean_delta_r2_training_selected = mean((sse_base-sse_selected)/denominator),
    n_ever_training_selected = sum(selected_folds > 0),
    mean_training_selected_fraction = mean(selected_folds/cv_folds))
}
cv <- rbindlist(cv_rows)
cv_summary <- cv[, .(
  repeats = .N, mean_r2_rna_age = mean(r2_rna_age), mean_r2_apa_rna_age = mean(r2_apa_rna_age),
  mean_delta_r2 = mean(delta_r2_fixed_candidate),
  split_low95 = quantile(delta_r2_fixed_candidate, .025),
  split_high95 = quantile(delta_r2_fixed_candidate, .975),
  positive_repeat_fraction = mean(delta_r2_fixed_candidate > 0),
  mean_delta_r2_training_selected = mean(delta_r2_training_selected),
  training_selection_fraction = mean(training_selection_fraction)), by = .(site_id, gene_symbol)]
cv_summary[, frozen_primary_candidate := site_id %chin% ranked$site_id]
fwrite(cv_summary, file.path(mod, "04_cv_all_pairs.tsv.gz"), sep = "\t")
fwrite(cv_summary[match(ranked$site_id, site_id)], file.path(mod, "04_cv_ranked_pairs.tsv"), sep = "\t")
fwrite(cv[site_id %chin% ranked$site_id], file.path(mod, "04_cv_ranked_repeats.tsv"), sep = "\t")
fwrite(rbindlist(cv_selection_rows), file.path(mod, "04_cv_training_selection.tsv"), sep = "\t")
