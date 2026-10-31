# Donor-level APA matrix and covariates for lifespan trajectories.
suppressPackageStartupMessages(library(data.table))

# --- paths ---

root <- Sys.getenv("APA_ROOT")
APA_IN <- file.path(root, "data/01_quant/apa_after_combat.tsv")
META_IN <- file.path(root, "data/04_qtl/CHB.all.meta")
OUT <- file.path(root, "data/02_traj")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- sample metadata ---

apa_samples <- names(fread(APA_IN, nrows = 0L, check.names = FALSE))[-1L]
meta <- fread(META_IN)
setnames(meta, 1L, "sample")
meta[, sample := as.character(sample)]
stopifnot(!anyDuplicated(meta$sample), !length(setdiff(apa_samples, meta$sample)))
samples <- meta[sample %chin% apa_samples, sample]
meta <- meta[match(samples, sample)]
meta[, individual_id := sub("_rep$", "", sample)]

# --- technical replicates ---

individual_qc <- meta[, .(
  n_profiles = .N,
  n_age_values = uniqueN(age),
  n_sex_values = uniqueN(tolower(trimws(sex))),
  n_rin_values = uniqueN(RIN)
), by = individual_id]
stopifnot(
  all(individual_qc$n_profiles %in% 1:2),
  all(individual_qc$n_age_values == 1L),
  all(individual_qc$n_sex_values == 1L)
)
# sex must be coded as Female/Male, F/M or 1/0
stopifnot(all(tolower(trimws(as.character(meta$sex))) %in% c("female", "f", "1", "male", "m", "0")))

# --- profile and donor covariates ---

meta[, profile_weight := 1 / .N, by = individual_id]
profile_cov <- meta[, .(
  sample,
  individual_id,
  profile_weight,
  age = as.numeric(age),
  sex_binary = as.integer(tolower(trimws(sex)) %in% c("male", "m", "0")),
  RIN_z = as.numeric(scale(RIN))
)]
donor_cov <- meta[, .(
  sample = first(individual_id),
  profile_weight = 1,
  age = first(as.numeric(age)),
  sex_binary = as.integer(first(tolower(trimws(sex))) %in% c("male", "m", "0")),
  RIN = mean(as.numeric(RIN))
), by = individual_id]
donor_cov[, RIN_z := as.numeric(scale(RIN))]
donor_cov[, RIN := NULL]
setcolorder(donor_cov, c("sample", "individual_id", "profile_weight", "age", "sex_binary", "RIN_z"))

# --- donor-level APA matrix ---

apa <- fread(APA_IN, check.names = FALSE, showProgress = FALSE)
setnames(apa, 1L, "site_id")
profile_matrix <- as.matrix(apa[, ..samples])
storage.mode(profile_matrix) <- "numeric"
stopifnot(all(is.finite(profile_matrix)))

donor_counts <- meta[, .N, by = individual_id]
donor_by_site <- rowsum(t(profile_matrix), meta$individual_id, reorder = FALSE)
donor_by_site <- donor_by_site / donor_counts$N
donor_matrix <- t(donor_by_site)
colnames(donor_matrix) <- donor_counts$individual_id
donor_apa <- cbind(data.table(site_id = as.character(apa$site_id)), as.data.table(donor_matrix))

# --- QC summary and outputs ---

qc_summary <- data.table(
  metric = c(
    "n_profiles", "n_individuals", "n_analysis_rows", "n_replicate_pairs",
    "n_age_discordant_individuals", "n_sex_discordant_individuals",
    "n_rin_discordant_pairs"
  ),
  value = c(
    nrow(meta), nrow(donor_cov), nrow(donor_cov),
    sum(individual_qc$n_profiles == 2L),
    sum(individual_qc$n_age_values > 1L),
    sum(individual_qc$n_sex_values > 1L),
    sum(individual_qc$n_profiles == 2L & individual_qc$n_rin_values > 1L)
  )
)
fwrite(donor_apa, file.path(OUT, "apa_after_combat_donor_mean.tsv"), sep = "\t")
fwrite(donor_cov, file.path(OUT, "00_cov.tsv"), sep = "\t")
fwrite(profile_cov, file.path(OUT, "00_cov_profiles.tsv"), sep = "\t")
fwrite(qc_summary, file.path(OUT, "01_individual_qc.tsv"), sep = "\t")
