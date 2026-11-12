#!/usr/bin/env Rscript
# Merge frozen-GWAS sLDSC, MAGMA results into the module 06 tables and compare with the previous run.
suppressMessages(library(data.table))

# --- config ---

root <- Sys.getenv("APA_ROOT")
run_dir <- file.path(root, "data", "06_disease")
targets <- c("EAS_AD", "EAS_SCZ", "EUR_ANX", "EUR_MDD")

# --- helper functions ---

read_gsa <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fread(text = paste(lines[!startsWith(lines, "#")], collapse = "\n"))
}

num_equal <- function(x, y, tolerance = 1e-12) {
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) <= tolerance * pmax(1, abs(x), abs(y)))
}

# --- binary sldsc ---

binary_new <- rbindlist(lapply(targets, function(gwas_id) {
  files <- list.files(file.path(run_dir, "binary_sldsc", gwas_id), pattern = "\\.results$", full.names = TRUE)
  rbindlist(lapply(files, function(path) {
    x <- fread(path)[.N]
    x[, `:=`(
      Trait = sub("^(EAS|EUR)_", "", gwas_id),
      Population = sub("_.*$", "", gwas_id),
      Annotation = sub(paste0("^", gwas_id, "_"), "", sub("\\.results$", "", basename(path)))
    )]
    x
  }))
}))

binary_old <- fread(file.path(run_dir, "sldsc_enrichment_summary.tsv"))
binary_old[, gwas_id := paste(Population, Trait, sep = "_")]
binary_keep <- binary_old[!gwas_id %in% targets]
binary_replace <- binary_new[, .(
  Trait, Population, Annotation, Coefficient, Coefficient_std_error,
  `Coefficient_z-score`, Enrichment, Enrichment_std_error, Enrichment_p
)]
binary_full <- rbindlist(list(binary_keep[, names(binary_replace), with = FALSE], binary_replace), use.names = TRUE)
binary_full[, FDR_p := p.adjust(Enrichment_p, method = "BH"), by = .(Population, Trait)]
binary_full[, FDR_p_global := p.adjust(Enrichment_p, method = "BH")]
binary_full[, n_tests_family := .N, by = .(Population, Trait)]
setorder(binary_full, Population, Trait, Annotation)
fwrite(binary_full, file.path(run_dir, "binary_sldsc_frozen_full.tsv"), sep = "\t")

# --- joint maxcpp sldsc ---

joint_old <- fread(file.path(root, "data", "04_qtl", "03_sldsc_tau_star.tsv"))
sd_dt <- unique(
  fread(file.path(root, "data", "04_qtl", "annot_sd.tsv"))[
    annot_set == "maxCPP", .(population, stage, annotation, sd_annot = sd, M_common = n_snps)
  ],
  by = c("population", "stage", "annotation")
)

joint_new <- rbindlist(lapply(targets, function(gwas_id) {
  rbindlist(lapply(c("prenatal", "postnatal"), function(stage_name) {
    result_path <- file.path(run_dir, "joint_sldsc", gwas_id, stage_name, paste0(gwas_id, ".results"))
    log_path <- file.path(run_dir, "joint_sldsc", gwas_id, stage_name, paste0(gwas_id, ".log"))
    x <- fread(result_path)[grepl("^(apaQTL|eQTL)L2_[0-9]+$", Category)]
    h2_line <- grep("^Total Observed scale h2:", readLines(log_path, warn = FALSE), value = TRUE)[1]
    h2_match <- regmatches(h2_line, regexec("h2:\\s*([0-9.eE+-]+)\\s*\\(([0-9.eE+-]+)\\)", h2_line))[[1]]
    stopifnot(length(h2_match) == 3L)
    x[, `:=`(
      population = sub("_.*$", "", gwas_id),
      stage = stage_name,
      annotation = fifelse(startsWith(Category, "apaQTL"), "apaQTL", "eQTL"),
      trait = gwas_id,
      model = "maxcpp",
      h2 = as.numeric(h2_match[2]),
      h2_se = as.numeric(h2_match[3])
    )]
    x
  }))
}))
setnames(joint_new, "Coefficient_z-score", "coef_z")
joint_new[, coef_p_1s := pnorm(coef_z, lower.tail = FALSE)]
joint_new[, coef_fdr := p.adjust(coef_p_1s, method = "BH"), by = .(population, stage, trait)]
joint_new[, n_tests_family := .N, by = .(population, stage, trait)]
joint_new[, well_powered := (h2 / h2_se) >= 7]
joint_new <- merge(joint_new, sd_dt, by = c("population", "stage", "annotation"), all.x = TRUE)
stopifnot(!anyNA(joint_new[, .(sd_annot, M_common)]))
joint_new[, tau_star := fifelse(h2 > 0, Coefficient * M_common * sd_annot / h2, NA_real_)]
joint_new[, tau_star_se := fifelse(h2 > 0, Coefficient_std_error * M_common * sd_annot / h2, NA_real_)]
joint_new[, `:=`(tau_lo = tau_star - 1.96 * tau_star_se, tau_hi = tau_star + 1.96 * tau_star_se)]

joint_columns <- names(joint_old)
joint_full <- rbindlist(list(
  joint_old[!trait %in% targets],
  joint_new[, joint_columns, with = FALSE]
), use.names = TRUE)
setorder(joint_full, population, stage, annotation, trait)
stopifnot(!anyDuplicated(joint_full[, .(population, stage, annotation, trait)]))
fwrite(joint_full, file.path(run_dir, "joint_sldsc_frozen_full.tsv"), sep = "\t")

# --- magma ---

magma_new <- rbindlist(lapply(targets, function(gwas_id) {
  x <- read_gsa(file.path(run_dir, "magma", gwas_id, "apa_qtl_sets.gsa.out"))
  setnames(x, tolower(names(x)))
  x[, `:=`(gwas_id = gwas_id, population = sub("_.*$", "", gwas_id))]
  x
}))

magma_old <- fread(file.path(run_dir, "magma_gsa_summary.tsv"))
magma_replace <- magma_new[, .(
  gwas_id, population, set_name = variable, type, ngenes, beta, beta_std, se, p,
  source_gsa = file.path(run_dir, "magma", gwas_id, "apa_qtl_sets.gsa.out")
)]
magma_full <- rbindlist(list(
  magma_old[!gwas_id %in% targets, .(gwas_id, population, set_name, type, ngenes, beta, beta_std, se, p, source_gsa)],
  magma_replace
), use.names = TRUE)
magma_full[, fdr_within_gwas := p.adjust(p, method = "BH"), by = gwas_id]
magma_full[, significant_within_gwas := fdr_within_gwas < 0.05]
magma_full[, `:=`(
  trait_id = sub("^(EAS|EUR)_", "", gwas_id),
  trait_axis = gwas_id,
  trait_label = sub("^(EAS|EUR)_", "", gwas_id),
  neg_log10_p = -log10(p),
  hit_label = fifelse(significant_within_gwas, gwas_id, "")
)]
setorder(magma_full, population, trait_id, set_name)
fwrite(magma_full, file.path(run_dir, "magma_frozen_full.tsv"), sep = "\t")

magma_conditional <- rbindlist(lapply(targets, function(gwas_id) {
  x <- read_gsa(file.path(run_dir, "magma", gwas_id, "prenatal_condition_postnatal.gsa.out"))
  setnames(x, tolower(names(x)))
  x[, `:=`(gwas_id = gwas_id, population = sub("_.*$", "", gwas_id))]
  x
}))
setcolorder(magma_conditional, c("gwas_id", "population", setdiff(names(magma_conditional), c("gwas_id", "population"))))
fwrite(magma_conditional, file.path(run_dir, "magma_conditional_frozen_affected.tsv"), sep = "\t")

# --- old versus frozen ---

binary_compare <- merge(
  binary_old[gwas_id %in% targets, .(
    branch = "binary_sldsc", gwas_id, test = Annotation,
    old_estimate = Enrichment, old_p = Enrichment_p, old_fdr = FDR_p
  )],
  binary_full[paste(Population, Trait, sep = "_") %in% targets, .(
    gwas_id = paste(Population, Trait, sep = "_"), test = Annotation,
    new_estimate = Enrichment, new_p = Enrichment_p, new_fdr = FDR_p
  )],
  by = c("gwas_id", "test")
)
binary_compare[, branch := "binary_sldsc"]

joint_compare <- merge(
  joint_old[trait %in% targets, .(
    gwas_id = trait, test = paste(stage, annotation, sep = ":"),
    old_estimate = tau_star, old_p = coef_p_1s, old_fdr = coef_fdr
  )],
  joint_full[trait %in% targets, .(
    gwas_id = trait, test = paste(stage, annotation, sep = ":"),
    new_estimate = tau_star, new_p = coef_p_1s, new_fdr = coef_fdr
  )],
  by = c("gwas_id", "test")
)
joint_compare[, branch := "joint_sldsc"]

magma_compare <- merge(
  magma_old[gwas_id %in% targets, .(
    gwas_id, test = set_name, old_estimate = beta_std, old_p = p, old_fdr = fdr_within_gwas
  )],
  magma_full[gwas_id %in% targets, .(
    gwas_id, test = as.character(set_name), new_estimate = beta_std, new_p = p, new_fdr = fdr_within_gwas
  )],
  by = c("gwas_id", "test")
)
magma_compare[, branch := "magma"]

comparison <- rbindlist(list(binary_compare, joint_compare, magma_compare), use.names = TRUE)
comparison[, `:=`(
  estimate_delta = new_estimate - old_estimate,
  p_delta = new_p - old_p,
  fdr_delta = new_fdr - old_fdr,
  value_changed = !(num_equal(old_estimate, new_estimate) & num_equal(old_p, new_p)),
  significance_changed = (old_fdr < 0.05) != (new_fdr < 0.05)
)]
setcolorder(comparison, c(
  "branch", "gwas_id", "test", "old_estimate", "new_estimate", "estimate_delta",
  "old_p", "new_p", "p_delta", "old_fdr", "new_fdr", "fdr_delta",
  "value_changed", "significance_changed"
))
setorder(comparison, branch, gwas_id, test)
fwrite(comparison, file.path(run_dir, "old_vs_frozen_results.tsv"), sep = "\t")

cat(sprintf("binary rows: %d\n", nrow(binary_full)))
cat(sprintf("joint rows: %d\n", nrow(joint_full)))
cat(sprintf("magma rows: %d\n", nrow(magma_full)))
cat(sprintf("changed comparisons: %d/%d\n", comparison[value_changed == TRUE, .N], nrow(comparison)))
cat(sprintf("significance changes: %d\n", comparison[significance_changed == TRUE, .N]))
