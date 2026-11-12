# Signal-level master table: module-07 mechanism signals joined to coloc, mediation, mash posteriors,
# robust core and GWAS disease chains. Writes data/04_qtl/06_master.tsv and 06_join.tsv.
suppressMessages(library(data.table))

APA_ROOT <- Sys.getenv("APA_ROOT")
D <- file.path(APA_ROOT, "data", "04_qtl")
M7 <- file.path(APA_ROOT, "data", "07_coloc")
GWAS_PP4_MIN <- 0.5

# --- signal table (mechanism layer, module 07) ---
sig <- fread(cmd = sprintf("zcat %s", shQuote(file.path(M7, "02_signal_mechanisms.tsv.gz"))))
sig_cols <- c("signal_id", "candidate_id", "stage", "phenotype_id", "gene_name", "variant_id", "chr", "pos",
              "ref", "alt", "af", "ma_count", "pval_nominal", "slope", "slope_se", "beta_alt_PDUI", "apa_direction",
              "start_distance", "end_distance", "pip", "finemap_tier", "credible_set_member", "signal_assignment_method",
              "direct_pas_evidence", "proximal_pas_core_change", "distal_pas_core_change", "mechanism_direction_concordant",
              "local_rbp_motif_change", "local_clip_support", "local_structure_change", "postar_variant_overlap_count",
              "variant_phastcons", "APARENT_delta_proximal", "APARENT_delta_distal", "APARENT_predicted_PDUI_delta",
              "aparent_beta_direction_consistent", "snp_to_apa_evidence_label", "mechanism_star_eligible",
              "gwas_coloc_max_pp4", "has_gwas_coloc_result", "interaction_group_id", "gxstage_class")
master <- sig[, intersect(sig_cols, names(sig)), with = FALSE]
for (cl in c("af", "pip", "beta_alt_PDUI", "variant_phastcons", "gwas_coloc_max_pp4",
             "APARENT_predicted_PDUI_delta", "pval_nominal", "slope")) {
  if (cl %in% names(master)) set(master, j = cl, value = suppressWarnings(as.numeric(master[[cl]])))
}
master[, maf := pmin(af, 1 - af, na.rm = TRUE)]
master[af > 0.5 & !is.na(af), maf := 1 - af]

# --- coloc (pair level: stage + apa phenotype) ---
coloc <- fread(file.path(D, "coloc_all_stages.tsv"))
coloc_sub <- coloc[, .(stage, phenotype_id = apa_pheno, ensg0, coloc_class,
                       PP_H2, PP_H3, PP_H4, coloc_direction = direction,
                       h4_snp, h4_snp_pp, h4_slope_eqtl, h4_slope_apa, is_egene)]
master <- merge(master, coloc_sub, by = c("stage", "phenotype_id"), all.x = TRUE)
master[is.na(coloc_class), coloc_class := "not_tested"]

# --- mediation (stage + gene; the table has no APA phenotype id) ---
med <- fread(file.path(D, "04_mediation_results.tsv"))
stopifnot(uniqueN(med, by = c("stage", "ensg0")) == nrow(med))
med_sub <- med[, .(stage, ensg0, mediation_call, ab_fwd, fdr_ab_fwd, prop_med, apa_effect_on_expr)]
master <- merge(master, med_sub, by = c("stage", "ensg0"), all.x = TRUE)

# --- mash posteriors (variant level) ---
mash <- fread(file.path(D, "05_mash_posterior_apa.tsv"))
mash_sub <- mash[, .(phenotype_id, variant_id, lfsr_prenatal, lfsr_postnatal, mash_class, shared_by_magnitude)]
master <- merge(master, mash_sub, by = c("phenotype_id", "variant_id"), all.x = TRUE)

# --- robust core: prenatal event replication across downsampled postnatal seeds ---
# The core denominator is the prenatal significant set, so "not core" does not mean unstable across stages.
core <- fread(file.path(D, "05_controlled_core_groups.tsv"))
core_sub <- core[, .(phenotype_id = prenatal_phenotype_id, core_n_seeds = n_replicating_seeds,
                     core_replicated_any = replicated_any, core_replicated_majority = replicated_majority,
                     core_replicated_all = replicated_all)]
master <- merge(master, core_sub, by = "phenotype_id", all.x = TRUE)
master[, core_applicable := !is.na(core_n_seeds)]

# --- GWAS disease chains (module 07), aggregated from locus to pair ---
gwas <- fread(cmd = sprintf("zcat %s", shQuote(file.path(M7, "04_disease_chains.tsv.gz"))))
gwas <- gwas[stage %in% c("prenatal", "postnatal")]
gwas[, gwas_coloc_pp4 := suppressWarnings(as.numeric(gwas_coloc_pp4))]
gwas_agg <- gwas[, {
  hit <- !is.na(gwas_coloc_pp4)
  best_i <- if (any(hit)) which.max(ifelse(hit, gwas_coloc_pp4, -Inf)) else NA_integer_
  .(
    gwas_best_pp4 = if (any(hit)) gwas_coloc_pp4[best_i] else NA_real_,
    gwas_best_trait = if (any(hit)) gwas_trait[best_i] else NA_character_,
    gwas_n_loci_ge_min = sum(hit & gwas_coloc_pp4 >= GWAS_PP4_MIN),
    gwas_traits_ge_min = paste(unique(gwas_trait[hit & gwas_coloc_pp4 >= GWAS_PP4_MIN]), collapse = ";"),
    gwas_exploratory = any(disease_chain_evidence_label == "exploratory_disease_chain", na.rm = TRUE),
    gwas_prior_qc_pass = any(coloc_prior_sensitivity_acceptable == "True", na.rm = TRUE),
    gwas_proposed_variants = paste(unique(proposed_functional_variant[proposed_functional_variant != ""]), collapse = ";")
  )
}, by = .(stage, phenotype_id)]
master <- merge(master, gwas_agg, by = c("stage", "phenotype_id"), all.x = TRUE)
master[is.na(gwas_n_loci_ge_min), gwas_n_loci_ge_min := 0L]
master[is.na(gwas_exploratory), gwas_exploratory := FALSE]

# --- join diagnostics ---
diag <- rbindlist(list(
  data.table(layer = "signals_total", n = nrow(master)),
  data.table(layer = "distinct_phenotypes", n = uniqueN(master$phenotype_id)),
  data.table(layer = "coloc_joined", n = sum(master$coloc_class != "not_tested")),
  data.table(layer = "coloc_not_tested", n = sum(master$coloc_class == "not_tested")),
  data.table(layer = "mediation_joined", n = sum(!is.na(master$fdr_ab_fwd))),
  data.table(layer = "mash_joined", n = sum(!is.na(master$mash_class))),
  data.table(layer = "core_joined(prenatal)", n = sum(!is.na(master$core_n_seeds))),
  data.table(layer = "gwas_any_locus", n = sum(!is.na(master$gwas_best_pp4))),
  data.table(layer = sprintf("gwas_ge_%.1f", GWAS_PP4_MIN), n = sum(master$gwas_n_loci_ge_min > 0)),
  data.table(layer = "gwas_exploratory", n = sum(master$gwas_exploratory))
))
fwrite(master, file.path(D, "06_master.tsv"), sep = "\t")
fwrite(diag, file.path(D, "06_join.tsv"), sep = "\t")
