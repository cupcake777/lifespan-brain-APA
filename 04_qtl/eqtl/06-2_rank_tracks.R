# Rank APA events for mRNA expression validation: mechanism, coloc, mediation, disease and core
# labels per event, then priority classes. Reads data/04_qtl/06_master.tsv; writes 06_candidate.tsv, 06_summary.tsv.
suppressMessages(library(data.table))

APA_ROOT <- Sys.getenv("APA_ROOT")
D <- file.path(APA_ROOT, "data", "04_qtl")
gtf_path <- file.path(APA_ROOT, "ref", "gencode.v40.annotation.gtf.gz")

PP_H4_MIN <- 0.8
FDR_CUT <- 0.05
GWAS_PP4_MIN <- 0.5
MAF_MIN <- 0.05
DIRECT_APA_DISTANCE_BP <- 1000
DELTA_UTR_MIN <- 300
DELTA_UTR_MAX <- 2000

# --- helpers ---
as_flag <- function(x) {
  y <- as.logical(x)
  fifelse(is.na(y), FALSE, y)
}
first_nonempty <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0L) "" else paste(x, collapse = ";")
}
safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}
safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

# --- inputs ---
master <- fread(file.path(D, "06_master.tsv"))
gene_annot <- fread(
  cmd = sprintf(
    "zcat %s | awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"} $3==\"gene\" {if(match($9,/gene_id \\\"([^\\\"]+)\\\"/,a) && match($9,/gene_name \\\"([^\\\"]+)\\\"/,b) && match($9,/gene_type \\\"([^\\\"]+)\\\"/,c)) print a[1],b[1],c[1]}'",
    shQuote(gtf_path)),
  col.names = c("gene_id", "gene_symbol", "gene_type"))
gene_annot[, ensg0 := sub("\\..*$", "", gene_id)]
gene_annot <- unique(gene_annot[, .(ensg0, gene_symbol, gene_type)])
apa_state <- unique(fread(file.path(D, "05_gene_state_apa.tsv"))[, .(gene_name, development_state_05 = gene_state)])

num_cols <- c("af", "maf", "pval_nominal", "slope", "beta_alt_PDUI", "start_distance", "end_distance",
              "pip", "variant_phastcons", "PP_H2", "PP_H3", "PP_H4", "h4_snp_pp", "ab_fwd", "fdr_ab_fwd",
              "prop_med", "lfsr_prenatal", "lfsr_postnatal", "gwas_best_pp4")
for (cl in intersect(num_cols, names(master))) set(master, j = cl, value = suppressWarnings(as.numeric(master[[cl]])))

# --- per-variant flags ---
master[, pas_mechanism := as_flag(direct_pas_evidence)]
master[, rbp_mechanism := as_flag(local_rbp_motif_change) | as_flag(local_clip_support)]
master[, aparent_mechanism := as_flag(aparent_beta_direction_consistent)]
master[, local_mechanism := pas_mechanism | rbp_mechanism | aparent_mechanism]
master[, h4_shared := coloc_class == "H4 shared causal variant" & !is.na(PP_H4) & PP_H4 >= PP_H4_MIN]
master[, apa_to_expr := !is.na(fdr_ab_fwd) & fdr_ab_fwd < FDR_CUT & mediation_call %in% c("APA -> EXP", "Partial mediation")]
master[, disease_coloc := !is.na(gwas_best_pp4) & gwas_best_pp4 >= GWAS_PP4_MIN]
master[, abs_beta := abs(beta_alt_PDUI)]
master[, event_distance_bp := pmin(abs(start_distance), abs(end_distance), na.rm = TRUE)]
master[is.infinite(event_distance_bp), event_distance_bp := NA_real_]
master[, direct_apa_variant_near_event := !is.na(event_distance_bp) & event_distance_bp <= DIRECT_APA_DISTANCE_BP]

# --- one row per APA event (representative variant = best by coloc, mediation, disease, mechanism, PIP) ---
candidate <- master[, {
  rep_i <- order(-as.integer(h4_shared), -as.integer(apa_to_expr), -as.integer(disease_coloc),
                 -as.integer(local_mechanism), -fifelse(is.na(pip), 0, pip),
                 -fifelse(is.na(abs_beta), 0, abs_beta), na.last = TRUE)[1L]
  mechanism_labels <- c("PAS", "RBP_or_CLIP", "APARENT")
  mechanism_flags <- c(any(pas_mechanism), any(rbp_mechanism), any(aparent_mechanism))
  best_gwas_i <- which.max(fifelse(is.na(gwas_best_pp4), -Inf, gwas_best_pp4))
  best_gwas_trait <- if (length(best_gwas_i) == 0L) "" else first_nonempty(gwas_best_trait[best_gwas_i])
  .(
    gene_name = gene_name[rep_i],
    ensg0 = ensg0[rep_i],
    n_signal_variants = .N,
    representative_variant = variant_id[rep_i],
    representative_maf = maf[rep_i],
    representative_pip = pip[rep_i],
    representative_beta_alt_PDUI = beta_alt_PDUI[rep_i],
    representative_apa_direction = apa_direction[rep_i],
    min_event_distance_bp = safe_min(event_distance_bp),
    n_variants_within_1kb_event = sum(direct_apa_variant_near_event, na.rm = TRUE),
    local_mechanism_bonus = any(local_mechanism),
    local_mechanism_type = { out <- mechanism_labels[mechanism_flags %in% TRUE]; if (length(out) == 0L) "none" else paste(unique(out), collapse = ";") },
    h4_shared = any(h4_shared),
    max_PP_H4 = safe_max(PP_H4),
    coloc_class = coloc_class[rep_i],
    h4_snp = first_nonempty(h4_snp),
    h4_snp_pp = safe_max(h4_snp_pp),
    apa_to_expr = any(apa_to_expr),
    mediation_call = first_nonempty(mediation_call),
    min_fdr_ab_fwd = safe_min(fdr_ab_fwd),
    max_abs_prop_med = safe_max(abs(prop_med)),
    apa_effect_on_expr = first_nonempty(apa_effect_on_expr),
    disease_coloc = any(disease_coloc),
    gwas_best_trait = best_gwas_trait,
    gwas_best_pp4 = safe_max(gwas_best_pp4),
    gwas_traits_ge_min = first_nonempty(gwas_traits_ge_min),
    gwas_n_loci_ge_min = safe_max(gwas_n_loci_ge_min),
    gwas_exploratory = any(as_flag(gwas_exploratory)),
    gwas_prior_qc_pass = any(as_flag(gwas_prior_qc_pass)),
    gwas_proposed_variants = first_nonempty(gwas_proposed_variants),
    mash_class = first_nonempty(mash_class),
    min_lfsr = safe_min(c(lfsr_prenatal, lfsr_postnatal)),
    core_n_seeds = safe_max(core_n_seeds),
    core_replicated_majority = any(as_flag(core_replicated_majority)),
    core_replicated_all = any(as_flag(core_replicated_all))
  )
}, by = .(stage, phenotype_id)]

# --- gene symbol and biotype from GENCODE (by ENSG, then by symbol) ---
candidate <- merge(candidate, gene_annot[, .(ensg0, gene_symbol_annot = gene_symbol, gene_type_ensg = gene_type)],
                   by = "ensg0", all.x = TRUE)
candidate[, gene_symbol := fifelse(is.na(gene_symbol_annot) | gene_symbol_annot == "", gene_name, gene_symbol_annot)]
candidate <- merge(candidate, unique(gene_annot[, .(gene_symbol, gene_type_by_symbol = gene_type)]), by = "gene_symbol", all.x = TRUE)
candidate[, gene_type := fifelse(!is.na(gene_type_ensg), gene_type_ensg, gene_type_by_symbol)]
candidate[, c("gene_symbol_annot", "gene_type_ensg", "gene_type_by_symbol") := NULL]
candidate <- unique(candidate)
candidate <- merge(candidate, apa_state, by = "gene_name", all.x = TRUE)
candidate[is.na(development_state_05), development_state_05 := "Not in 05 APA state table"]

phen_parts <- tstrsplit(candidate$phenotype_id, "\\|")
event_span <- tstrsplit(phen_parts[[3]], "[:-]")
candidate[, delta_utr_length_bp := as.integer(event_span[[3]]) - as.integer(event_span[[2]]) + 1L]
candidate[, delta_utr_preferred := !is.na(delta_utr_length_bp) & delta_utr_length_bp >= DELTA_UTR_MIN & delta_utr_length_bp <= DELTA_UTR_MAX]

for (cl in c("min_event_distance_bp", "max_PP_H4", "h4_snp_pp", "min_fdr_ab_fwd", "max_abs_prop_med",
             "gwas_best_pp4", "gwas_n_loci_ge_min", "min_lfsr", "core_n_seeds")) {
  candidate[is.infinite(get(cl)), (cl) := NA_real_]
}
candidate[, validation_level := fifelse(n_variants_within_1kb_event > 0 |
                                          (!is.na(min_event_distance_bp) & min_event_distance_bp <= DIRECT_APA_DISTANCE_BP),
                                        "APA-level or mRNA-level", "mRNA-level")]
candidate[, representative_variant_feasible := !is.na(representative_maf) & representative_maf >= MAF_MIN]
candidate[, core_status_05 := fcase(core_replicated_all, "Strict core", core_replicated_majority, "Consensus core",
                                    !is.na(core_n_seeds) & core_n_seeds > 0, "Any-seed core", default = "Not core")]

# --- priority classes and ranking ---
candidate[, priority_class := fcase(
  h4_shared & apa_to_expr & disease_coloc, "APA_expression_disease",
  h4_shared & apa_to_expr, "APA_expression",
  h4_shared & disease_coloc, "APA_disease_no_mediation",
  disease_coloc, "disease_relevance_only",
  h4_shared, "APA_eQTL_coloc_only",
  default = "lower_priority")]
candidate[, priority_rank := fcase(
  priority_class == "APA_expression_disease", 1L,
  priority_class == "APA_expression", 2L,
  priority_class == "APA_disease_no_mediation", 3L,
  priority_class == "disease_relevance_only", 5L,
  priority_class == "APA_eQTL_coloc_only", 6L,
  default = 7L)]
setorder(candidate, priority_rank, -representative_variant_feasible, -gwas_prior_qc_pass, -gwas_best_pp4,
         min_fdr_ab_fwd, -max_abs_prop_med, -max_PP_H4, -local_mechanism_bonus, min_event_distance_bp, na.last = TRUE)
candidate[, candidate_rank := .I]
candidate[, report_eas_disease_followup := grepl("^EAS_", gwas_best_trait) & disease_coloc &
            priority_class %in% c("APA_expression_disease", "APA_disease_no_mediation", "disease_relevance_only")]
candidate[, report_prenatal_discovery := stage == "prenatal" &
            priority_class %in% c("APA_expression", "APA_disease_no_mediation", "disease_relevance_only")]

summary <- candidate[, .(
  n_candidates = .N,
  n_genes = uniqueN(gene_name),
  n_h4_shared = sum(h4_shared),
  n_apa_to_expr = sum(apa_to_expr),
  n_disease_coloc = sum(disease_coloc),
  n_local_mechanism_bonus = sum(local_mechanism_bonus),
  n_representative_variant_feasible = sum(representative_variant_feasible)
), by = .(priority_rank, priority_class)]
setorder(summary, priority_rank)

fwrite(candidate, file.path(D, "06_candidate.tsv"), sep = "\t")
fwrite(summary, file.path(D, "06_summary.tsv"), sep = "\t")
