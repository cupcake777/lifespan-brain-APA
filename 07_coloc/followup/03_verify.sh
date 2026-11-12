#!/usr/bin/env bash
# Checks the followup outputs (R1 in data/03_omics, R2-R3 in data/07_coloc, R4 in data/04_qtl) against expected row counts, keys and gzip integrity.
# Exits non-zero on the first failed check; prints PASS when all pass.
: "${APA_ROOT:?set APA_ROOT (see release/config/paths.sh)}"
set -euo pipefail

out_dir="${APA_ROOT}/data/07_coloc"
omics_dir="${APA_ROOT}/data/03_omics"
qtl_dir="${APA_ROOT}/data/04_qtl"

# --- R1 ---
gzip -t "$omics_dir/R1_expr_TPM.tsv.gz" "$omics_dir/R1_expr_junction.tsv.gz"
gzip -t "$omics_dir/R1b_cis_prot_models.tsv.gz" "$omics_dir/R1c_cis_apa_mrna_models.tsv.gz"
test "$(awk 'END {print NR-1}' "$omics_dir/R1a_mediation_expr_TPM.tsv")" -eq 772
test "$(awk 'END {print NR-1}' "$omics_dir/R1a_mediation_junction.tsv")" -eq 772
test "$(awk -F '\t' 'NR>1 && $53=="tested" {n++} END {print n+0}' "$omics_dir/R1a_mediation_expr_TPM.tsv")" -eq 772
test "$(awk -F '\t' 'NR>1 && $53=="tested" {n++} END {print n+0}' "$omics_dir/R1a_mediation_junction.tsv")" -eq 685
test "$(awk -F '\t' 'NR>1 && $1=="QTL" && $2==679 && $3==679 {n++} END {print n+0}' "$omics_dir/R1_SJ_sample_coverage_sum.tsv")" -eq 1
test "$(awk -F '\t' 'NR>1 && $1=="protein" && $2==371 && $3==371 {n++} END {print n+0}' "$omics_dir/R1_SJ_sample_coverage_sum.tsv")" -eq 1
test "$(gzip -cd "$omics_dir/R1b_cis_prot_models.tsv.gz" | awk -F '\t' 'NR>1 && $7!="" {n++} END {print n+0}')" -gt 40000
test "$(gzip -cd "$omics_dir/R1c_cis_apa_mrna_models.tsv.gz" | awk -F '\t' 'NR>1 && $7!="" {n++} END {print n+0}')" -gt 60000
test "$(awk 'END {print NR-1}' "$omics_dir/R1bc_missing_exact_contract.tsv")" -eq 2

# --- R2 ---
test "$(awk 'END {print NR-1}' "$out_dir/R2a_PRDX6_sample_data.tsv")" -eq 697
test "$(awk 'END {print NR-1}' "$out_dir/R2a_HID1_sample_data.tsv")" -eq 697
test "$(awk 'END {print NR-1}' "$out_dir/R2b_selected_bam_samples.tsv")" -eq 12
test "$(gzip -cd "$out_dir/R2b_depth_cpm.tsv.gz" | awk 'END {print NR-1}')" -eq 52044
# --- R3 ---
jq -e '.meta.requested_index_snp == "rs11136325" and .meta.requested_index_present == true and ([.points[].snp] | index("rs11136325") != null)' "$out_dir/Fig5e_NRBP2_locuscompare.json" >/dev/null
jq -e '.meta.locus == "rs7251" and ([.points[].snp] | index("rs7047") != null) and ([.points[].snp] | index("rs7251") != null)' "$out_dir/Fig5d_SCAF1_locuscompare.json" >/dev/null
test "$(awk 'END {print NR-1}' "$out_dir/R3_EUR_AD_provenance.tsv")" -eq 1
# --- R4 ---
test "$(awk 'END {print NR-1}' "$qtl_dir/R4_nan2026_source_contract.tsv")" -eq 8
test "$(awk -F '\t' 'NR > 1 && $1 ~ /Supplementary Table 11 workbook/ && $5 ~ /6706 valid significant lead g-APA rows/ {n++} END {print n+0}' "$qtl_dir/R4_nan2026_source_contract.tsv")" -eq 1
test "$(awk -F '\t' 'NR > 1 && $3 ~ /sensitivity only/ {n++} END {print n+0}' "$qtl_dir/R4_nan2026_model_sens.tsv")" -gt 0
gzip -t "$qtl_dir/R4_nan2026_nominal_matches.tsv.gz"
test "$(awk 'END {print NR-1}' "$qtl_dir/R4_intx_nan_exact_pairs.tsv")" -eq 1158
test "$(awk -F '\t' 'NR>1 && $15=="TRUE" {n++} END {print n+0}' "$qtl_dir/R4_intx_nan_exact_pairs.tsv")" -eq 193
test "$(awk -F '\t' 'NR>1 && $1=="Amplified" && $2==40 && $3==20 {n++} END {print n+0}' "$qtl_dir/R4_intx_nan_replication_rates.tsv")" -eq 1
test "$(awk -F '\t' 'NR>1 && $1=="Attenuated" && $2==33 && $3==12 {n++} END {print n+0}' "$qtl_dir/R4_intx_nan_replication_rates.tsv")" -eq 1
test "$(awk -F '\t' 'NR>1 && $1=="Opposing" && $2==120 && $3==17 {n++} END {print n+0}' "$qtl_dir/R4_intx_nan_replication_rates.tsv")" -eq 1
test "$(awk -F '\t' 'NR>1 && $1=="Prenatal-specific" && $2==106 && $3==23 {n++} END {print n+0}' "$qtl_dir/R4_stage_nan_replication_rates.tsv")" -eq 1
test "$(awk -F '\t' 'NR>1 && $1=="Postnatal-specific" && $2==1686 && $3==479 {n++} END {print n+0}' "$qtl_dir/R4_stage_nan_replication_rates.tsv")" -eq 1
test "$(awk -F '\t' 'NR>1 && $1=="Shared" && $2==563 && $3==194 {n++} END {print n+0}' "$qtl_dir/R4_stage_nan_replication_rates.tsv")" -eq 1
test "$(awk -F '\t' 'NR>1 && $2=="likelihood-ratio" && $7 != "" {n++} END {print n+0}' "$qtl_dir/R4_nan2026_exact_model_outputs.tsv")" -eq 2
test "$(awk -F '\t' 'NR>1 && $4 !~ /not downloaded/ {n++} END {print n+0}' "$qtl_dir/R4_nan2026_source_contract.tsv")" -eq 8

# --- summary ---
printf 'PASS\tR1_mediation=772_each\tR1_SJ=679_QTL+371_protein\tR1_cis=finite\tR2a_rows=697_each\tR2b_cells=12\tdepth_rows=52044\tR3_json_indices=present\tEUR_AD_signals=1\tR4_exact_models=accepted\n'
