#!/usr/bin/env bash
# EUR AD2: apaQTL vs eQTL S-LDSC enrichment per stage, printed as one table
: "${APA_ROOT:?source config/paths.sh}"
: "${APA_HPC_HOME:?source config/paths.sh}"
set -euo pipefail

# --- paths ---
ref_dir=${APA_GROUP_SHARE}/.ref/sLDSC/EUR_hg38
ldsc_python=${APA_HPC_HOME}/.local/share/mamba/envs/ldsc3/bin/python
ldsc_bin=${APA_HPC_HOME}/.local/share/mamba/envs/ldsc3/bin/ldsc.py
sumstats="${APA_ROOT}/data/06_disease/gwas/EUR_AD2.sumstats.gz"
ldscore_dir="${APA_ROOT}/data/06_disease/ld_scores/EUR"
apa_dir="${APA_ROOT}/data/06_disease/regression"

tmp_dir=$(mktemp -d /dev/shm/fig5a_eur_ad2_eqtl.XXXXXXXX)
trap 'rm -rf "$tmp_dir"' EXIT

printf 'GWAS\tPopulation\tTrait\tStage\tapaQTL_enrichment\tapaQTL_enrichment_se\tapaQTL_enrichment_p\teQTL_enrichment\teQTL_enrichment_se\teQTL_enrichment_p\tapaQTL_annotation\teQTL_annotation\n'
for stage in prenatal postnatal; do
  annotation="${stage}_eqtl"
  output="${tmp_dir}/EUR_AD2_${annotation}"
  "${ldsc_python}" "${ldsc_bin}" \
    --h2 "${sumstats}" \
    --ref-ld-chr "${ref_dir}/baselineLD_v2.2/baselineLD.,${ldscore_dir}/${annotation}." \
    --w-ld-chr "${ref_dir}/weights/weights.hm3_noMHC." \
    --frqfile-chr "${ref_dir}/plink_files/1000G.EUR.hg38." \
    --overlap-annot --print-coefficients --out "${output}" \
    > "${output}.stdout.log" 2>&1
  apa=$(awk -F '\t' 'END {print $5 "\t" $6 "\t" $7}' "${apa_dir}/EUR_AD2_full_${stage}_apa_qtl.results")
  eqtl=$(awk -F '\t' 'END {print $5 "\t" $6 "\t" $7}' "${output}.results")
  printf 'EUR_AD2\tEUR\tAD2\t%s\t%s\t%s\t%s\t%s\n' "$stage" "$apa" "$eqtl" "${stage}_apa_qtl" "$annotation"
done
