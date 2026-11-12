#!/usr/bin/env bash
# Frozen-GWAS downstream runs: binary sLDSC, joint maxCPP sLDSC, MAGMA, SMR/HEIDI, S-PrediXcan TWAS.
# Submit from $APA_ROOT so that logs/ resolves to $APA_ROOT/logs.
#SBATCH --job-name=apa_frozen
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --partition=cu

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# --- config ---

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="${APA_ROOT}"
run_dir="${root}/data/06_disease"
manifest="${root}/04_coloc/work/followup_staging/gwas_source_checksum_manifest_20261028.tsv"
targets=(EAS_AD EAS_SCZ EUR_ANX EUR_MDD)

ldsc_python="${APA_HPC_HOME}/.local/share/mamba/envs/ldsc/bin/python"
ldsc_dir="${root}/../ref_repo/ldsc"
sldsc_ref="${APA_HPC_HOME}/project/.ref/sLDSC"
magma="${MAGMA:-magma}"

mode="${MODE:?Set MODE to binary_sldsc, joint_sldsc, magma, smr, or twas}"
gwas_id="${GWAS_ID:?Set GWAS_ID}"

if [[ ! " ${targets[*]} " =~ " ${gwas_id} " ]]; then
  echo "Unsupported GWAS_ID: ${gwas_id}" >&2
  exit 2
fi

if [[ "${mode}" == "twas" ]]; then
  twas_stage="${TWAS_STAGE:?Set TWAS_STAGE to prenatal or postnatal}"
  [[ "${twas_stage}" == "prenatal" || "${twas_stage}" == "postnatal" ]]
  out_dir="${run_dir}/twas/${gwas_id}/${twas_stage}"
else
  out_dir="${run_dir}/${mode}/${gwas_id}"
fi
mkdir -p "${root}/logs" "${out_dir}"

# --- frozen input identity ---

read -r frozen_path expected_sha < <(awk -F '\t' -v id="${gwas_id}" '$1 == id {print $2, $5}' "${manifest}")
[[ -s "${frozen_path}" && -n "${expected_sha}" ]]
observed_sha="$(sha256sum "${frozen_path}" | awk '{print $1}')"
[[ "${observed_sha}" == "${expected_sha}" ]]
printf 'gwas_id\tfrozen_path\texpected_sha256\tobserved_sha256\tproducer_mode\tslurm_job_id\n' > "${out_dir}/input_provenance.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${gwas_id}" "${frozen_path}" "${expected_sha}" "${observed_sha}" "${mode}" "${SLURM_JOB_ID}" \
  >> "${out_dir}/input_provenance.tsv"

population="${gwas_id%%_*}"

# --- binary-window sldsc ---

if [[ "${mode}" == "binary_sldsc" ]]; then
  sumstats="${out_dir}/${gwas_id}.sumstats.gz"
  python3 "${script_dir}/../sldsc/05_format_gwas.py" "${frozen_path}" "${sumstats}"

  for annot in prenatal_apa_qtl postnatal_apa_qtl interaction_apa_qtl lifespan_apa_qtl; do
    "${ldsc_python}" "${ldsc_dir}/ldsc.py" \
      --h2 "${sumstats}" \
      --ref-ld-chr "${sldsc_ref}/${population}_hg38/baselineLD_v2.2/baselineLD.,${run_dir}/ld_scores/${population}/${annot}." \
      --w-ld-chr "${sldsc_ref}/${population}_hg38/weights/weights.hm3_noMHC." \
      --frqfile-chr "${sldsc_ref}/${population}_hg38/plink_files/1000G.${population}.hg38." \
      --overlap-annot --print-coefficients \
      --out "${out_dir}/${gwas_id}_${annot}"
  done
  [[ "$(find "${out_dir}" -maxdepth 1 -name '*.results' | wc -l)" -eq 4 ]]
fi

# --- joint apa-qtl and eqtl maxcpp sldsc ---

if [[ "${mode}" == "joint_sldsc" ]]; then
  munged="${out_dir}/${gwas_id}"
  "${ldsc_python}" "${ldsc_dir}/munge_sumstats.py" \
    --sumstats "${frozen_path}" --snp SNPID --N-col N --a1 EA --a2 NEA --p P \
    --signed-sumstats BETA,0 \
    --merge-alleles "${run_dir}/joint_maxcpp/merge_alleles/${population}.hm3_merge_alleles.tsv" \
    --out "${munged}"

  for stage in prenatal postnatal; do
    stage_dir="${out_dir}/${stage}"
    mkdir -p "${stage_dir}"
    "${ldsc_python}" "${ldsc_dir}/ldsc.py" \
      --h2 "${munged}.sumstats.gz" \
      --ref-ld-chr "${sldsc_ref}/${population}_hg38/baselineLD_v2.2/baselineLD.,${run_dir}/joint_maxcpp/ldscores/${population}/${stage}_maxCPP." \
      --w-ld-chr "${sldsc_ref}/${population}_hg38/weights/weights.hm3_noMHC." \
      --frqfile-chr "${sldsc_ref}/${population}_hg38/plink_files/1000G.${population}.hg38." \
      --overlap-annot --print-coefficients \
      --out "${stage_dir}/${gwas_id}"
  done
  [[ "$(find "${out_dir}" -name '*.results' | wc -l)" -eq 2 ]]
fi

# --- magma default and prenatal conditioned on postnatal ---

if [[ "${mode}" == "magma" ]]; then
  pval_table="${out_dir}/${gwas_id}.pval.tsv"
  gzip -cd "${frozen_path}" | awk -F '\t' -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      if (!("SNPID" in col) || !("P" in col) || !("N" in col)) exit 2
      print "SNPID", "P", "N"
      next
    }
    $col["SNPID"] != "" && $col["P"] != "" && $col["N"] != "" {
      print $col["SNPID"], $col["P"], $col["N"]
    }
  ' > "${pval_table}"
  [[ "$(wc -l < "${pval_table}")" -gt 1 ]]

  baseline="${out_dir}/baseline"
  "${magma}" \
    --bfile "${run_dir}/magma/work/${population}/${population}_hg38_genome" \
    --pval "${pval_table}" use=SNPID,P ncol=N \
    --gene-annot "${run_dir}/annotations/${population}/baseline.genes.annot" \
    --out "${baseline}"

  "${magma}" \
    --gene-results "${baseline}.genes.raw" \
    --set-annot "${run_dir}/apa_qtl_gene_sets.sets" \
    --out "${out_dir}/apa_qtl_sets"

  "${magma}" \
    --gene-results "${baseline}.genes.raw" \
    --set-annot "${root}/data/07_coloc/eur_scz_pre_post.sets" \
    --model condition=postnatal \
    --out "${out_dir}/prenatal_condition_postnatal"

  [[ -s "${out_dir}/apa_qtl_sets.gsa.out" && -s "${out_dir}/prenatal_condition_postnatal.gsa.out" ]]
fi

# --- smr and heidi ---

if [[ "${mode}" == "smr" ]]; then
  smr="${SMR:-smr}"
  ld_ref="${run_dir}/refs/all_hg38_bed_uniqueid"
  keep_file="${APA_HPC_HOME}/project/.ref/1KG/hg38/samples/${population}.keep"
  ma_dir="${out_dir}/gwas_ma"
  neff_dir="${out_dir}/gwas_neff"
  smr_dir="${out_dir}/smr_output"
  mkdir -p "${ma_dir}" "${neff_dir}" "${smr_dir}"
  smr_keep_file="${out_dir}/smr_keep.tsv"
  awk '{print "0\t" $1}' "${keep_file}" > "${smr_keep_file}"
  [[ -s "${smr_keep_file}" ]]

  "${smr}" --version > "${out_dir}/smr_version.txt" 2>&1 || true
  python3 "${script_dir}/../../07_coloc/smr/02_harmony_to_ma.py" \
    --harmony-dir "${APA_HPC_HOME}/project/.ref/GWAS" \
    --out-dir "${ma_dir}" --include-prefix "${gwas_id}"
  ma_file="${ma_dir}/${gwas_id}.ma"
  [[ -s "${ma_file}" ]]

  awk 'NR > 1 {print $1}' "${ma_file}" > "${neff_dir}/${gwas_id}.snplist.tmp"
  plink2 --bfile "${ld_ref}" --keep "${smr_keep_file}" \
    --extract "${neff_dir}/${gwas_id}.snplist.tmp" \
    --indep-pairwise 50 5 0.2 --threads "${SLURM_CPUS_PER_TASK}" --memory 90000 \
    --out "${neff_dir}/${gwas_id}"
  n_gwas_rows="$(wc -l < "${neff_dir}/${gwas_id}.snplist.tmp")"
  n_effective="$(wc -l < "${neff_dir}/${gwas_id}.prune.in")"
  n_removed="$(wc -l < "${neff_dir}/${gwas_id}.prune.out")"
  n_ld_overlap=$((n_effective + n_removed))
  gwas_threshold="$(awk -v n="${n_effective}" 'BEGIN {printf "%.12g", 0.05 / n}')"
  printf 'gwas_id\tpopulation\tn_gwas_rows\tn_ld_overlap\tn_effective_snp\tgwas_alpha\tgwas_threshold\tprune_window_variants\tprune_step_variants\tprune_r2\n' \
    > "${neff_dir}/${gwas_id}.metrics.tsv"
  printf '%s\t%s\t%s\t%s\t%s\t0.05\t%s\t50\t5\t0.2\n' \
    "${gwas_id}" "${population}" "${n_gwas_rows}" "${n_ld_overlap}" "${n_effective}" "${gwas_threshold}" \
    >> "${neff_dir}/${gwas_id}.metrics.tsv"

  for qtl_type in prenatal postnatal interaction; do
    cmd=(
      "${smr}"
      --bfile "${ld_ref}"
      --keep "${smr_keep_file}"
      --gwas-summary "${ma_file}"
      --beqtl-summary "${root}/06_smr/results/besd/${qtl_type}"
      --smr-multi --heidi-mtd 1 --heidi-min-m 3 --heidi-max-m 20
      --cis-wind 2000 --peqtl-smr 1e-3 --peqtl-heidi 1.5e-3
      --out "${smr_dir}/${qtl_type}_${gwas_id}"
      --thread-num "${SLURM_CPUS_PER_TASK}"
    )
    if [[ "${population}" == "EUR" ]]; then
      cmd+=(--diff-freq-prop 1.0)
    fi
    "${cmd[@]}"
  done
  [[ "$(find "${smr_dir}" -maxdepth 1 -name '*.msmr' | wc -l)" -eq 3 ]]
fi

# --- twas ---

if [[ "${mode}" == "twas" ]]; then
  tmp_gwas="/tmp/apa_twas_${SLURM_JOB_ID}_${gwas_id}_${twas_stage}.tsv"
  trap 'rm -f "${tmp_gwas}"' EXIT
  gzip -cd "${frozen_path}" > "${tmp_gwas}"
  python3 "${script_dir}/03_prepare_gwas_and_run_spredixcan.py" \
    --stage "${twas_stage}" \
    --gwas-id "${gwas_id}" \
    --gwas-source "${tmp_gwas}" \
    --model-dir "${run_dir}/twas/models" \
    --out-dir "${run_dir}/twas/spredixcan" \
    --format-gwas --run-spredixcan
  [[ -s "${run_dir}/twas/spredixcan/${twas_stage}/${gwas_id}.csv.gz" ]]
  [[ -s "${run_dir}/twas/qc/by_run/${twas_stage}/${gwas_id}.gwas_alignment_audit.tsv" ]]
fi

printf 'COMPLETE\t%s\t%s\t%s\n' "${mode}" "${gwas_id}" "$(date -Is)" > "${out_dir}/COMPLETE.tsv"
