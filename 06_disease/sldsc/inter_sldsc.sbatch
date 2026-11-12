#!/usr/bin/env bash
#SBATCH --job-name=interaction_sldsc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=CU
#SBATCH --output=interaction_sldsc_%j.out
#SBATCH --error=interaction_sldsc_%j.err
#SBATCH --mail-type=FAIL,END
: "${APA_ROOT:?source config/paths.sh}" "${APA_HPC_HOME:?source config/paths.sh}"
# sLDSC for the APA-QTL interaction annotation: BED, annot, LD scores, regression.
# Submit from $APA_ROOT (the #SBATCH log paths are relative to the submit directory).

set -euo pipefail

# --- config ---
source_dir="${APA_ROOT}/release/06_disease/sldsc"
python_bin="${APA_HPC_HOME}/.local/share/mamba/envs/ldsc3/bin/python"
ldsc_bin="${APA_HPC_HOME}/.local/share/mamba/envs/ldsc3/bin/ldsc.py"
ref_dir="${APA_GROUP_SHARE}/.ref/sLDSC"
eur_bfile="${ref_dir}/EUR_hg38/plink_files/1000G.EUR.hg38"
eas_bfile="${ref_dir}/EAS_hg38/plink_files/1000G.EAS.hg38"
hm3="${ref_dir}/hm3_no_MHC.list.txt"
eur_baseline="${ref_dir}/EUR_hg38/baselineLD_v2.2/baselineLD."
eas_baseline="${ref_dir}/EAS_hg38/baselineLD_v2.2/baselineLD."
eur_weights="${ref_dir}/EUR_hg38/weights/weights.hm3_noMHC."
eas_weights="${ref_dir}/EAS_hg38/weights/weights.hm3_noMHC."
eur_frqfile="${ref_dir}/EUR_hg38/plink_files/1000G.EUR.hg38."
eas_frqfile="${ref_dir}/EAS_hg38/plink_files/1000G.EAS.hg38."
gwas_dir="${APA_ROOT}/data/06_disease/gwas"
pair_file="${APA_ROOT}/data/07_coloc/intx.signif_pairs.txt.gz"
stage_dir="${APA_ROOT}/data/06_disease/interaction_sldsc"
bed_file="${stage_dir}/qtl_annot/interaction_apa_qtl.bed"
annot_dir="${stage_dir}/annot"
ldscore_dir="${stage_dir}/ld_scores"
reg_dir="${stage_dir}/regression"
eur_traits=(MDD SCZ BIP ADHD PTSD AD ASD ANX INS INT PD EP)
eas_traits=(SCZ BIP AD PD AIS ALS EP)
workers="${SLURM_CPUS_PER_TASK:-4}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# --- interaction annotation ---
mkdir -p "${stage_dir}/qtl_annot" "${annot_dir}/EUR" "${annot_dir}/EAS" \
  "${ldscore_dir}/EUR" "${ldscore_dir}/EAS" "${reg_dir}"
"${python_bin}" "${source_dir}/02_annot_bed.py" "${pair_file}" "${bed_file}" 500

run_population_chr() {
  local population="$1" chr="$2" bfile annot ldscore
  if [[ "${population}" == EUR ]]; then
    bfile="${eur_bfile}"
  else
    bfile="${eas_bfile}"
  fi
  annot="${annot_dir}/${population}/interaction_apa_qtl.${chr}.annot.gz"
  ldscore="${ldscore_dir}/${population}/interaction_apa_qtl.${chr}"
  if [[ ! -s "${annot}" ]]; then
    "${python_bin}" "${source_dir}/03_annot.py" "${bed_file}" "${bfile}.${chr}.bim" interaction_apa_qtl "${annot}"
  fi
  ln -sf "${annot}" "${ldscore}.annot.gz"
  if [[ ! -s "${ldscore}.l2.ldscore.gz" ]]; then
    "${python_bin}" "${ldsc_bin}" --l2 --bfile "${bfile}.${chr}" \
      --ld-wind-cm 1.0 --annot "${annot}" --out "${ldscore}" --print-snps "${hm3}"
  fi
}

export -f run_population_chr
export python_bin source_dir eur_bfile eas_bfile hm3 annot_dir ldscore_dir bed_file
{
  for population in EUR EAS; do
    for chr in $(seq 1 22); do
      printf '%s\t%s\n' "${population}" "${chr}"
    done
  done
} | xargs -r -P "${workers}" -n 2 bash -c 'set -euo pipefail; run_population_chr "$@"' _

# --- interaction sLDSC regressions ---
run_sldsc_regression() {
  local population="$1" trait="$2" output baseline weights frqfile
  if [[ "${population}" == EUR ]]; then
    baseline="${eur_baseline}"
    weights="${eur_weights}"
    frqfile="${eur_frqfile}"
  else
    baseline="${eas_baseline}"
    weights="${eas_weights}"
    frqfile="${eas_frqfile}"
  fi
  output="${reg_dir}/${population}_${trait}_interaction_apa_qtl"
  if [[ -s "${output}.results" ]]; then
    return 0
  fi
  "${python_bin}" "${ldsc_bin}" --h2 "${gwas_dir}/${population}_${trait}.sumstats.gz" \
    --ref-ld-chr "${baseline},${ldscore_dir}/${population}/interaction_apa_qtl." \
    --w-ld-chr "${weights}" --frqfile-chr "${frqfile}" \
    --overlap-annot --out "${output}" --print-coefficients
}

export -f run_sldsc_regression
export python_bin ldsc_bin gwas_dir ldscore_dir reg_dir eur_baseline eas_baseline \
  eur_weights eas_weights eur_frqfile eas_frqfile
{
  for trait in "${eur_traits[@]}"; do
    printf 'EUR\t%s\n' "${trait}"
  done
  for trait in "${eas_traits[@]}"; do
    printf 'EAS\t%s\n' "${trait}"
  done
} | xargs -r -P "${workers}" -n 2 bash -c 'set -euo pipefail; run_sldsc_regression "$@"' _
