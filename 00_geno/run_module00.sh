#!/usr/bin/env bash
# Submit the module 00 Slurm chain: joint genotyping, merge, VQSR, pre-QC, imputation prep.
set -euo pipefail

# --- submit chain ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

joint_job=$(sbatch --parsable "${SCRIPT_DIR}/01_joint_geno.sh")
merge_job=$(sbatch --parsable --dependency="afterok:${joint_job}" "${SCRIPT_DIR}/02_merge_clean.sh")
vqsr_job=$(sbatch --parsable --dependency="afterok:${merge_job}" "${SCRIPT_DIR}/03_VQSR.sh")
qc_job=$(sbatch --parsable --dependency="afterok:${vqsr_job}" "${SCRIPT_DIR}/04_preQC.sh")
upload_job=$(sbatch --parsable --dependency="afterok:${qc_job}" "${SCRIPT_DIR}/05_prep_impute.sh")

printf 'joint=%s\nmerge=%s\nvqsr=%s\npre_qc=%s\nupload_prep=%s\n' \
    "${joint_job}" "${merge_job}" "${vqsr_job}" "${qc_job}" "${upload_job}"
echo 'ChinaMAP imputation is a web checkpoint after upload_prep.'
echo 'After the imputed VCF is downloaded, run 00_4_stage_filter.sh with RAW_VCF set to it.'
