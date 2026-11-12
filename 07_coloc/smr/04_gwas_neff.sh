#!/usr/bin/env bash
# Effective SNP count per GWAS: ancestry-matched 1KG LD pruning, one Slurm array task per GWAS.
# Submit from $APA_ROOT: sbatch release/07_coloc/smr/04_gwas_neff.sh
#SBATCH --job-name=smr_neff
#SBATCH --output=logs/neff_%A_%a.out
#SBATCH --error=logs/neff_%A_%a.err
#SBATCH --array=1-19%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=CU

set -euo pipefail
: "${APA_ROOT:?source release/config/paths.local.sh}" "${APA_HPC_HOME:?source release/config/paths.local.sh}"

# --- config ---

GWAS_DIR="${APA_ROOT}/data/07_coloc"
LD_REF="${GWAS_DIR}/refs/plink_hg38_uniqueid_runtime/1kgp_hg38_chrpos_ref_alt"
KEEP_DIR="${APA_HPC_HOME}/project/work/apa/04_coloc/refs/plink_keep_runtime"
OUT_DIR="${APA_ROOT}/data/07_coloc"

WINDOW_VARIANTS=50
STEP_VARIANTS=5
R2_THRESHOLD=0.2
ALPHA=0.05

GWAS_IDS=(
    EAS_AD EAS_AIS EAS_ALS EAS_BIP EAS_EP EAS_PD EAS_SCZ
    EUR_AD EUR_ADHD EUR_ANX EUR_ASD EUR_BIP EUR_EP EUR_INS EUR_INT
    EUR_MDD EUR_PD EUR_PTSD EUR_SCZ
)

# --- select GWAS ---

TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
GWAS_ID="${GWAS_IDS[${TASK_INDEX}]}"
POPULATION="${GWAS_ID%%_*}"
GWAS_MA="${GWAS_DIR}/${GWAS_ID}.ma"
KEEP_FILE="${KEEP_DIR}/${POPULATION}.keep"
PREFIX="${OUT_DIR}/${GWAS_ID}"
SNPLIST="${PREFIX}.snplist.tmp"

mkdir -p "${OUT_DIR}" logs

echo "[$(date)] GWAS=${GWAS_ID} population=${POPULATION}"

# --- GWAS SNP set ---

awk 'NR > 1 {print $1}' "${GWAS_MA}" > "${SNPLIST}"
N_GWAS_ROWS=$(wc -l < "${SNPLIST}")

# --- ancestry-matched LD pruning ---

plink2 \
    --bfile "${LD_REF}" \
    --keep "${KEEP_FILE}" \
    --extract "${SNPLIST}" \
    --indep-pairwise "${WINDOW_VARIANTS}" "${STEP_VARIANTS}" "${R2_THRESHOLD}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --out "${PREFIX}"

N_EFFECTIVE=$(wc -l < "${PREFIX}.prune.in")
N_REMOVED=$(wc -l < "${PREFIX}.prune.out")
N_LD_OVERLAP=$((N_EFFECTIVE + N_REMOVED))
GWAS_THRESHOLD=$(awk -v alpha="${ALPHA}" -v n="${N_EFFECTIVE}" 'BEGIN {printf "%.12g", alpha / n}')

# --- task metrics ---

printf 'gwas_id\tpopulation\tn_gwas_rows\tn_ld_overlap\tn_effective_snp\tgwas_alpha\tgwas_threshold\tprune_window_variants\tprune_step_variants\tprune_r2\n' > "${PREFIX}.metrics.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${GWAS_ID}" "${POPULATION}" "${N_GWAS_ROWS}" "${N_LD_OVERLAP}" \
    "${N_EFFECTIVE}" "${ALPHA}" "${GWAS_THRESHOLD}" "${WINDOW_VARIANTS}" \
    "${STEP_VARIANTS}" "${R2_THRESHOLD}" >> "${PREFIX}.metrics.tsv"

rm -f "${SNPLIST}" "${PREFIX}.prune.out"

echo "[$(date)] N_GWAS=${N_GWAS_ROWS} N_LD_OVERLAP=${N_LD_OVERLAP} N_EFF=${N_EFFECTIVE} GWAS_threshold=${GWAS_THRESHOLD}"
