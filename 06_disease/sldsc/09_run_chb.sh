#!/usr/bin/env bash
# Stage-specific CHB eQTL sLDSC: annotations, LD scores and regressions, then paired with APA-QTL results
# Independent chromosome and population tasks run in parallel, bounded by SLURM_CPUS_PER_TASK
: "${APA_ROOT:?source config/paths.sh}"
: "${APA_HPC_HOME:?source config/paths.sh}"
set -euo pipefail

# --- paths ---
SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EQTL_DIR="${APA_HPC_HOME}/project/work/effect/EXP/stage/data"
LDSC_BIN="${APA_HPC_HOME}/.local/share/mamba/envs/ldsc3/bin/ldsc.py"
REF_ROOT="${APA_GROUP_SHARE}/.ref/sLDSC"
EUR_BFILE="${REF_ROOT}/EUR_hg38/plink_files/1000G.EUR.hg38"
EAS_BFILE="${REF_ROOT}/EAS_hg38/plink_files/1000G.EAS.hg38"
EUR_BASELINE="${REF_ROOT}/EUR_hg38/baselineLD_v2.2/baselineLD."
EAS_BASELINE="${REF_ROOT}/EAS_hg38/baselineLD_v2.2/baselineLD."
EUR_WEIGHTS="${REF_ROOT}/EUR_hg38/weights/weights.hm3_noMHC."
EAS_WEIGHTS="${REF_ROOT}/EAS_hg38/weights/weights.hm3_noMHC."
EUR_FRQFILE="${REF_ROOT}/EUR_hg38/plink_files/1000G.EUR.hg38."
EAS_FRQFILE="${REF_ROOT}/EAS_hg38/plink_files/1000G.EAS.hg38."
HM3="${REF_ROOT}/hm3_no_MHC.list.txt"
MAX_JOBS="${SLURM_CPUS_PER_TASK:-8}"
WINDOW_BP=500
EUR_TRAITS=(MDD SCZ BIP ADHD PTSD AD ASD ANX INS INT PD EP)
EAS_TRAITS=(SCZ BIP AD PD AIS ALS EP)
STAGES=(prenatal postnatal)

BED_DIR="${APA_ROOT}/data/06_disease/qtl_annot"
ANNOT_DIR="${APA_ROOT}/data/06_disease/annot"
LDSCORE_DIR="${APA_ROOT}/data/06_disease/ld_scores"
GWAS_DIR="${APA_ROOT}/data/06_disease/gwas"
REG_DIR="${APA_ROOT}/data/06_disease/regression"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# --- bounded parallel pool ---
PIDS=()
POOL_FAILED=0

wait_oldest() {
    local pid="${PIDS[0]}"
    if ! wait "${pid}"; then
        POOL_FAILED=1
    fi
    PIDS=("${PIDS[@]:1}")
}

launch_bounded() {
    "$@" &
    PIDS+=("$!")
    if (( ${#PIDS[@]} >= MAX_JOBS )); then
        wait_oldest
    fi
}

finish_pool() {
    while (( ${#PIDS[@]} > 0 )); do
        wait_oldest
    done
    if (( POOL_FAILED != 0 )); then
        echo "[ERROR] At least one parallel task failed" >&2
        exit 1
    fi
}

# --- per-job functions ---
make_annotation() {
    local population="$1" annotation="$2" chromosome="$3" bfile="$4"
    local output="${ANNOT_DIR}/${population}/${annotation}.${chromosome}.annot.gz"
    if [[ -s "${output}" ]]; then
        return
    fi
    python "${SDIR}/03_annot.py" \
        --bed "${BED_DIR}/${annotation}.bed" \
        --bim "${bfile}.${chromosome}.bim" \
        --annot-name "${annotation}" \
        --out "${output}.tmp"
    mv "${output}.tmp" "${output}"
}

compute_ld_score() {
    local population="$1" annotation="$2" chromosome="$3" bfile="$4"
    local prefix="${LDSCORE_DIR}/${population}/${annotation}.${chromosome}"
    ln -sfn "${ANNOT_DIR}/${population}/${annotation}.${chromosome}.annot.gz" "${prefix}.annot.gz"
    if [[ -s "${prefix}.l2.ldscore.gz" && -s "${prefix}.l2.M" && -s "${prefix}.l2.M_5_50" ]]; then
        return
    fi
    python "${LDSC_BIN}" \
        --l2 \
        --bfile "${bfile}.${chromosome}" \
        --ld-wind-cm 1.0 \
        --annot "${ANNOT_DIR}/${population}/${annotation}.${chromosome}.annot.gz" \
        --out "${prefix}" \
        --print-snps "${HM3}"
}

run_regression() {
    local population="$1" trait="$2" annotation="$3" baseline="$4" weights="$5" frqfile="$6"
    local output="${REG_DIR}/${population}_${trait}_${annotation}"
    if [[ -s "${output}.results" ]]; then
        return
    fi
    python "${LDSC_BIN}" \
        --h2 "${GWAS_DIR}/${population}_${trait}.sumstats.gz" \
        --ref-ld-chr "${baseline},${LDSCORE_DIR}/${population}/${annotation}." \
        --w-ld-chr "${weights}" \
        --frqfile-chr "${frqfile}" \
        --overlap-annot \
        --out "${output}" \
        --print-coefficients
}

mkdir -p "${BED_DIR}" "${ANNOT_DIR}/EUR" "${ANNOT_DIR}/EAS" \
    "${LDSCORE_DIR}/EUR" "${LDSCORE_DIR}/EAS" "${REG_DIR}"

# --- eQTL windows as BED ---
for stage in "${STAGES[@]}"; do
    input="${EQTL_DIR}/${stage}.sig_QTL.txt.gz"
    python "${SDIR}/02_annot_bed.py" \
        "${input}" "${BED_DIR}/${stage}_eqtl.bed.tmp" "${WINDOW_BP}"
    mv "${BED_DIR}/${stage}_eqtl.bed.tmp" "${BED_DIR}/${stage}_eqtl.bed"
done

# --- annotation files ---
for stage in "${STAGES[@]}"; do
    annotation="${stage}_eqtl"
    for chromosome in $(seq 1 22); do
        launch_bounded make_annotation EUR "${annotation}" "${chromosome}" "${EUR_BFILE}"
        launch_bounded make_annotation EAS "${annotation}" "${chromosome}" "${EAS_BFILE}"
    done
done
finish_pool

# --- LD scores ---
for stage in "${STAGES[@]}"; do
    annotation="${stage}_eqtl"
    for chromosome in $(seq 1 22); do
        launch_bounded compute_ld_score EUR "${annotation}" "${chromosome}" "${EUR_BFILE}"
        launch_bounded compute_ld_score EAS "${annotation}" "${chromosome}" "${EAS_BFILE}"
    done
done
finish_pool

# --- sLDSC regressions ---
for stage in "${STAGES[@]}"; do
    annotation="${stage}_eqtl"
    for trait in "${EUR_TRAITS[@]}"; do
        launch_bounded run_regression EUR "${trait}" "${annotation}" "${EUR_BASELINE}" "${EUR_WEIGHTS}" "${EUR_FRQFILE}"
    done
    for trait in "${EAS_TRAITS[@]}"; do
        launch_bounded run_regression EAS "${trait}" "${annotation}" "${EAS_BASELINE}" "${EAS_WEIGHTS}" "${EAS_FRQFILE}"
    done
done
finish_pool

# --- pair with APA-QTL results ---
python "${SDIR}/07_apa_eqtl_cmp.py"
