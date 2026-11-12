#!/bin/bash
#SBATCH --job-name=sldsc_pipeline
#SBATCH --output=logs/sldsc_pipeline_%j.out
#SBATCH --error=logs/sldsc_pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=cu
#SBATCH --mail-type=FAIL,END
# sLDSC phases 3-7: annotation, LD scores, GWAS sumstats, regression, summary.
# Submit from $APA_ROOT (the #SBATCH log paths are relative to the submit directory).

set -euo pipefail
: "${APA_ROOT:?source config/paths.sh}" "${APA_HPC_HOME:?source config/paths.sh}"

# --- config ---
SCRIPT_DIR="${SCRIPT_DIR:-${APA_ROOT}/release/06_disease/sldsc}"
LDSC_ENV="ldsc3"
LDSC_BIN="${LDSC_BIN:-${APA_HPC_HOME}/.local/share/mamba/envs/ldsc3/bin/ldsc.py}"
REF="${APA_GROUP_SHARE}/.ref/sLDSC"
EUR_BFILE="${EUR_BFILE:-${REF}/EUR_hg38/plink_files/1000G.EUR.hg38}"
EAS_BFILE="${EAS_BFILE:-${REF}/EAS_hg38/plink_files/1000G.EAS.hg38}"
HM3="${HM3:-${REF}/hm3_no_MHC.list.txt}"
EUR_BASELINE="${EUR_BASELINE:-${REF}/EUR_hg38/baselineLD_v2.2/baselineLD.}"
EAS_BASELINE="${EAS_BASELINE:-${REF}/EAS_hg38/baselineLD_v2.2/baselineLD.}"
EUR_WEIGHTS="${EUR_WEIGHTS:-${REF}/EUR_hg38/weights/weights.hm3_noMHC.}"
EAS_WEIGHTS="${EAS_WEIGHTS:-${REF}/EAS_hg38/weights/weights.hm3_noMHC.}"
EUR_FRQFILE="${EUR_FRQFILE:-${REF}/EUR_hg38/plink_files/1000G.EUR.hg38.}"
EAS_FRQFILE="${EAS_FRQFILE:-${REF}/EAS_hg38/plink_files/1000G.EAS.hg38.}"
GWAS_HARM_DIR="${GWAS_HARM_DIR:-${APA_GROUP_SHARE}/.ref/GWAS/harmony}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RUN_PLOTS="${RUN_PLOTS:-1}"
CHROMOSOMES="${CHROMOSOMES:-$(seq 1 22)}"
GWAS_DIR="${APA_ROOT}/data/06_disease/gwas"
BED_DIR="${APA_ROOT}/data/06_disease/qtl_annot"
ANNOT_DIR="${APA_ROOT}/data/06_disease/annot"
LDSCORE_DIR="${APA_ROOT}/data/06_disease/ld_scores"
REG_DIR="${APA_ROOT}/data/06_disease/regression"
SUMMARY_DIR="${APA_ROOT}/data/06_disease"

# Annotation categories
ANNOT_NAMES=("prenatal_apa_qtl" "postnatal_apa_qtl" "interaction_apa_qtl" "lifespan_apa_qtl")

# EUR and EAS traits
read -r -a EUR_TRAITS <<< "${EUR_TRAITS_OVERRIDE:-MDD SCZ BIP ADHD PTSD AD ASD ANX INS INT PD EP}"
read -r -a EAS_TRAITS <<< "${EAS_TRAITS_OVERRIDE:-SCZ BIP AD PD AIS ALS EP}"

if [[ -n "${LDSC_PYTHON:-}" ]]; then
    LDSC_CMD=("${LDSC_PYTHON}")
else
    LDSC_CMD=(micromamba run -n "${LDSC_ENV}" python)
fi

mkdir -p "${ANNOT_DIR}/EUR" "${ANNOT_DIR}/EAS"
mkdir -p "${LDSCORE_DIR}/EUR" "${LDSCORE_DIR}/EAS"
mkdir -p "${REG_DIR}" "${SUMMARY_DIR}" "${GWAS_DIR}" "${APA_ROOT}/logs"

# --- phase 3: annotation files (.annot.gz) ---
for ANNOT in "${ANNOT_NAMES[@]}"; do
    for CHR in ${CHROMOSOMES}; do
        EUR_OUT="${ANNOT_DIR}/EUR/${ANNOT}.${CHR}.annot.gz"
        if [[ ! -f "${EUR_OUT}" ]]; then
            "${PYTHON_BIN}" "${SCRIPT_DIR}/03_annot.py" "${BED_DIR}/${ANNOT}.bed" "${EUR_BFILE}.${CHR}.bim" "${ANNOT}" "${EUR_OUT}"
        fi
        EAS_OUT="${ANNOT_DIR}/EAS/${ANNOT}.${CHR}.annot.gz"
        if [[ ! -f "${EAS_OUT}" ]]; then
            "${PYTHON_BIN}" "${SCRIPT_DIR}/03_annot.py" "${BED_DIR}/${ANNOT}.bed" "${EAS_BFILE}.${CHR}.bim" "${ANNOT}" "${EAS_OUT}"
        fi
    done
done

# --- phase 4: stratified LD scores ---
for ANNOT in "${ANNOT_NAMES[@]}"; do
    for CHR in ${CHROMOSOMES}; do
        EUR_OUT="${LDSCORE_DIR}/EUR/${ANNOT}.${CHR}"
        ln -sf "${ANNOT_DIR}/EUR/${ANNOT}.${CHR}.annot.gz" "${LDSCORE_DIR}/EUR/${ANNOT}.${CHR}.annot.gz"
        if [[ ! -f "${EUR_OUT}.l2.ldscore.gz" ]]; then
            "${LDSC_CMD[@]}" "${LDSC_BIN}" \
                --l2 \
                --bfile "${EUR_BFILE}.${CHR}" \
                --ld-wind-cm 1.0 \
                --annot "${ANNOT_DIR}/EUR/${ANNOT}.${CHR}.annot.gz" \
                --out "${EUR_OUT}" \
                --print-snps "${HM3}"
        fi

        EAS_OUT="${LDSCORE_DIR}/EAS/${ANNOT}.${CHR}"
        ln -sf "${ANNOT_DIR}/EAS/${ANNOT}.${CHR}.annot.gz" "${LDSCORE_DIR}/EAS/${ANNOT}.${CHR}.annot.gz"
        if [[ ! -f "${EAS_OUT}.l2.ldscore.gz" ]]; then
            "${LDSC_CMD[@]}" "${LDSC_BIN}" \
                --l2 \
                --bfile "${EAS_BFILE}.${CHR}" \
                --ld-wind-cm 1.0 \
                --annot "${ANNOT_DIR}/EAS/${ANNOT}.${CHR}.annot.gz" \
                --out "${EAS_OUT}" \
                --print-snps "${HM3}"
        fi
    done
done

# --- phase 5: GWAS sumstats ---
for TRAIT in "${EUR_TRAITS[@]}"; do
    OUTPUT_FILE="${GWAS_DIR}/EUR_${TRAIT}.sumstats.gz"
    if [[ -f "${OUTPUT_FILE}" ]]; then
        continue
    fi
    INPUT_FILE="${GWAS_HARM_DIR}/EUR_${TRAIT}_b19to38_harmonized.tsv"
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "WARNING: EUR_${TRAIT} input not found at ${INPUT_FILE}, skipping"
        continue
    fi
    "${PYTHON_BIN}" "${SCRIPT_DIR}/05_format_gwas.py" "${INPUT_FILE}" "${OUTPUT_FILE}"
done

for TRAIT in "${EAS_TRAITS[@]}"; do
    OUTPUT_FILE="${GWAS_DIR}/EAS_${TRAIT}.sumstats.gz"
    if [[ -f "${OUTPUT_FILE}" ]]; then
        continue
    fi
    if [[ "${TRAIT}" == "AD" ]]; then
        INPUT_FILE="${GWAS_HARM_DIR}/EAS_AD_b38to38_harmonized.tsv"
    else
        INPUT_FILE="${GWAS_HARM_DIR}/EAS_${TRAIT}_b19to38_harmonized.tsv"
    fi
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "WARNING: EAS_${TRAIT} input not found at ${INPUT_FILE}, skipping"
        continue
    fi
    "${PYTHON_BIN}" "${SCRIPT_DIR}/05_format_gwas.py" "${INPUT_FILE}" "${OUTPUT_FILE}"
done

# --- phase 6: stratified LD score regression ---
for ANNOT in "${ANNOT_NAMES[@]}"; do
    for TRAIT in "${EUR_TRAITS[@]}"; do
        OUTPUT="${REG_DIR}/EUR_${TRAIT}_${ANNOT}"
        if [[ ! -f "${GWAS_DIR}/EUR_${TRAIT}.sumstats.gz" ]]; then
            echo "WARNING: EUR_${TRAIT} sumstats not found, skipping"
            continue
        fi
        if [[ -s "${OUTPUT}.results" ]]; then
            continue
        fi
        "${LDSC_CMD[@]}" "${LDSC_BIN}" \
            --h2 "${GWAS_DIR}/EUR_${TRAIT}.sumstats.gz" \
            --ref-ld-chr "${EUR_BASELINE},${LDSCORE_DIR}/EUR/${ANNOT}." \
            --w-ld-chr "${EUR_WEIGHTS}" \
            --frqfile-chr "${EUR_FRQFILE}" \
            --overlap-annot \
            --out "${OUTPUT}" \
            --print-coefficients
    done

    for TRAIT in "${EAS_TRAITS[@]}"; do
        OUTPUT="${REG_DIR}/EAS_${TRAIT}_${ANNOT}"
        if [[ ! -f "${GWAS_DIR}/EAS_${TRAIT}.sumstats.gz" ]]; then
            echo "WARNING: EAS_${TRAIT} sumstats not found, skipping"
            continue
        fi
        if [[ -s "${OUTPUT}.results" ]]; then
            continue
        fi
        "${LDSC_CMD[@]}" "${LDSC_BIN}" \
            --h2 "${GWAS_DIR}/EAS_${TRAIT}.sumstats.gz" \
            --ref-ld-chr "${EAS_BASELINE},${LDSCORE_DIR}/EAS/${ANNOT}." \
            --w-ld-chr "${EAS_WEIGHTS}" \
            --frqfile-chr "${EAS_FRQFILE}" \
            --overlap-annot \
            --out "${OUTPUT}" \
            --print-coefficients
    done
done

# --- phase 7: collect and plot ---
"${LDSC_CMD[@]}" "${SCRIPT_DIR}/07_collect.py" "${REG_DIR}" "${SUMMARY_DIR}/sldsc_enrichment_summary.tsv"
if [[ "${RUN_PLOTS}" == "1" ]]; then
    "${LDSC_CMD[@]}" "${SCRIPT_DIR}/08_plot_sldsc.py"
fi
