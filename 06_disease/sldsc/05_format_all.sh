#!/bin/bash
# Format the harmonized EUR and EAS GWAS tables to sLDSC sumstats
: "${APA_ROOT:?source config/paths.local.sh}"
set -euo pipefail

# --- paths ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${APA_GROUP_SHARE}/.ref/GWAS/harmony"
OUTPUT_DIR="${APA_ROOT}/data/06_disease"

# --- traits ---
EUR_TRAITS=("MDD" "SCZ" "BIP" "ADHD" "PTSD" "AD" "ASD" "ANX" "INS" "INT" "PD" "EP")
EAS_TRAITS=("SCZ" "BIP" "AD" "PD" "AIS" "ALS" "EP")
mkdir -p "${OUTPUT_DIR}"

# --- EUR traits ---
for TRAIT in "${EUR_TRAITS[@]}"; do
    INPUT_FILE="${INPUT_DIR}/EUR_${TRAIT}_b19to38_harmonized.tsv"
    OUTPUT_FILE="${OUTPUT_DIR}/EUR_${TRAIT}.sumstats.gz"
    # optional traits: skip the ones whose harmonized file is absent
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "WARNING: Input file not found: ${INPUT_FILE}; skipping EUR_${TRAIT}"
        continue
    fi
    python "${SCRIPT_DIR}/05_format_gwas.py" "${INPUT_FILE}" "${OUTPUT_FILE}"
done

# --- EAS traits (AD is already build 38) ---
for TRAIT in "${EAS_TRAITS[@]}"; do
    if [[ "${TRAIT}" == "AD" ]]; then
        INPUT_FILE="${INPUT_DIR}/EAS_AD_b38to38_harmonized.tsv"
    else
        INPUT_FILE="${INPUT_DIR}/EAS_${TRAIT}_b19to38_harmonized.tsv"
    fi
    OUTPUT_FILE="${OUTPUT_DIR}/EAS_${TRAIT}.sumstats.gz"
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "WARNING: Input file not found: ${INPUT_FILE}; skipping EAS_${TRAIT}"
        continue
    fi
    python "${SCRIPT_DIR}/05_format_gwas.py" "${INPUT_FILE}" "${OUTPUT_FILE}"
done
