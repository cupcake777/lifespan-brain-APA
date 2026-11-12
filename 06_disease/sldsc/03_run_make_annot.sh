#!/bin/bash
# Write per-chromosome .annot.gz files for the four QTL window sets, EUR and EAS panels
: "${APA_ROOT:?source config/paths.local.sh}"
set -euo pipefail

# --- paths ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EUR_BFILE="${APA_GROUP_SHARE}/.ref/sLDSC/EUR_hg38/plink_files/1000G.EUR.hg38"
EAS_BFILE="${APA_GROUP_SHARE}/.ref/sLDSC/EAS_hg38/plink_files/1000G.EAS.hg38"
BED_DIR="${APA_ROOT}/data/06_disease"
OUT_DIR="${APA_ROOT}/data/06_disease/annot"

# --- annotation sets ---
ANNOT_NAMES=("prenatal_apa_qtl" "postnatal_apa_qtl" "interaction_apa_qtl" "lifespan_apa_qtl")
mkdir -p "${OUT_DIR}/EUR" "${OUT_DIR}/EAS"

# --- one .annot.gz per annotation and chromosome ---
for ANNOT in "${ANNOT_NAMES[@]}"; do
    for CHR in $(seq 1 22); do
        python3 "${SCRIPT_DIR}/03_annot.py" --bed "${BED_DIR}/${ANNOT}.bed" --bim "${EUR_BFILE}.${CHR}.bim" --annot-name "${ANNOT}" --out "${OUT_DIR}/EUR/${ANNOT}.${CHR}.annot.gz"
        python3 "${SCRIPT_DIR}/03_annot.py" --bed "${BED_DIR}/${ANNOT}.bed" --bim "${EAS_BFILE}.${CHR}.bim" --annot-name "${ANNOT}" --out "${OUT_DIR}/EAS/${ANNOT}.${CHR}.annot.gz"
    done
done
