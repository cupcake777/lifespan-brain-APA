#!/bin/bash
# Compute stratified LD scores for each QTL annotation (ldsc3 on the cluster)
: "${APA_ROOT:?source config/paths.local.sh}"
set -euo pipefail

# --- paths ---
TOOL_DIR="${APA_HPC_HOME}/.local/share/mamba/envs/ldsc3/bin"
LDSC_ENV="ldsc3"
EUR_BFILE="${APA_GROUP_SHARE}/.ref/sLDSC/EUR_hg38/plink_files/1000G.EUR.hg38"
EAS_BFILE="${APA_GROUP_SHARE}/.ref/sLDSC/EAS_hg38/plink_files/1000G.EAS.hg38"
HM3="${APA_GROUP_SHARE}/.ref/sLDSC/hm3_no_MHC.list.txt"
ANNOT_DIR="${APA_ROOT}/data/06_disease/annot"
LDSCORE_DIR="${APA_ROOT}/data/06_disease/ld_scores"

# --- annotation sets ---
ANNOT_NAMES=("prenatal_apa_qtl" "postnatal_apa_qtl" "interaction_apa_qtl" "lifespan_apa_qtl")
mkdir -p "${LDSCORE_DIR}/EUR" "${LDSCORE_DIR}/EAS"

for ANNOT in "${ANNOT_NAMES[@]}"; do
    for CHR in $(seq 1 22); do
        # EUR
        ln -sf "${ANNOT_DIR}/EUR/${ANNOT}.${CHR}.annot.gz" "${LDSCORE_DIR}/EUR/${ANNOT}.${CHR}.annot.gz"
        micromamba run -n "${LDSC_ENV}" python "${TOOL_DIR}/ldsc.py" \
            --l2 \
            --bfile "${EUR_BFILE}.${CHR}" \
            --ld-wind-cm 1.0 \
            --annot "${ANNOT_DIR}/EUR/${ANNOT}.${CHR}.annot.gz" \
            --out "${LDSCORE_DIR}/EUR/${ANNOT}.${CHR}" \
            --print-snps "${HM3}"

        # EAS
        ln -sf "${ANNOT_DIR}/EAS/${ANNOT}.${CHR}.annot.gz" "${LDSCORE_DIR}/EAS/${ANNOT}.${CHR}.annot.gz"
        micromamba run -n "${LDSC_ENV}" python "${TOOL_DIR}/ldsc.py" \
            --l2 \
            --bfile "${EAS_BFILE}.${CHR}" \
            --ld-wind-cm 1.0 \
            --annot "${ANNOT_DIR}/EAS/${ANNOT}.${CHR}.annot.gz" \
            --out "${LDSCORE_DIR}/EAS/${ANNOT}.${CHR}" \
            --print-snps "${HM3}"
    done
done
