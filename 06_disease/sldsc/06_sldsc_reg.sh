#!/bin/bash
# Stratified LD score regression for each trait x QTL annotation (skips finished runs)
: "${APA_ROOT:?source config/paths.local.sh}"
: "${APA_HPC_HOME:?source config/paths.local.sh}"
set -euo pipefail

# --- paths ---
TOOL_DIR="${APA_HPC_HOME}/.local/share/mamba/envs/ldsc3/bin"
LDSC_ENV="ldsc3"
REF="${APA_GROUP_SHARE}/.ref/sLDSC"
EUR_BASELINE="${REF}/EUR_hg38/baselineLD_v2.2/baselineLD."
EAS_BASELINE="${REF}/EAS_hg38/baselineLD_v2.2/baselineLD."
EUR_WEIGHTS="${REF}/EUR_hg38/weights/weights.hm3_noMHC."
EAS_WEIGHTS="${REF}/EAS_hg38/weights/weights.hm3_noMHC."
EUR_FRQFILE="${REF}/EUR_hg38/plink_files/1000G.EUR.hg38."
EAS_FRQFILE="${REF}/EAS_hg38/plink_files/1000G.EAS.hg38."

GWAS_DIR="${APA_ROOT}/data/06_disease/gwas"
LDSCORE_DIR="${APA_ROOT}/data/06_disease/ld_scores"
OUT_DIR="${APA_ROOT}/data/06_disease/regression"

# --- traits and annotations ---
EUR_TRAITS=("MDD" "SCZ" "BIP" "ADHD" "PTSD" "AD" "ASD" "ANX" "INS" "INT" "PD" "EP")
EAS_TRAITS=("SCZ" "BIP" "AD" "PD" "AIS" "ALS" "EP")
ANNOT_NAMES=("prenatal_apa_qtl" "postnatal_apa_qtl" "interaction_apa_qtl" "lifespan_apa_qtl")

mkdir -p "${OUT_DIR}"

for ANNOT in "${ANNOT_NAMES[@]}"; do
    # --- EUR traits ---
    for TRAIT in "${EUR_TRAITS[@]}"; do
        if [[ ! -f "${GWAS_DIR}/EUR_${TRAIT}.sumstats.gz" ]]; then
            echo "WARNING: Sumstats file not found: EUR_${TRAIT}.sumstats.gz; skipping"
            continue
        fi
        if [[ -s "${OUT_DIR}/EUR_${TRAIT}_${ANNOT}.results" ]]; then
            continue
        fi
        micromamba run -n "${LDSC_ENV}" python "${TOOL_DIR}/ldsc.py" \
            --h2 "${GWAS_DIR}/EUR_${TRAIT}.sumstats.gz" \
            --ref-ld-chr "${EUR_BASELINE},${LDSCORE_DIR}/EUR/${ANNOT}." \
            --w-ld-chr "${EUR_WEIGHTS}" \
            --frqfile-chr "${EUR_FRQFILE}" \
            --overlap-annot \
            --out "${OUT_DIR}/EUR_${TRAIT}_${ANNOT}" \
            --print-coefficients
    done

    # --- EAS traits ---
    for TRAIT in "${EAS_TRAITS[@]}"; do
        if [[ ! -f "${GWAS_DIR}/EAS_${TRAIT}.sumstats.gz" ]]; then
            echo "WARNING: Sumstats file not found: EAS_${TRAIT}.sumstats.gz; skipping"
            continue
        fi
        if [[ -s "${OUT_DIR}/EAS_${TRAIT}_${ANNOT}.results" ]]; then
            continue
        fi
        micromamba run -n "${LDSC_ENV}" python "${TOOL_DIR}/ldsc.py" \
            --h2 "${GWAS_DIR}/EAS_${TRAIT}.sumstats.gz" \
            --ref-ld-chr "${EAS_BASELINE},${LDSCORE_DIR}/EAS/${ANNOT}." \
            --w-ld-chr "${EAS_WEIGHTS}" \
            --frqfile-chr "${EAS_FRQFILE}" \
            --overlap-annot \
            --out "${OUT_DIR}/EAS_${TRAIT}_${ANNOT}" \
            --print-coefficients
    done
done
