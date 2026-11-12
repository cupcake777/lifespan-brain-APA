#!/usr/bin/env bash
# Run tensorQTL for every stage with a phenotype BED in pre_input (module 04).
# Outputs go to data/04_qtl/{nominal,permutation,susie}/ through qtl_pipe.sh.

set -euo pipefail

: "${APA_ROOT:?set APA_ROOT (source config/paths.sh)}"

# --- settings ---

PRE_DIR="${PRE_DIR:-${APA_ROOT}/data/04_qtl/pre_input}"
GENO_DIR="${GENO_DIR:-${PRE_DIR}}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORE_QTL_PIPE="${CORE_QTL_PIPE:-${HERE}/qtl_pipe.sh}"

# --- run each stage ---

shopt -s nullglob
phenotype_files=("${PRE_DIR}"/phenotype.*.GeneRange.bed.gz)
shopt -u nullglob
if [[ ${#phenotype_files[@]} -eq 0 ]]; then
  echo "[ERROR] no phenotype.*.GeneRange.bed.gz in ${PRE_DIR}; run preprocess.R first" >&2
  exit 1
fi

for phenotype_file in "${phenotype_files[@]}"; do
  file_name="$(basename "${phenotype_file}")"
  stage="${file_name#phenotype.}"
  stage="${stage%.GeneRange.bed.gz}"
  stage_lower="$(printf '%s' "${stage}" | tr '[:upper:]' '[:lower:]')"

  export GROUP_MAPPING=true
  export GENO_PREFIX="${GENO_DIR}/genotype.${stage_lower}"
  export PHENO_BED="${PRE_DIR}/phenotype.${stage}.GeneRange.bed.gz"
  export COV_FILE="${PRE_DIR}/covariates_for_qtl.${stage}.txt"
  export GROUP_FILE="${PRE_DIR}/phenotype.${stage}.group"
  export RUN_PREFIX="${stage_lower}"
  export NOMINAL_DIR="${NOMINAL_DIR:-${APA_ROOT}/data/04_qtl/nominal}"
  export PERMUTATION_DIR="${PERMUTATION_DIR:-${APA_ROOT}/data/04_qtl/permutation}"
  export SUSIE_DIR="${SUSIE_DIR:-${APA_ROOT}/data/04_qtl/susie}"
  export SEED="${SEED:-2026}"

  echo "[INFO] full cohort stage=${stage}"
  bash "${CORE_QTL_PIPE}"
done
