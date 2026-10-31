#!/usr/bin/env bash
# Run qtl_pipe.sh for each downsample seed and stage (module 04).
# Outputs go to data/04_qtl/downsample/seed_<seed>/qtl/<stage>/.
# Completed outputs are skipped when SKIP_COMPLETED=1.

set -euo pipefail

: "${APA_ROOT:?set APA_ROOT (source config/paths.sh)}"

# --- settings ---

SEEDS="${SEEDS:-2026 2027 2028 2029 2030 2031 2032 2033 2034 2035 2036 2037 2038 2039 2040 2041 2042 2043 2044 2045}"
SEEDS="${SEEDS//,/ }"
BASE_DIR="${BASE_DIR:-${APA_ROOT}/data/04_qtl/downsample}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORE_QTL_PIPE="${CORE_QTL_PIPE:-${HERE}/qtl_pipe.sh}"
RUN_NOMINAL="${RUN_NOMINAL:-1}"
RUN_PERMUTATION="${RUN_PERMUTATION:-1}"
RUN_SUSIE="${RUN_SUSIE:-1}"
SKIP_COMPLETED="${SKIP_COMPLETED:-1}"
STAGES="${STAGES:-}"
STAGES="${STAGES//,/ }"

# --- run each seed and stage ---

for seed in ${SEEDS}; do
  pre_dir="${BASE_DIR}/seed_${seed}/pre_input"
  if [[ ! -d "${pre_dir}" ]]; then
    echo "[WARN] skip seed ${seed}: missing ${pre_dir}" >&2
    continue
  fi

  shopt -s nullglob
  phenotype_files=("${pre_dir}"/phenotype.*.GeneRange.bed.gz)
  shopt -u nullglob
  if [[ ${#phenotype_files[@]} -eq 0 ]]; then
    echo "[WARN] skip seed ${seed}: no phenotype.*.GeneRange.bed.gz in ${pre_dir}" >&2
    continue
  fi

  for phenotype_file in "${phenotype_files[@]}"; do
    file_name="$(basename "${phenotype_file}")"
    stage="${file_name#phenotype.}"
    stage="${stage%.GeneRange.bed.gz}"

    stage_enabled=1
    if [[ -n "${STAGES}" ]]; then
      stage_enabled=0
      for requested in ${STAGES}; do
        if [[ "${requested}" == "${stage}" ]]; then
          stage_enabled=1
        fi
      done
    fi
    if [[ "${stage_enabled}" == "0" ]]; then
      echo "[INFO] skip stage ${stage} for seed ${seed}: STAGES=${STAGES}"
      continue
    fi

    stage_lower="$(printf '%s' "${stage}" | tr '[:upper:]' '[:lower:]')"
    qtl_dir="${BASE_DIR}/seed_${seed}/qtl/${stage}"

    export GROUP_MAPPING=true
    export GENO_PREFIX="${pre_dir}/genotype.${stage_lower}"
    export PHENO_BED="${pre_dir}/phenotype.${stage}.GeneRange.bed.gz"
    export COV_FILE="${pre_dir}/covariates_for_qtl.${stage}.txt"
    export GROUP_FILE="${pre_dir}/phenotype.${stage}.group"
    export RUN_PREFIX="seed_${seed}_${stage_lower}"
    export NOMINAL_DIR="${qtl_dir}/nominal"
    export PERMUTATION_DIR="${qtl_dir}/permutation"
    export SUSIE_DIR="${qtl_dir}/susie"
    export SEED="${seed}"

    echo "[INFO] seed=${seed} stage=${stage}"

    if [[ "${RUN_NOMINAL}" == "1" ]]; then
      shopt -s nullglob
      nominal_files=("${NOMINAL_DIR}/${RUN_PREFIX}".cis_qtl_pairs.*.parquet)
      shopt -u nullglob
      if [[ "${SKIP_COMPLETED}" == "1" && ${#nominal_files[@]} -ge 22 ]]; then
        echo "[INFO] skip completed nominal: ${RUN_PREFIX}"
      else
        export RUN_QTL_MODE=nominal
        bash "${CORE_QTL_PIPE}"
      fi
    fi

    if [[ "${RUN_PERMUTATION}" == "1" ]]; then
      if [[ "${SKIP_COMPLETED}" == "1" && -s "${PERMUTATION_DIR}/${RUN_PREFIX}.cis_qtl.txt.gz" ]]; then
        echo "[INFO] skip completed permutation: ${RUN_PREFIX}"
      else
        export RUN_QTL_MODE=permutation
        bash "${CORE_QTL_PIPE}"
      fi
    fi

    if [[ "${RUN_SUSIE}" == "1" ]]; then
      if [[ "${SKIP_COMPLETED}" == "1" && -s "${SUSIE_DIR}/${RUN_PREFIX}.SuSiE_summary.parquet" ]]; then
        echo "[INFO] skip completed susie: ${RUN_PREFIX}"
      else
        export RUN_QTL_MODE=susie
        bash "${CORE_QTL_PIPE}"
      fi
    fi
  done
done
