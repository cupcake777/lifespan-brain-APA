#!/usr/bin/env bash
# Run tensorQTL for one dataset: cis nominal, cis permutation, and cis SuSiE.
# Required env vars: GENO_PREFIX (PLINK2 pgen prefix), PHENO_BED, COV_FILE.
# Optional env vars: GROUP_MAPPING, GROUP_FILE, RUN_PREFIX, RUN_QTL_MODE, SEED, MAF_THRESHOLD,
# PYTHON_BIN, TENSORQTL_MODULE, and the NOMINAL_DIR, PERMUTATION_DIR, SUSIE_DIR output dirs.

set -euo pipefail

: "${APA_ROOT:?set APA_ROOT (source config/paths.sh)}"

# --- settings ---

GENO_PREFIX="${GENO_PREFIX:-}"
PHENO_BED="${PHENO_BED:-}"
COV_FILE="${COV_FILE:-}"
GROUP_MAPPING="${GROUP_MAPPING:-false}"
GROUP_FILE="${GROUP_FILE:-}"
RUN_PREFIX="${RUN_PREFIX:-}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
TENSORQTL_MODULE="${TENSORQTL_MODULE:-tensorqtl}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.05}"
SEED="${SEED:-2026}"
NOMINAL_DIR="${NOMINAL_DIR:-${APA_ROOT}/data/04_qtl/nominal}"
PERMUTATION_DIR="${PERMUTATION_DIR:-${APA_ROOT}/data/04_qtl/permutation}"
SUSIE_DIR="${SUSIE_DIR:-${APA_ROOT}/data/04_qtl/susie}"
RUN_QTL_MODE="${RUN_QTL_MODE:-all}"

if [[ -z "${GENO_PREFIX}" || -z "${PHENO_BED}" || -z "${COV_FILE}" ]]; then
  echo "[ERROR] GENO_PREFIX, PHENO_BED and COV_FILE are required" >&2
  exit 1
fi

case "${GROUP_MAPPING,,}" in
  true|1|yes|y) GROUP_MAPPING=true ;;
  false|0|no|n) GROUP_MAPPING=false ;;
  *)
    echo "[ERROR] GROUP_MAPPING must be true or false, got: ${GROUP_MAPPING}" >&2
    exit 1
    ;;
esac

RUN_PREFIX="${RUN_PREFIX:-$(basename "${GENO_PREFIX}")}"

case "${RUN_QTL_MODE}" in
  all|nominal|permutation|susie) ;;
  *)
    echo "[ERROR] RUN_QTL_MODE must be all, nominal, permutation, or susie; got: ${RUN_QTL_MODE}" >&2
    exit 1
    ;;
esac

COMMON_ARGS=(
  --covariates "${COV_FILE}"
  --maf_threshold "${MAF_THRESHOLD}"
  --seed "${SEED}"
)

if [[ "${GROUP_MAPPING}" == "true" ]]; then
  if [[ -z "${GROUP_FILE}" ]]; then
    echo "[ERROR] GROUP_FILE is required when GROUP_MAPPING=true" >&2
    exit 1
  fi
  COMMON_ARGS+=(--phenotype_groups "${GROUP_FILE}")
fi

mkdir -p "${NOMINAL_DIR}" "${PERMUTATION_DIR}" "${SUSIE_DIR}"

echo "[INFO] prefix=${RUN_PREFIX} mode=${RUN_QTL_MODE} group_mapping=${GROUP_MAPPING}"

# --- tensorQTL passes ---

if [[ "${RUN_QTL_MODE}" == "all" || "${RUN_QTL_MODE}" == "nominal" ]]; then
  "${PYTHON_BIN}" -m "${TENSORQTL_MODULE}" "${GENO_PREFIX}" "${PHENO_BED}" "${RUN_PREFIX}" \
    --mode cis_nominal \
    "${COMMON_ARGS[@]}" \
    -o "${NOMINAL_DIR}"
fi

if [[ "${RUN_QTL_MODE}" == "all" || "${RUN_QTL_MODE}" == "permutation" ]]; then
  "${PYTHON_BIN}" -m "${TENSORQTL_MODULE}" "${GENO_PREFIX}" "${PHENO_BED}" "${RUN_PREFIX}" \
    --mode cis \
    "${COMMON_ARGS[@]}" \
    -o "${PERMUTATION_DIR}"
fi

if [[ "${RUN_QTL_MODE}" == "all" || "${RUN_QTL_MODE}" == "susie" ]]; then
  "${PYTHON_BIN}" -m "${TENSORQTL_MODULE}" "${GENO_PREFIX}" "${PHENO_BED}" "${RUN_PREFIX}" \
    --mode cis_susie \
    --cis_output "${PERMUTATION_DIR}/${RUN_PREFIX}.cis_qtl.txt.gz" \
    "${COMMON_ARGS[@]}" \
    -o "${SUSIE_DIR}"
fi

echo "[INFO] end: $(date '+%F %T')"
