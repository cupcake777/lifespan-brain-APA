#!/usr/bin/env bash
# tensorQTL cis mapping per life stage: nominal, permutation and SuSiE.
# Reads the eQTL pre_input from 01_preprocess.R; writes tensorQTL outputs to data/04_qtl/eqtl/{nominal,permutation,susie}.
set -euo pipefail
: "${APA_ROOT:?source config/paths.sh}"
: "${APA_HPC_HOME:?source config/paths.sh}"

D="${APA_ROOT}/data/04_qtl"
EQ="${D}/eqtl"
PRE="${EQ}/pre_input"
NOM="${EQ}/nominal"
PERM="${EQ}/permutation"
SUSIE="${EQ}/susie"
# absolute interpreter: a batch environment does not inherit the login PATH
PYTHON_BIN="${PYTHON_BIN:-${APA_HPC_HOME}/.local/share/mamba/bin/python}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.05}"
SEED="${SEED:-2026}"
RUN_QTL_MODE="${RUN_QTL_MODE:-all}"   # all | nominal | permutation | susie

# --- helpers ---
die() { echo "[ERROR] $*" >&2; exit 1; }

# --- check mode ---
case "${RUN_QTL_MODE}" in
  all|nominal|permutation|susie) ;;
  *) die "RUN_QTL_MODE must be all, nominal, permutation or susie; got: ${RUN_QTL_MODE}" ;;
esac

# --- run per stage ---
mkdir -p "${NOM}" "${PERM}" "${SUSIE}"
for phenotype_file in "${PRE}"/phenotype.*.bed.gz; do
  stage="$(basename "${phenotype_file}" .bed.gz)"
  stage="${stage#phenotype.}"
  stage_lower="$(printf '%s' "${stage}" | tr '[:upper:]' '[:lower:]')"
  geno_prefix="${PRE}/genotype.${stage_lower}"
  cov_file="${PRE}/covariates_for_qtl.${stage}.txt"

  # read the chromosome style from data rows, not ##contig headers (plink2 re-codes chr1 to 1)
  # a mismatch is silent: tensorQTL would report no cis variants
  if grep -m1 -v '^#' "${geno_prefix}.pvar" | cut -f1 | grep -q '^chr'; then geno_style=chr; else geno_style=plain; fi
  if zcat "${phenotype_file}" | sed -n '2p' | cut -f1 | grep -q '^chr'; then pheno_style=chr; else pheno_style=plain; fi
  [[ "${geno_style}" == "${pheno_style}" ]] || die "[${stage}] chromosome naming mismatch: genotype=${geno_style}-style, phenotype=${pheno_style}-style"

  common_args=(--covariates "${cov_file}" --maf_threshold "${MAF_THRESHOLD}" --seed "${SEED}")

  if [[ "${RUN_QTL_MODE}" == "all" || "${RUN_QTL_MODE}" == "nominal" ]]; then
    "${PYTHON_BIN}" -m tensorqtl "${geno_prefix}" "${phenotype_file}" "${stage_lower}" \
      --mode cis_nominal "${common_args[@]}" -o "${NOM}"
  fi

  if [[ "${RUN_QTL_MODE}" == "all" || "${RUN_QTL_MODE}" == "permutation" ]]; then
    "${PYTHON_BIN}" -m tensorqtl "${geno_prefix}" "${phenotype_file}" "${stage_lower}" \
      --mode cis "${common_args[@]}" -o "${PERM}"
  fi

  if [[ "${RUN_QTL_MODE}" == "all" || "${RUN_QTL_MODE}" == "susie" ]]; then
    "${PYTHON_BIN}" -m tensorqtl "${geno_prefix}" "${phenotype_file}" "${stage_lower}" \
      --mode cis_susie --cis_output "${PERM}/${stage_lower}.cis_qtl.txt.gz" "${common_args[@]}" -o "${SUSIE}"
  fi
done
