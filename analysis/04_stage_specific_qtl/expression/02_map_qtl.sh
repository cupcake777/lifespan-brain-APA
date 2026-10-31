#!/usr/bin/env bash
# tensorQTL nominal, permutation and SuSiE cis mapping per life stage.
# Inputs: data/04_qtl/pre_input files from 01_preprocess.R. Outputs: data/04_qtl/<stage>.*.
set -euo pipefail
: "${APA_ROOT:?source config/paths.sh}"
: "${APA_HPC_HOME:?source config/paths.sh}"

D="$APA_ROOT/data/04_qtl"
PYTHON_BIN="${PYTHON_BIN:-$APA_HPC_HOME/.local/share/mamba/bin/python}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.05}"
SEED="${SEED:-2026}"
RUN_QTL_MODE="${RUN_QTL_MODE:-all}"   # all | nominal | permutation | susie

# --- helpers ---
die() { echo "[ERROR] $*" >&2; exit 1; }

# chromosome naming must match between genotype and phenotype, otherwise tensorQTL
# silently reports no cis variants
chrom_style() {
  if [[ "$1" == *.gz ]]; then
    zcat "$1" | grep -m1 -v '^#' | cut -f1 | grep -q '^chr' && echo chr || echo plain
  else
    grep -m1 -v '^#' "$1" | cut -f1 | grep -q '^chr' && echo chr || echo plain
  fi
}

run_stage() {
  local Stage="$1" stage_lower="$2"
  local geno="$D/genotype.$stage_lower"
  local pheno="$D/phenotype.$Stage.bed.gz"
  local cov="$D/covariates_for_qtl.$Stage.txt"
  local gs ps
  gs="$(chrom_style "$geno.pvar")"
  ps="$(chrom_style "$pheno")"
  [[ "$gs" == "$ps" ]] || die "[$Stage] chromosome naming mismatch: genotype=$gs, phenotype=$ps"

  local common=(--covariates "$cov" --maf_threshold "$MAF_THRESHOLD" --seed "$SEED")
  echo "[INFO] [$Stage] start $(date '+%F %T') | mode=$RUN_QTL_MODE maf=$MAF_THRESHOLD seed=$SEED"
  if [[ "$RUN_QTL_MODE" == all || "$RUN_QTL_MODE" == nominal ]]; then
    "$PYTHON_BIN" -m tensorqtl "$geno" "$pheno" "$stage_lower" --mode cis_nominal "${common[@]}" -o "$D"
  fi
  if [[ "$RUN_QTL_MODE" == all || "$RUN_QTL_MODE" == permutation ]]; then
    "$PYTHON_BIN" -m tensorqtl "$geno" "$pheno" "$stage_lower" --mode cis "${common[@]}" -o "$D"
  fi
  if [[ "$RUN_QTL_MODE" == all || "$RUN_QTL_MODE" == susie ]]; then
    "$PYTHON_BIN" -m tensorqtl "$geno" "$pheno" "$stage_lower" --mode cis_susie \
      --cis_output "$D/$stage_lower.cis_qtl.txt.gz" "${common[@]}" -o "$D"
  fi
  echo "[INFO] [$Stage] end $(date '+%F %T')"
}

# --- run ---
case "$RUN_QTL_MODE" in
  all|nominal|permutation|susie) ;;
  *) die "RUN_QTL_MODE must be all, nominal, permutation or susie; got: $RUN_QTL_MODE" ;;
esac
run_stage Prenatal prenatal
run_stage Postnatal postnatal
