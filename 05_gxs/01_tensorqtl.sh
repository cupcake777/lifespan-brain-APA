#!/usr/bin/env bash
# Nominal cis interaction mapping (tensorQTL) for G x Stage.
set -euo pipefail
: "${APA_ROOT:?source config/paths.sh}"
: "${APA_HPC_HOME:?source config/paths.sh}"
D="$APA_ROOT/data/05_gxs"
mkdir -p "$D"

# --- nominal interaction mapping ---
python -m tensorqtl "$APA_HPC_HOME/project/work/apa/00_raw_data/genotype" "$D/phenotype.GeneRange.bed.gz" interaction \
  --covariates "$D/covariates_for_qtl.txt" \
  --interaction "$D/LifeStage.txt" \
  --mode cis_nominal \
  --output_dir "$D"
