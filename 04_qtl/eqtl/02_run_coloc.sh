#!/usr/bin/env bash
# apaQTL x eQTL colocalization: extraction, coloc, summary and locus plots.
#SBATCH -J apa_eqtl_coloc
#SBATCH -p cu
#SBATCH -c 16
#SBATCH --mem=96G
: "${APA_HPC_HOME:?source config/paths.sh}"
: "${APA_ROOT:?source config/paths.sh}"

set -euo pipefail

SCRIPT_DIR=$(cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" && pwd)

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  mkdir -p "${APA_ROOT}/logs"
  sbatch --chdir="${SCRIPT_DIR}" -o "${APA_ROOT}/logs/02_coloc_%j.out" -e "${APA_ROOT}/logs/02_coloc_%j.err" "$0"
  exit 0
fi

cd "${SCRIPT_DIR}"

echo "[$(date '+%F %T')] stage 1/4: extraction"
python 02-1_extract.py
echo "[$(date '+%F %T')] stage 2/4: coloc"
Rscript --vanilla 02-2_coloc.R
echo "[$(date '+%F %T')] stage 3/4: summary"
Rscript --vanilla 02-3_summary.R
echo "[$(date '+%F %T')] stage 4/4: locus"
Rscript --vanilla 02-4_locus.R
echo "[$(date '+%F %T')] coloc workflow complete"
