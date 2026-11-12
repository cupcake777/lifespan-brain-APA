#!/bin/bash
# Collect sLDSC regression results into the enrichment summary, then plot them
: "${APA_ROOT:?source config/paths.sh}"
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${APA_ROOT}/data/06_disease/regression"
OUTPUT_FILE="${APA_ROOT}/data/06_disease/sldsc_enrichment_summary.tsv"
LDSC_ENV="ldsc3"

micromamba run -n "${LDSC_ENV}" python "${SCRIPT_DIR}/07_collect.py" "${INPUT_DIR}" "${OUTPUT_FILE}"
micromamba run -n "${LDSC_ENV}" python "${SCRIPT_DIR}/08_plot_sldsc.py"
