#!/usr/bin/env bash
set -euo pipefail

# --- config ---

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${APA_ROOT:?source config/paths.local.sh}"

cd "${project_dir}"
python 00_preflight.py
python 01_prepare.py
python 02_mechanism.py
Rscript 03_deconvolution.R
python 04_consequence.py
Rscript 05_coloc.R
python 06_disease.py
Rscript 07_enrichment.R
python 08_plot_cases.py
Rscript 09_plot.R
python 10_finalize.py
