#!/usr/bin/env bash
# Refresh the Figure 6d SMR source and rebuild the handoff manifest in data/07_coloc.
: "${APA_HPC_HOME:?source config/paths.sh}"
: "${APA_ROOT:?source config/paths.sh}"
set -euo pipefail

cd ${APA_ROOT}/data/07_coloc
awk -F '\t' 'NR == 1 || $3 == "ARL3"' ${APA_HPC_HOME}/project/work/apa/06_smr/results/smr_all_results.tsv > Figure6d_ARL3_smr.tsv
python assemble_main_figures_handoff2.py 5
python assemble_main_figures_handoff.py 6
python ${APA_ROOT}/release/07_coloc/handoff/build_handoff.py
