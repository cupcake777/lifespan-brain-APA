#!/usr/bin/env bash
# Runs the module 07 colocalization steps in order; run from the func directory.
set -euo pipefail

python 00_preflight.py
python 01_qtl_ref.py
python 01b_inter_qtl.py
python 02_func.py
python 03_run_plan.py
python 02b_protein_ev.py
python 02d_gwas_coloc.py
python 02e_cross_stage.py
if [[ "${APA05_RUN_EQTL_COLOC:-0}" == "1" ]]; then
  Rscript 02c_eqtl_coloc.R
  export APA05_USE_EQTL_COLOC=1
else
  printf 'Skipping optional eQTL coloc; set APA05_RUN_EQTL_COLOC=1 to enable.\n'
  export APA05_USE_EQTL_COLOC=0
fi
python 05b_flip_protein.py
python 03_run_plan.py
python 04_ranked_cands.py
python 05_report.py
echo 'Finish!'
