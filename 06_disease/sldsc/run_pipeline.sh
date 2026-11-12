#!/bin/bash
# Run the S-LDSC steps: bash run_pipeline.sh [phase]   (phase 1-7, default: all)
# Phases 1, 2 and 5 (QTL sumstats, annotation BEDs, GWAS sumstats) run before 3-7 in "all".
set -euo pipefail

SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHASE=${1:-"all"}

run_phase() {
    echo "Phase $1: $2"
    cd "${SDIR}"
    bash "$3"
}

case ${PHASE} in
    1) run_phase 1 "Format APA-QTL as sLDSC sumstats" 01_format_qtl_all.sh ;;
    2) run_phase 2 "Generate custom annotation BED files" 02_annot_bed_all.sh ;;
    3) run_phase 3 "Create custom annotation files (.annot.gz)" 03_run_make_annot.sh ;;
    4) run_phase 4 "Compute stratified LD scores" 04_ld_scores.sh ;;
    5) run_phase 5 "Format GWAS for sLDSC input" 05_format_all.sh ;;
    6) run_phase 6 "Run sLDSC regression" 06_sldsc_reg.sh ;;
    7) run_phase 7 "Collect and summarize results" 07_collect.sh ;;
    all)
        echo "Phases 1, 2 and 5 must have run first (QTL sumstats, annotation BEDs, GWAS sumstats)."
        run_phase 3 "Create custom annotation files (.annot.gz)" 03_run_make_annot.sh
        run_phase 4 "Compute stratified LD scores" 04_ld_scores.sh
        run_phase 6 "Run sLDSC regression" 06_sldsc_reg.sh
        run_phase 7 "Collect and summarize results" 07_collect.sh
        ;;
    *)
        echo "Error: Unknown phase number '${PHASE}'"
        echo "Valid phases: 1, 2, 3, 4, 5, 6, 7, all"
        exit 1
        ;;
esac
