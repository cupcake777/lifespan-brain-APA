#!/usr/bin/env bash
# DaPars2 APA quantification: coverage, inputs, per-chromosome runs, merge.
set -euo pipefail

: "${APA_ROOT:?source config/paths.local.sh}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${APA_ROOT}/data/01_quant/apa"
CHR_LIST_FILE="${CHR_LIST_FILE:-${APA_ROOT}/ref/chrList.txt}"
N_JOBS="${N_JOBS:-5}"
export SAMPLE_LIST="${SAMPLE_LIST:-${APA_ROOT}/raw/bam/sample_list.txt}"

bash "${SCRIPT_DIR}/00_bam_to_wig.sh"
bash "${SCRIPT_DIR}/prep_apa_inputs.sh"

(cd "${OUT}" && WORKDIR="${OUT}" CHR_LIST_FILE="${CHR_LIST_FILE}" \
    CONFIG_FILE="${OUT}/dapars2_config.txt" DAPARS_SCRIPT="${SCRIPT_DIR}/Dapars2_Multi_Sample.py" \
    N_JOBS="${N_JOBS}" bash "${SCRIPT_DIR}/dapars2_by_chr.sh")

WORKDIR="${OUT}" FINAL_OUTPUT="${OUT}/all.apa.count" bash "${SCRIPT_DIR}/dapars2_merge.sh"

echo "DaPars2 workflow complete: ${OUT}/all.apa.count"
