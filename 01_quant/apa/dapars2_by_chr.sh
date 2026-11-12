#!/usr/bin/env bash
# Run DaPars2 once per chromosome in parallel (GNU parallel).
set -euo pipefail

WORKDIR="${WORKDIR:?set WORKDIR (run_dapars2.sh sets it)}"
CHR_LIST_FILE="${CHR_LIST_FILE:?}"
CONFIG_FILE="${CONFIG_FILE:?}"
DAPARS_SCRIPT="${DAPARS_SCRIPT:?}"
N_JOBS="${N_JOBS:-5}"
LOGDIR="${APA_ROOT}/logs/01_quant_dapars2"
JOBLOG="${LOGDIR}/parallel_joblog.txt"
mkdir -p "${LOGDIR}"

grep -v '^#' "${CHR_LIST_FILE}" | parallel -j "${N_JOBS}" --bar --joblog "${JOBLOG}" \
    "python3 ${DAPARS_SCRIPT} ${CONFIG_FILE} {} > ${LOGDIR}/dapars2_{}.log 2>&1"
echo "DaPars2 per-chromosome runs finished; joblog ${JOBLOG}"