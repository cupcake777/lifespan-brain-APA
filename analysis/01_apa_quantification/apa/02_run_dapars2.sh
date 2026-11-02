#!/usr/bin/env bash
# Run DaPars2 independently for the chromosomes listed in CHR_LIST_FILE.
set -euo pipefail

: "${APA_ROOT:?source config/paths.local.sh}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${APA_ROOT}/data/01_quant/apa"
CHR_LIST_FILE="${CHR_LIST_FILE:-${APA_ROOT}/ref/chrList.txt}"
N_JOBS="${N_JOBS:-5}"
PYTHON_BIN="${PYTHON_BIN:-python}"
CONFIG_FILE="${OUT}/dapars2_config.txt"
DAPARS_SCRIPT="${SCRIPT_DIR}/DaPars2_Multi_Sample.py"

[[ -s "${CONFIG_FILE}" ]] || { echo "[ERROR] missing ${CONFIG_FILE}; run 01_prepare_inputs.sh" >&2; exit 1; }
[[ -s "${CHR_LIST_FILE}" ]] || { echo "[ERROR] missing chromosome list: ${CHR_LIST_FILE}" >&2; exit 1; }
[[ "${N_JOBS}" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] N_JOBS must be a positive integer" >&2; exit 1; }

mapfile -t chromosomes < <(awk 'NF && $1 !~ /^#/ {print $1}' "${CHR_LIST_FILE}")
[[ ${#chromosomes[@]} -gt 0 ]] || { echo "[ERROR] no chromosomes in ${CHR_LIST_FILE}" >&2; exit 1; }

printf '%s\n' "${chromosomes[@]}" | xargs -r -n1 -P "${N_JOBS}" \
    bash -c "\"\${1}\" \"\${2}\" \"\${3}\" \"\${4}\"" _ "${PYTHON_BIN}" "${DAPARS_SCRIPT}" "${CONFIG_FILE}"

for chromosome in "${chromosomes[@]}"; do
    result="${OUT}/Dapars2_out_${chromosome}/Dapars2_result_temp.${chromosome}.txt"
    [[ -s "${result}" ]] || { echo "[ERROR] missing DaPars2 result: ${result}" >&2; exit 1; }
done

echo "DaPars2 chromosome runs complete: ${#chromosomes[@]} chromosomes"
