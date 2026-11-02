#!/usr/bin/env bash
# DaPars2 inputs: 3'UTR reference, wig and read-depth list, configure file.
set -euo pipefail

: "${APA_ROOT:?source config/paths.local.sh}"
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENE_ANNO="${GENE_ANNO:-${APA_ROOT}/ref/refseq/refseq.annotation.bed}"
REF_ID_MAP="${REF_ID_MAP:-${APA_ROOT}/ref/refseq/refseq_to_gene.txt}"
OUT="${APA_ROOT}/data/01_quant/apa"
WIG_DIR="${OUT}/wig"
FLAGSTAT_DIR="${OUT}/flagstat"
COVERAGE_THRESHOLD="${COVERAGE_THRESHOLD:-15}"
THREADS="${THREADS:-8}"
CONFIG="${OUT}/dapars2_config.txt"
DEPTH_FILE="${OUT}/wig_and_read_depth.txt"
ANNOTATION_FILE="${OUT}/refseq_3utr_annotation.bed"

mkdir -p "${OUT}"

python "${SRC_DIR}/DaPars_Extract_Anno.py" -b "${GENE_ANNO}" -s "${REF_ID_MAP}" -o "${ANNOTATION_FILE}"

python "${SRC_DIR}/read_depth.py" --sample_list "${SAMPLE_LIST:-${APA_ROOT}/raw/bam/sample_list.txt}" \
    --path_wig "${WIG_DIR}" --path_flagstat "${FLAGSTAT_DIR}" --output "${DEPTH_FILE}"

wig_files=()
while IFS=$'\t' read -r wig_file read_depth; do
    [[ -s "${wig_file}" ]] || { echo "[ERROR] missing or empty wig: ${wig_file}" >&2; exit 1; }
    [[ "${read_depth}" =~ ^[0-9]+$ && "${read_depth}" -gt 0 ]] || {
        echo "[ERROR] invalid mapped-read count for ${wig_file}: ${read_depth}" >&2
        exit 1
    }
    wig_files+=("${wig_file}")
done < "${DEPTH_FILE}"
[[ ${#wig_files[@]} -gt 0 ]] || { echo "[ERROR] no samples in ${DEPTH_FILE}" >&2; exit 1; }

aligned_wigs="$(IFS=,; echo "${wig_files[*]}")"
{
    printf 'Annotated_3UTR=%s\n' "${ANNOTATION_FILE}"
    printf 'Aligned_Wig_files=%s\n' "${aligned_wigs}"
    printf 'Output_directory=%s\n' "${OUT}/Dapars2_out"
    printf 'Output_result_file=Dapars2\n'
    printf 'Coverage_threshold=%s\n' "${COVERAGE_THRESHOLD}"
    printf 'Num_Threads=%s\n' "${THREADS}"
    printf 'sequencing_depth_file=%s\n' "${DEPTH_FILE}"
} > "${CONFIG}"

echo "DaPars2 inputs written: ${CONFIG}"
