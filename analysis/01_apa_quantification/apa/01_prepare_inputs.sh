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

mkdir -p "${OUT}"

python "${SRC_DIR}/DaPars_Extract_Anno.py" -b "${GENE_ANNO}" -s "${REF_ID_MAP}" -o "${OUT}/refseq_3utr_annotation.bed"

python "${SRC_DIR}/read_depth.py" --sample_list "${SAMPLE_LIST:-${APA_ROOT}/raw/bam/sample_list.txt}" \
    --path_wig "${WIG_DIR}" --path_flagstat "${FLAGSTAT_DIR}" --output "${OUT}/wig_and_read_depth.txt"

python "${SRC_DIR}/dapars2_config.py" --annotation_3utr "${OUT}/refseq_3utr_annotation.bed" \
    --wigFile_depth "${OUT}/wig_and_read_depth.txt" --coverage_threshold "${COVERAGE_THRESHOLD}" \
    --threads "${THREADS}" --out_config_name "${CONFIG}"
