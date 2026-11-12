#!/usr/bin/env bash
# DaPars2 inputs: 3'UTR reference, wig and read-depth list, configure file.
set -euo pipefail

# --- config ---
: "${APA_ROOT:?source config/paths.local.sh}"
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_LIST="${SAMPLE_LIST:-${APA_ROOT}/raw/bam/sample_list.txt}"
GENE_ANNO="${GENE_ANNO:-${APA_ROOT}/ref/refseq/refseq.annotation.bed}"
REF_ID_MAP="${REF_ID_MAP:-${APA_ROOT}/ref/refseq/refseq_to_gene.txt}"
OUT="${APA_ROOT}/data/01_quant/apa"
WIG_DIR="${WIG_DIR:-${OUT}/wig}"
FLAGSTAT_DIR="${FLAGSTAT_DIR:-${OUT}/flagstat}"
COVERAGE_THRESHOLD="${COVERAGE_THRESHOLD:-15}"
THREADS="${THREADS:-8}"
CONFIG="${OUT}/dapars2_config.txt"

mkdir -p "${OUT}"

# --- 3'UTR reference ---
python "${SRC_DIR}/DaPars_Extract_Anno.py" -b "${GENE_ANNO}" -s "${REF_ID_MAP}" -o "${OUT}/refseq_3utr_annotation.bed"

# --- wig and read depth ---
python "${SRC_DIR}/read_depth.py" --sample_list "${SAMPLE_LIST}" \
    --path_wig "${WIG_DIR}" --path_flagstat "${FLAGSTAT_DIR}" --output "${OUT}/wigFile_and_readDepth.txt"

# --- configure file ---
python "${SRC_DIR}/dapars2_config.py" --annotation_3utr "${OUT}/refseq_3utr_annotation.bed" \
    --wigFile_depth "${OUT}/wigFile_and_readDepth.txt" --coverage_threshold "${COVERAGE_THRESHOLD}" \
    --threads "${THREADS}" --out_config_name "${CONFIG}"
