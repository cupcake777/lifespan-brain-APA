#!/usr/bin/env bash
# Per-sample coverage (bedtools) and flagstat for DaPars2.
set -euo pipefail

: "${APA_ROOT:?source config/paths.local.sh}"
SAMPLE_LIST="${SAMPLE_LIST:-${APA_ROOT}/raw/bam/sample_list.txt}"
OUT="${APA_ROOT}/data/01_quant/apa"
WIG_DIR="${OUT}/wig"
FLAGSTAT_DIR="${OUT}/flagstat"
THREADS="${THREADS:-8}"

mkdir -p "${WIG_DIR}" "${FLAGSTAT_DIR}"
while IFS=$'\t' read -r sample bam_path _; do
    [[ -n "${sample}" ]] || continue
    samtools quickcheck "${bam_path}"
    samtools flagstat -@ "${THREADS}" "${bam_path}" > "${FLAGSTAT_DIR}/${sample}.flagstat"
    bedtools genomecov -ibam "${bam_path}" -bg -split > "${WIG_DIR}/${sample}.wig"
done < "${SAMPLE_LIST}"
