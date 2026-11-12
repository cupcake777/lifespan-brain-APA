#!/usr/bin/env bash
# Per-sample flagstat and bedGraph coverage for DaPars2, one row per BAM in the sample list.
set -euo pipefail

# --- config ---
: "${APA_ROOT:?source config/paths.local.sh}"
SAMPLE_LIST="${SAMPLE_LIST:-${APA_ROOT}/raw/bam/sample_list.txt}"
OUT="${APA_ROOT}/data/01_quant/apa"
WIG_DIR="${WIG_DIR:-${OUT}/wig}"
FLAGSTAT_DIR="${FLAGSTAT_DIR:-${OUT}/flagstat}"
THREADS="${THREADS:-8}"

mkdir -p "${WIG_DIR}" "${FLAGSTAT_DIR}"

# --- per sample ---
while IFS=$'\t' read -r sample bam_path _; do
    [[ -n "${sample}" ]] || continue
    [[ -s "${bam_path}" ]] || { echo "Missing BAM for ${sample}: ${bam_path}" >&2; exit 1; }

    samtools quickcheck "${bam_path}"
    samtools flagstat -@ "${THREADS}" "${bam_path}" > "${FLAGSTAT_DIR}/${sample}.flagstat.tmp"
    bedtools genomecov -ibam "${bam_path}" -bg -split > "${WIG_DIR}/${sample}.wig.tmp"
    mv "${FLAGSTAT_DIR}/${sample}.flagstat.tmp" "${FLAGSTAT_DIR}/${sample}.flagstat"
    mv "${WIG_DIR}/${sample}.wig.tmp" "${WIG_DIR}/${sample}.wig"
done < "${SAMPLE_LIST}"
