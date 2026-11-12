#!/usr/bin/env bash
# RSEM/STAR quantification per sample, then merged gene and isoform count/TPM matrices.
set -euo pipefail

: "${APA_ROOT:?source config/paths.local.sh}"

# --- config ---
THREADS="${THREADS:-8}"
STAR_DIR="${STAR_DIR:-/path/to/STAR-2.7.9a/bin/Linux_x86_64}"
RSEM_REF="${RSEM_REF:-/path/to/reference/gencodev40}"

DATA_DIR="${DATA_DIR:-${APA_ROOT}/raw/fastq}"
SAMPLE_LIST="${SAMPLE_LIST:-${DATA_DIR}/sample.v40.list}"

OUT_DIR="${APA_ROOT}/data/01_quant"
QC_DIR="${OUT_DIR}/qc_fastqc"
RSEM_DIR="${OUT_DIR}/rsem_v40"
STATUS_DIR="${APA_ROOT}/logs/01_quant_status"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

mkdir -p "${QC_DIR}" "${RSEM_DIR}" "${STATUS_DIR}"
log "[INFO] Pipeline started for samples in ${SAMPLE_LIST}."

PROCESSED_SAMPLES=()

while IFS= read -r raw_sample || [[ -n "${raw_sample}" ]]; do
    sample=$(echo "${raw_sample}" | tr -d '\r[:space:]')
    [[ -n "${sample}" ]] || continue

    fq1="${DATA_DIR}/${sample}_1.fq.gz"
    fq2="${DATA_DIR}/${sample}_2.fq.gz"

    sample_done="${STATUS_DIR}/${sample}_v40.done"
    if [[ -f "${sample_done}" ]]; then
        log "[SKIP] ${sample} already processed."
        PROCESSED_SAMPLES+=("${sample}")
        continue
    fi

    log "----- Processing Sample: ${sample} -----"

    # 1. FastQC
    fastqc -t "${THREADS}" -o "${QC_DIR}" "${fq1}" "${fq2}" || { log "[ERROR] FastQC failed: ${sample}"; exit 1; }

    # 2. RSEM + STAR
    if rsem-calculate-expression --paired-end --star --star-path "${STAR_DIR}" \
        --star-gzipped-read-file -p "${THREADS}" \
        "${fq1}" "${fq2}" "${RSEM_REF}" \
        "${RSEM_DIR}/${sample}"; then
        touch "${sample_done}"
        PROCESSED_SAMPLES+=("${sample}")
        log "[SUCCESS] Finished ${sample}."
    else
        log "[ERROR] RSEM failed: ${sample}"
        exit 1
    fi
done < "${SAMPLE_LIST}"

(( ${#PROCESSED_SAMPLES[@]} > 0 )) || { log "[ERROR] No samples completed successfully."; exit 1; }

# 3. MultiQC
if command -v multiqc >/dev/null 2>&1; then
    log "[INFO] Aggregating QC reports..."
    multiqc -f "${QC_DIR}" "${RSEM_DIR}" -o "${OUT_DIR}/multiqc_report"
fi

# 4. Merge gene and transcript matrices
GENE_COUNT_FILES=()
GENE_TPM_FILES=()
ISOFORM_COUNT_FILES=()
ISOFORM_TPM_FILES=()
SOURCE_SAMPLES=()
MERGED_SAMPLES=()
declare -A SEEN_MATRIX_SAMPLES=()
for s in "${PROCESSED_SAMPLES[@]}"; do
    gr="${RSEM_DIR}/${s}.genes.results"
    ir="${RSEM_DIR}/${s}.isoforms.results"
    if [[ ! -s "${gr}" || ! -s "${ir}" ]]; then
        log "[WARN] Incomplete RSEM results for ${s}; excluding it from all four matrices."
        continue
    fi
    if [[ ${#MERGED_SAMPLES[@]} -eq 0 ]]; then
        awk '{print $1}' "${gr}" > "${RSEM_DIR}/gene_id.tmp"
        awk '{print $1}' "${ir}" > "${RSEM_DIR}/transcript_id.tmp"
    elif ! awk '{print $1}' "${gr}" | cmp -s - "${RSEM_DIR}/gene_id.tmp" ||
         ! awk '{print $1}' "${ir}" | cmp -s - "${RSEM_DIR}/transcript_id.tmp"; then
        log "[ERROR] RSEM feature order differs for ${s}; refusing positional merge."
        exit 1
    fi
    # These aliases were verified by exact gene-count versus summed-transcript-count matching.
    case "${s}" in
        XYA20180516B) matrix_sample=XYA20180516 ;;
        XYA20210813B) matrix_sample=XYA20210813 ;;
        XYA20230626B) matrix_sample=XYA20230626 ;;
        XYA20231026B) matrix_sample=XYA20231026 ;;
        XYA20240412A) matrix_sample=XYA20240412 ;;
        *) matrix_sample="${s}" ;;
    esac
    if [[ -n "${SEEN_MATRIX_SAMPLES[${matrix_sample}]:-}" ]]; then
        log "[ERROR] Duplicate matrix sample ID after alias normalization: ${matrix_sample}."
        exit 1
    fi
    SEEN_MATRIX_SAMPLES[${matrix_sample}]=1
    awk -v id="${matrix_sample}" 'NR==1{print id} NR>1{print $5}' "${gr}" > "${RSEM_DIR}/${s}.cnt.tmp"
    awk -v id="${matrix_sample}" 'NR==1{print id} NR>1{print $6}' "${gr}" > "${RSEM_DIR}/${s}.tpm.tmp"
    awk -v id="${matrix_sample}" 'NR==1{print id} NR>1{print $5}' "${ir}" > "${RSEM_DIR}/${s}.iso.cnt.tmp"
    awk -v id="${matrix_sample}" 'NR==1{print id} NR>1{print $6}' "${ir}" > "${RSEM_DIR}/${s}.iso.tpm.tmp"
    GENE_COUNT_FILES+=("${RSEM_DIR}/${s}.cnt.tmp")
    GENE_TPM_FILES+=("${RSEM_DIR}/${s}.tpm.tmp")
    ISOFORM_COUNT_FILES+=("${RSEM_DIR}/${s}.iso.cnt.tmp")
    ISOFORM_TPM_FILES+=("${RSEM_DIR}/${s}.iso.tpm.tmp")
    SOURCE_SAMPLES+=("${s}")
    MERGED_SAMPLES+=("${matrix_sample}")
done

if [[ ${#MERGED_SAMPLES[@]} -eq 0 ]]; then
    log "[ERROR] No complete gene/isoform RSEM result pairs were found."
    exit 1
fi
matrix_prefix="${OUT_DIR}/all.${#MERGED_SAMPLES[@]}"
paste "${RSEM_DIR}/gene_id.tmp" "${GENE_COUNT_FILES[@]}" > "${matrix_prefix}.expr.count"
paste "${RSEM_DIR}/gene_id.tmp" "${GENE_TPM_FILES[@]}" > "${matrix_prefix}.expr.tpm"
paste "${RSEM_DIR}/transcript_id.tmp" "${ISOFORM_COUNT_FILES[@]}" > "${matrix_prefix}.isoform.count"
paste "${RSEM_DIR}/transcript_id.tmp" "${ISOFORM_TPM_FILES[@]}" > "${matrix_prefix}.isoform.tpm"
printf 'rsem_sample\tmatrix_sample\n' > "${matrix_prefix}.sample_manifest.tsv"
for i in "${!MERGED_SAMPLES[@]}"; do
    printf '%s\t%s\n' "${SOURCE_SAMPLES[${i}]}" "${MERGED_SAMPLES[${i}]}" >> "${matrix_prefix}.sample_manifest.tsv"
done
rm -f "${RSEM_DIR}/gene_id.tmp" "${RSEM_DIR}/transcript_id.tmp" \
    "${GENE_COUNT_FILES[@]}" "${GENE_TPM_FILES[@]}" \
    "${ISOFORM_COUNT_FILES[@]}" "${ISOFORM_TPM_FILES[@]}"
log "[MERGE] Four matrices saved with prefix ${matrix_prefix}."

log "[INFO] Pipeline finished successfully."
