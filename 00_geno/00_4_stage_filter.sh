#!/usr/bin/env bash
# Stage-specific MAF/MAC/HWE filter on the imputed VCF; keeps variants passing both stages.
# RAW_VCF: imputed VCF from ChinaMAP. META_FILE: sample table with a Prenatal/Postnatal column.
set -euo pipefail
: "${APA_ROOT:?source config/paths.local.sh}"

# --- paths ---
RAW_VCF="${RAW_VCF:-${APA_ROOT}/raw/CHB_all_anno.vcf.gz}"
META_FILE="${META_FILE:-${APA_ROOT}/data/04_qtl/CHB.all.meta}"
OUT="${APA_ROOT}/data/00_geno"
WORK="${TMPDIR:-/tmp}/stage_qc_$$"
THREADS="${THREADS:-8}"
mkdir -p "${WORK}" "${OUT}"

# --- sample lists ---
awk '$0 ~ /Prenatal/ {print $1}' "${META_FILE}" > "${WORK}/pre.sample"
awk '$0 ~ /Postnatal/ {print $1}' "${META_FILE}" > "${WORK}/post.sample"
cat "${WORK}/pre.sample" "${WORK}/post.sample" > "${WORK}/all.sample"

# --- stage-specific variant QC ---
plink2 --vcf "${RAW_VCF}" --keep "${WORK}/pre.sample" \
    --maf 0.05 --mac 20 --hwe 1e-6 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 500 \
    --write-snplist --threads "${THREADS}" --out "${WORK}/pre_qc"

plink2 --vcf "${RAW_VCF}" --keep "${WORK}/post.sample" \
    --maf 0.05 --mac 20 --hwe 1e-6 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 500 \
    --write-snplist --threads "${THREADS}" --out "${WORK}/post_qc"

# --- shared variant set ---
awk 'FNR == NR { keep[$1] = 1; next } ($1 in keep) { print $1 }' \
    "${WORK}/pre_qc.snplist" "${WORK}/post_qc.snplist" \
    | sed 's/^/chr/' > "${WORK}/common_variants.txt"

# --- final cohort PGEN ---
plink2 --vcf "${RAW_VCF}" --keep "${WORK}/all.sample" \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 500 \
    --output-chr chrM --make-pgen --threads "${THREADS}" \
    --out "${WORK}/genotype_intermediate"

plink2 --pfile "${WORK}/genotype_intermediate" \
    --extract "${WORK}/common_variants.txt" \
    --make-pgen --output-chr chrM --threads "${THREADS}" \
    --out "${OUT}/genotype_clean"

echo "Stage variant QC complete: ${OUT}/genotype_clean"
