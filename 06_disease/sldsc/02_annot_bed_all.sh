#!/bin/bash
# Build the prenatal, postnatal, interaction and lifespan APA-QTL window BED files
: "${APA_ROOT:?source config/paths.local.sh}"
set -euo pipefail

# --- paths ---
SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IN="${APA_ROOT}/data/07_coloc"
OUT="${APA_ROOT}/data/06_disease"
WINDOW=500
mkdir -p "${OUT}"

# --- one window BED per QTL set ---
python "${SDIR}/02_annot_bed.py" "${IN}/prenatal.signif_pairs.txt.gz" "${OUT}/prenatal_apa_qtl.bed" "${WINDOW}"
python "${SDIR}/02_annot_bed.py" "${IN}/postnatal.signif_pairs.txt.gz" "${OUT}/postnatal_apa_qtl.bed" "${WINDOW}"
python "${SDIR}/02_annot_bed.py" "${IN}/intx.signif_pairs.txt.gz" "${OUT}/interaction_apa_qtl.bed" "${WINDOW}"

# --- lifespan = union of prenatal and postnatal windows ---
sort -u "${OUT}/prenatal_apa_qtl.bed" "${OUT}/postnatal_apa_qtl.bed" > "${OUT}/lifespan_apa_qtl.bed"
