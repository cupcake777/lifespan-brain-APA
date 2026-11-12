#!/bin/bash
# Format the prenatal, postnatal and interaction APA-QTL signif pairs as sLDSC sumstats
: "${APA_ROOT:?source config/paths.local.sh}"
set -euo pipefail

# --- paths ---
SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IN="${APA_ROOT}/data/07_coloc"
OUT="${APA_ROOT}/data/06_disease"
mkdir -p "${OUT}"

# --- one sumstats file per QTL set ---
python "${SDIR}/01_format_qtl.py" "${IN}/prenatal.signif_pairs.txt.gz" "${OUT}/prenatal.sumstats.gz"
python "${SDIR}/01_format_qtl.py" "${IN}/postnatal.signif_pairs.txt.gz" "${OUT}/postnatal.sumstats.gz"
python "${SDIR}/01_format_qtl.py" "${IN}/intx.signif_pairs.txt.gz" "${OUT}/interaction.sumstats.gz"
