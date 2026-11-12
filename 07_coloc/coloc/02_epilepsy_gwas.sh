#!/usr/bin/env bash
# Normalize the ILAE3 epilepsy GWAS (whitespace-delimited .tbl) to a gzipped TSV in data/07_coloc.
: "${APA_HPC_HOME:?set APA_HPC_HOME (see release/config/paths.sh)}"
: "${APA_ROOT:?set APA_ROOT (see release/config/paths.sh)}"
set -euo pipefail

# --- config ---
input="${APA_HPC_HOME}/share_group_folder/.ref/GWAS/raw/ILAE3/ILAE3_Caucasian_all_epilepsy_final.tbl"
output="${APA_ROOT}/data/07_coloc/ILAE3_Caucasian_all_epilepsy_final.tsv.gz"
tmp="${output}.tmp"

# --- skip if already normalized ---
if [[ -s "${output}" ]]; then
  echo "Existing normalized ILAE3 file found: ${output}"
  exit 0
fi

# --- normalize whitespace-delimited tail columns ---
LC_ALL=C awk 'BEGIN { OFS="\t" }
  NR == 1 {
    print "CHR", "BP", "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "Effective_N", "Z-score", "P-value", "Direction", "Beta", "SE"
    next
  }
  { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13 }
' "${input}" | gzip -c > "${tmp}"

mv "${tmp}" "${output}"
echo "Normalized ILAE3 file: ${output}"
