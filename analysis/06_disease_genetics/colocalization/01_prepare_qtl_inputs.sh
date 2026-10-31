#!/usr/bin/env bash
# Build bgzipped, tabix-indexed QTL pair tables for SMR and EasyColoc in data/07_coloc.
# Needs tabix and bgzip on PATH.
: "${APA_ROOT:?set APA_ROOT (see config/paths.local.sh)}"
set -euo pipefail

# --- config ---
OUT_DIR="${APA_ROOT}/data/07_coloc"
QTL_DIR="${QTL_DIR:-${APA_ROOT}/data/04_qtl}"
INTERACTION_DIR="${INTERACTION_DIR:-${APA_ROOT}/data/05_gxs}"
SORT_PARALLEL="${SORT_PARALLEL:-4}"
SORT_MEM="${SORT_MEM:-8G}"
cd "${OUT_DIR}"

# --- helpers ---
# stdin: table with header; sort by chr/pos, bgzip, tabix. Header sorts first (chr < chr1).
write_tabix() {
    LC_ALL=C sort --parallel="${SORT_PARALLEL}" -S "${SORT_MEM}" -t$'\t' -k1,1 -k2,2n | bgzip -c > "$1"
    tabix -s 1 -b 2 -e 2 -S 1 -f "$1"
    echo "  Done: $1"
}

# --- nominal and significant pairs (prenatal, postnatal) ---
for stage in prenatal postnatal; do
    echo "[$(date +%H:%M:%S)] Processing ${stage}.all_pairs.txt.gz"
    src="${QTL_DIR}/nominal/${stage}.cis_qtl_pairs.txt.gz"
    {
        echo -e "chr\tpos\tref\talt\t$(zcat "$src" | head -1)"
        zcat "$src" | tail -n +2 | awk -F'\t' 'BEGIN{OFS="\t"}{split($2,p,":"); print p[1],p[2],p[3],p[4],$0}'
    } | write_tabix "${stage}.all_pairs.txt.gz"

    echo "[$(date +%H:%M:%S)] Processing ${stage}.signif_pairs.txt.gz"
    src="${QTL_DIR}/${stage}.sig_QTL.txt.gz"
    {
        echo -e "chr\tpos\tref\talt\tphenotype_id\tvariant_id\tstart_distance\tend_distance\taf\tma_samples\tma_count\tpval_nominal\tslope\tslope_se"
        zcat "$src" | tail -n +2 | awk -F'\t' 'BEGIN{OFS="\t"}{split($2,p,":"); print p[1],p[2],p[3],p[4],$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'
    } | write_tabix "${stage}.signif_pairs.txt.gz"
done

# --- interaction (all pairs; significant at p < 0.05) ---
echo "[$(date +%H:%M:%S)] Processing interaction.all_pairs.txt.gz"
{
    echo -e "chr\tpos\tref\talt\tphenotype_id\tvariant_id\tstart_distance\tend_distance\taf\tma_samples\tma_count\tpval_nominal\tslope\tslope_se"
    zcat "${INTERACTION_DIR}/interaction.cis_qtl_pairs.txt.gz" | tail -n +2 | awk -F'\t' 'BEGIN{OFS="\t"}{split($2,p,":"); print p[1],p[2],p[3],p[4],$1,$2,$3,$4,$5,$6,$7,$14,$15,$16}'
} | write_tabix interaction.all_pairs.txt.gz

echo "[$(date +%H:%M:%S)] Processing interaction.signif_pairs.txt.gz"
{
    echo -e "chr\tpos\tref\talt\tphenotype_id\tvariant_id\tstart_distance\tend_distance\taf\tma_samples\tma_count\tpval_nominal\tslope\tslope_se"
    zcat "${INTERACTION_DIR}/interaction.cis_qtl_top_assoc.txt.gz" | tail -n +2 | awk -F'\t' '$19<0.05{OFS="\t"; split($2,p,":"); print p[1],p[2],p[3],p[4],$1,$2,$3,$4,$5,$6,$7,$14,$15,$16}'
} | write_tabix interaction.signif_pairs.txt.gz

echo "=== QTL input tables updated in ${OUT_DIR} ==="
ls -lh ./*.gz
