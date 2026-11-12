#!/usr/bin/env bash
: "${APA_HPC_HOME:?Set APA_HPC_HOME (source config/paths.sh)}"
: "${APA_ROOT:?Set APA_ROOT (source config/paths.sh)}"
# eQTL nominal statistics at the SCAF1 and NRBP2 locus windows (Fig 5 tracks), summary stats only.
# Input: cluster nominal eQTL cache, falling back to the full tensorQTL nominal file. Run: bash fig5_eqtl_tracks.sh
set -euo pipefail
CACHE=${APA_HPC_HOME}/project/work/apa/11_eqtl/cache
NOM=${APA_ROOT}/data/04_qtl/eqtl/nominal
OUT=${APA_ROOT}/data/07_coloc
mkdir -p "$OUT" "${APA_ROOT}/logs"
TSV="$OUT/fig5_v4_eqtl_tracks.tsv"; LOG="${APA_ROOT}/logs/fig5_v4_eqtl_tracks.log"
printf "stage\tgene\tgene_id\tvariant_id\taf\tma_count\tpval_nominal\tslope\tslope_se\tsource\n" > "$TSV"; : > "$LOG"
# gene  ENSG(no version)  chrom  window_start  window_end   (windows = Fig5d/e locus JSON windows)
GENES="SCAF1 ENSG00000126461 chr19 49160219 50159469
NRBP2 ENSG00000185189 chr8 143288857 144288036"
for stage in prenatal postnatal; do
  while read -r g id chr a b; do
    n=0
    if [ -s "$CACHE/nominal_eqtl_${stage}.tsv.gz" ]; then
      # cache keeps the full tensorQTL header -> header-driven column lookup
      zcat "$CACHE/nominal_eqtl_${stage}.tsv.gz" | awk -F'\t' -v OFS='\t' -v s=$stage -v g=$g -v id=$id '
        NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
        {k=$c["phenotype_id"]; sub(/\..*/,"",k); if(k==id) print s,g,$c["phenotype_id"],$c["variant_id"],$c["af"],$c["ma_count"],$c["pval_nominal"],$c["slope"],$c["slope_se"],"cache_oct15"}' > "$OUT/.tmp"
      n=$(wc -l < "$OUT/.tmp")
    fi
    if [ "$n" -eq 0 ] && [ -s "$NOM/${stage}.cis_qtl_pairs.txt.gz" ]; then
      # header-driven column lookup for the full tensorQTL nominal file
      zcat "$NOM/${stage}.cis_qtl_pairs.txt.gz" | awk -F'\t' -v OFS='\t' -v s=$stage -v g=$g -v id=$id '
        NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
        {k=$c["phenotype_id"]; sub(/\..*/,"",k); if(k==id) print s,g,$c["phenotype_id"],$c["variant_id"],$c["af"],$c["ma_count"],$c["pval_nominal"],$c["slope"],$c["slope_se"],"nominal_full"}' > "$OUT/.tmp"
      n=$(wc -l < "$OUT/.tmp")
    fi
    cat "$OUT/.tmp" >> "$TSV"
    echo -e "$stage\t$g\t$id\tn_rows=$n\t$(cut -f10 "$OUT/.tmp" | sort -u | tr '\n' ',')" | tee -a "$LOG"
  done <<< "$GENES"
done
rm -f "$OUT/.tmp"; gzip -f "$TSV"
echo "done -> $TSV.gz"
