#!/usr/bin/env bash
# Read depth over the R2b windows for each selected BAM, CPM-normalized by mapped reads.
# Reads data/07_coloc/R2b_selected_bam_samples.tsv; writes data/07_coloc/R2b_depth_cpm.tsv.gz.
: "${APA_ROOT:?set APA_ROOT (see release/config/paths.sh)}"
set -euo pipefail

out_dir="${APA_ROOT}/data/07_coloc"
selection="$out_dir/R2b_selected_bam_samples.tsv"
output="$out_dir/R2b_depth_cpm.tsv.gz"

# --- depth per selected BAM ---
printf 'gene\tstage\tgenotype_group\tsample\tchrom\tposition\tdepth\tmapped_reads\tdepth_cpm\n' | gzip -c > "$output"

tail -n +2 "$selection" | while IFS=$'\t' read -r gene event variant stage genotype_group sample dosage phenotype group_mean age batch tpm median_tpm chrom start end bam change; do
  mapped=$(samtools idxstats "$bam" | awk '{sum += $3} END {print sum + 0}')
  samtools depth -a -r "${chrom}:${start}-${end}" "$bam" |
    awk -v gene="$gene" -v stage="$stage" -v group="$genotype_group" -v sample="$sample" -v mapped="$mapped" \
      'BEGIN {OFS="\t"} {print gene,stage,group,sample,$1,$2,$3,mapped,$3/(mapped/1000000)}' |
    gzip -c >> "$output"
done

# --- coverage guard ---
zcat "$output" | awk 'NR > 1 {n[$1 FS $4]++} END {for (k in n) if (n[k] < 4000) exit 1}'
