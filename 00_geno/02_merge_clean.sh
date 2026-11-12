#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output=merge_%j.out
#SBATCH --error=merge_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --partition=TEST
# Concatenate the 24 joint-genotyped chromosome VCFs into one indexed VCF.
set -euo pipefail
: "${APA_ROOT:?source config/paths.local.sh}"

# --- paths ---
OUTDIR=${APA_ROOT}/data/00_geno/joint_calling
MERGED=${APA_ROOT}/data/00_geno/joint_all_chr.vcf.gz
cd "${OUTDIR}"

# --- chromosome list ---
CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
: > vcf_list.txt
for chr in "${CHRS[@]}"; do
  vcf="joint_genotyped.${chr}.vcf.gz"
  [[ -s "${vcf}" ]] || { echo "Missing joint-genotyped VCF: ${OUTDIR}/${vcf}" >&2; exit 1; }
  printf '%s\n' "${vcf}" >> vcf_list.txt
done

# --- index and concatenate ---
while IFS= read -r f; do
  [[ -f "${f}.tbi" ]] || tabix -p vcf "${f}"
done

bcftools concat -Oz -o "${MERGED}" -f vcf_list.txt
tabix -p vcf "${MERGED}"

echo 'finish joint genotyping and merging!'
