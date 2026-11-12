#!/bin/bash
#SBATCH --job-name=chinamsp
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=TEST
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
# Prepare per-chromosome VCFs for ChinaMAP imputation (upload in small batches).
set -euo pipefail
: "${APA_ROOT:?source config/paths.local.sh}"

# --- paths ---
DATA=${APA_ROOT}/data/00_geno
RES=${APA_ROOT}/ref/chinamap
WORK=${DATA}/imputation_input
PYTHON_BIN=${PYTHON_BIN:-python}
mkdir -p "${WORK}"

echo "prepare for China-map imputation!"

# --- drop INFO, sort, rename chromosomes ---
bcftools annotate -x INFO -Oz -o ${WORK}/concise.vcf.gz ${DATA}/qc.final.vcf.gz
bcftools sort --output-type z --output ${WORK}/sort.vcf.gz ${WORK}/concise.vcf.gz
bcftools index -t ${WORK}/sort.vcf.gz
bcftools annotate --rename-chrs ${RES}/chrname.txt --output-type z --output ${WORK}/sort.rename.vcf.gz ${WORK}/sort.vcf.gz
bcftools index -t ${WORK}/sort.rename.vcf.gz

# --- ChinaMAP reference filter ---
"${PYTHON_BIN}" "${RES}/ChinaMAP_filterVCF.py" --input "${WORK}/sort.rename.vcf.gz" \
    --reference ${RES}/mbiobank_ChinaMAP.phase1.vcf.gz \
    --exclude_monomorphic \
    --outputDir ${WORK}
tabix -p vcf ${WORK}/sort.rename_filtered.vcf.gz

# --- split by chromosome ---
for chr in {1..22}
do
  bcftools view -r $chr --output-type z --output ${WORK}/chr${chr}.vcf.gz ${WORK}/sort.rename_filtered.vcf.gz
done

echo 'Please upload the data in small batches, chromosome by chromosome, during imputation.'
