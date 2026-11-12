#!/bin/bash
# GATK VQSR (SNP then INDEL) on the merged joint-genotyped VCF.
#SBATCH --job-name=vqsr
#SBATCH --output=vqsr_%j.out
#SBATCH --error=vqsr_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --partition=CU
: "${APA_ROOT:?source config/paths.local.sh}"
: "${APA_HPC_HOME:?source config/paths.local.sh}"
: "${APA_GROUP_SHARE:?source config/paths.local.sh}"

set -euo pipefail

# --- paths ---
unset LD_LIBRARY_PATH
export JAVA_HOME=${APA_HPC_HOME}/bin/jdk-17.0.11+9
export PATH=$JAVA_HOME/bin:$PATH

GATK_PATH=${APA_HPC_HOME}/tools/gatk-4.5.0.0/gatk
REFERENCE=${APA_ROOT}/ref/hg38.fa
BUNDLE=${APA_ROOT}/ref/gatk_bundle
INPUT_VCF=${APA_ROOT}/data/00_geno/joint_all_chr.vcf.gz
OUTPUT_DIR=${APA_ROOT}/data/00_geno/vqsr
TMP_DIR=${APA_GROUP_SHARE}/data/CHFB/process/processing/tmp
HAPMAP=${BUNDLE}/hapmap_3.3.hg38.vcf.gz
OMNI=${BUNDLE}/1000G_omni2.5.hg38.vcf.gz
PHASE1_1000G=${BUNDLE}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
DBSNP=${BUNDLE}/Homo_sapiens_assembly38.dbsnp138.vcf.gz
MILLS=${BUNDLE}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

mkdir -p "${TMP_DIR}" "${OUTPUT_DIR}"
tabix -p vcf "${INPUT_VCF}"

# --- SNP recalibration ---
${GATK_PATH} --java-options "-Xmx60G -Djava.io.tmpdir=${TMP_DIR}" VariantRecalibrator \
    -R ${REFERENCE} \
    -V ${INPUT_VCF} \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 ${OMNI} \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${PHASE1_1000G} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode SNP \
    -O ${OUTPUT_DIR}/recal_SNP.recal \
    --tranches-file ${OUTPUT_DIR}/recal_SNP.tranches \
    --rscript-file ${OUTPUT_DIR}/recal_SNP_plots.R \
    --tranche 90.0 --tranche 93.0 --tranche 95.0 --tranche 97.0 --tranche 99.0 --tranche 99.9 --tranche 100.0

if [ -f "${OUTPUT_DIR}/recal_SNP_plots.R" ]; then
    sed -i 's/space = "rgb"/space = "Lab"/g' "${OUTPUT_DIR}/recal_SNP_plots.R"
    Rscript "${OUTPUT_DIR}/recal_SNP_plots.R" || true
fi

# --- SNP apply ---
${GATK_PATH} --java-options "-Xmx60G -Djava.io.tmpdir=${TMP_DIR}" ApplyVQSR \
    -R ${REFERENCE} \
    -V ${INPUT_VCF} \
    -mode SNP \
    --recal-file ${OUTPUT_DIR}/recal_SNP.recal \
    --tranches-file ${OUTPUT_DIR}/recal_SNP.tranches \
    -O ${OUTPUT_DIR}/recal_snp.vcf.gz

# --- INDEL recalibration ---
${GATK_PATH} --java-options "-Xmx60G -Djava.io.tmpdir=${TMP_DIR}" VariantRecalibrator \
    -R ${REFERENCE} \
    -V ${OUTPUT_DIR}/recal_snp.vcf.gz \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 ${MILLS} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode INDEL \
    --max-gaussians 4 \
    -O ${OUTPUT_DIR}/recal_INDEL.recal \
    --tranches-file ${OUTPUT_DIR}/recal_INDEL.tranches \
    --rscript-file ${OUTPUT_DIR}/recal_INDEL_plots.R \
    --tranche 90.0 --tranche 93.0 --tranche 95.0 --tranche 97.0 --tranche 99.0 --tranche 99.9 --tranche 100.0

if [ -f "${OUTPUT_DIR}/recal_INDEL_plots.R" ]; then
    sed -i 's/space = "rgb"/space = "Lab"/g' "${OUTPUT_DIR}/recal_INDEL_plots.R"
    Rscript "${OUTPUT_DIR}/recal_INDEL_plots.R" || true
fi

# --- INDEL apply ---
${GATK_PATH} --java-options "-Xmx60G -Djava.io.tmpdir=${TMP_DIR}" ApplyVQSR \
    -R ${REFERENCE} \
    -V ${OUTPUT_DIR}/recal_snp.vcf.gz \
    -mode INDEL \
    --recal-file ${OUTPUT_DIR}/recal_INDEL.recal \
    --tranches-file ${OUTPUT_DIR}/recal_INDEL.tranches \
    -O ${OUTPUT_DIR}/vqsr.vcf.gz

rm ${OUTPUT_DIR}/recal_snp.vcf.gz

echo "VQSR completed successfully."
