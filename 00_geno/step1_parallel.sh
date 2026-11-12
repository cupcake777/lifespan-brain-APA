#!/bin/bash
# Per-sample Sentieon BWA alignment, dedup, recalibration and gVCF calling (array job).
#SBATCH --job-name=sentieon_array
#SBATCH --output=logs/sentieon_%A_%a.out
#SBATCH --error=logs/sentieon_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=CU
#SBATCH --array=0-1
: "${APA_ROOT:?source config/paths.local.sh}"
: "${APA_GROUP_SHARE:?source config/paths.local.sh}"
: "${SENTIEON_HOME:?source config/paths.local.sh}"
: "${SENTIEON_LICENSE:?source config/paths.local.sh}"

set -euxo pipefail

# --- paths ---
export SENTIEON_INSTALL_DIR=${SENTIEON_HOME}/sentieon-genomics-202503
export SENTIEON_LICENSE=${SENTIEON_LICENSE}
export PATH=$SENTIEON_INSTALL_DIR/bin:$PATH

REFERENCE=${APA_ROOT}/ref/gatk_bundle/Homo_sapiens_assembly38.fasta
DBSNP=${APA_ROOT}/ref/gatk_bundle/Homo_sapiens_assembly38.dbsnp138.vcf.gz
KNOWN_1000G_INDELS=${APA_ROOT}/ref/gatk_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz
KNOWN_MILLS_INDELS=${APA_ROOT}/ref/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

NT=${SLURM_CPUS_PER_TASK:-8}
SAMPLE_PARENT_DIR=${APA_GROUP_SHARE}/new_chb/adult_fastq
RESULT_DIR=${APA_GROUP_SHARE}/new_chb/adult_gvcf_sentieon_alt_aware

# --- select sample by array index ---
SAMPLES=(${SAMPLE_PARENT_DIR}/PTB522)

SAMPLE_PATH=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
if [ -z "$SAMPLE_PATH" ]; then
    echo "Error: No sample found for ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

SAMPLE_NAME=$(basename $SAMPLE_PATH)
FASTQ_1=${SAMPLE_PATH}/PTB522_1.fq.gz
FASTQ_2=${SAMPLE_PATH}/PTB522_2.fq.gz
WORKDIR=$RESULT_DIR/$SAMPLE_NAME
mkdir -p $WORKDIR
cd $WORKDIR

# --- S1 alignment and sort ---
if [ ! -f sorted.bam ]; then
    sentieon bwa mem -M -R "@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME\tPL:DNBSEQ" -t $NT -K 10000000 $REFERENCE $FASTQ_1 $FASTQ_2 \
    | sentieon util sort --bam_compression 1 -r $REFERENCE -o sorted.bam -t $NT --sam2bam -i -
fi

# --- S2 duplicate removal ---
if [ ! -f deduped.bam ]; then
    sentieon driver -t $NT -i sorted.bam --algo LocusCollector --fun score_info score.txt
    sentieon driver -t $NT -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt --bam_compression 1 deduped.bam
fi

# --- S3 recalibration ---
if [ ! -f recal_data.table ]; then
    sentieon driver -r $REFERENCE -t $NT -i deduped.bam \
    --algo QualCal -k $DBSNP -k $KNOWN_MILLS_INDELS -k $KNOWN_1000G_INDELS recal_data.table
fi

# --- S4 variant calling (gVCF) ---
if [ ! -f ${SAMPLE_NAME}.gvcf.gz.tbi ]; then
sentieon driver -r $REFERENCE -t $NT -i deduped.bam -q recal_data.table \
    --algo Haplotyper \
    -d $DBSNP \
    --emit_conf=10 \
    --call_conf=10 \
    --emit_mode gvcf \
    ${SAMPLE_NAME}.gvcf.gz
fi

echo "Sample $SAMPLE_NAME finished."
