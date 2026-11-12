#!/bin/bash
# Per-chromosome GenomicsDB import and GATK GenotypeGVCFs (one array task per chromosome).
#SBATCH --job-name=joint
#SBATCH --output=joint_%A_%a.out
#SBATCH --error=joint_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=96:00:00
#SBATCH --partition=CU
#SBATCH --array=0-23
: "${APA_ROOT:?source config/paths.local.sh}"
: "${APA_HPC_HOME:?source config/paths.local.sh}"
: "${APA_GROUP_SHARE:?source config/paths.local.sh}"

# --- paths ---
unset LD_LIBRARY_PATH
export JAVA_HOME=${APA_HPC_HOME}/bin/jdk-17.0.11+9
export PATH=$JAVA_HOME/bin:$PATH

GATK=${APA_HPC_HOME}/tools/gatk-4.5.0.0/gatk
REFERENCE=${APA_ROOT}/ref/hg38.fa
SAMPLEDIR=${APA_GROUP_SHARE}/data/CHFB/process/processing/sample
TMPDIR_BASE=${APA_GROUP_SHARE}/data/CHFB/process/processing/tmp
OUTDIR=${APA_ROOT}/data/00_geno/joint_calling
LOGDIR=${APA_ROOT}/logs
mkdir -p ${OUTDIR} ${LOGDIR}

# --- one chromosome per array task ---
CHRS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
TMPDIR=$TMPDIR_BASE/${CHR}_${SLURM_JOB_ID}
mkdir -p $TMPDIR
GVCFLIST=$OUTDIR/gvcf_${CHR}.list
ls $SAMPLEDIR/*.g.vcf.gz > $GVCFLIST

# --- GenomicsDB import ---
GENOMICSDB_DIR=$TMPDIR/genomicsdb_${CHR}
rm -rf $GENOMICSDB_DIR

$GATK --java-options "-Xmx48G -Djava.io.tmpdir=$TMPDIR" GenomicsDBImport \
  --genomicsdb-workspace-path $GENOMICSDB_DIR \
  --batch-size 50 \
  -L $CHR \
  $(cat $GVCFLIST | awk '{print "-V "$0}') \
  &> ${LOGDIR}/00_geno_GenomicsDBImport_${CHR}.log

if [ $? -ne 0 ]; then
  echo "GenomicsDBImport failed for $CHR" >&2
  exit 1
fi

# --- genotype and index ---
$GATK --java-options "-Xmx48G -Djava.io.tmpdir=$TMPDIR" GenotypeGVCFs \
  -R $REFERENCE \
  -V gendb://$GENOMICSDB_DIR \
  -O $OUTDIR/joint_genotyped.${CHR}.vcf.gz \
  &> ${LOGDIR}/00_geno_GenotypeGVCFs_${CHR}.log

if [ $? -ne 0 ]; then
  echo "GenotypeGVCFs failed for $CHR" >&2
  exit 1
fi

$GATK IndexFeatureFile \
  -I $OUTDIR/joint_genotyped.${CHR}.vcf.gz \
  &> ${LOGDIR}/00_geno_IndexFeatureFile_${CHR}.log

rm -rf $GENOMICSDB_DIR $TMPDIR
echo "$CHR finished!"
