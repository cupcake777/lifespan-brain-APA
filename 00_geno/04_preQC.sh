#!/bin/bash
#SBATCH --job-name=QC
#SBATCH --partition=CU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/pre_imputation%j.log
#SBATCH --error=logs/pre_imputation%j.err
# Pre-imputation QC: PASS sites, het and IBD sample filters, SNP filters; writes qc.final.vcf.gz.
set -euo pipefail
: "${APA_ROOT:?source config/paths.local.sh}"
export APA_ROOT

# --- paths ---
INPUT_VCF="${APA_ROOT}/data/00_geno/vqsr/vqsr.vcf.gz"
OUT="${APA_ROOT}/data/00_geno"
CHRNAME="${APA_ROOT}/ref/resource/Chrname.txt"
WORK="${TMPDIR:-/tmp}/preqc_${SLURM_JOB_ID:-$$}"
mkdir -p "${WORK}" "${OUT}" "${APA_ROOT}/figures/00_geno"
cd "${WORK}"

# --- het filter script ---
cat > het_filter.R << 'EOF'
het <- read.table("qc2.het", header=TRUE, check.names=FALSE)
mean_f <- mean(het$F, na.rm=TRUE)
sd_f <- sd(het$F, na.rm=TRUE)
png(file.path(Sys.getenv("APA_ROOT"), "figures/00_geno/het_F_distribution.png"))
hist(het$F, xlab="F coefficient", main="F Distribution")
dev.off()
summary(het$F)
outliers <- subset(het, F < (mean_f - 3*sd_f) | F > (mean_f + 3*sd_f))
write.table(outliers[,c("FID","IID")], file.path(Sys.getenv("APA_ROOT"), "data/00_geno/het_outliers.ID"), row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("samples:", nrow(het), "\n")
cat("outliers:", nrow(outliers), "\n")
EOF

# --- PASS sites and LD pruning ---
bcftools view -f PASS -Oz -o raw.vcf.gz "${INPUT_VCF}"
bcftools index -t raw.vcf.gz
bcftools annotate --rename-chrs "${CHRNAME}" raw.vcf.gz -Oz -o raw.nochr.vcf.gz
bcftools index -t raw.nochr.vcf.gz
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' raw.nochr.vcf.gz -Oz -o raw.nochr.setid.vcf.gz
plink --vcf raw.nochr.setid.vcf.gz --make-bed --out qc0
plink --bfile qc0 --mind 0.1 --make-bed --out qc1
plink --bfile qc1 --indep-pairwise 50 5 0.2 --out qc2
plink --bfile qc1 --extract qc2.prune.in --het --out qc2

# --- heterozygosity outliers ---
Rscript het_filter.R
OUTLIERS=$(wc -l < "${OUT}/het_outliers.ID")
TOTAL=$(wc -l < qc1.fam)
echo "Heterozygosity outliers: $OUTLIERS of $TOTAL samples"
if [ "$OUTLIERS" -gt $((TOTAL / 10)) ]; then
  echo "Warning: Too many heterozygosity outliers detected! Please check the F coefficient distribution and outlier threshold."
fi
plink --bfile qc1 --remove "${OUT}/het_outliers.ID" --make-bed --out qc1_clean

# --- relatedness (IBD) filter ---
plink2 --bfile qc1_clean --king-cutoff 0.0884 --out king_filtered
plink --bfile qc1_clean --keep king_filtered.king.cutoff.in.id --make-bed --out qc3
TOTAL_SAMPLES=$(wc -l < qc1_clean.fam)
KEPT_SAMPLES=$(wc -l < qc3.fam)
REMOVED=$(($TOTAL_SAMPLES - $KEPT_SAMPLES))
echo "IBD Filtering: Total samples = ${TOTAL_SAMPLES}, Remaining = ${KEPT_SAMPLES}, Removed = ${REMOVED}"

# --- SNP filters ---
plink --bfile qc3 --geno 0.05 --make-bed --out qc4
plink --bfile qc4 --set-hh-missing --make-bed --out qc5
plink --bfile qc5 --hwe 1e-6 --make-bed --out qc5.hwe
plink --bfile qc5.hwe --test-mishap --out qc5.mis
awk 'BEGIN{OFS="\t"}{if($8 < 1e-9 && $9=="HETERO"){print $1}}' qc5.mis.missing.hap | sort | uniq > mishap.txt
plink --bfile qc5.hwe --exclude mishap.txt --make-bed --out qc6.mis
plink --bfile qc6.mis --maf 0.01 --make-bed --out qc7
plink --bfile qc7 --freq --out qc7

# --- palindromic SNPs with MAF > 0.4 ---
awk 'NR>1 && (($3=="A" && $4=="T") || ($3=="T" && $4=="A") || ($3=="G" && $4=="C") || ($3=="C" && $4=="G")) && $5>0.4 {print $2}' qc7.frq > qc7.ambig.remove.txt
plink --bfile qc7 --exclude qc7.ambig.remove.txt --make-bed --out qc8
plink --bfile qc8 --recode vcf-iid bgz --out "${OUT}/qc.final"

echo 'pre-imputation QC finished!'
