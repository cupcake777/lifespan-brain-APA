# Module 00: genotype QC and imputation

Paths come from `config/paths.local.sh` (`APA_ROOT`, `APA_HPC_HOME`, `APA_GROUP_SHARE`).
Outputs go to `data/00_geno/`, figures to `figures/00_geno/`, logs to `logs/`.

Chain (Slurm, submitted by `run_module00.sh`):

1. `01_joint_geno.sh` GATK GenomicsDB import and GenotypeGVCFs, one array task per chromosome.
2. `02_merge_clean.sh` concatenates the chromosome VCFs into `data/00_geno/joint_all_chr.vcf.gz`.
3. `03_VQSR.sh` SNP and INDEL VQSR, output `data/00_geno/vqsr/vqsr.vcf.gz`.
4. `04_preQC.sh` PASS sites, heterozygosity and IBD sample filters, SNP filters, output `data/00_geno/qc.final.vcf.gz`.
5. `05_prep_impute.sh` per-chromosome VCFs for ChinaMAP, output `data/00_geno/imputation_input/`.

ChinaMAP imputation is a web step. After the imputed VCF is downloaded, run
`00_4_stage_filter.sh` with `RAW_VCF` set to it. It applies MAF 0.05, MAC 20 and HWE 1e-6
within each life stage, keeps variants passing both stages, and writes
`data/00_geno/genotype_clean.{pgen,pvar,psam}`.

`step1_parallel.sh` is a Sentieon test workflow (MHC/alt-aware). It is not the main-cohort
producer; the main cohort uses GATK HaplotypeCaller.

The DaPars2 steps (formerly `01_sep_chr.sh` and `02_combine.sh`) are in `../01_quant/apa/`
as `dapars2_by_chr.sh` and `dapars2_merge.sh`.

`gatk_vqsr_env.yml` is the conda environment for the VQSR step.
