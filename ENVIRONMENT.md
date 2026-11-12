# Software environment (as reported in the manuscript Methods, v22)

This is the environment **recorded** for the accepted analyses, collected from the manuscript Methods.
It is not a tested lockfile or container (see `docs/frozen_1015/RELEASE.md`, gate "Locked software environment").

| Layer | Software (version) |
|---|---|
| Sequencing QC / alignment | SOAPnuke v1.5.6; STAR v2.7.9a; RSEM v1.3.0; SAMtools v1.2; Picard v1.119 |
| Genotypes | GATK v4.5.0.0; PLINK v1.9 (clumping) and PLINK 2; ChinaMAP imputation server (external) |
| APA | DaPars2 v2.1; GENCODE v40 3′UTR annotation; PolyASite 2.0 (reference) |
| Expression / batch | edgeR v3.42.4; limma v3.56.2; ruv v0.9.7.1; variancePartition v1.30.2; WGCNA v1.71; irlba v2.3.5.1; softImpute v1.4-1 |
| Trajectories / clustering | mgcv v1.9-3; cluster v2.1.6; bluster v1.10.0 |
| Cell composition | BisqueRNA v1.0.5; limSolve v2.0.3 |
| Proteomics | MaxQuant v2.6.4.0 (DIA); sva (surrogate variables) |
| QTL | tensorQTL v1.0.10; mashr v0.2.79; flashier v1.0.7; ebnm v1.1-42 |
| Disease genetics | ldsc v1.0.1 (baseline-LD v2.2; EAS baseline-LD v2.2); MAGMA v1.10; MetaXcan v0.7.5 (S-PrediXcan); glmnet v5.0; SMR v1.3.2; coloc v5.2.3; susieR v0.14.2; EasyColoc (vendored, commit 850178b5) |
| Annotation | org.Hs.eg.db v3.17.0; MSigDB v2024.1; gnomAD v4.1.1 (LOEUF); GTEx v11 |
| Languages | R v4.3.3 (up to disease analyses) and v4.5.3 (expression, protein, mash); Python v3.11.11 |
| R / Python libraries | data.table v1.17.8; splines; ggplot2 v4.0.3; patchwork v1.3.2; pandas v3.0.3; numpy v1.26.4; scipy v1.16.0; matplotlib v3.10.5 |

Figures (`figures/`) were rendered locally with R (conda env `r`, ragg, Arial from mscorefonts) and Python matplotlib.
