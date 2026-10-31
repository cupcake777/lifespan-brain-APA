# Disease genetics

## Purpose

This module evaluates whether APA-regulatory annotations and genes are enriched for complex-trait association and tests colocalization between molecular QTLs and GWAS signals.

## Directory structure

```text
sldsc/           stratified LD score regression input and result processing
magma/           gene-set definitions and MAGMA enrichment
colocalization/  QTL input preparation and study-specific colocalization summaries
```

## S-LDSC

1. `sldsc/01_format_qtl.py` — format molecular-QTL summary statistics.
2. `sldsc/02_make_annotations.py` — construct merged variant windows.
3. `sldsc/03_make_ldsc_annotations.py` — write chromosome-level annotations.
4. `sldsc/04_format_gwas.py` — harmonize GWAS inputs for regression.
5. run S-LDSC using the upstream software installation.
6. `sldsc/05_collect_results.py` and `sldsc/06_plot_results.py` — collect, correct and summarize results.

## MAGMA

1. `magma/01_define_gene_sets.R` — define APA-QTL gene sets.
2. `magma/02_run_magma.sh` — run gene analysis and set enrichment.
3. `magma/03_collect_enrichment.R` — collect and summarize enrichment results.

## Colocalization

1. `colocalization/01_prepare_qtl_inputs.sh` — create indexed molecular-QTL inputs.
2. `colocalization/02_run_coloc.sh` — call an external EasyColoc installation through an explicit path.
3. `colocalization/03_plot_summary.R` and `colocalization/04_plot_matrix.R` — summarize gene-level and disease-level results.

EasyColoc, S-LDSC and MAGMA are installed separately. Their versions are recorded with the corresponding analysis run.

## Inputs

- stage-specific molecular-QTL summaries;
- harmonized GWAS summary statistics;
- LD reference panels matched to the analyzed ancestry;
- gene annotations and baseline S-LDSC resources.

## Principal outputs

- annotation-level S-LDSC enrichment;
- MAGMA gene-set associations;
- molecular-QTL/GWAS colocalization summaries;
- disease-by-gene and representative-locus tables.

These outputs support Figures 5–6 and Extended Data Figure 9.
