# Genotype-by-stage interaction analysis

## Purpose

This module tests whether genetic effects on APA differ between prenatal and postnatal brain and evaluates the robustness, external support and functional context of the detected interactions.

## Execution order

1. `01_preprocess.R` — prepare APA phenotypes, covariates and genotype inputs.
2. `02_tensorqtl.sh` — run nominal cis genotype-by-stage interaction mapping.
3. `03_group_eigenmt.py` — correct grouped tests with LD-aware eigenMT and control FDR.
4. `04_eigenmt_diagnostics.py` — audit effective test counts and group-level calls.
5. `05_hc3_refit.py` — refit lead associations with HC3 standard errors and evaluate boundary/genotype support.
6. `06_external_comparison.py` — compare prespecified leads with an external adult brain molecular-QTL map.
7. `07_functional_enrichment.R` — test enrichment using declared foregrounds and backgrounds.
8. `08_representative_loci.R` — prepare genotype-stratified locus summaries.

## Inputs

- harmonized APA phenotypes;
- genotypes, stage labels and model covariates;
- nominal interaction results;
- external adult brain QTL summaries;
- annotation and gene-set resources.

## Principal outputs

- group-level interaction calls;
- effective-test and FDR diagnostics;
- HC3-refitted effects and sensitivity summaries;
- descriptive reversed, amplified and attenuated effect classes;
- external comparison and enrichment results;
- representative locus tables.

## Interpretation boundary

The reversed, amplified and attenuated labels summarize estimated stage-specific effects rather than establishing a molecular mechanism on their own.

These outputs support Figure 4, parts of Figure 6 and Extended Data Figure 8.
