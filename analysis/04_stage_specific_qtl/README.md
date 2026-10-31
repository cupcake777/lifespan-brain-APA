# Stage-specific molecular QTL analysis

## Purpose

This module maps APA and expression QTLs separately in prenatal and postnatal brain, compares their effects across developmental stages, and evaluates overlap, colocalization, mediation and external support.

## Directory structure

```text
apa/           APA phenotype preparation and tensorQTL mapping
expression/    expression phenotype preparation and eQTL mapping
integration/   apaQTL/eQTL overlap, colocalization, mediation and mashr
validation/    comparison with an external single-cell 3′aQTL resource
```

## Execution order

### APA mapping

1. `apa/01_preprocess.R`
2. `apa/02_map_qtl.sh`
3. `apa/03_tensorqtl_driver.sh`
4. `apa/04_downsampling.sh`
5. `apa/05_screen_stage_pairs.R`

### Expression mapping

1. `expression/01_preprocess.R`
2. `expression/02_map_qtl.sh`

### Cross-phenotype integration

1. `integration/01_overlap.R`
2. `integration/02_extract_coloc_inputs.py`
3. `integration/03_colocalization.R`
4. `integration/04_coloc_summary.R`
5. `integration/05_mediation.R`
6. `integration/07_mashr.sh`
7. `integration/06_mashr.R`

The shell step prepares the mashr inputs; the R step fits and summarizes the model.

### External validation

- `validation/01_nan2026_exact.R` performs exact transcript/position/allele matching to the external resource.
- `validation/02_nan2026_model_sensitivity.R` records source limitations and sensitivity analyses when the lead table does not retain the inputs required for the fully adjusted pair-level model.

## Inputs

- harmonized APA and expression phenotypes;
- genotypes and genotype principal components;
- stage labels and covariates;
- GRCh38/GENCODE annotation;
- external 3′aQTL summary data.

## Principal outputs

- prenatal and postnatal apaQTL/eQTL maps;
- nominal, permutation and fine-mapping summaries;
- cross-stage sharing estimates;
- apaQTL–eQTL overlap, colocalization and mediation tables;
- exact external matches and sensitivity summaries.

These outputs support Figure 3 and Extended Data Figures 6–7.
