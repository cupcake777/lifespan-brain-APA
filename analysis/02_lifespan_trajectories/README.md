# Lifespan APA trajectories

## Purpose

This module identifies APA events associated with age, defines the age intervals supported by the observed cohort, groups directional trajectories, and evaluates their robustness and external support.

## Execution order

1. `01_prepare_covariates.R` — align donor-level APA measurements and covariates.
2. `02_fit_gam.R` — fit per-site age models using generalized additive models.
3. `03_define_supported_age_grid.R` — define age-grid regions supported by unique donors using leave-one-out and subsampling diagnostics.
4. `04_call_dynamic_events.R` — call directional events against individual-block permutation nulls.
5. `05_cluster_trajectories.R` — cluster fitted trajectories and quantify seed/bootstrap stability.
6. `06_functional_enrichment.R` — test functional annotations and gene constraint using the prespecified background.
7. `07_brainspan_validation.R` — apply the discovery trajectory template to external developmental brain data.
8. `08_brainspan_heatmap.R` — visualize matched external events.
9. `09_external_pattern_validation.R` — test external age-pattern replication with discovery labels fixed.

## Inputs

- harmonized donor-level PDUI matrix from module 01;
- donor age, sex, RIN and analysis covariates;
- gene annotation and enrichment resources;
- external BrainSpan/BrainGVEX APA profiles.

## Principal outputs

- fitted age curves and the supported age grid;
- age-associated and directional APA-event tables;
- trajectory-cluster assignments and stability summaries;
- functional-enrichment tables;
- external validation and sensitivity summaries.

## Statistical design visible to reviewers

- inference is performed at the donor level;
- model support is evaluated over the observed age distribution;
- dynamic-event thresholds are calibrated with individual-block permutations;
- trajectory labels are fixed before external validation;
- enrichment analyses retain an explicit tested background.

These outputs support Figure 1 and Extended Data Figures 1–4.
