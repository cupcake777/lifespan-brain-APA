# APA–RNA–protein integration

## Purpose

This module tests whether lifespan APA trajectories track changes in host-gene RNA and protein abundance, defines matched molecular-layer universes, and evaluates robustness to model and expression definitions.

## Execution order

1. `01_define_universe.R` — define feature and donor sets with matched APA, RNA and protein evidence.
2. `02_fit_layer_trajectories.R` — fit comparable trajectories for all three molecular layers.
3. `03_integrate_layers.R` — estimate cross-layer coupling and rank informative APA–RNA–protein relationships.
4. `04_enrichment.R` — test functional enrichment against the measured analysis background.
5. `05_gam_sensitivity.R` — compare fixed-spline and adaptive-GAM trajectory shapes.
6. `06_prepare_matched_matrices.R` — generate explicitly matched matrices used by integration analyses.
7. `07_utr_insensitive_expression.R` — repeat key comparisons with expression definitions less sensitive to 3′UTR usage.

## Inputs

- lifespan APA trajectories from module 02;
- harmonized RNA abundance;
- protein abundance and its cohort metadata;
- transcript annotation and external functional resources.

## Principal outputs

- a declared matched analysis universe;
- layer-specific fitted trajectories;
- APA–RNA and APA–protein coupling summaries;
- ranked cross-layer relationships;
- model and UTR-sensitivity analyses;
- functional-enrichment results with measured backgrounds.

These outputs support Figure 2 and Extended Data Figure 5.
