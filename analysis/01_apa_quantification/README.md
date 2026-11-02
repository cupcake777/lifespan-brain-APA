# APA quantification and cohort harmonization

## Purpose

This module converts alignment-derived coverage into a cohort APA matrix, evaluates batch-correction strategies using bridge samples, and produces the harmonized donor-level inputs used by the lifespan and QTL modules.

## Execution order

### APA quantification

1. `apa/00_bam_to_wig.sh` — generate strand-aware coverage and mapping summaries.
2. `apa/01_prepare_inputs.sh` — construct the 3′UTR annotation, read-depth list and DaPars2 configuration.
3. `apa/02_run_dapars2.sh` — run chromosome-level DaPars2 jobs.
4. `apa/03_merge_dapars2.sh` — merge autosomal results into a single APA matrix.

The Python components called by these wrappers are:

- `apa/DaPars2_Multi_Sample.py` — multi-sample DaPars2 inference;
- `apa/DaPars_Extract_Anno.py` — extract the DaPars-compatible 3′UTR annotation;
- `apa/extract_3utr.py` — derive transcript-level 3′UTR intervals;
- `apa/read_depth.py` — collect mapped-read depth for DaPars2 normalization.

The first two files are derived from DaPars2 v2.1 commit `fb81c6c` and remain under GNU GPL v2; see `apa/LICENSE.DaPars2`. The numbered shell scripts are study-specific workflow code.

### Cohort harmonization

5. `01_batch_correction.R` — compare ComBat, limma and RUV correction strategies.
6. `02_bridge_selection.R` — evaluate agreement on bridge donors and assess method selection with paired bootstrap resampling.
7. `03_cohort_overview.R` — summarize donor-weighted age and assay coverage.
8. `04_rna_collapse_check.R` — check within-batch expression-matrix collapse where required.

## Inputs

- coordinate-sorted RNA-seq alignments;
- sample sheet and mapped-read summaries;
- GRCh38/GENCODE transcript annotation;
- bridge-sample metadata and cohort covariates.

## Principal outputs

- merged per-sample PDUI matrix;
- harmonized APA matrix;
- donor-level APA matrix;
- batch-method diagnostics and bridge-pair summaries;
- cohort summary tables used by downstream modules and Figure 1.

## Review notes

The batch method is selected from bridge-donor agreement rather than from visual batch removal alone. Technical replicates can inform harmonization, but downstream lifespan inference uses donor-level observations.
