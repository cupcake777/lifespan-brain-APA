# Code accompanying the lifespan human brain APA study

This repository contains the analysis code for a study of alternative polyadenylation (APA) across the human lifespan and its genetic regulation in prenatal and postnatal brain.

The repository is organized for scientific review: analytical decisions, execution order, model definitions and links between analysis stages are documented alongside the code that implements them.

## Start here

A reviewer can follow the study in three steps:

1. Read the workflow below for the scientific order of the analyses.
2. Open the README in an analysis module for its inputs, methods, script order and outputs.
3. Use [`docs/figure_map.md`](docs/figure_map.md) to trace a manuscript figure back to the analysis modules that support it.

## Analysis workflow

```text
01 APA quantification and harmonization
          |
02 Lifespan trajectory modelling
          |
03 APA–RNA–protein integration
          |
04 Prenatal and postnatal apaQTL/eQTL mapping
          |
05 Genotype-by-stage interaction analysis
          |
06 Disease enrichment and colocalization
          |
   Manuscript figure generation
```

| Step | Scientific purpose | Documentation |
|---|---|---|
| 01 | Quantify APA and harmonize batches | [`analysis/01_apa_quantification/`](analysis/01_apa_quantification/) |
| 02 | Identify and validate lifespan APA trajectories | [`analysis/02_lifespan_trajectories/`](analysis/02_lifespan_trajectories/) |
| 03 | Compare APA trajectories with RNA and protein | [`analysis/03_multiomics_integration/`](analysis/03_multiomics_integration/) |
| 04 | Map stage-specific apaQTLs and eQTLs | [`analysis/04_stage_specific_qtl/`](analysis/04_stage_specific_qtl/) |
| 05 | Test genotype-by-stage effects | [`analysis/05_genotype_by_stage/`](analysis/05_genotype_by_stage/) |
| 06 | Relate APA regulation to complex-trait genetics | [`analysis/06_disease_genetics/`](analysis/06_disease_genetics/) |
| Figures | Render main and Extended Data figures from finalized analysis outputs | [`figure_generation/`](figure_generation/) |

Numeric prefixes describe execution order. Within each module, scripts are likewise numbered in the order in which they are normally run.

## Repository layout

```text
analysis/                  study-specific analysis code and module READMEs
figure_generation/         final manuscript figure code
config/                    local path template
scripts/                   repository checks and figure entry point
docs/
  figure_map.md            manuscript figure to analysis-code map
  input_layout.md          expected input and intermediate directory layout
```

## Inputs and outputs

The code expects a local data root defined by `APA_ROOT`. Each module README documents the required inputs and the outputs passed to the next module. Paths are expressed relative to that data root rather than to an institutional server.
See [`docs/input_layout.md`](docs/input_layout.md) for the common directory contract.

## Configuration

```bash
cp config/paths.example.sh config/paths.local.sh
source config/paths.local.sh
```

The principal variables are:

- `APA_CODE_ROOT`: this repository;
- `APA_ROOT`: study data and reference root;
- `APA_FIGURE_DIR`: output directory for rendered figures;
- `APA_FONT_DIR`: optional publication-font directory.

Cluster-scale scripts also document the external tools they call. Third-party packages and pipelines are linked rather than copied, except for the two DaPars2-derived source files documented in `analysis/01_apa_quantification/README.md`.

## Figure generation

Analysis and plotting code are maintained separately. Figure scripts consume finalized outputs from the analysis modules, keeping statistical inference distinct from presentation.

```bash
bash scripts/render_figures.sh figure_3
bash scripts/render_figures.sh extended_data_7
```

See [`figure_generation/README.md`](figure_generation/README.md) and the panel-level [`docs/figure_map.md`](docs/figure_map.md).

## Software environment

The project uses R, Python and established genomics tools, including DaPars2, tensorQTL, mashr, coloc, susieR, S-LDSC and MAGMA. Module READMEs identify method-specific requirements. `environment.yml` describes the portable figure-generation environment; it is not presented as a container for the full HPC workflow.
