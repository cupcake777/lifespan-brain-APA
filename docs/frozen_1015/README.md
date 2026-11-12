# APA publication analysis repository

This repository contains the analysis code and evidence registers for a study
of alternative polyadenylation (APA) across the human lifespan, its genetic
regulation, multi-omic relationships, and links to disease genetics.

Evidence status: **Partial reproducibility package**. The current numerical
rerun chain is closed, but the repository must not be released publicly until
the open gates in [`RELEASE.md`](RELEASE.md) are resolved.

## Analysis modules

The publication-facing scripts are under [`1011/`](1011/). This tree is a
materialized, symlink-free snapshot that can be transferred independently.

| Module | Scientific scope | Publication entrance | Current evidence status |
|---|---|---|---|
| `00` | WGS joint calling, QC, ChinaMAP imputation handoff, stage genotypes and PCs | `1011/00_genotype_qc_imputation/` | Partial; ChinaMAP remains an external web checkpoint |
| `01` | RNA-seq/isoform quantification, APA quantification and batch correction | `1011/01_rnaseq_apa_quant/` | Accepted/promoted corrected APA output |
| `02` | Lifespan trajectory fitting, clustering and external validation | `1011/02_lifespan_trajectories/` | Accepted donor-level trajectory output |
| `03` | APA-RNA-protein integration and robustness analyses | `1011/03_multiomics_uncoupling/` | Current corrected tail accepted |
| `04` | Stage-specific apaQTL/eQTL mapping, mashr and mediation | `1011/04_stage_specific_apaqtl/` | Numerical analyses complete; release provenance remains open |
| `05` | GxStage interaction mapping and mechanism follow-up | `1011/05_gxstage_interaction/` | Accepted interaction groups and downstream summaries |
| `06` | sLDSC, MAGMA and TWAS | `1011/06_disease_heritability_geneset/` | Numerical summaries current; public environment incomplete |
| `07` | Colocalization, SMR and functional follow-up | `1011/07_locus_colocalization_smr/` | Results complete; public configs and dependency restoration remain open |

## Repository layout

| Path | Role | Release handling |
|---|---|---|
| `1010/` | Retained analysis producers in their original working layout | Keep as local provenance; it is not required when transferring `1011/` |
| `1011/` | Symlink-free publication scripts, smoke tests and provisional assembled figures | Transfer this tree after final source/status review |
| `1012/` | Exploratory notes | Internal review; not a canonical analysis entrance |
| `1013/` | Rerun and focused follow-up workspace | Select code/results only; local BrainGVEX inputs are ignored |
| `1014/` | RNA refresh and figure-source staging | Reconcile selected sources to canonical locations before release |
| `raw/` | Local human molecular/genotype inputs | Never publish or stage in Git |
| `raw_manifest.tsv` | Local input inventory without data contents | Include only if approved by the data-governance review |

`1010/` preserves the working producer history. `1011/` contains materialized
copies selected for publication; use `1011/SOURCE_MANIFEST.tsv` to trace each
copied path back to its source.

## Documentation

| File | Purpose |
|---|---|
| [`framework.md`](framework.md) | Scientific workflow and method-to-script map |
| [`audit.md`](audit.md) | Evidence-backed correctness and release audit |
| [`results.md`](results.md) | Accepted result and denominator register |
| [`figure.md`](figure.md) | Figure claims, sources and synchronization status |
| [`plan.md`](plan.md) | Execution order, dependencies and remaining work |
| [`RELEASE.md`](RELEASE.md) | Public-release contents, exclusions and acceptance gates |

## Data and software boundaries

- Individual-level molecular, metadata and genotype inputs are local and
  Git-ignored. Public release requires approved accessions, controlled-access
  instructions, or de-identified derived inputs.
- Module 07 includes a tracked-file snapshot of EasyColoc from
  `https://github.com/cupcake777/EasyColoc.git`, commit
  `850178b5c99899d3e00b8655ffda716f705aa09c`, including its MIT license and
  environment file.
- Some canonical scripts still contain original HPC paths. These record the
  execution environment but are not yet a portable fresh-clone configuration.
- No unified locked environment is currently available. Do not describe this
  repository as end-to-end reproducible until the environment gate in
  `RELEASE.md` is closed.

## Bounded validation

The included smoke tests validate selected handoff contracts without using
human data or running the full analyses:

```bash
bash 1011/tests/smoke_module00.sh
bash 1011/tests/smoke_module06.sh
git diff --check
```

These tests do not establish full pipeline or scientific-result
reproducibility. See `audit.md` for the tested and untested boundaries.

## Citation and reuse

Release requires an explicit `LICENSE` and `CITATION.cff`. They are not created
until the repository owner supplies the license choice, author list, citation
title/version, and DOI or repository URL.
