# Public release preparation

This file defines the boundary between the working analysis repository and a
reviewer-facing code release. It does not replace the scientific evidence audit
in `audit.md`.

## Intended release contents

| Content | Decision | Notes |
|---|---|---|
| Root README and evidence registers | Include | `README.md`, `framework.md`, `audit.md`, `results.md`, `figure.md`, `plan.md` |
| Retained working producers | Local provenance only | `1010/` is not required in the transferred publication package |
| Publication scripts and smoke tests | Include after review | `1011/` is materialized and contains no filesystem symlinks |
| Provisional assembled figures | Include only when synchronized | Current files remain provisional until figure/manuscript freeze |
| Exploratory notes | Exclude from the canonical workflow | `1012/` may be archived separately if needed for history |
| Rerun/follow-up workspace | Selective | Review `1013/` file by file; local BrainGVEX matrix and sample list are ignored |
| Figure/RNA staging workspace | Do not publish as a second pipeline | Reconcile selected `1014/` sources to canonical producers first |
| Human molecular/genotype inputs | Exclude | `raw/` and protected `1013/` inputs must remain outside Git |
| Logs, caches and temporary files | Exclude | Covered by `.gitignore` patterns |

## Release gates

| Gate | Status | Acceptance evidence or required action |
|---|---|---|
| Numerical rerun closure | Pass | Current Module 01/02/03/05 and interaction maxCPP acceptance is recorded in `results.md` and `audit.md` |
| Private-data Git boundary | Pass locally | `raw/`, `1013/after.impute.all.raw`, and `1013/sample.Adult.final` are ignored; recheck before staging |
| Canonical module map | Pass | `README.md`, `framework.md`, and `1011/README.md` identify Modules 00-07 |
| EasyColoc dependency identity | Pass | Module 07 vendors the 98 tracked files from pinned commit `850178b5c99899d3e00b8655ffda716f705aa09c`, including its license/environment |
| Public input/config contract | Open | Replace or parameterize remaining HPC paths and provide de-identified config templates plus input/reference manifests |
| Locked software environment | Open | Capture tested R, Python and CLI versions in a machine-readable lockfile or container; do not infer unrecorded historical versions |
| Cell-QTL QC delivery | Open | Locate and validate the accepted `cell_qtl_covariate_qc_summary.tsv` artifact |
| Figure/manuscript synchronization | Open | Regenerate dependent panels and synchronize accepted denominators and claims before freeze |
| License | Owner input required | Select explicit code and documentation reuse terms |
| Citation metadata | Owner input required | Supply author list, title, version and DOI/repository URL for `CITATION.cff` |
| Data access statement | Owner input required | Supply public/controlled accessions and consent-compatible access language |

## Dependency boundary

Observed analysis dependencies include R/Bioconductor packages, Python
scientific packages, GATK 4.5.0.0, PLINK 2, bcftools, samtools, tensorQTL,
LDSC, MAGMA, SMR and EasyColoc. SOAPnuke v1.5.6 and ChinaMAP are external
checkpoints. The repository currently records these across producers and
method documents rather than in one tested lockfile.

The release environment must describe the environment that generated the
accepted outputs. A newly solved environment that merely installs similar
packages is not equivalent evidence.

## Pre-release checks

Run these checks from the repository root after the reviewed release set is
selected:

```bash
git status --short
git check-ignore -v raw/chb.full.vcf.gz \
  1013/after.impute.all.raw 1013/sample.Adult.final
bash 1011/tests/smoke_module00.sh
bash 1011/tests/smoke_module06.sh
git diff --check
```

Before tagging, also verify that no staged file contains personal contact
details, credentials, individual-level identifiers, unapproved data, or
filesystem symlinks, and that the manuscript numbers match the accepted result
and figure registers.
