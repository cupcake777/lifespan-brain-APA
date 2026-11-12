# Genetic regulation of alternative polyadenylation across the human lifespan — analysis code

Code for the manuscript *Genetic regulation of alternative polyadenylation across the human lifespan in prefrontal cortex*.
Analyses are organized as modules `00`–`07` in execution order, followed by `figures/`.
No individual-level data are included (see *Data*).

## Modules (run in this order)

| Module | Scope | Main entry point(s) | Manuscript |
|---|---|---|---|
| `00_geno` | WGS joint calling (GATK), variant QC, ChinaMAP imputation hand-off (web step), stage genotypes and PCs | `run_module00.sh`; `00_4_stage_filter.sh` after imputation | Methods (genotyping) |
| `01_quant` | RNA-seq quantification (STAR/RSEM), APA quantification (DaPars2), QC and batch correction (ComBat_REF; bridge-method selection) | `rna_quant.sh` → `rna_pre_step01.r` → `rna_pre_step02.r` → `batch_correction.R`, `bridge_select.R`; `apa/run_dapars2.sh` | ED Fig. 1–2 |
| `02_traj` | Donor-level PDUI trajectories (GAM), directional events, C1/C2 clustering, gene classes, enrichment, LOEUF, BrainSpan replication | `run.sh` or `run.slurm` → `01_cov.R` … `05_func.R`, `08_brainspan.R`; `followup/` | Fig. 1, ED Fig. 3–4 |
| `03_omics` | APA–mRNA–protein integration in 371 donors; trans-gene nulls; regulatory elements; UTR-insensitive expression | `run_protein.sh`; `01_universe.R` … `10_gam_shape.R`; `s2_element.py`; `r1_run.sh` | Fig. 2, ED Fig. 5 |
| `04_qtl` | Stage-specific apaQTL (tensorQTL), downsampling, eQTL, APA–eQTL coloc, mediation, mashr, adult-map (BrainGVEX, Nan 2026) comparisons | `apaqtl/qtl_pipe.sh`, `apaqtl/qtl_DS.sh`; `eqtl/01_preprocess.R` → `02_qtl_mapping.sh` → `02-2_coloc.R` → `04-1_mediate.R` → `05_mashr.R`; `followup/` | Fig. 3, ED Fig. 6–7 |
| `05_gxs` | Genotype-by-stage interaction mapping (tensorQTL nominal, eigenMT, HC3 refit), classes, gating and saturation | `00_preprocess.R` → `01_tensorqtl.sh` → `02_group_eigenmt.py` → `04_hc3_refit.py` → `06`–`09`; `reconstructed/recon_s4_gating.py` | Fig. 4, ED Fig. 8 |
| `06_disease` | S-LDSC (binary and joint apaQTL/eQTL maxCPP), MAGMA gene-set analysis, S-PrediXcan TWAS | `run_module06.sh`; `sldsc/01`–`08`; `magma/run_enrichment.R`; `twas/01`–`08`; `frozen/` | Fig. 5a,b, ED Fig. 9 |
| `07_coloc` | GWAS colocalization (coloc.abf, coloc.susie via vendored EasyColoc), SMR/HEIDI, APA-only classification, representative loci | `coloc/run_coloc.sh`; `smr/run_smr.sh`; `func/run_analysis.sh`; `handoff/run_handoff.sh` | Fig. 5c–e, Fig. 6 |
| `figures` | Final composites for Fig. 1–6 and ED Fig. 1–9, theme, locus helpers, data-prep producers | see `figures/README.md` | all figures |

`followup/` folders hold analyses run after the 2026-10-15 code freeze (Codex hand-offs and local analyses).

## Provenance
- `MANIFEST.tsv`: one row per file currently under `release/`, with its pre-restyle path (`old_path`, empty for new files), purpose, manuscript items and status (`final_verified`, `final_not_rerun`, `server_only_outputs`, `superseded_kept_for_provenance`). `results_final_path` points to analysis-tree folders that are not part of this release.
- `docs/number_registry.tsv`: historical v22 number audit with producer paths and verification status. It has not been reconciled with the current manuscript.
- `_manifests/` and `docs/frozen_1015/` are historical build records; they are kept unchanged and still use the old paths.
- `docs/frozen_1015/`: framework, audit, results and release notes of the 2026-10-15 freeze (`RELEASE.md` lists the release gates).

## Data
No individual-level data, sample metadata, genotypes or third-party reference data are included.
The raw sequencing data accession will be added upon release. Required external resources (GWAS summary statistics, LD reference panels, GENCODE v40,
PolyASite 2.0, GTEx v11, BrainGVEX, Nan et al. 2026 single-cell 3′aQTL statistics, gene-set libraries) are listed in
`_manifests/*EXTERNAL_INPUTS.tsv` with their sources.

Analysis root (`APA_ROOT`) layout:
- `data/<module>/` (`data/00_geno` … `data/07_coloc`) holds module inputs and outputs.
- `ref/` holds reference tables; `raw/` holds raw inputs; `figures/` holds figure outputs.
- The scripts in `release/` are read from `$APA_ROOT/release/`.

## Configuration
Paths come from environment variables. Copy `config/paths.sh` to `config/paths.local.sh`, set the values, and
`source config/paths.local.sh` before running any module. R scripts use `Sys.getenv()`, Python scripts use `os.environ`,
shell scripts read variables such as `${APA_HPC_HOME}`. `07_coloc/func/config.yml` values are expanded when loaded.

Slurm logs go to `$APA_ROOT/logs`. Create it once with `mkdir -p $APA_ROOT/logs` before submitting.
Some `#SBATCH --output` paths are relative (`logs/...`), so submit those jobs from `$APA_ROOT`.
No script uses `#SBATCH --chdir`.

## Environment
Software versions recorded for the accepted analyses are in `ENVIRONMENT.md` (not a tested lockfile).

## Known limitations
- Some outputs exist only on the analysis server (`server_only_outputs` in `MANIFEST.tsv`; e.g. TWAS model training, MAGMA BD rerun).
- Some inputs are not in this release, including the supplementary workbook used by `04_qtl/followup/nan2026_exact.R`. `05_gxs/00_preprocess.R` generates `LifeStage.txt` after sample QC and genotype alignment.
- `07_coloc/func/config.yml` records the complete server configuration, including evidence thresholds and model locations. Its absolute input and reference paths must be adapted to another environment before running `func/`.
- The stage-specific GO/KEGG producer `06_disease/magma/run_enrichment.R` filters the prenatal and postnatal gene sets separately. Enrichment tables from the earlier producer with the unfiltered 2,024-gene input must be regenerated before use.
- Fig. 2 trans-gene null summaries (S1): `figures/producers/s22_supp.R` uses a 2-d.f. age spline (dfres = n - 5). This reproduces the figure inputs `data/03_omics/S1_*_371.tsv` exactly (pi1 0.160 vs 0.118 for APA-protein, 0.67 vs 0.64 for APA-RNA). The original header's 3-d.f. model does not.
- Open items that need an author decision are marked `open` in the `notes` column of `MANIFEST.tsv`. The full list is in the fix-pass table (`restyle_bugs_1111.tsv`), which is not part of this release.
- Module `07` vendors EasyColoc (commit 850178b5) under its own MIT license.

## Tests
Six scripts. All use synthetic inputs; none needs human data.
- `tests/disease/run.sh`, `tests/omics/run.sh`, `tests/qtl/run.sh`: synthetic tests. They need `Rscript` on PATH or `RSCRIPT` set, and Python with numpy and pandas (`PYTHON`). The disease tests also need coloc 5.x and susieR.
- `tests/smoke_module00.sh`, `tests/smoke_module06.sh`: smoke tests. They create their own scratch `APA_ROOT`.
- `tests/check_portable_tree.sh`: checks required files, symlinks, shell syntax and the MANIFEST column count.

## License
The vendored tools retain their own licenses. Licensing of the original analysis code has not yet been specified by the authors.
