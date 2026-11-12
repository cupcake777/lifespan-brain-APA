# NG Results Register

Audit date: 2026-10-15

This file separates verified current model outputs from release-final results.
The authoritative location is branch-specific: Furong
(`${APA_ROOT}/data/07_coloc`) is used for remote production
branches, while local `${APA_ROOT}/data/07_coloc/11_eqtl` is used for the
completed main maxCPP rerun. A model output can be numerically verified while
still requiring a summary refresh or an upstream-dependent rerun. The current
Module 01 -> 02 -> 03 rerun chain is now accepted; remaining holds are release
provenance and figure/text synchronization.

## CONFIRMED: current model outputs and corrected interpretation

### 1. Binary-window APA-QTL sLDSC (`09_sLDSC`)

Producer/run evidence:

- Producer: `${APA_ROOT}/data/07_coloc/09_sLDSC/src/run_sldsc_pipeline.sbatch`
- Run log: `${APA_ROOT}/data/07_coloc/09_sLDSC/slurm_logs/sldsc_pipeline_83320.out`
- Completion: the 57 unaffected prenatal/postnatal/lifespan regressions are
  retained; corrected job 46274987 completed 44 interaction annotation files,
  44 interaction LD-score files and 19/19 interaction regressions. The 19
  corrected results replaced the superseded interaction rows before resummary.
- Summary: `${APA_ROOT}/data/07_coloc/09_sLDSC/output/summary/sldsc_enrichment_summary.tsv`
- Estimand: binary APA-QTL variant-window annotations (`prenatal`, `postnatal`,
  `lifespan` union, and `interaction`) conditioned on ancestry-specific LDSC
  baseline annotations. This is not the joint APA/eQTL maxCPP model.

The refreshed summary applies BH within each population-trait GWAS across its
four annotations, giving 13/76 current-snapshot hits. The eight core figures
and four annotation figures were regenerated from this summary; their layout
was visually reviewed. The refreshed table and figures are current for the
existing binary-window model outputs:

| Population | Trait | Annotation | Enrichment | SE | P | within-GWAS FDR |
|---|---|---|---:|---:|---:|---:|
| EUR | BIP | lifespan APA-QTL | 1.5849 | 0.1620 | 4.03e-4 | 8.93e-4 |
| EUR | BIP | postnatal APA-QTL | 1.5913 | 0.1652 | 4.47e-4 | 8.93e-4 |
| EAS | BIP | lifespan APA-QTL | 2.6682 | 0.7112 | 1.47e-3 | 0.00366 |
| EAS | BIP | postnatal APA-QTL | 2.6547 | 0.7128 | 1.83e-3 | 0.00366 |
| EUR | SCZ | postnatal APA-QTL | 1.4128 | 0.1328 | 2.05e-3 | 0.00500 |
| EUR | SCZ | lifespan APA-QTL | 1.3945 | 0.1295 | 2.50e-3 | 0.00500 |
| EAS | ALS | prenatal APA-QTL | -3.1431 | 2.2893 | 4.15e-3 | 0.01661 |
| EUR | BIP | prenatal APA-QTL | 1.8313 | 0.3302 | 1.34e-2 | 0.01790 |
| EAS | BIP | prenatal APA-QTL | 3.5951 | 1.3289 | 2.78e-2 | 0.03712 |
| EAS | PD | lifespan APA-QTL | 3.8089 | 2.0603 | 2.20e-2 | 0.04730 |
| EAS | PD | postnatal APA-QTL | 3.7941 | 2.0492 | 2.37e-2 | 0.04730 |
| EUR | PD | lifespan APA-QTL | 1.9054 | 0.4130 | 1.98e-2 | 0.04577 |
| EUR | PD | postnatal APA-QTL | 1.8901 | 0.4144 | 2.29e-2 | 0.04577 |

Safe interpretation:

- The BIP and EUR-SCZ rows are positive heritability enrichment signals under
  this binary-window estimand.
- The EAS ALS prenatal row is a negative, highly uncertain estimate
  (`-3.143 +/- 2.289`, EAS ALS GWAS N = 4,084). Report it as an unstable
  depletion-like signal, not as a robust biological depletion claim.
- EUR and EAS use different baseline annotation sets (EUR baselineLD v2.2,
  97 annotations; EAS local 83-annotation baseline). Do not interpret the
  cross-population effect sizes as a like-for-like comparison.
- The exploratory 19-test postnatal-versus-prenatal stage comparison has zero
  nominal P < 0.05 and zero BH FDR < 0.05.
- None of the 19 corrected interaction rows passes within-GWAS BH; the minimum
  interaction P is 0.06763 (EUR SCZ; within-GWAS FDR 0.09017).

Figure note: all 12 PNG/PDF figure pairs were regenerated. The plotting log
reports zero clipped labels and zero same-axes text overlaps for every panel;
the previously flagged panels were also visually inspected.

### 2. Joint APA-QTL/eQTL maxCPP sLDSC (`11_eqtl`)

Producer/run evidence:

- Producer: `${APA_ROOT}/data/07_coloc/11_eqtl/03_run_sldsc_v2.sh`
- Per-regression logs: `${APA_ROOT}/data/07_coloc/11_eqtl/sldsc/logs/h2.maxcpp.*.log`
- Summary: `${APA_ROOT}/data/06_disease/03_sldsc_tau_star.tsv`
- Scope: 38 population-stage-trait regressions (19 trait-population combinations
  x 2 stages), each with APA-QTL and eQTL coefficients, giving 76 coefficient
  rows; `well_powered` is the explicit `h2 / h2_se >= 7` gate.

The summary and four figure pairs have been regenerated from the existing
coefficient p-values. Correction is within each population-stage-trait
regression across its two annotation coefficients:

- 76 coefficient rows were produced.
- 32 rows pass the well-powered gate.
- Each population-stage-trait regression is corrected separately across its two
  annotation coefficients. A bounded recomputation gives 9/76 coefficients at
  within-family BH FDR < 0.05; 8 pass the well-powered gate and one additional
  EUR INS postnatal eQTL coefficient does not.
- These significant coefficients do not establish an APA-QTL-versus-
  eQTL coefficient difference: this model estimates each annotation's
  coefficient separately and is not a direct coefficient-contrast test.

Safe manuscript wording: "In the joint maxCPP sLDSC model, 32/76 coefficient
tests were well-powered. Using BH correction within each GWAS-stage regression
across the APA-QTL and eQTL coefficients, 9 coefficients passed FDR < 0.05;
eight were well-powered and one EUR INS postnatal eQTL coefficient did not pass
the power gate."
Do not write that APA-QTL and eQTL contributions differ without a direct contrast.
The refreshed `03_sldsc_tau_star.tsv` has 76 rows, 38 two-coefficient families,
and no cross-GWAS/global FDR gate. The 38 regressions were not repeated. Four
figure pairs were regenerated; the heatmap caption was shortened after visual
review identified clipping.

### 3. MAGMA competitive gene-set enrichment (`10_magma`)

The current table at
`${APA_ROOT}/data/07_coloc/10_magma/enrichment/magma_gsa_summary.tsv`
contains all 57 expected rows (19 GWAS x 3 gene sets) with no missing P values.
Five current-snapshot rows pass within-GWAS BH:

| GWAS | Gene set | Genes | P | within-GWAS FDR |
|---|---|---:|---:|---:|
| EAS BIP | prenatal | 501 | 0.020748 | 0.031122 |
| EAS BIP | postnatal | 1,764 | 0.013585 | 0.031122 |
| EUR BIP | postnatal | 1,821 | 0.0066987 | 0.0200961 |
| EUR SCZ | prenatal | 526 | 0.0012054 | 0.0036162 |
| EUR SCZ | postnatal | 1,824 | 0.019127 | 0.0286905 |

The refreshed plotting run regenerated two figure pairs and a source summary
using BH within each GWAS across the three gene sets; cross-GWAS `fdr_global` is
not a primary gate. The five rows describe the current 57-row snapshot only.
The interaction set now uses the accepted 193-group Module 05 export (154 mapped
Entrez genes). None of its 19 rows is significant; its minimum P is 0.10599.

### 4. Genotype-by-stage interaction (`03_interaction`)

The label-aligned rerun passed all 12 eigenMT diagnostic checks and reproduced
198 significant groups. HC3 refitting retained 193 core groups and excluded 5.
The cohort-LD SuSiE branch produced 386 stage-group summaries, including 99
usable fits, and 2,397 credible-set/mechanism rows across 79 groups. Canonical
consumer tables and the refreshed Figure 4 now use these accepted outputs.

### 5. Protein measured-background ORA (`12_protein`)

The live local and Furong tables are byte-identical. One GO BP term passes the
database-specific BH result already emitted by `clusterProfiler`:

| Database | Term | P | within-database FDR |
|---|---|---:|---:|
| GO BP | phosphatidylinositol-3-phosphate biosynthetic process | 9.33e-5 | 0.03814 |

The current rerun table has 641 terms and one GO BP term at within-database FDR
< 0.05; the across-database gate is not the primary correction family. The table
was regenerated after accepted Modules 01/02 inputs; figure and manuscript text
still need synchronization from this current output.

### 6. Lifespan functional enrichment (`02_pattern`)

The current Furong producer already applies BH within each annotation library and
gene architecture. It does not use the stale local global-FDR workflow. In the
current Furong snapshot, Shared HPO Autism has OR 3.62, library FDR 0.00894,
adjusted OR 2.18 and adjusted FDR 0.03607. These replace the older global-FDR
wording, but the complete Module 02 result remains held until the donor-level
trajectory rerun.

### 7. Branches without a cross-trait correction defect

- TWAS corrects each GWAS-stage feature family separately.
- SMR corrects within `qtl_type x gwas_id`.
- Interaction functional enrichment corrects within annotation library.
- Module 04 coloc, mediation and mashr results are not invalidated by the
  cross-trait-FDR issue; their existing provenance and estimand gates remain.

## HOLD: current tables exist, but not final NG-safe results

These values are recorded for traceability only. They must not be copied into
the main Results section as final numbers until the listed release closure is complete.

| Branch | Current observed table value | Why held |
|---|---|---|
| Module 01 RNA/APA | accepted ComBat_REF matrix: 63,097 loci x 815 samples | Figure and manuscript counts still need synchronization. |
| Module 02 trajectories | accepted donor-level output: 63,097 sites, 750 donors, 39,580 dynamic/clustered sites | Current trajectory figures and enrichment panels still need refresh/freeze. |
| Module 03 protein | current rerun: 20,450 primary pairs, 46 beta-gate pairs, 641 ORA terms | Current outputs exist; publication figure/text synchronization remains. |
| Interaction maxCPP | 19-row `interaction_maxcpp_tau_star.tsv`, category `L2_1` | Collector is accepted; interaction-dependent figure/text freeze remains. |
| Module 04 coloc | Corrected eQTL mapping: 4,735 prenatal and 13,929 postnatal eGenes; PP4 >= 0.7: prenatal 186/665, postnatal 655/2,229; classification gate PP_H4 >= 0.8 | Refreshed Oct 15 from corrected ComBat eQTL inputs; keep stage eGene counts separate from APA-locus coloc denominators and PP_H4 classification. |
| Module 04 mediation | 772 pairs total (164 prenatal, 608 postnatal); 2,000 donor resamples | Refreshed Oct 15. P-values use delta-method/normal approximation; CIs are percentile bootstrap. |
| Module 04 mashr | Strong-set sharing: apaQTL 0.84957 (12,918-row strong set), eQTL 0.63294 (54,880-row strong set); detected counts 2,792 and 20,901 | Refreshed Oct 15; detected-pair counts and strong-set sharing denominators are distinct. |
| Module 07 coloc/SMR/functional | Corrected internal eQTL contract: 298,563 rows across 4,524 candidates; 1,375 same-gene supported candidates. Disease integration: 3,084 comparisons and 321,097 evidence tiers (Furong job 46325111, Oct 16). | New integration outputs and producer summary are current; public reconstruction remains held because wrappers depend on external inputs and executable roots. |

## Execution classification

- **Rerun:** Module 01 -> Module 02 -> Module 03 is accepted through job 84024.
  Interaction maxCPP chromosomes/regressions were already complete; collector
  job 84026 wrote and validated the 19-row tau-star summary.
- **Recheck:** Modules 00, 04, 07 and 13 plus TWAS/SMR provenance and output
  contracts are complete. TWAS now has a unique 38-row master summary; external
  input manifests remain a release-packaging issue, not a model-rerun reason.
- **Resummary:** binary sLDSC, main maxCPP, MAGMA primary wording/figures,
  protein ORA and all manuscript statements that currently use a global null.

## Release decision

The repository is not yet a final-results package. Binary sLDSC, main maxCPP,
corrected Module 05, corrected MAGMA, the Module 01 -> 02 -> 03 chain, and
interaction maxCPP collection are accepted. Release-final packaging still waits
on provenance, cell-QTL delivery, and synchronized figures/text.

## Current closure record (2026-10-15)

- Module 01 acceptance: `accepted_and_promoted`, ComBat_REF, 63,097 loci x 815 samples, zero nonfinite values.
- Module 02 acceptance: 63,097 sites, 750 donors, 39,580 dynamic/clustered sites.
- Module 03 corrected tail: Slurm job 84024 completed `0:0` with `MODULE03_TAIL_COMPLETE`; the current bootstrap, triad, validation, and figure outputs are present.
- Interaction maxCPP collection: Slurm job 84026 completed `0:0`; `interaction_maxcpp_tau_star.tsv` contains 19 rows, 19 unique population-trait keys, and the producer category `L2_1`.
- Release holds: current interaction-dependent figures and manuscript prose still need synchronization; the accepted cell-QTL QC summary and public input/environment provenance are not yet attached to this publication register.

This register supersedes older narrative claims that conflict with the current
tables, including the historical statement that EUR SCZ postnatal APA-QTL was
FDR-significant in the maxCPP branch.
