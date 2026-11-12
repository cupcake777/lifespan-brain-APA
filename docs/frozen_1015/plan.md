# APA Study Analysis Plan

## Scope and execution boundary

This plan records the current study workflow before any new flow edits or reruns.
The publication repository under `pub/` is a local mirror and organization layer.
The authoritative producer workspace and large inputs may reside on Furong at
`${APA_HPC_HOME}/project/work/apa`; local absence is therefore a
remote-execution boundary, not evidence that an analysis or dataset is missing.

`14_context` on Furong is a dedicated checking/context directory and is not part
of the primary analysis workflow.

## Main workflow

```text
WGS and genotype QC
  -> ChinaMAP imputation checkpoint
  -> stage PGEN and genotype PCs

RNA-seq, isoform and APA inputs
  -> sample/feature QC and alias harmonization
  -> LifeStage-protected ComBat_REF
  -> donor-level APA/PDUI matrix

donor-level APA
  -> lifespan GAM trajectories and C1/C2 clusters
  -> functional enrichment, constraint and external validation

stage APA phenotypes + stage genotype
  -> prenatal/postnatal apaQTL
  -> eQTL overlap, coloc, mediation, mashr and power control
  -> GxStage interaction, HC3 refitting and saturation/gating analyses

apaQTL outputs
  -> binary-window sLDSC (`09_sLDSC`)
  -> joint APA/eQTL maxCPP sLDSC (`11_eqtl`)
  -> MAGMA gene-set enrichment (`10_magma`)
  -> transcript-level PDUI TWAS (`08_twas`)
  -> coloc/SMR and locus follow-up

APA/RNA/protein matched data
  -> cross-layer association, ranking, robustness and functional enrichment

external scAPA resources
  -> cell-type and tissue context for selected APA signals
```

## Module order and purpose

### 00. Genotype preparation

GATK joint genotyping, VQSR and sample/variant QC are followed by the external
ChinaMAP imputation checkpoint. The downloaded imputed VCF is filtered for the
stage cohorts, converted to stage-specific PGEN, and used to compute genotype PCs.

### 01. RNA, isoform and APA preprocessing

Use the merged matrices `all.912.expr.count/tpm`, `all.912.isoform.count/tpm`,
`all.apa.count`, and the sample metadata. No per-sample RSEM files are required
for the retained-matrix workflow. Apply verified aliases, RIN/sex/QC filters,
DaPars2/PDUI processing where needed, bridge evaluation and LifeStage-protected
ComBat_REF. The corrected v2 outputs are the inputs for all downstream reruns.

### 02. Lifespan trajectories

Average technical replicates at donor level, fit age GAMs with sex and donor RIN,
define supported age ranges and directional events, cluster trajectories into C1/C2,
and run BrainSpan, functional-enrichment and constraint analyses. This module
produces the Figure 1 evidence.

### 03. Multi-omics uncoupling and protein validation

Residualize APA and RNA, preprocess matched DIA proteomics, fit protein models
with APA/RNA/age terms, rank the tested site-protein pairs, and run robustness,
bootstrap, cross-validation and functional enrichment analyses. The current
cross-layer universe is about 20,447 pairs with 13 ranked pairs retained for the
primary story; this remains association evidence, not causal mediation.
Because the matched APA/RNA inputs and the dynamic-site universe come from
Modules 01 and 02, the protein branch must be rerun after the corrected
preprocessing and donor-level trajectory outputs are frozen.

### 04. Stage-specific APA-QTL and eQTL integration

Map prenatal and postnatal cis-apaQTLs with stage-specific genotype and covariates;
map eQTLs in the same samples; run coloc, PP4-gated mediation, mashr effect
consistency, matched-power/downsampling, BrainGVEX replication and sequence-level
PAS/RBP/miRNA follow-up. This is the main APA-vs-eQTL evidence chain.

### 05. Genotype-by-stage interaction

Fit `APA ~ G + Stage + G:Stage + covariates` in the combined cohort, apply eigenMT,
HC3 robust refitting and expression-gating tests, then inspect representative
gating loci and ATP2B1 saturation. Corrected stage-mask results must propagate to
all dependent tables and Figure 4 panels.

### 06/08/09/10/11. Disease association branches

- `09_sLDSC`: binary APA-QTL annotation enrichment, including prenatal,
  postnatal, lifespan and interaction annotations. This is a separate binary-window
  estimand. BH is applied within each population-trait GWAS across its four
  annotation tests; different GWAS traits are not one correction family.
- `11_eqtl`: joint APA-QTL/eQTL maxCPP with ancestry-matched baseline annotations;
  current scope is 19 trait-population combinations x 2 stages = 38 regressions
  and 76 coefficient rows. Each regression is its own two-coefficient BH family.
  It is not a direct APA-vs-eQTL coefficient-difference test.
- `10_magma`: competitive enrichment of prenatal, postnatal and interaction
  APA-QTL gene sets across GWAS traits. BH is applied across the three gene sets
  within each GWAS; a cross-GWAS global FDR is not a primary decision rule.
- `08_twas`: transcript-level PDUI PrediXcan-style model training and GWAS
  association, with each GWAS-stage analysis corrected separately.

### 07 and supporting locus modules

Use coloc/coloc.susie, SMR/HEIDI, GTEx/CHB eQTL context and selected locus plots
to move from genome-wide disease enrichment to SCAF1, NRBP2, ATP2B1, 12q24.31,
ARL3 and other report-facing candidates. `04_coloc`, `06_smr`, `07_integrate` and
`05_func` are producer/support branches for this locus and mechanism layer.

### 13. scAPA context

Use external scAPA, RBP and miRNA resources to annotate tissue/cell-type context
for selected APA events. This is contextual support and does not replace the bulk
APA-QTL or stage-interaction analyses.

## Figure mapping

- Figure 1: Module 02 lifespan APA trajectories and functional classes.
- Figure 2: Module 03 / `12_protein` APA-RNA-protein uncoupling.
- Figure 3: Module 04 / `11_eqtl` stage apaQTL, coloc, mediation and sharing.
- Figure 4: Module 05 GxStage, gating and saturation.
- Figure 5: Modules 06, 08, 09, 10 and 11 disease-level evidence and candidates.
- Figure 6: locus-level coloc/SMR/functional follow-up and selected mechanisms.

## Current workflow position

The scientific design and main output branches exist, but the workflow is not
release-closed. Two independent rerun chains can start in parallel:

```text
Module 01 corrected preprocessing
  -> Module 02 donor-level trajectories and enrichment
  -> Module 03 matched APA/RNA/protein analysis

Module 05 label-aligned interaction analysis
  -> interaction-dependent MAGMA, sLDSC and maxCPP inputs
```

The completed binary-sLDSC, main maxCPP, MAGMA and protein-ORA tables also need
correction-family resummaries. Those resummaries reuse existing model p-values;
they are not a reason to rerun unaffected regressions. New `1012` exploratory
analyses and final figure synchronization wait for both rerun chains.

## Execution status (2026-10-13 18:35 CST)

- Module 01 APA correction is running as Slurm job `83903`. Three parametric
  ComBat branches are cached; exact nonparametric `ComBat_REF_NP` is active.
- Job `83918` is dependency-held after Module 01. It validates 63,097 unique
  loci, 815 aligned retained samples, finite values, method identity and all four
  ComBat caches before atomic promotion, then runs and validates Module 02.
- Module 03 job `83919` is dependency-held after accepted Module 02. Its v2
  wrapper includes matched-input regeneration, primary models, robustness,
  donor bootstrap, cross-validation, enrichment and Figure 2 producers.
- Corrected Module 05, Figure 4, interaction MAGMA and binary-window interaction
  sLDSC are complete. The full interaction maxCPP chromosome array `83904` is
  active; dependent annotation/regression job `83905` is queued.
- Group 2 scientific-run rechecks are complete: retained Module 04, TWAS, SMR,
  functional and Module 13 jobs have successful scheduler/log evidence. Their
  remaining external-input and revision-pinning gaps are release packaging, not
  reasons to rerun statistically complete models.
- Group 3 resummaries are complete for binary sLDSC (13/76 within-family hits),
  main maxCPP (8/76 coefficient rows), corrected MAGMA (5/57 within-GWAS hits),
  and the 38-row TWAS master table. Protein ORA refresh waits for Module 03.

## Parallel task division

Window ownership is exclusive while the three groups run: `rerun` owns producer
scripts and regenerated result tables; `recheck` is read-only and records
provenance findings; `resummary` owns correction-family summaries, figures and
manuscript-facing text, and consumes only outputs accepted by `rerun`. Only the
`resummary` window updates `results.md` during parallel execution.

### Group 1: rerun

1. Run Module 01 v2 RNA/isoform/APA preprocessing and validate the aligned matrices.
2. Run Module 02 from donor aggregation through trajectories, enrichment and figures.
3. Run Module 03 after Modules 01/02 are accepted.
4. In parallel with steps 1-3, run the corrected Module 05 eigenMT, HC3,
   fine-mapping and mechanism chain.
5. Propagate the accepted Module 05 interaction set into MAGMA and any
   interaction-dependent sLDSC/maxCPP branches.

### Group 2: recheck

1. Check Module 00 handoff documentation and retained genotype contracts.
2. Recheck Module 04 coloc, mediation and mashr producer-to-result provenance;
   their trait-specific correction families do not require rerun for the
   cross-trait-FDR issue.
3. Recheck Module 07 locus/SMR/functional manifests and successful-run evidence.
4. Recheck TWAS, SMR and Module 13 contextual outputs without changing their
   already trait-scoped multiple-testing families.

### Group 3: resummary

1. Recompute binary-sLDSC BH within each population-trait GWAS.
2. Regenerate main maxCPP summaries/figures using the two coefficients within
   each population-stage-trait regression as the family.
3. Remove global-across-GWAS MAGMA and global-across-database protein ORA gates
   from primary claims and refresh their tables/figures.
4. Synchronize `results.md`, figure notes and manuscript prose only from accepted
   result tables; do not retain global-null wording as a downgrade.

## Safety rules before changes

1. Inspect authoritative Furong scripts/results before rerunning a remote branch.
2. Use the login node only for shell inspection, network transfer and Slurm
   control. Run R/Python analysis on an allocated compute node (`cu01` when the
   allocation is verified) or through the project Slurm wrapper.
3. Match process/thread counts to the actual allocation and cap nested BLAS
   threads when process-level parallelism is used.
4. Keep `raw/` outside Git; `raw_manifest.tsv` records the current local data boundary.
5. Validate schemas, sample IDs, row/column keys and output manifests after each
   dependency-level rerun.
6. Do not promote a figure or manuscript number until its producer and current
   result table agree.

## Current workflow position (2026-10-15)

- Module 01 is accepted/promoted at 63,097 loci x 815 samples with zero nonfinite values.
- Module 02 is accepted at 63,097 sites, 750 donors, and 39,580 dynamic/clustered sites.
- Module 03 corrected tail job 84024 completed successfully; current bootstrap, triad, validation, and figure outputs are present.
- Module 05 accepted outputs remain promoted, and interaction maxCPP collection job 84026 completed a 19-row summary from the existing regression outputs.
- Remaining work is release synchronization and provenance closure: refresh dependent figures/text, deliver the accepted cell-QTL QC artifact, and record public input/environment contracts. No model rerun is pending in the current numerical chain.
