# Final Figure Framework

This is the NG-facing figure plan for the APA lifespan study. The narrative
order follows the biological argument:

```text
lifespan remodeling -> post-transcriptional layer -> inherited APA regulation
-> stage-dependent effects -> disease relevance -> locus-level mechanisms
```

The figure titles are declarative, but the wording is calibrated to the current
computational evidence. A number marked `TBD after rerun` must not be frozen in
the manuscript before the corrected producer chain finishes.

## Evidence labels

- **Confirmed**: current producer, result table and scope agree; safe for direct
  manuscript use with the stated estimand.
- **Provisional**: the signal is usable for figure planning, but producer
  provenance, a dependency rerun or a sensitivity check remains open.
- **Exploratory**: hypothesis-generating or underpowered; retain in a secondary
  panel, Extended Data or clearly qualified text.

Main figures may contain provisional panels while the analysis is being closed,
but only confirmed values are frozen in the final Results text.

## Main figures

### Figure 1 | Lifespan APA remodeling follows two coordinated developmental programs

**Message:** 3′UTR usage changes across the human lifespan in reproducible
prenatal and postnatal trajectory programs, with distinct functional contexts.

- **A. Cohort and life-course design**: prenatal to late-adult sampling,
  donor/technical-replicate structure, stage bridge and supported age ranges.
- **B. Age coverage and sampling density**: show the observed age distribution
  and the intervals used for inference; do not imply uniform coverage.
- **C. C1/C2 trajectory map**: donor-level APA/PDUI heatmap and fitted curves for
  the two directional trajectory classes.
- **D. Event and gene architecture**: C1-specific, C2-specific and shared
  events/genes, with counts updated after the corrected Module 01/02 rerun.
- **E. Functional separation**: pathway enrichment for C1 and C2 using the
  declared background and multiple-testing scope.
- **F. Constraint and disease context**: LOEUF/autism enrichment only if the
  corrected trajectory universe reproduces the current result; otherwise move
  this panel to Extended Data.

**Status:** the corrected donor-level producer is accepted at 63,097 sites and
750 donors. A-C and the general trajectory claim are eligible for refresh; event
counts and dependent enrichment panels must be regenerated from the accepted
current tables before final figure freeze. BrainSpan replication goes to Extended
Data.

### Figure 2 | APA variation forms a post-transcriptional layer that explains protein variation beyond RNA

**Message:** individual-level APA variation contains information not reducible
to bulk RNA abundance and is associated with protein variation in matched data.

- **A. RNA/APA decoupling**: residualized APA versus RNA relationships after the
  declared covariates; report effect sizes and sample scope rather than a target
  correlation value.
- **B. Matched multi-omics background**: RNA-protein and APA-protein dependence
  in the matched 366-sample set, with uncertainty shown.
- **C. Incremental protein variance**: APA contribution conditional on RNA,
  using the 20,447-pair universe and the prespecified 13-pair beta gate as
  separate denominators.
- **D. Representative cross-layer trajectories**: APA, RNA and protein for
  selected genes, with fitted uncertainty and no causal arrows.
- **E. INPP4A mechanism candidate**: 3′UTR extension and miRNA/RBP sequence
  context as a mechanistic hypothesis, not proof of post-transcriptional
  causality.

**Status:** the corrected matched-input chain and Module 03 tail are complete;
the panel remains `Provisional` until the current producer manifest is attached
to the publication package. Use “incremental association” or “additional
variance explained”, not “direct regulation” or “causal mediation”.

### Figure 3 | Stage-specific apaQTLs are widespread, developmentally conserved, and only partly shared with eQTLs

**Message:** inherited variation affects APA before and after birth, with strong
cross-stage conservation and a substantial component that is not represented by
the tested eQTL signal.

- **A. Prenatal and postnatal apaQTL landscape**: discovery counts, genomic
  distribution and PAS-proximal enrichment, with counts updated after the
  dependency closure.
- **B. APA/eQTL relationship**: coloc proportions with `PP4 >= 0.7` and the
  separate `PP_H4 >= 0.8` mediation subset; show denominators explicitly.
- **C. PAS sequence example**: an AAUAAA/PAS-disrupting allele with stage-wise
  PDUI effects, labelled as sequence-level support.
- **D. Cross-stage effect conservation**: matched prenatal/postnatal effect
  scatter, correlation and direction concordance with the tested pair universe.

**Status:** the claim is `Provisional` until the Module 04 producer-to-result
record is closed. BrainGVEX replication and downsampling belong in Extended
Data. Do not call this a direct APA-versus-eQTL superiority test; the current
maxCPP branch estimates separate coefficients.

### Figure 4 | Stage-dependent apaQTL effects arise from expression gating and phenotype saturation

**Message:** developmental interaction effects are explained by two interpretable
boundaries of the molecular phenotype: expression must be present, and the APA
measurement must retain dynamic range.

- **A. Interaction landscape**: stage-dependent effect classes, with all counts
  regenerated from the label-aligned stage mask.
- **B. Prenatal versus postnatal effect map**: effect-size comparison with
  uncertainty and the interaction estimand made explicit.
- **C. Expression-gating test**: relationship between stage expression change
  and interaction effect change, including the declared covariates and test.
- **D. TESPA1/SLCO1A2 gating example**: low-expression stage versus activated
  stage, with sample-level distributions.
- **E. ATP2B1 saturation example**: stage-specific dynamic range and the
  genotype-phenotype relationship; label the plateau as saturation evidence,
  not a solved disease mechanism.

**Status:** `Confirmed` from the accepted label-aligned rerun: 198 significant
groups, 193 HC3 core groups and 5 exclusions. The corrected SuSiE/mechanism
outputs and refreshed Figure 4 are promoted; the superseded sorted-half outputs
are no longer used by canonical consumers.

### Figure 5 | APA-QTL annotations reveal disease-relevant genetic signal beyond the tested eQTL context

**Message:** APA-QTL annotations carry disease-relevant heritability and identify
candidate loci that are not consistently supported by the tested brain eQTL
resources.

- **A. Binary-window sLDSC**: confirmed APA-QTL enrichment across traits and
  stages. The current summary contains 13/76 within-population-trait
  BH-significant rows, including positive EUR BIP and EUR SCZ signals; the
  unstable EAS ALS negative estimate is shown with its uncertainty and not
  framed as biological depletion. The corrected interaction subset has 0/19
  hits (minimum P=0.06763; within-GWAS FDR=0.09017).
- **B. Joint APA/eQTL maxCPP**: show the two annotation coefficient families
  separately, including the 32/76 well-powered rows and 8/76 coefficients
  passing within-regression two-coefficient BH. This panel is descriptive and
  is not a coefficient-contrast claim.
- **C. Disease-linked gene matrix**: coloc/TWAS/locus summary with explicit
  resource coverage, tested-stage scope and candidate status. The historical
  72-gene count remains provisional until producer provenance is closed.
- **D. SCAF1 locus**: prenatal apaQTL-SCZ colocalization and the tested eQTL
  resources, described as “no detected eQTL support in the tested resources”.
- **E. NRBP2 locus**: same evidence layout and wording as SCAF1; do not label
  either gene causal solely from PP4.

**Status:** A and B are `Confirmed` under their separate estimands. C-E are
`Provisional` pending coloc/input-manifest closure. Corrected MAGMA and the
interaction maxCPP collector are confirmed from accepted inputs; protein ORA is
numerically refreshed, but its final figure/text freeze remains pending. The
title deliberately says “beyond the tested eQTL context”, not “uncaptured
heritability” or “APA is superior to eQTL”.

### Figure 6 | Developmental stage and local regulatory context resolve selected disease-associated APA loci

**Message:** selected loci illustrate how stage, local gene context and multiple
regulatory layers can explain disease-associated APA signals.

- **A. ATP2B1**: stage-specific APA signal, saturation pattern and disease
  colocalization, with association-level wording.
- **B. 12q24.31**: neighbouring-gene signal separation across traits/stages,
  showing LD and annotation boundaries rather than asserting causal gene
  assignment.
- **C. ARL3**: APA/eQTL cross-layer evidence and the tested disease context,
  described as coordinated evidence rather than mediation proof.

**Status:** `Provisional`. These are mechanism-oriented case studies and should
remain in the main figure only after the coloc/SMR/functional producer logs and
input manifests are closed. Otherwise retain the same panels in Extended Data.

## Extended Data figures

| Figure | Content | Role |
|---|---|---|
| ED Fig. 1 | Sample QC, sex, PMI, RIN, donor and stage composition | Establish cohort quality and exclusions. |
| ED Fig. 2 | Batch structure and LifeStage-protected correction | Document preprocessing choices and sensitivity. |
| ED Fig. 3 | Trajectory robustness, leave-out/bridge checks and BrainSpan replication | Support Figure 1 without crowding the main narrative. |
| ED Fig. 4 | Full functional enrichment and background/sensitivity analyses | Preserve all tested pathways and null comparisons. |
| ED Fig. 5 | Protein QC, missingness, SVA and matched-input checks | Close the Figure 2 reproducibility boundary. |
| ED Fig. 6 | Cell-type/tissue context and BrainGVEX replication | Contextual support, not replacement for bulk evidence. |
| ED Fig. 7 | Cross-stage effect sensitivity, mashr and downsampling | Show alternate denominators and robustness analyses. |
| ED Fig. 8 | Full disease-trait panels, MAGMA/TWAS details and locus supplements | Keep exploratory disease results visible without overloading Figure 5. |

## Figure-freeze rules

1. Freeze final numbers only from the corrected producer result table, not from
   old manuscript prose or retained superseded outputs.
2. Keep `PP4 >= 0.7`, `PP_H4 >= 0.8`, strong-set sharing and detected-pair counts
   as separate estimands and denominators.
3. Use “associated with”, “supports”, “consistent with” and “candidate” for
   observational computational evidence; reserve “causal”, “demonstrates” and
   “solves” for evidence that actually identifies those claims.
4. Every final panel must show the tested sample/pair/gene denominator in the
   caption or an accompanying source table.
5. Main-figure layout follows the current order even while provisional panels
   are being rerun; only the panel data and captions change at freeze time.
