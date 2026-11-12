# Publication Workflow Audit

## Current repository verdict

Audit date: 2026-10-15

Evidence status: **Partial**. The repository contains code for most major
analysis branches described in `framework.md`, but it does not yet provide a
closed `documented input -> executable driver -> validated output` chain for
all eight modules. It is **not ready for public upload in its current state**.

| Module | Framework coverage | Current release status |
|---|---|---|
| 00 genotype | Core GATK/QC, ChinaMAP handoff, stage filtering, stage PGEN and PCA code are present | Partial: the post-imputation handoff is documented; ChinaMAP remains an external checkpoint and public configuration is incomplete |
| 01 RNA/APA | Major quantification and correction producers are present | Accepted/promoted: corrected ComBat_REF output is 63,097 loci x 815 samples with zero nonfinite values |
| 02 trajectories | Canonical donor-level driver and figure workflow are present | Accepted: 63,097 sites, 750 donors, 39,580 dynamic/clustered sites |
| 03 multiomics | Five-step protein validation driver is present | Accepted current rerun/tail: job 84024 completed with current matched inputs and refreshed validation outputs |
| 04 stage QTL | Mapping, coloc, mediation, mashr and downsampling code are present | Scientific runs complete; Slurm accounting closes coloc/mediation/mashr, but exact current-source revision pinning remains open |
| 05 GxStage | Interaction, eigenMT, HC3 and follow-up code are present | Accepted: label-aligned eigenMT, HC3, cohort-LD SuSiE/mechanism outputs and Figure 4 are promoted |
| 06 disease association | Two sLDSC estimands are connected | Numerical sLDSC/MAGMA and interaction maxCPP outputs are current; public-driver and environment gaps remain |
| 07 locus follow-up | Coloc, SMR and functional producers are represented | Scientific runs complete; public reconstruction remains partial because inputs and executable roots are external |

## Module 00: genotype QC and imputation

Current evidence status: **core methods represented, end-to-end orchestration
partial**. Company/manual WGS processing and web-based ChinaMAP imputation are
explicit checkpoints rather than locally reproducible stages.

### M00-01. The publication driver stops before the stage-specific QTL genotype handoff

- Classification: **Presentation only — Major workflow-coverage gap**
- Evidence: `1010/preprocess/run_module00.sh:8-16` submits joint genotyping, chromosome merge, VQSR, pre-imputation QC and upload preparation, then stops at the ChinaMAP web checkpoint. It only prints an instruction to run `00_4_stage_variant_filtering.sh`; it does not connect the post-imputation VCF to the stage-specific PGEN/PCA producer.
- Downstream producer: `1010/QTL/preprocess.R:262-309` extracts prenatal/postnatal PGEN files, performs `--indep-pairwise 50 5 0.2`, and computes 20 genotype PCs. These steps correspond to framework items 00_5 and 00_6, but `framework.md:39-40` and the Module 00 README do not identify this producer.
- Consequence: the scientific steps exist, but a public user cannot follow one documented chain from the downloaded ChinaMAP VCF to the stage QTL inputs.
- Resolution: `1011/00_genotype_qc_imputation/README.md` now records the exact
  `downloaded VCF -> 00_4_stage_variant_filtering.sh -> genotype_clean PGEN ->
  QTL/preprocess.R -> stage PGEN/20 PCs` contract, including MAF 0.05, MAC 20,
  HWE 1e-6 and LD-pruning `50 5 0.2`. ChinaMAP remains an explicit web checkpoint.

### M00-02. Module 00 configuration is not sufficient for a fresh public checkout

- Classification: **User/workspace change — Major reproducibility gap**
- Evidence: the submitted scripts default `SHARE_ROOT` to `/path/to/share_group_folder` (`01_joint_genotyping.sh:17`, `02_merge_clean.sh:13`, `03_VQSR.sh:18`, `04_preQC.sh:14`, and `05_prepare_imputation.sh:12`). `run_module00.sh` does not set or validate this root, and `1011/00_genotype_qc_imputation/README.md` does not list the required environment variables or reference resources.
- Test boundary: `1011/tests/smoke_module00.sh` tests chromosome ordering and stage-variant intersection with command stubs; it does not exercise GATK, VQSR, sample QC, ChinaMAP preparation, or stage PGEN/PCA generation.
- Consequence: the smoke test establishes two handoff contracts, but it is not evidence that a fresh checkout can execute the complete Module 00 workflow.
- Required action: provide one example config or environment file listing the required reference/resource inputs and document the post-imputation continuation command.

## Module 01: RNA-seq and APA quantification

Audit date: 2026-10-11

### Scope and evidence status

This audit evaluates the scientific workflow rather than server-specific paths. It focuses on data identity, transformations, statistical choices, correction scope, denominators, and whether each downstream matrix can be traced to a producer. Absolute paths are relevant only when they change the input data or break that trace.

Current evidence status: **accepted/promoted corrected rerun**. The selected ComBat_REF matrix and RNA/isoform outputs passed the Module 01 schema, sample, key, finite-value, and method gates.

### Findings

#### M01-01. RNA ComBat does not protect LifeStage

- Classification: **Wrong result — Major**
- Original evidence: `1010/preprocess/rna_pre_step01.r:198` calls `ComBat(..., mod = NULL)`.
- Scientific issue: batch and LifeStage are strongly aliased in the inspected metadata: batches 1–4 are postnatal, batches 5–8 are prenatal, and batch9 contains both stages. Without a biological design matrix, ComBat can remove genuine prenatal/postnatal signal together with batch effects.
- Cross-check: the active producer `1010/preprocess/rna_pre_step01_v2.r:240-245` uses `mod = ~ LifeStage` and `ref.batch = batch9` when that batch is present.
- Consequence: the legacy RNA preprocessing script is scientifically inconsistent with the implementation supporting the retained analysis and should not be presented as the publication producer.
- Required action: identify the authoritative RNA correction producer and retain LifeStage in the ComBat model. State the protected covariate and reference batch in Methods.
- Resolution implemented in code: `1010/preprocess/rna_pre_step01_v2.r` processes the combined cohort, performs stage-stratified outlier detection, protects LifeStage in the ComBat design, and uses batch9 as reference when present. A metadata file without `LifeStage` is supported by deriving the two stages from the sign of `age`, matching the supplied cohort convention. The Module 01 `rna_pre_step01.r` entrance now resolves to this v2 producer.

#### M01-02. RSEM matrix merge has no retained identity contract

- Classification: **Wrong result — Major for the retained isoform branch; producer and reader fixes implemented but not rerun**
- Original evidence: `1010/preprocess/rna_quant.sh:68-85` merges only RSEM `.genes.results` files.
- Consumer: the active Module 01 preprocessing producer reads the retained gene and isoform matrices before correction; Module 02 consumes the resulting `data/process/expr_after_combat.tsv` and `data/process/isoform_after_combat.tsv` files.
- Data identity check: the supplied gene matrices have 912 columns but duplicate `XYA20171125`; metadata contains 911 unique sample IDs. The retained isoform matrices contain 912 unique sample IDs but differ from the gene/metadata naming at six columns: `20171125C`, five A/B-suffixed XYA IDs, and their unsuffixed counterparts.
- Direct verification: summing transcript expected counts by GENCODE v40 gene produces Pearson correlation 1.0 against each corresponding gene-count column. The two duplicated gene columns map separately to `XYA20171125` and `20171125C`; the other five suffixed isoform IDs map to the corresponding unsuffixed metadata IDs, with mean absolute discrepancies of 0.00114 or less from RSEM rounding.
- Consequence: the retained isoform consumer silently excluded five metadata-matched samples because it intersected unharmonized column names; `20171125C` is separately excluded because no metadata row was supplied. Among metadata profiles with RIN >= 5, the raw name intersection increases from 822 to 827 after the verified aliases are applied. The retained 805-profile `isoform_after_combat.tsv` predates this fix and contains none of the five recovered IDs.
- Resolution implemented in code: `1010/preprocess/rna_quant_v2.sh` follows the supplied BrainGvex merge contract, exports paired gene/isoform count and TPM matrices from the same complete RSEM result set, canonicalizes the five verified aliases, preserves `20171125C` as a distinct ID, writes a source-to-matrix sample manifest, and refuses duplicated IDs or positional merges with inconsistent feature order. `1010/preprocess/rna_pre_step01_v2.r` also disambiguates the verified second `XYA20171125` gene-matrix column as `20171125C`; the retained large matrices and downstream isoform output have not been rewritten.

The gene and isoform sex-QC implementations are now aligned: both use `set.seed(2026)`, two-center k-means on `log1p` XIST/UTY expression, thresholds `expm1(mean(cluster centers))`, and the same known/unknown/mismatch/ambiguous decision logic with UTY-positive precedence. The isoform branch sums GENCODE v40 XIST/UTY transcripts before applying that rule. The existing isoform output was generated before this unification and is not evidence that the two branches retain identical samples.

#### M01-03. APA correction-method selection uses the downstream age outcome

- Classification: **Wrong result — Major**
- Original evidence: `1010/batch/00_batch_correction.R:1883-1897` ranks correction methods using downstream age-GAM significance, partial deviance, age-effect magnitude, and k-index failures; age-signal recovery receives weight 2.
- Scientific issue: the same cohort and corrected matrix are then used for lifespan trajectory discovery. Selecting preprocessing partly by recovery of the target age association is outcome-aware and can make downstream age effects optimistic.
- Consequence: the selected `ComBat_REF` matrix may be technically well behaved, but the strength and number of age associations are not independent of the preprocessing-selection criterion.
- Required action: disclose the selection rule explicitly in Methods and select the correction without using the downstream age outcome.
- Resolution implemented in code: `1010/batch/00_batch_correction_v2.R` retains age diagnostics for reporting but excludes them from selection. Corrected candidates are ranked across four equally weighted domains: bridge concordance (Pearson and CCC), bridge error (MAE and RMSE), batch removal (residual batch coefficients, batch variance, and batch ASW), and biological preservation (LifeStage variance and ASW). Correlated metrics are first averaged within their domain or subcomponent, and `Raw` is retained only as an uncorrected baseline. The Module 01 `batch_correction.R` entrance now resolves to this v2 producer.
- Direct comparison: `1010/batch/01_bridge_method_selection.R` evaluates 65 donor pairs and 10,000 paired-bootstrap resamples over the common finite locus set. ComBat_REF exceeds limma for Pearson by 0.01432 (95% CI 0.01209–0.01678), CCC by 0.01458 (0.01208–0.01725), and MAE by 0.07552 (0.05312–0.09889); RMSE favors limma by 0.03773 (0.01782–0.05687). The mixed error result is retained rather than presented as uniform superiority. Across all corrected candidates, the four-domain composite ranks ComBat_REF first (2.0625), followed by ComBat_REF_NP (2.7083) and limma (3.4167).

#### M01-04. The selected APA matrix is not connected to the downstream TSV by a recorded export step

- Classification: **User/workspace change — Major reproducibility gap**
- Original evidence: `1010/batch/00_batch_correction.R` saves candidate method matrices and selects `ComBat_REF`, but does not export `00_raw_data/batch_correction/rds/method_matrix_ComBat_REF.rds` to `02_pattern/data/raw/apa_after_combat.tsv`.
- Observed contracts: the selected RDS manifest reports 63,097 loci by 815 profiles. The downstream TSV contains 815 profiles, comprising 750 unique donors and 65 replicate profiles.
- Scientific issue: matching dimensions do not establish value-level identity, ordering, or whether any additional filtering/transformation occurred during handoff.
- Consequence: downstream analyses cannot yet be traced deterministically to the selected batch-corrected object.
- Required action: recover the original export command or add one explicit handoff step, then verify feature IDs, sample IDs, ordering, dimensions, and numeric equality or document any intended transformation.
- Resolution implemented in code: `1010/batch/00_batch_correction_v2.R` exports the selected in-memory matrix with an explicit `site_id` column to `02_pattern/data/raw/apa_after_combat.tsv`. The v2 producer has not yet been run, so the existing TSV has not been replaced.

#### M01-05. DaPars2 filtering is incompletely described

- Classification: **Presentation only — Minor Methods omission**
- Evidence: `Dapars2_Multi_Sample.py:223-241` averages Python slice `[0:99]` (99 positions), applies strict coverage `>15`, requires more than 30% of profiles to pass, and requires a 3' UTR span of at least 150 bp.
- Scientific issue: the previous framework reported only `>=15`, omitting the exact window, prevalence rule, and minimum 3' UTR span. The implementation uses strict `>15`.
- Consequence: the implemented site-retention rule cannot be reconstructed accurately from the current Methods description.
- Required action: report all three filters and the exact inequality used by the implementation.

#### M01-06. DaPars2 source provenance was not linked

- Classification: **Presentation only — resolved provenance**
- Evidence: the authoritative checkout at `${APA_HPC_HOME}/project/work/ref_repo/DaPars2` is clean at commit `fb81c6c`, reports `v2.1-14-gfb81c6c`, and its README labels the software DaPars v2.1. The collected Module 01 Python source differs only by a multiprocessing spawn/main guard.
- Resolution: framework provenance now records the exact commit/describe value and the local compatibility patch rather than relying on an unsupported version claim.

### Confirmed workflow elements

- `all.apa.count` contains 827 profiles, including 66 replicate-labelled profiles.
- The selected `ComBat_REF` output has 63,097 loci and 815 profiles.
- The downstream APA matrix has 815 profiles: 750 unique donors and 65 replicate-labelled profiles.
- The current ranking table selects `ComBat_REF`.
- The retained status report records zero non-finite values in the selected matrix.
- The APA correction workflow implements LifeStage-protected ComBat, batch9 as reference, bridge-pair checks, and comparison of eight correction methods.
- Reviewed shell scripts pass `bash -n` syntax checking.

### Publication gate for Module 01

Module 01 now defines ComBat_REF as the primary correction under the no-age, four-domain selection framework. The retained bridge analysis provides a direct comparison with limma, including the opposing RMSE result; a full duplicate limma downstream pipeline is therefore not part of the selection requirement. Publication release still requires execution of the corrected RNA/isoform readers and an end-to-end v2 batch-correction run. Regeneration of the four aligned matrices from the original RSEM files remains the preferred provenance path when those files become available; the current readers provide an explicit, verified compatibility path for the retained matrices.

### Not verified in this audit

- Value-level identity between the selected `ComBat_REF` RDS and downstream `apa_after_combat.tsv`.
- Metadata provenance for the distinct `20171125C` RSEM sample; it remains excluded because no metadata row was supplied.
- A retained execution log tying each published RNA/APA matrix to a specific command and software environment.
- Re-execution of RNA quantification, RNA preprocessing, DaPars2, batch correction v2, or downstream trajectory analysis.
- Re-execution of the corrected isoform reader; the retained `isoform_after_combat.tsv` still reflects the pre-fix sample intersection.

## Module 02: lifespan trajectories

Current evidence status: **accepted donor-level rerun**. The refreshed trajectory
branch uses 750 donors and 63,097 accepted sites; downstream counts and figures
must still be synchronized before release freeze.

### M02-01. Technical replicates were treated as separate GAM observations

- Classification: **Wrong result — Major; producer fix implemented but not rerun**
- Evidence: the retained covariate table contains 815 profiles from 750 donors, including 65 technical-replicate pairs. The superseded profile-level producer assigned each pair weights of 0.5, while `02_gam.R` still supplied all 815 rows to `mgcv::gam`.
- Scientific issue: prior weights equalize each donor's contribution but do not model within-donor dependence; mgcv residual degrees of freedom remain based on the number of observations. The age-effect p-values and BH-FDR therefore did not use 750 independent analysis units.
- Consequence: the retained 52,265 age-associated and 39,555 directional events require recomputation before publication use.
- Resolution implemented in code: `1011/02_lifespan_trajectories/pipeline/01_cov.R` averages ComBat-corrected logit-PDUI values within donor, averages RIN within donor, and writes a 750-donor matrix plus one covariate row per donor. GAM, prediction, directional-event, functional, and discovery cross-fit consumers now use this donor-level matrix. Selected RNA and isoform rows are averaged with the same donor mapping before cross-layer calculations. The profile-level sensitivity branch is not part of the publication workflow.

## Module 03: multiomics uncoupling and protein validation

Current evidence status: **accepted current rerun/tail**. The corrected Module 01/02
inputs were consumed by the v4 protein chain; job 84024 completed the tail and
refreshed bootstrap, model, triad, and figure outputs. Historical retained-output
path gaps remain provenance notes for the superseded snapshot.

### M03-01. Protein validation entrypoint is disconnected from the retained matched inputs

- Classification: **User/workspace change — Major reproducibility gap**
- Evidence: `1010/protein/scripts/run_protein_validation.sh:7-15` calls the five producer scripts. The producer configuration uses `${APA_HPC_HOME}/project/work/apa/fix/protein` and `${APA_HPC_HOME}/project/work/apa/02_pattern`, for example `01_build_analysis_universe.R:9-21` and `03_integrate_cross_layer_candidates_v4.R:14-30`. Those the former cluster-scratch `fix/protein` directories are absent at audit time.
- Input contract: `${APA_HPC_HOME}/project/work/apa/12_protein/raw/input_manifest.tsv:2-6` records the five matched artifacts with 366 samples, 63,097 APA features, 28,035 RNA features, and 8,522 protein features, but every matched artifact path in those rows points to the absent the former cluster-scratch `fix/protein` tree.
- Retained result scope: `${APA_HPC_HOME}/project/work/apa/12_protein/results/03_site_protein_pairs.tsv` has 20,447 data rows; `${APA_HPC_HOME}/project/work/apa/12_protein/results/03_vif_diagnostics.tsv:2-3` records 20,447 primary pairs and 13 beta-gate pairs. These outputs are not evidence that the current publication wrapper can reconstruct the same inputs.
- Consequence: the five-step protein producer cannot be replayed from the current checked-out paths, and the value-level relationship between the retained outputs and the `${APA_HPC_HOME}/project/work/apa/12_protein` source tree is not established by the wrapper.
- Required action: recover the exact matched artifacts or repoint the producer only after verifying sample order, feature order, and the manifest MD5 values; then record a successful producer log and output manifest before using the protein claims as reproducible publication outputs.

### Not verified in Module 03

- A retained successful execution log for the exact five-script wrapper and its current source revision is **未记录**.
- Whether the retained 20,447-pair and 13-pair tables were generated from the current `${APA_HPC_HOME}/project/work/apa/12_protein` files or from the missing the former cluster-scratch `fix/protein` tree is **未记录**.

## Module 04: stage-specific apaQTL, colocalization, mediation, and mashr

Current evidence status: **results and producer definitions are available; reporting estimands must remain separated**. The coloc summary explicitly uses a different reporting threshold from the classification threshold, while the mediation and mashr prose still needs provenance repair.

### M04-01. Coloc classification and PP4 reporting use different thresholds and denominators

- Classification: **Presentation only — estimand wording gate**
- Evidence: `${APA_HPC_HOME}/project/work/apa/11_eqtl/02-2_coloc.R:19-21,229-234` classifies loci at `PP_H4 >= 0.8`. `${APA_HPC_HOME}/project/work/apa/11_eqtl/02-3_summary.R:19-20,167-180` reports PP4 at `0.7` over significant APA loci whose gene was eQTL-tested, retaining loci with no coloc result in the denominator.
- Authoritative summary: `${APA_HPC_HOME}/project/work/apa/11_eqtl/results/coloc_pp4_definition_summary.tsv:2-3` gives prenatal `165/666 = 24.7748%` and postnatal `576/2231 = 25.8180%`; each stage has four denominator entries without a coloc result.
- Required action: keep `PP_H4 >= 0.8` for classification and `PP4 >= 0.7` for this reporting summary, and state the numerator/denominator and four missing-result entries in any manuscript or figure text. These estimands must not be merged into one percentage.

### M04-02. Mediation p-values and bootstrap intervals have different statistical provenance

- Classification: **Presentation only — Major methods-provenance gap**
- Evidence: `${APA_HPC_HOME}/project/work/apa/11_eqtl/04-1_mediate.R:21-27` selects `PP_H4 >= 0.8` pairs and sets `N_BOOT = 2000`. The indirect-effect p-values at `:156-160` use a delta-method standard-error formula followed by a normal approximation; the bootstrap at `:162-183` supplies percentile confidence intervals.
- Result scope: `${APA_HPC_HOME}/project/work/apa/11_eqtl/results/04_mediation_summary.tsv:2-7` contains 137 prenatal pairs and 521 postnatal pairs, 658 pairs total.
- Consequence: describing the result only as “2,000-bootstrap mediation” hides that the reported indirect-effect p-values are not bootstrap p-values.
- Required action: Methods and captions must state the PP4 gate, 2,000 donor resamples, delta-method/normal-approximation p-values, and percentile bootstrap CIs. Do not label the current p-values as bootstrap-derived unless the producer is changed and rerun.

### M04-03. mashr sharing values are reported with the detected-pair denominator

- Classification: **Wrong result — reporting estimand/denominator**
- Evidence: `${APA_HPC_HOME}/project/work/apa/11_eqtl/05_mashr.R:24-27,165-175` defines `share_factor = 0.5` and fits mashr on the `strong` set. The posterior outputs contain 12,918 APA rows and 55,459 eQTL rows after excluding headers. The sharing metric at `:333-337` is returned by `get_pairwise_sharing()` on the fitted strong-set posterior, whereas detected counts at `:358-371` are calculated only after excluding `Not significant` classes.
- Current results: `${APA_HPC_HOME}/project/work/apa/11_eqtl/results/05_mash_summary.tsv:10-11` reports pairwise sharing `0.8485827054` for apaQTL and `0.6566184352` for eQTL. `${APA_HPC_HOME}/project/work/apa/11_eqtl/results/05_effect_consistency_direction_summary.tsv:2-5` separately reports 2,787 and 14,928 detected pairs.
- Reporting conflict: `${APA_HPC_HOME}/project/work/apa/11_eqtl/framework.md:112-119` describes the sharing metric inside the “at least one stage detected” subset and retains stale counts 2,792/15,017 with rounded sharing 0.850/0.654. `${APA_HPC_HOME}/project/work/apa/draft/parts/results.md:35` also presents 2,787 detected pairs together with the sharing value without stating the strong-set scope.
- Consequence: the current prose makes the sharing value appear to be estimated over 2,787/14,928 detected pairs, while the producer calculates it over the larger strong posterior sets. The reported values are not comparable until the scope is labeled correctly.
- Required action: update the framework and Results prose to distinguish strong-set denominators (12,918 and 55,459) from detected-class counts (2,787 and 14,928), and replace the stale counts/rounded values with the current authoritative result values.

### Not verified in Module 04

- A single retained execution log linking all current coloc, mediation, and mashr tables to exact source revisions and environment versions is **未记录**.

## Module 05: genotype-by-stage interaction

Current evidence status: **accepted label-aligned rerun**. Corrected eigenMT,
HC3, cohort-LD SuSiE/mechanism outputs and Figure 4 are promoted to their
canonical consumer paths.

### M05-01. eigenMT uses a sorted-half mask instead of the LifeStage labels

- Classification: **Wrong result — Major**
- Evidence: `${APA_HPC_HOME}/project/work/apa/03_interaction/LifeStage.txt` contains 189 stage-0 samples and 508 stage-1 samples. In `${APA_HPC_HOME}/project/work/apa/03_interaction/scripts/02_group_eigenmt.py:105-112`, the code reads the interaction table, sorts the first interaction column, and sets the first `interaction_df.shape[0] // 2` rows to `False`; it never constructs the mask from the `LifeStage` identities.
- Consequence: for the 697-sample interaction table, the current code creates an approximately 348/349 split rather than the observed 189/508 stage split. That mask is passed into `filter_maf_interaction()` at `:118-126` and changes the variants retained for the eigenMT count at `:129-132`.
- Resolution: `LifeStage.txt` is joined to phenotype sample IDs with the 189/508
  label gate. The rerun retained 198 significant groups, 193 HC3 core groups
  and 5 exclusions; all 12 eigenMT diagnostic checks passed. Cohort-LD SuSiE
  produced 386 stage-group summaries, including 99 usable fits, and the accepted
  credible-set/mechanism table contains 2,397 rows across 79 groups.

### Remaining Module 05 boundary

- Interaction maxCPP is a downstream disease branch rather than an acceptance
  condition for the corrected interaction-group result. Chromosome array 83904
  is active and dependent regression job 83905 is queued.

## Module 06: disease heritability and gene-set analysis

Current evidence status: **binary-window sLDSC, main joint maxCPP and corrected
MAGMA numerical outputs are present, but the framework module is not one closed
public workflow and its completion gate remains too weak**.

### M06-01. The Module 06 completion gate accepts partial result sets

- Classification: **Wrong result — Major completion-gate defect**
- Evidence: `1011/06_disease_heritability_geneset/run_module06.sh:46-50` requires only one APA-only `.results` file and one joint maxCPP `.results` file. The APA-only producer warns and skips missing GWAS inputs at `1010/sLDSC/run_sldsc_pipeline.sbatch:210-246`, so a run with most traits absent can still reach the manifest-writing step.
- Test contrast: `1011/tests/smoke_module06.sh:172-176` checks exact mock counts (8 APA-only and 4 joint results), but the production wrapper does not derive or enforce the expected counts from its annotation, stage and trait sets.
- Consequence: an incomplete trait-by-annotation matrix can be labelled as a completed Module 06 run and then summarized as if the requested scope were closed.
- Required action: calculate the expected unique result keys from configured traits, annotations and stages; fail unless every expected non-skipped result and summary row is present, or record an explicit skipped-trait manifest.

### M06-02. The public Module 06 driver covers sLDSC only

- Classification: **Presentation only — Major framework-to-code mismatch**
- Evidence: `framework.md:104-111` defines GWAS harmonization, two sLDSC estimands, MAGMA and TWAS. `1011/06_disease_heritability_geneset/README.md:19-20` explicitly states that MAGMA and TWAS are not invoked by `run_module06.sh`.
- Upstream gap: the APA-only sLDSC producer begins at Phase 3 and requires pre-existing annotation BEDs (`run_sldsc_pipeline.sbatch:93-101`); the smoke test supplies these directly at `smoke_module06.sh:113-116`. The joint producer and TWAS scripts consume files named `*_harmonized.tsv`, but no GWAS liftOver/allele-harmonization producer is exposed by the Module 06 entrance.
- Consequence: `run_module06.sh` is accurately an sLDSC driver, not the complete framework Module 06 workflow.
- Required action: either rename/document it as the sLDSC-only driver and add explicit commands for 06_1, 06_3 and 06_4, or create one thin Module 06 driver that calls the existing producers in dependency order.

### M06-03. The framework claims a MAGMA conditional analysis that the producer does not run

- Classification: **Presentation only — Major Methods mismatch**
- Evidence: `framework.md:110` states that the competitive gene-set test uses conditional analysis to exclude background-gene effects. `1010/MAGMA/run_magma_enrichment.sh:97-102` calls MAGMA with `--gene-results`, `--set-annot`, and `--out`; it does not specify a `condition=` model or another conditional covariate/set.
- Consequence: the retained producer supports a standard competitive gene-set analysis, but it does not support the additional conditional-analysis statement in the framework.
- Required action: remove the conditional-analysis wording if it was not part of the reported result, or implement the intended condition and rerun before retaining that claim.

### M06-04. The copied TWAS submission wrapper points outside the copied TWAS directory

- Classification: **User/workspace change — Major reproducibility gap**
- Evidence: `1010/TWAS/submit_full_twas_training.sh:23-31,47,61` resolves `ROOT` to `1010/` and calls `1010/scripts/02_train_full_twas_models.py`; that path does not exist. The copied worker is `1010/TWAS/02_train_full_twas_models.py`. The same flattened-layout mismatch affects `01_prepare_full_twas_inputs.py:31-33`, which resolves data/output roots one directory above `TWAS/`.
- Scope gap: the submission wrapper defaults to prenatal training only and does not call input preparation, postnatal training, GWAS alignment/SPrediXcan, result summarization, or final plotting.
- Consequence: the current public TWAS copy cannot execute framework item 06_4 from its documented location.
- Required action: restore the original TWAS directory layout or update the copied paths consistently, then add one documented command sequence covering preparation, both stage models, SPrediXcan and summary generation.

### Remaining Module 06 boundary

- Corrected interaction binary-window sLDSC job 46274987 completed 44/44
  annotations, 44/44 LD-score files and 19/19 regressions. The combined 76-row
  resummary has 13 within-family hits and 0/19 interaction hits. Corrected MAGMA
  contains 57 rows and five within-GWAS hits, with 0/19 interaction hits.
- Interaction maxCPP chromosome/regression outputs and the 19-row tau-star
  collector summary are complete. TWAS did not require a model rerun: the
  19-trait overview jobs completed successfully, and the corrected master
  summary contains all 38 GWAS-stage rows. Exact training-environment capture
  remains incomplete for public release.

## Module 07: locus colocalization, SMR, and functional follow-up

Current evidence status: **authoritative coloc, SMR and functional runs have
successful scheduler/log closure, but a fresh public clone cannot reconstruct
their external inputs and executable roots**.

### M07-01. Coloc depends on an external authoritative configuration

- Classification: **User/workspace change — Major provenance gap**
- Current evidence: `1010/coloc/run_coloc.sh:16-24` defaults to the authoritative remote `04_coloc` config/result root and `:17,35-43` executes the separate EasyColoc checkout. Neither the real GWAS/QTL config nor its referenced inputs are included in this publication repository.
- Consequence: the wrapper selects the intended remote analysis on the original system, but a fresh public clone cannot identify or reconstruct the exact configured datasets.
- Required action: publish de-identified config templates plus an input manifest with stable accessions/checksums, and pin the EasyColoc dependency before citing the wrapper as reproducible.

### M07-02. SMR has successful output-closing logs but no public input manifest

- Classification: **User/workspace change — Major provenance gap**
- Current evidence: `1010/SMR/run_smr.sh:23-65` now defaults to the authoritative remote SMR, harmony, QTL and 1KG roots. Those inputs are outside this repository. The authoritative retained output tree contains 57 `.msmr` files, 19 effective-test records, 57 gene-summary records, and 8 gene-level Bonferroni rows after excluding headers.
- Log evidence: Furong job `44649810` completed `0:0` after 2:43:49, and August summary jobs `44907045` and `44907046` also completed `0:0`. The latter log records 1,367,247 SMR tests and the retained aggregate outputs. The older `smr_42454069.err` failure is superseded, not the terminal workflow state.
- Consequence: the retained SMR result chain is complete on the authoritative system, but its external input/reference roots are not reconstructable from this repository alone.
- Required action: publish the input/reference manifest and pin the authoritative wrapper/environment; do not rerun the statistically complete SMR models solely to create provenance.

### M07-03. Functional follow-up delegates to an external project absent from the public checkout

- Classification: **User/workspace change — Major provenance gap**
- Current evidence: `1010/func/run.sh:6-9` delegates directly to `${APA_HPC_HOME}/project/work/apa/05_func/run_analysis.sh` unless `FUNCTIONAL_PROJECT_DIR` is overridden. The copied `1010/func/` scripts do not include that authoritative `run_analysis.sh`, its config, raw inputs or reference tree.
- Consequence: the wrapper preserves the original authoritative execution boundary but is not executable from a fresh public clone.
- Required action: include or externally pin the authoritative functional producer and publish its input/config manifest. Until then, functional follow-up claims remain provenance-blocked.

### Not verified in Module 07

- A successful invocation log for the current publication-facing coloc wrapper is **未记录**, but the authoritative EasyColoc run is complete: `04_coloc/results/runtime/monitor_snapshot.json` records `pipeline_complete`, 19 GWAS and 3,084 completed tasks; `output_manifest_summary.json` records 2,070 files; the run log at lines 22845-22847 records 3,084 merged rows and 151 PP4 >= 0.7 rows.
- The current PP4 >= 0.7 CSV reproduces 104 unique GWAS/QTL/phenotype groups and 72 unique gene tokens (151 locus rows). The separate 40% APA-only component remains **未记录** because no retained producer/result mapping for the GTEx/CHB classification was found.
- EAS_ALS was processed and produced no significant loci at the configured threshold (run log line 9644); its absence is not a failed pipeline run.
- SMR correction is within `qtl_type x gwas_id` (`1010/SMR/run_smr.sh:347-350`); jobs `44649810`, `44907045` and `44907046` provide successful producer/aggregation closure for the retained `.msmr` and summary outputs.
- Functional job `44899287` completed `0:0`, validated the canonical package and wrote the manifest; later mechanism jobs `45874499`, `45874500` and HID1 job `46256367` also completed `0:0`. The current authoritative manifest has 125/125 passing rows. Public reproduction still requires the external producer/input tree.

## Recheck update: cross-module contracts

- Module 00 retains 707 samples and 5,014,001 variants; stage sets contain 189 prenatal and 508 postnatal samples. The PGEN/PCA handoff is now documented; the external ChinaMAP execution record and public resource configuration remain open.
- Module 04 keeps distinct estimands: coloc reports PP4 >= 0.7 (165/666 prenatal, 576/2,231 postnatal; four missing results per stage), while classification and mediation use PP_H4 >= 0.8. Mediation covers 658 pairs with 2,000 bootstraps. mashr sharing uses strong-set denominators (12,918 APA; 55,459 eQTL), not detected subsets (2,787; 14,928). Slurm jobs `83118`, `83120`, `83121`, `83124`-`83127`, `83190` and `83191` completed `0:0`; exact current-source revision pinning remains open because several scripts postdate those runs.
- TWAS applies BH per GWAS-stage family. The canonical table has 17 GWAS/34 stage pairs and supplemental EAS_EP/EUR_EP rows bring the contract to 19 GWAS/38 pairs. `06_prepare_multi_gwas_visualization_v2.py` now merges both sources and hard-checks the complete 38-row master summary.
- Module 13 job 83697 closed successfully: 242 candidates were checked; exact overlaps were 0 Hippocampus, 15 Retina and 9 Eye. Retina/Eye matches are cross-tissue context, not same-tissue validation.

## Resolution update: remote execution boundary and minimum code repairs

The findings above describe the pre-patch state. The following resolution
updates include both static repairs and later accepted numerical reruns:

- Module 03 protein: `1010/protein/scripts/run_protein_validation.sh` now defaults to the recorded remote protein and pattern roots, exports `PROTEIN_PROJECT_DIR`/`PATTERN_PROJECT_DIR`, and keeps `SCRIPT_DIR` plus the existing `PROJECT_DIR` override. The five R producers now retain the same remote defaults. The remote matched-input manifest remains authoritative; no local artifact paths were substituted.
- Module 04 reporting: `11_eqtl/framework.md`, `draft/parts/results.md`, and `draft/parts/method.md` now distinguish the mash strong-set denominators (12,918 apaQTL; 55,459 eQTL) from detected subsets (2,787; 14,928), and state the PP_H4 >= 0.8 mediation gate, delta-method/normal-approximation indirect-effect p-values, and 2,000-resample percentile CIs.
- Module 05 interaction: both copies of `02_group_eigenmt.py` align `LifeStage.txt` to phenotype sample IDs, require binary labels with counts 189/508, and pass `stage == 1` to TensorQTL MAF filtering. The corrected eigenMT/HC3/SuSiE/mechanism outputs have now been validated and promoted.
- Module 07: `1010/coloc/run_coloc.sh` now defaults to `${APA_HPC_HOME}/project/work/apa/04_coloc/config` and its remote results root; `1010/SMR/run_smr.sh` now defaults to the recorded remote SMR, harmony, and 1KG roots while accepting an absolute `RESULTS_DIR`; `1010/func/run.sh` now dispatches to the authoritative remote `05_func/run_analysis.sh` with `FUNCTIONAL_PROJECT_DIR` override.
- Verification: wrapper shell syntax, both interaction Python syntax checks,
  `git diff --check`, `LifeStage.txt` counts and equality of the two interaction
  scripts passed. Later label-aligned Module 05, interaction sLDSC and MAGMA
  numerical reruns closed their stated output contracts.

## Public repository release gate

### R01. Human molecular and genotype inputs are excluded from Git

- Classification: **User/workspace change — Resolved local release-safety gate**
- Evidence: `.gitignore` excludes `raw/`, `1013/after.impute.all.raw`, and `1013/sample.Adult.final`. `raw_manifest.tsv` is a valid four-column TSV describing the eight current local inputs without embedding their contents.
- Consequence: the currently known individual-level molecular, metadata and genotype inputs are protected from a broad `git add .`.
- Release check: re-run `git check-ignore` on the protected paths and inspect all staged files before every public tag; publish only consent-compatible accessions, instructions, or approved derived inputs.

### R02. Corrected outputs still need final figure/manuscript synchronization

- Classification: **Presentation only — Major synchronization gap**
- Evidence: Module 01 and Module 02 acceptance gates are closed, and Module 03 job 84024 completed after consuming those accepted inputs. `results.md` records the accepted closure, while dependent figures and manuscript text remain provisional.
- Consequence: a code release made before figure/text synchronization could pair current producers with stale presentation artifacts.
- Required action: regenerate dependent panels and synchronize manuscript numbers and claims from the accepted current outputs before tagging a release.

### R03. Release metadata and environment capture remain incomplete

- Classification: **Presentation only — Major release blocker**
- Evidence: `README.md` now identifies the study scope, module map, data boundary and dependency limitations, and `RELEASE.md` defines the intended release set and gates. The root still has no `LICENSE`, `CITATION.cff`, package manifest or lockfile.
- Consequence: users can navigate the repository, but they still cannot determine reuse terms, cite a fixed release, or reconstruct the accepted software environments.
- Required action: obtain the owner-supplied license and citation metadata; capture the tested R/Python/CLI environments in a machine-readable lockfile or container; stage only the reviewed release set.

### R04. The coloc dependency is included in the publication tree

- Classification: **User/workspace change — Major reproducibility gap**
- Evidence: `1011/07_locus_colocalization_smr/EasyColoc/` contains the 98 tracked files from upstream commit `850178b5c99899d3e00b8655ffda716f705aa09c`, including `LICENSE` and `environment.yml`; the publication coloc wrapper resolves this local copy by default.
- Consequence: transferring `1011/` preserves the coloc executable and its declared environment without relying on a parent-repository symlink or nested `.git` directory.
- Release check: retain the upstream commit in `SOURCE_MANIFEST.tsv` and preserve the vendored MIT license.

### R05. Public configuration and test coverage are incomplete

- Classification: **Presentation only — Major release-readiness gap**
- Evidence: absolute local/HPC paths occur in scripts across every module; no root example configuration maps them to required inputs. Only Module 00 and the two sLDSC branches have smoke tests under `1011/tests/`. Personal Slurm mail addresses were removed from the two current sLDSC submission scripts.
- Consequence: paths are recoverable for the original environment but do not define a portable public input contract, and the critical Module 01-05/07 entrypoints have no bounded fresh-checkout checks.
- Required action: replace personal contact details, document each required external input/reference and override variable, and add small contract checks for the canonical drivers that currently carry correctness-critical handoffs.

## Audit completion update

Rechecked on 2026-10-12 before the pre-rerun baseline commit:

- The audit covers Modules 00-07, the public-release gates, remote execution
  boundaries, and the currently observed producer/result contracts.
- The `raw/` protection finding is resolved for this local repository state:
  `raw/` is now Git-ignored and `raw_manifest.tsv` records the local file names,
  sizes and modification times without staging human molecular or genotype data.
- The remaining findings are intentionally open blockers, not unfinished audit
  sections: corrected Module 01/02 outputs still need reruns; Module 03/07
  producer provenance remains remote; Module 06 orchestration and release
  metadata remain incomplete.
- Therefore `audit.md` is complete as an evidence-backed audit of the current
  stage, while the repository itself remains not ready for public release.

### 2026-10-15 closure update

- Module 01, Module 02, Module 03, and Module 05 numerical rerun gates are closed.
- Interaction maxCPP collection is closed by `04-4_collect_interaction_v2.R`; the accepted summary has 19 rows, 19 unique population-trait keys, and category `L2_1`.
- Release packaging now includes a root module/data/dependency map and an explicit release checklist; private inputs under `raw/` and `1013/` are protected from Git.
- Remaining blockers are owner-supplied license/citation/data-access metadata, public input/environment provenance, missing cell-QTL QC delivery, and final figure/manuscript synchronization.

### Release verdict

The repository **covers the manuscript's major analytical topics and now has a symlink-free reviewer-facing script tree**, but it does not yet provide the complete fresh-clone data/environment workflow promised by `framework.md`. The numerical rerun order, local private-data boundary, script materialization, and EasyColoc source delivery are closed. Public release remains blocked by license/citation/data-access metadata, Module 03/07 input and environment provenance, cell-QTL delivery, and synchronized figures/manuscript text.
