# Figure-to-code map

Each main figure carries one step of the scientific argument. Extended Data figures contain validation, sensitivity analyses and supporting measurements.

## Main figures

### Figure 1 — Lifespan APA landscape

**Entry point:** `figure_generation/main_figures/figure_1.R`

- cohort age, sex and assay coverage;
- lifespan-associated APA events;
- directional APA trajectories;
- trajectory-associated gene sets, functional enrichment and constraint.

**Analysis source:** `analysis/01_apa_quantification/`, `analysis/02_lifespan_trajectories/`.

### Figure 2 — Relationship between APA and other molecular layers

**Entry point:** `figure_generation/main_figures/figure_2.R`

- APA–RNA trajectory coupling;
- cis versus trans protein association evidence;
- sequence-motif and brain CLIP-RBP enrichment;
- representative PIK3CA transcript context;
- matched APA, RNA and protein trajectories.

**Panel producer:** `figure_generation/panels/figure_2d_gene_model.py`.

**Analysis source:** `analysis/03_multiomics_integration/`.

### Figure 3 — Stage-specific genetic regulation of APA

**Entry point:** `figure_generation/main_figures/figure_3.R`

- prenatal/postnatal apaQTL and eQTL overlap;
- fine-mapping signal around cleavage sites;
- PAS-hexamer-changing variants;
- CEP41 locus and genotype-stratified molecular phenotypes;
- apaQTL–eQTL colocalization and mediation;
- ACSF3 example.

**Panel producer:** `figure_generation/panels/figure_3d_cep41_locus.py`.

**Analysis source:** `analysis/04_stage_specific_qtl/`.

### Figure 4 — Genotype-by-stage effects

**Entry point:** `figure_generation/main_figures/figure_4.R`

- reversed, amplified and attenuated stage effects;
- representative genotype-stratified loci;
- detection in adult molecular-QTL maps;
- conceptual gating and saturation models.

**Panel producer:** `figure_generation/panels/figure_4d_mechanism.py`.

**Analysis source:** `analysis/05_genotype_by_stage/`.

### Figure 5 — Disease enrichment and colocalization landscape

**Entry point:** `figure_generation/main_figures/figure_5.R`

- joint S-LDSC comparison of apaQTL and eQTL annotations;
- MAGMA gene-set associations;
- apaQTL–GWAS colocalization overview;
- representative SCAF1 and KLC1 loci.

**Helpers:** `figure_5_overview.R`, `figure_5_loci.R`.

**Analysis source:** `analysis/06_disease_genetics/`.

### Figure 6 — Representative disease-linked regulatory loci

**Entry point:** `figure_generation/main_figures/figure_6.R`

- ATP2B1 genotype-by-stage signal and ADHD colocalization;
- MTRFR and OGFOD2 locus comparisons;
- ancestry-dependent ARL3 colocalization and molecular effects.

**Analysis source:** `analysis/04_stage_specific_qtl/`, `analysis/05_genotype_by_stage/`, `analysis/06_disease_genetics/`.

## Extended Data

- **ED Figure 1:** cohort flow, age coverage and APA-site validation — [`extended_data_1.R`](../figure_generation/extended_data_figures/extended_data_1.R)
- **ED Figure 2:** batch structure and correction diagnostics — [`extended_data_2.R`](../figure_generation/extended_data_figures/extended_data_2.R)
- **ED Figure 3:** trajectory thresholds, bootstrap stability and BrainSpan comparison — [`extended_data_3.R`](../figure_generation/extended_data_figures/extended_data_3.R)
- **ED Figure 4:** enrichment-background and gene-constraint sensitivity — [`extended_data_4.R`](../figure_generation/extended_data_figures/extended_data_4.R)
- **ED Figure 5:** multi-omic validation and neurite/soma enrichment — [`extended_data_5.R`](../figure_generation/extended_data_figures/extended_data_5.R)
- **ED Figure 6:** external apaQTL replication and cell-composition checks — [`extended_data_6.R`](../figure_generation/extended_data_figures/extended_data_6.R)
- **ED Figure 7:** cross-stage QTL sharing, external replication and lead-variant distance — [`extended_data_7.R`](../figure_generation/extended_data_figures/extended_data_7.R)
- **ED Figure 8:** genotype-by-stage sensitivity and proposed regulatory routes — [`extended_data_8.R`](../figure_generation/extended_data_figures/extended_data_8.R)
- **ED Figure 9:** disease-enrichment and background sensitivity — [`extended_data_9.R`](../figure_generation/extended_data_figures/extended_data_9.R)
