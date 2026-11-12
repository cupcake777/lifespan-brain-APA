# Module 01: RNA-seq and APA quantification

Paths come from `config/paths.local.sh` (`APA_ROOT`, `APA_HPC_HOME`).
Outputs go to `data/01_quant/`, figures to `figures/01_quant/`, logs to `logs/`.

RNA chain
- `rna_quant.sh`: RSEM/STAR per sample (GENCODE v40), then merged gene and isoform
  count/TPM matrices `data/01_quant/all.<N>.*` with a sample manifest.
- `rna_pre_step01.r [gene|transcript]`: alias fix, XIST/UTY sex check, expression filter,
  TMM and VST, outlier removal, mean-only ComBat (LifeStage protected), INT, and for genes
  the QTL covariates. Reads `raw/all.912.*`.
- `rna_pre_step02.r`: aligns genotype, expression and covariate samples; writes the PGEN,
  tensorQTL BED and covariate table to `data/01_quant/qtl_input/`.
- `batch_correction.R`: compares Raw, ComBat variants, limma, RUV-III and RUVs on bridge pairs;
  writes the diagnostics bundle flat into `data/01_quant/` and selects one method into
  `data/04_qtl/apa_after_combat.tsv`.
- `bridge_select.R`: bridge agreement and donor-paired bootstrap of ComBat_REF vs limma.
- `rna_collapse_check.R`: within-batch collapse check (no caller; kept as a manual QC).
- `cohort_overview.R`: cohort age overview table and figure.
- `bridge_select.R` reads the flat diagnostics in `data/01_quant/` (see `batch_correction.R`).

APA chain (`apa/`)
- `run_dapars2.sh` drives `00_bam_to_wig.sh`, `prep_apa_inputs.sh`, `dapars2_by_chr.sh`, `dapars2_merge.sh`.
- `DaPars_Extract_Anno.py` and `Dapars2_Multi_Sample.py` are vendored DaPars2 code (unedited).
- `read_depth.py`, `dapars2_config.py`, `extract_3utr.py` are helpers for the APA inputs.
- Output: `data/01_quant/apa/all.apa.count`; the batch step reads the consolidated copy `raw/all.apa.count`.

SOAPnuke v1.5.6 filtering was done by the sequencing company (Methods checkpoint; FASTQ not retained).
