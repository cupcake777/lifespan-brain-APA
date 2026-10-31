# Input layout

The analysis uses a common data root defined by `APA_ROOT`. Scripts resolve study inputs, intermediate products and reference resources relative to this directory.

```text
${APA_ROOT}/
  raw/                    cohort source files
  ref/                    GRCh38, GENCODE and external reference resources
  data/
    01_quant/             APA quantification and harmonization
    02_traj/              lifespan trajectory models
    03_omics/             APA–RNA–protein integration
    04_qtl/               stage-specific molecular QTL analyses
    05_gxs/               genotype-by-stage analyses
    06_disease/           disease-enrichment analyses
    07_coloc/             colocalization analyses
    figures/              finalized inputs consumed by figure scripts
```

Each analysis-module README lists its required inputs and principal outputs. Figure-specific aggregate tables are organized beneath `data/figures/<figure_id>/`.

External resources include GRCh38, GENCODE, PolyASite, GTEx, GWAS summary statistics and the Nan et al. 3′aQTL resource. Exact resource versions and checksums should be recorded with the corresponding analysis run.
