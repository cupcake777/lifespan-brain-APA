# Figures (manuscript v27, 2026-10-28)

Final code for main Figs 1–6 and Extended Data Figs 1–9. Every script was re-run from this folder on 2026-10-26 (Figs 1–4, ED1–ED9; byte- or pixel-identical to manuscript v25 except the seeded ED4/ED7 jitter) on 2026-10-28 (fig5_abc.R panel c, fig6.R) and on 2026-10-29 (fig5.R).

## Setup
`source release/config/paths.local.sh` (sets `APA_ROOT`, the analysis root that contains `release/`, `data/`, `ref/` and `raw/`, and optional `APA_FONT_DIR` with Arial TTFs). R scripts run with `Rscript <script>` in the conda environment `r` (see ENVIRONMENT.md); Python panel scripts with `python <script>`. Outputs (PNG 600 dpi + PDF) are written next to the scripts and are git-ignored.

## Run order
1. Python panels: `fig2_panel_d.py`, `fig3_panel_c.py`, `fig4_panel_c.py`.
2. Composites (any order) as listed below.
3. Helpers sourced automatically: `theme.R` (locked palette/theme), `locus.R`. `producers/` prepares summary inputs (see v22 notes in `_manifests/FIGURES_MANIFEST.tsv`).

| Figure | Composite | Also needs |
|---|---|---|
| Fig 1 | `fig1.R` | — |
| Fig 2 | `fig2.R` | fig2_panel_d.py (panel d PNG) |
| Fig 3 | `fig3.R` | fig3_panel_c.py (panel d PNG) |
| Fig 4 | `fig4.R` | fig4_panel_c.py (panel d PNG; icons in assets/) |
| Fig 5 | `fig5.R` (sources `fig5_abc.R` a–c and `fig5_de.R` d,e) | c reads the repaired coloc run (see Inputs); d = SCAF1 prenatal, e = KLC1 postnatal; d,e tracks from `producers/prep_fig5_tracks.py` (reads the individual-level VCF locally for r², writes summary tables only) |
| Fig 6 | `fig6.R` | sources `locus.R`; all PP4 and regional data from the repaired coloc run |
| ED Fig 1 | `ed1.R` | — |
| ED Fig 2 | `ed2.R` | — |
| ED Fig 3 | `ed3.R` | — |
| ED Fig 4 | `ed4.R` | jitter seed fixed (seed = 20261026) |
| ED Fig 5 | `ed5.R` | — |
| ED Fig 6 | `ed6.R` | — |
| ED Fig 7 | `ed7.R` | writes ED7_v6_sources.tsv; jitter seed fixed (seed = 20261026) |
| ED Fig 8 | `ed8.R` | — |
| ED Fig 9 | `ed9.R` | — |

## Not included
Per-individual inputs (sample lists, per-sample covariates, BAM coverage used for Fig 3d) are listed in `_manifests/*WITHHELD.tsv`. Superseded figure versions were moved out of the release (`figure_archive_20261026/` in the analysis tree).

## Colocalization inputs (manuscript v27)
Fig 5c and Fig 6 read the accepted repaired 19-GWAS coloc run under `APA_ROOT`:
`data/07_coloc/coloc_run/all_colocalization.csv` (accepted 10-28 run), `data/figures/fig6/{cohort,gtex}_eqtl_gwas_coloc_long.tsv`,
`data/figures/fig6/fig6_reg_tables_all.tsv`, plus the 193-gene G x Stage list
(`data/05_gxs/intx_hc3_core_groups.tsv`); interaction-term colocalizations are used only for those 193 genes.
The manuscript coloc numbers are recomputed with `07_coloc/ms_numbers.py`.
