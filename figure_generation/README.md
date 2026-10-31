# Manuscript figure generation

This directory contains the final plotting layer. It is downstream of the analysis modules and is kept separate so that a reviewer can distinguish statistical inference from presentation code.

## Layout

```text
main_figures/              Figure 1 through Figure 6
extended_data_figures/     Extended Data Figure 1 through Figure 9
panels/                    Python producers for selected locus/schematic panels
shared/                    common theme and locus plotting functions
```

## Running a figure

After configuring `APA_CODE_ROOT`, `APA_ROOT` and `APA_FIGURE_DIR`:

```bash
bash scripts/render_figures.sh figure_1
bash scripts/render_figures.sh extended_data_4
```

The entry point runs any required Python panel producer before assembling its corresponding R figure.

## Input contract

Figure scripts read finalized aggregate outputs from the analysis modules. They may check expected columns, classes or row relationships, but do not refit the underlying models. Expected paths follow [`../docs/input_layout.md`](../docs/input_layout.md).

## Traceability

[`../docs/figure_map.md`](../docs/figure_map.md) maps every main and Extended Data figure to its plotting entry point and supporting analysis module.

Generated PDF/PNG files are written to `APA_FIGURE_DIR`.
