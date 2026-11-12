# Module 07: locus colocalization and SMR

`coloc/`, `smr/` and `func/` contain the retained locus-level, SMR and
functional follow-up producers used for the manuscript figures and candidate
loci. `easycoloc/` is the tracked-file snapshot of upstream commit
`850178b5c99899d3e00b8655ffda716f705aa09c`; its `LICENSE` and
`environment.yml` are included.

The functional branch includes `run_analysis.sh`, `run.sh`, and the original
server configuration in `func/config.yml`. Its absolute input and reference
paths need adjustment on other systems. `run.sh` executes the configured
analysis; `run_analysis.sh` retains the earlier follow-up workflow.
