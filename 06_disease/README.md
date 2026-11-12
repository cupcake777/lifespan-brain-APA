# Module 06: disease heritability and gene sets

`run_module06.sh` runs the two sLDSC estimands and keeps their outputs separate:

- APA-only: `sldsc/run_sldsc.sbatch`, with prenatal, postnatal,
  interaction and lifespan APA annotations.
- Joint apaQTL+eQTL maxCPP: `../04_qtl/eqtl/03_run_sldsc.sh`, with separate
  apaQTL and eQTL annotation columns for prenatal and postnatal stages.

Set `MODULE06_RUN_ROOT` for the output root. The joint branch accepts
`JOINT_APA_DIR`, `JOINT_EQTL_DIR`, `JOINT_SLDSC_REF`, `JOINT_GWAS_DIR`,
`JOINT_LDSC_DIR`, `JOINT_MERGE_ALLELES_DIR`, and the corresponding executable
overrides. Set `RUN_LOCAL=1` for a direct local or mock-data run; without it,
the module entrypoint submits itself to Slurm.

The APA-only summary is written below `apa_only/output/summary/`; joint
maxCPP results are written below `joint_maxcpp/sldsc/results/maxcpp/`.

`magma/` and `twas/` remain exposed as retained producers and are not invoked by
the sLDSC entrypoint.
