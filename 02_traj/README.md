# Module 02: lifespan trajectories

The trajectory and external-validation producers are in this folder (`01_cov.R` … `12_ext_pattern.R`),
and the figure scripts are the `fig_*.R` files. Run the complete workflow
from this directory with `bash run.sh`, or submit it with `sbatch run.slurm`.
Inputs are read from `$APA_ROOT/data/02_traj/`, `$APA_ROOT/data/01_quant/` and `$APA_ROOT/ref/`;
intermediate and result tables are written to `$APA_ROOT/data/02_traj/`.

The workflow starts with `01_cov.R`. It averages the
ComBat-corrected logit-PDUI profiles within donor, writes one matrix column and
one covariate row per donor, and then fits all lifespan models on independent
donors. Technical-replicate profiles remain available to the upstream batch
and bridge checks but are not separate observations in Module 02 inference.
