# Module 02 input contract

Controlled and large analysis inputs are not included in the publication code
package. Populate the input folders before running `run.sh`.

Required input groups are defined directly in the scripts:

- `$APA_ROOT/data/04_qtl/`: corrected APA matrix (`apa_after_combat.tsv`) and cohort metadata (`CHB.all.meta`), phenotype groups, BrainSpan
  inputs, and batch9 holdout inputs.
- `$APA_ROOT/data/01_quant/`: corrected donor-level APA, RNA, isoform and transcript-to-gene
  tables (e.g. `apa_combat_donor_mean.tsv`). `01_cov.R` creates the donor-level APA table from
  `data/04_qtl/apa_after_combat.tsv` and `data/04_qtl/CHB.all.meta`.
- `$APA_ROOT/ref/`: gnomAD constraint, MSigDB GO/KEGG, GWAS Catalog and HPO reference
  tables used by `05_func.R`.

Use consent-compatible access instructions or approved derived inputs when
transferring this package. Do not place individual-level data in Git.
