# 05 reconstructed (G x Stage gating and HC3 tables)

Reconstructed 2026-10-19. Scripts here do not modify existing release files.

- `recon_s4_gating.py`: rebuilds `S4_gating_eqtl_fixed_RNA.tsv` (Fig 4c/4f-g, ED Fig 8b) and the Delta-expression statistics of registry R069-R076. Verification: partial (see `_manifests/RECONSTRUCTED_QTL.tsv`).
- `shared_set_ensg.txt`: membership of the 1,572 shared-apaQTL comparison genes. Taken from the existing table, not independently reproduced.

Run: `python3 recon_s4_gating.py --out-dir OUT [--compare REF.tsv]` (paths from APA_ROOT via config/paths.local.sh, or arguments).

## Bugs and broken links (existing files; not fixed here)

1. `../04_hc3_refit.py`: the pipeline runs at import time (no `if __name__ == "__main__":` guard). Proposed fix: move the top-level pipeline into `def main():` and call it under the guard. Effect: unit tests need ast extraction (done in tests).
2. `../04_hc3_refit.py` `point_estimate_class`: returns "Opposing" for sign reversals. The S4 table, Fig 4 and the manuscript use "Reversed" (122 rows). Proposed fix: return "Reversed" or map the label in the figure producer.
3. Registry `docs/number_registry.tsv`: `track_row_id` values R069-R078 are reused across m04_05 and m06_07 (e.g. R076 = N0001 and N0126). Proposed fix: make track_row_id unique, or key on reg_id.
4. Registry note and `results_final/RESULTS_MANIFEST.tsv` still say "no producer" / "NONE LOCATED" for S4_gating_eqtl_fixed_RNA.tsv. Proposed link: `release/05_gxs/reconstructed/recon_s4_gating.py`.
5. `results_final/05_gxstage_interaction/fig4/fig4g_shared_set_1579.tsv` (copy of 1014/figure/source_fig4): 8 rows flagged `matrix_all_stage_labelled_samples_approx` have no producer and do not reproduce (max |diff| 0.45). The text P = 6.4e-9 depends on them; the 1,572 S4 shared rows give 5.4e-9. Needs an author decision.
6. `../04_qtl/eqtl/04-1_mediate.R` lines 5-8 load clusterProfiler and org.Hs.eg.db at top level, so the mediation core cannot run without them. Proposed fix: guard the enrichment block (after line ~290) behind a flag.
7. The mediation count "162 of 772" (registry R028) has no documented rule in the script. Proposed: state in 04-1_mediate.R that mediated = "Partial mediation" or "APA -> EXP".
8. Sandbox: system Rscript (4.1.2) reports missing base packages (`rbinom` not found). Use a working R (conda R 4.5.3) for the synthetic test; `tests/synthetic/e2e/run_tests.sh` accepts RSCRIPT.
