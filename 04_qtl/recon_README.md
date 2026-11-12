# 04 reconstructed (ED7 seed recovery, apaQTL totals, mediation counts)

Reconstructed 2026-10-19. No existing release file is modified.

- `recon_ed7_counts.py`: ed7c_seed_recovery.tsv (exact), ed7_per_gene_recovery n_recovered_seeds (exact), apaQTL gene and pair totals (exact), mediation call counts and the 162/772 rule (exact recount of the producer's detail table).
- Not reconstructed: ED6 pi1 table (needs EUR nominal pairs, not local); ED7 per-gene columns other than n_recovered_seeds.

Run: `python3 recon_ed7_counts.py --out-dir OUT` (paths from APA_ROOT).

## Bugs and broken links (existing files; not fixed here)

1. `eqtl/04-1_mediate.R` lines 5-8: clusterProfiler and org.Hs.eg.db are loaded at top level. Mediation core cannot run without them. Proposed: guard the enrichment block behind a flag.
2. `eqtl/04-1_mediate.R`: the category rule behind the registry count 162 (= Partial mediation + APA -> EXP) is not written in the script. Proposed: document it in the script header.
3. Registry: R028/R029 "162 (21.0%)" lacks a stated denominator rule; see above.
