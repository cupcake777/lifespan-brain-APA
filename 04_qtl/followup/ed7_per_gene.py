#!/usr/bin/env python3
# ED Fig 7: per-gene recovery of prenatal apaQTL genes in 20 postnatal subsets (n = 189).
# Gene-level recovery = qval < 0.05 in seed_<k>_postnatal.cis_qtl.txt.gz; gene = 2nd field of phenotype_id.
# One row per prenatal apaQTL gene (qval < 0.05), strongest event.
import os, glob
import pandas as pd

D = os.environ["APA_ROOT"] + "/data/04_qtl/"

# --- prenatal apaQTL genes ---
pre = pd.read_csv(D + "prenatal.cis_qtl.txt.gz", sep="\t",
                  usecols=["phenotype_id", "variant_id", "pval_nominal", "slope", "slope_se", "af", "qval"])
pre["gene"] = pre.phenotype_id.str.split("|").str[1]
g = pre[pre.qval < 0.05].sort_values("pval_nominal").drop_duplicates("gene").copy()

# --- recovery in the 20 postnatal subsets ---
seeds = sorted(glob.glob(D + "seed_*_post.cis_qtl.txt.gz")); assert len(seeds) == 20
hits = []
for s in seeds:
    d = pd.read_csv(s, sep="\t", usecols=["phenotype_id", "qval"])
    hits.append(set(d.loc[d.qval < 0.05, "phenotype_id"].str.split("|").str[1]))
g["n_recovered_seeds"] = [sum(e in h for h in hits) for e in g.gene]

# --- full postnatal map and prenatal effect size ---
full = pd.read_csv(D + "postnatal.cis_qtl.txt.gz", sep="\t", usecols=["phenotype_id", "qval"])
g["sig_full_postnatal"] = g.gene.isin(set(full.loc[full.qval < 0.05, "phenotype_id"].str.split("|").str[1]))
g["abs_z_prenatal"] = (g.slope / g.slope_se).abs()

# --- check against the downsampling CDF ---
cdf = pd.read_csv(D + "05_downsample_replication_cdf.tsv", sep="\t").sort_values("k")
assert [int((g.n_recovered_seeds >= k).sum()) for k in range(21)] == cdf.n_genes.tolist(), "CDF check failed"

cols = ["gene", "phenotype_id", "variant_id", "af", "slope", "slope_se", "abs_z_prenatal", "pval_nominal", "qval",
        "n_recovered_seeds", "sig_full_postnatal"]
g[cols].sort_values(["n_recovered_seeds", "abs_z_prenatal"]).to_csv(D + "ed7_per_gene_recovery.tsv", sep="\t", index=False)
print("genes", len(g), "never", int((g.n_recovered_seeds == 0).sum()), "all20", int((g.n_recovered_seeds == 20).sum()))
