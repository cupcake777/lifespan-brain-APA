#!/usr/bin/env python3
# Look up our stage-specific apaQTL leads in Nan et al. 2026 cell-type 3'aQTL nominal sumstats.
# Usage: python nan2026_compare.py <celltype>; writes data/04_qtl/nan_lookup_<celltype>.tsv
import os, sys
import pandas as pd, numpy as np

ROOT = os.environ["APA_ROOT"]
D = ROOT + "/data/04_qtl/"
N = ROOT + "/ref/nan2026/"
ct = sys.argv[1]

# --- helpers ---
def leads(stage):
    d = pd.read_csv(D + f"{stage}.cis_qtl.txt.gz", sep="\t")
    d = d[d.qval < 0.05].copy()
    d["stage"] = stage
    return d

# --- our leads, keyed by gene and position ---
L = pd.concat([leads("prenatal"), leads("postnatal")])
v = L.variant_id.str.split(":", expand=True)
L["chr"], L["pos"], L["ref"], L["alt"] = v[0], v[1].astype(int), v[2], v[3]
L["gene"] = L.group_id
L["enst"] = L.phenotype_id.str.split("|").str[0].str.split(".").str[0]
L["key"] = L.gene + "@" + L.chr + ":" + L.pos.astype(str)
keys = set(L.key)

# --- stream Nan nominal sumstats and keep rows at our keys ---
out = []
for ch in pd.read_csv(N + f"{ct}.3aQTL.sumstats.txt.gz", sep="\t", chunksize=1_000_000,
                      usecols=["variant_id", "CHR", "position", "Reference_allele", "Effect_allele",
                               "phenotype_id", "maf", "pval_nominal", "slope", "slope_se"]):
    g = ch.phenotype_id.str.split("|").str[1]
    k = g + "@" + ch.CHR.astype(str) + ":" + ch.position.astype(str)
    m = k.isin(keys)
    if m.any():
        s = ch[m].copy(); s["key"] = k[m]
        out.append(s)
R = pd.concat(out) if out else pd.DataFrame()
R["nan_enst"] = R.phenotype_id.str.split("|").str[0].str.split(".").str[0]

# --- join and align alleles (our slope is per ALT allele) ---
M = L[["stage", "gene", "phenotype_id", "variant_id", "enst", "ref", "alt", "slope", "pval_nominal", "key"]].merge(
    R.rename(columns={"phenotype_id": "nan_phenotype", "slope": "nan_slope", "pval_nominal": "nan_p",
                      "variant_id": "rsid"}), on="key", how="inner")
same = (M.Effect_allele == M.alt) & (M.Reference_allele == M.ref)
flip = (M.Effect_allele == M.ref) & (M.Reference_allele == M.alt)
M = M[same | flip].copy()
M["nan_slope_aligned"] = np.where(same[same | flip], M.nan_slope, -M.nan_slope)
M["same_transcript"] = M.enst == M.nan_enst
M["celltype"] = ct
M.to_csv(D + f"nan_lookup_{ct}.tsv", sep="\t", index=False)
print(ct, len(L), "leads;", M.key.nunique(), "matched", len(M), "rows")
