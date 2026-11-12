#!/usr/bin/env python3
# Summarise the Nan 2026 cell-type lookups: replication of stage-specific apaQTL leads by cell type and in any cell type.
# Writes data/04_qtl/nan2026_comparison.xlsx (by_celltype, any_celltype, focal_genes, lookup_all).
import os, glob
import pandas as pd, numpy as np
from scipy.stats import fisher_exact

D = os.environ["APA_ROOT"] + "/data/04_qtl/"
L = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(D + "nan_lookup_*.tsv")])

# --- one row per our lead x cell type: prefer same transcript, else min p ---
L["o"] = (~L.same_transcript).astype(int)
L = L.sort_values(["o", "nan_p"]).drop_duplicates(["stage", "gene", "celltype"])
L["concordant"] = np.sign(L.slope) == np.sign(L.nan_slope_aligned)
L["rep"] = (L.nan_p < 0.05) & L.concordant

# --- gene groups ---
pre = pd.read_csv(D + "prenatal.cis_qtl.txt.gz", sep="\t"); preG = set(pre.loc[pre.qval < 0.05, "group_id"])
post = pd.read_csv(D + "postnatal.cis_qtl.txt.gz", sep="\t"); postG = set(post.loc[post.qval < 0.05, "group_id"])
rec = pd.read_excel(D + "Supplementary_Tables_draft.xlsx", sheet_name="S_downsampling")
never = set(rec.loc[rec.n_seeds_recovered == 0, "gene"]); absent = preG - postG

groups = []
for stage, gene in zip(L.stage, L.gene):
    if stage == "postnatal":
        groups.append("postnatal")
    elif gene in never:
        groups.append("prenatal: never recovered (144)")
    elif gene in absent:
        groups.append("prenatal: absent from postnatal map (106)")
    else:
        groups.append("prenatal: shared")
L["group"] = groups

# --- replication by group and cell type ---
rows = []
for (grp, ct), d in L.groupby(["group", "celltype"]):
    p = np.asarray(d.nan_p)
    pi1 = max(0.0, 1 - min(1, np.mean(p > 0.5) / 0.5)) if len(p) > 20 else np.nan
    rows.append(dict(group=grp, celltype=ct, n_matched=len(d), pi1=pi1, frac_p05=(d.nan_p < 0.05).mean(),
                     concord_p05=d.loc[d.nan_p < 0.05, "concordant"].mean(), replicated=d.rep.mean()))
T1 = pd.DataFrame(rows)

# --- replication in any cell type per gene ---
any_ct = L.groupby(["stage", "gene", "group"]).agg(n_ct=("celltype", "nunique"), replicated_any=("rep", "any")).reset_index()
T2 = any_ct.groupby("group").agg(genes_matched=("gene", "size"), replicated_any=("replicated_any", "mean")).reset_index()
nG = {"postnatal": len(postG), "prenatal: shared": len(preG - absent - never), "prenatal: absent from postnatal map (106)": len(absent - never),
      "prenatal: never recovered (144)": len(never)}
T2["genes_total"] = T2.group.map(nG); T2["matched_frac"] = T2.genes_matched / T2.genes_total

# --- Fisher: prenatal-specific (absent | never) vs postnatal, replicated in any cell type ---
a = any_ct[any_ct.group.str.contains("absent|never")]; b = any_ct[any_ct.group == "postnatal"]
tab = [[a.replicated_any.sum(), (~a.replicated_any).sum()], [b.replicated_any.sum(), (~b.replicated_any).sum()]]
OR, P = fisher_exact(tab)

# --- focal genes ---
FL = L[L.gene.isin(["SCAF1", "NRBP2", "ATP2B1", "MTRFR", "OGFOD2", "ARL3", "HID1", "TESPA1", "SLCO1A2"])][
    ["stage", "gene", "variant_id", "rsid", "celltype", "slope", "nan_slope_aligned", "nan_p", "same_transcript", "concordant"]].sort_values(["gene", "stage", "nan_p"])

with pd.ExcelWriter(D + "nan2026_comparison.xlsx") as w:
    T1.to_excel(w, sheet_name="by_celltype", index=False); T2.to_excel(w, sheet_name="any_celltype", index=False)
    FL.to_excel(w, sheet_name="focal_genes", index=False)
    L.drop(columns=["o", "key"]).to_excel(w, sheet_name="lookup_all", index=False)
print(T1.round(3).to_string()); print(T2.round(3).to_string())
print("prenatal-specific vs postnatal replicated_any: OR", round(OR, 2), "P", P, tab)
print(FL.groupby(["gene", "stage"]).apply(lambda d: f"{len(d)} ct; min p {d.nan_p.min():.2g}; concordant@min {d.sort_values('nan_p').concordant.iloc[0]}").to_string())
