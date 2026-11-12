#!/usr/bin/env python3
# PP4 >= 0.7 colocalization hit counts, per-gene max PP4 and coloc.susie summaries.
# Inputs are the merged coloc and coloc.susie tables.
import os

import numpy as np
import pandas as pd

GENES = ["SCAF1", "NRBP2", "ATP2B1", "MTRFR", "OGFOD2", "ARL3"]

# --- helpers ---
def summarize_coloc(ac, pp4_cut=0.7, genes=GENES):
    ac = ac.copy()
    ac["gene"] = ac["Phenotype"].str.split("|").str[1]
    hi = ac[ac["PP4"] >= pp4_cut]
    out = {
        "coloc_rows_total": len(ac),
        "coloc_hits_raw": len(hi),
        "coloc_hits_unique_gene_gwas_stage": len(hi.drop_duplicates(["gene", "GWAS_ID", "QTL_ID"])),
        "coloc_hit_genes": hi["gene"].nunique(),
        "coloc_hit_gwas": hi["GWAS_ID"].nunique(),
    }
    for g in genes:
        out[f"coloc_{g}_maxpp4"] = round(ac.loc[ac["gene"] == g, "PP4"].max(), 3)
        x = ac[(ac["gene"] == g) & (ac["PP4"] >= pp4_cut)]
        out[f"coloc_{g}_hit_qtl"] = ";".join(sorted(set(x["QTL_ID"]))) if len(x) else "none"
    sub = ac[(ac["gene"] == "OGFOD2") & (ac["GWAS_ID"] == "EUR_INT") & (ac["QTL_ID"] == "prenatal")]
    out["coloc_OGFOD2_EUR_INT_pre_pp4"] = round(sub["PP4"].max(), 3) if len(sub) else np.nan
    return out

def summarize_susie(sus, genes=GENES, pp4_cut=0.7):
    ps = sus.groupby("source_file")["PP.H4.abf"].max()
    out = {
        "susie_comparisons": len(ps),
        "susie_comparisons_pp4ge0.7": int((ps >= pp4_cut).sum()),
    }
    out["susie_SCAF1_max"] = round(ps[ps.index.str.contains("_SCAF1_susie")].max(), 3)
    out["susie_NRBP2_max"] = round(ps[ps.index.str.contains("_NRBP2_susie")].max(), 3)
    out["susie_ATP2B1_n"] = int(ps.index.str.contains("_ATP2B1_susie").sum())
    out["susie_OGFOD2_EUR_INT_n"] = int((ps.index.str.contains("_OGFOD2_susie") &
                                         ps.index.str.contains("EUR_INT")).sum())
    return out

if __name__ == "__main__":
    root = os.environ["APA_ROOT"]
    pp4_cut = 0.7
    ac = pd.read_csv(os.path.join(root, "data/07_coloc/all_colocalization.csv"))
    sus = pd.read_csv(os.path.join(root, "data/07_coloc/all_susie_results.csv"))
    vals = {}
    vals.update(summarize_coloc(ac, pp4_cut))
    vals.update(summarize_susie(sus, pp4_cut=pp4_cut))
    for k, v in vals.items():
        print(f"{k}\t{v}")

    # --- write ---
    pd.DataFrame([{"key": k, "value": v, "source": "coloc_pp4_summary.py"} for k, v in vals.items()]).to_csv(
        os.path.join(root, "data/07_coloc/coloc_susie_summary_keys.tsv"), sep="\t", index=False)
