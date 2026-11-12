# Reconstruct the G x Stage gating table (S4) and its Delta-expression statistics.
# Reads data/05_gxs, data/01_quant, data/04_qtl and ref; writes data/05_gxs by default.
import argparse
import gzip
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# --- paths and settings ---
HERE = Path(__file__).resolve().parent
APA_ROOT = Path(os.environ["APA_ROOT"])
D = APA_ROOT / "data" / "05_gxs"
np.random.seed(2026)

# --- helpers ---
def delta_expr(expr_path, stage_of):
    """Postnatal minus prenatal mean log2CPM per ENSG (version stripped, first duplicate kept)."""
    ex = pd.read_csv(expr_path, sep="\t", index_col=0)
    ex.index = pd.Index(ex.index.astype(str)).str.split(".").str[0]
    ex = ex[~ex.index.duplicated(keep="first")]
    pre = [c for c in ex.columns if stage_of.get(c) == "Prenatal"]
    post = [c for c in ex.columns if stage_of.get(c) == "Postnatal"]
    d = ex[post].mean(axis=1) - ex[pre].mean(axis=1)
    out = {"delta": d, "n_pre": len(pre), "n_post": len(post)}
    del ex
    return out


# --- arguments ---
p = argparse.ArgumentParser()
p.add_argument("--hc3", default=str(D / "intx_hc3_core_groups.tsv"))
p.add_argument("--expr-nonfixed", default=str(APA_ROOT / "data/01_quant/expr_after_combat.tsv"))
p.add_argument("--expr-fixed", default=str(D / "expr_combat_fixed_local.tsv"))
p.add_argument("--meta", default=str(APA_ROOT / "data/04_qtl/CHB.all.meta"))
p.add_argument("--gtf", default=str(APA_ROOT / "ref/gencode.v40.annotation.gtf.gz"))
p.add_argument("--shared-list", default=str(HERE / "shared_set_ensg.txt"))
p.add_argument("--out-dir", default=str(D))
p.add_argument("--compare", default=None, help="existing S4 table (read-only reference)")
a = p.parse_args()

# --- stage labels and GENCODE transcript-to-gene map ---
meta = pd.read_csv(a.meta, sep="\t", dtype=str, usecols=["sample", "LifeStage"])
stage_of = dict(zip(meta["sample"], meta["LifeStage"]))
tmap = {}
with gzip.open(a.gtf, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1).split(".")[0]
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1).split(".")[0]
        tmap[t] = g

# --- class rows from the HC3 lead pairs ---
hc3 = pd.read_csv(a.hc3, sep="\t")
lead_enst = hc3["lead_phenotype_id"].astype(str).str.split("|").str[0].str.split(".").str[0]
hc3_ensg = lead_enst.map(tmap)
if hc3_ensg.isna().any():
    raise ValueError("lead transcript without GENCODE gene: %d rows" % int(hc3_ensg.isna().sum()))

d_v3 = delta_expr(a.expr_nonfixed, stage_of)
d_v1 = delta_expr(a.expr_fixed, stage_of)

cls_rows = pd.DataFrame({
    "gene_id": hc3["gene_id"],
    "ensg": hc3_ensg,
    "cls": ["Reversed" if np.sign(bp) != np.sign(bq) else ("Amplified" if abs(bq) > abs(bp) else "Attenuated")
            for bp, bq in zip(hc3["beta_pre"], hc3["beta_post"])],
    "beta_pre": hc3["beta_pre"],
    "beta_post": hc3["beta_post"],
    "dabs": hc3["beta_post"].abs() - hc3["beta_pre"].abs(),
})
cls_rows["dexpr"] = cls_rows["ensg"].map(d_v1["delta"])
cls_rows["group"] = cls_rows["cls"]
cls_rows["dexpr_v3"] = cls_rows["ensg"].map(d_v3["delta"])

# --- shared apaQTL rows (membership from the input list) ---
shared = [s.strip() for s in open(a.shared_list) if s.strip() and not s.startswith("#")]
shared = [s.split(".")[0] for s in shared]
sh_rows = pd.DataFrame({
    "gene_id": "", "ensg": shared, "cls": "", "beta_pre": np.nan, "beta_post": np.nan,
    "dabs": np.nan, "dexpr": pd.Series(shared).map(d_v1["delta"]).values,
    "group": "Shared apaQTL", "dexpr_v3": pd.Series(shared).map(d_v3["delta"]).values,
})
table = pd.concat([cls_rows, sh_rows], ignore_index=True)[
    ["gene_id", "ensg", "cls", "beta_pre", "beta_post", "dabs", "dexpr", "group", "dexpr_v3"]]
info = {"n_pre_samples": d_v3["n_pre"], "n_post_samples": d_v3["n_post"],
        "n_class_rows": int(len(cls_rows)), "n_shared_rows": int(len(sh_rows)),
        "class_counts": cls_rows["cls"].value_counts().to_dict()}

# --- gating statistics quoted in the manuscript ---
ig = table[table.cls.isin(["Reversed", "Amplified", "Attenuated"])].dropna(subset=["dexpr_v3", "dabs"])
sh = table[table.cls.isna() | (table.cls == "")].dropna(subset=["dexpr_v3"])
amp = ig[ig.cls == "Amplified"].dexpr_v3
att = ig[ig.cls == "Attenuated"].dexpr_v3
rho = stats.spearmanr(ig.dabs, ig.dexpr_v3)
st = {
    "n_interaction_finite": int(len(ig)),
    "median_abs_dexpr_v3_interaction": float(ig.dexpr_v3.abs().median()),
    "median_abs_dexpr_v3_shared_S4": float(sh.dexpr_v3.abs().median()),
    "P_wilcoxon_interaction_vs_shared_S4": float(stats.mannwhitneyu(ig.dexpr_v3.abs(), sh.dexpr_v3.abs()).pvalue),
    "spearman_rho_dabs_dexpr_v3": float(rho.statistic),
    "spearman_P": float(rho.pvalue),
    "median_dexpr_v3_amplified": float(amp.median()),
    "median_dexpr_v3_attenuated": float(att.median()),
    "n_amplified": int(len(amp)), "n_attenuated": int(len(att)),
    "P_wilcoxon_amp_vs_att_dexpr_v3": float(stats.mannwhitneyu(amp, att).pvalue),
    "P_wilcoxon_amp_vs_att_dexpr_older": float(stats.mannwhitneyu(
        ig[ig.cls == "Amplified"].dexpr.dropna(), ig[ig.cls == "Attenuated"].dexpr.dropna()).pvalue),
}

# --- write outputs ---
out = Path(a.out_dir)
out.mkdir(parents=True, exist_ok=True)
table.to_csv(out / "S4_gating_eqtl_fixed_RNA.reconstructed.tsv", sep="\t", index=False)
pd.Series({**st, **{k: v for k, v in info.items() if k != "class_counts"}}).to_csv(
    out / "gating_stats.tsv", sep="\t", header=["value"])
summary = {"info": info, "stats": st}

if a.compare:
    ref = pd.read_csv(a.compare, sep="\t")
    m = table.merge(ref, on="ensg", how="outer", suffixes=("", "_ref"), indicator=True)
    res = {"rows_recon": int(len(table)), "rows_ref": int(len(ref)),
           "merge": m["_merge"].value_counts().to_dict()}
    both = m[m["_merge"] == "both"]
    for c in ["dexpr", "dexpr_v3", "beta_pre", "beta_post", "dabs"]:
        a_, b_ = pd.to_numeric(both[c], errors="coerce"), pd.to_numeric(both[c + "_ref"], errors="coerce")
        ok = a_.notna() & b_.notna()
        res[f"max_abs_diff_{c}"] = float((a_[ok] - b_[ok]).abs().max()) if ok.any() else None
        res[f"n_compared_{c}"] = int(ok.sum())
    res["cls_mismatch"] = int(((both["cls"].fillna("") != both["cls_ref"].fillna(""))).sum())
    res["group_mismatch"] = int(((both["group"].fillna("") != both["group_ref"].fillna(""))).sum())
    summary["compare"] = res

(out / "run_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(json.dumps(summary, indent=2, default=str))
