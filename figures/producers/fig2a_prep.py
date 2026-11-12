# Fig 2a input: per-event developmental change in PDUI and host-gene expression (same model as S5 coupling).
import os
import numpy as np
import pandas as pd
import patsy
from scipy.stats import spearmanr

ROOT = os.environ["APA_ROOT"]
D2 = f"{ROOT}/data/02_traj"; D4 = f"{ROOT}/data/04_qtl"; D1 = f"{ROOT}/data/01_quant"
OUT = f"{ROOT}/data/figures/fig2"
os.makedirs(OUT, exist_ok=True)
meta = pd.read_csv(f"{D4}/CHB.all.meta", sep="\t").set_index("sample")
sc = pd.read_csv(f"{D2}/02_sites_clustered.tsv", sep="\t", usecols=["site_id", "cluster_id", "gene_id"],
                 low_memory=False).dropna(subset=["cluster_id"]).set_index("site_id")
sup = pd.read_csv(f"{D2}/01_support.tsv", sep="\t")
ga = sup.grid_age.values; am = sup.age_model.values; ok = sup.is_supported.astype(bool).values
segs = np.split(np.arange(len(ok)), np.where(np.diff(ok.astype(int)) != 0)[0] + 1)
segs = [s for s in segs if ok[s[0]]]
e, l = segs[0], segs[1]

# APA matrix: chunked read, directional sites only, producer site order
want = set(sc.index); parts = []
for ch in pd.read_csv(f"{D4}/apa_after_combat.tsv", sep="\t", index_col=0, chunksize=4000):
    parts.append(ch.loc[ch.index[ch.index.isin(want)]])
Apa = pd.concat(parts)
Apa = Apa.loc[[s for s in sc.index if s in Apa.index]]
ma = meta.reindex(Apa.columns)
ka = (ma.age.notna() & ma.RIN.notna()).values
Apa = Apa.loc[:, ka]; ma = ma[ka]
xa = np.log10(np.maximum(ma.age.values + 1, 1e-6))
sa = (ma.sex.values == "Male").astype(float)
dmA = patsy.dmatrix("cr(x, df=5) + sexb + RIN", {"x": xa, "sexb": sa, "RIN": ma.RIN.values.astype(float)}, return_type="dataframe")
Y = Apa.values.T.astype(float)
bA = np.linalg.lstsq(dmA.values, Y, rcond=None)[0]
gA = np.asarray(patsy.build_design_matrices([dmA.design_info], {"x": am, "sexb": np.full(len(am), sa.mean()),
                "RIN": np.full(len(am), ma.RIN.astype(float).mean())})[0])
PA = gA @ bA
rA = (Y - dmA.values @ bA).std(axis=0, ddof=dmA.shape[1])
apa_early = (PA[e[-1]] - PA[e[0]]) / rA
apa_gap = (PA[l[0]] - PA[e[-1]]) / rA
apa_late = (PA[l[-1]] - PA[l[0]]) / rA

# host-gene expression
E = pd.read_csv(f"{D1}/expr_after_combat.tsv", sep="\t", index_col=0)
E.index = E.index.str.replace(r"\..*$", "", regex=True)
m = meta.reindex(E.columns)
keep = m.age.notna() & m.RIN.notna()
E = E.loc[:, keep]; m = m[keep]
am_s = np.log10(np.maximum(m.age.values + 1, 1e-6))
sexb = (m.sex.values == "Male").astype(float)
dm = patsy.dmatrix("cr(x, df=5) + sexb + RIN", {"x": am_s, "sexb": sexb, "RIN": m.RIN.values.astype(float)}, return_type="dataframe")
beta = np.linalg.lstsq(dm.values, E.values.T, rcond=None)[0]
grid = np.asarray(patsy.build_design_matrices([dm.design_info], {"x": am, "sexb": np.full(len(am), sexb.mean()),
                  "RIN": np.full(len(am), m.RIN.astype(float).mean())})[0])
P = grid @ beta
rsd = (E.values.T - dm.values @ beta).std(axis=0, ddof=dm.shape[1])

gagg = pd.read_csv(f"{D2}/05_gene_aggregation.tsv", sep="\t")
sym2ens = dict(zip(gagg.gene_symbol, gagg.gene_id_base))
ev = pd.DataFrame({"site_id": list(Apa.index), "cl": sc.loc[list(Apa.index), "cluster_id"].values,
                   "apa_early": apa_early, "apa_gap": apa_gap, "apa_late": apa_late})
ev["gene_id"] = ev.site_id.str.split("|").str[1].map(sym2ens)
ev["apa_dev"] = ev.apa_early + ev.apa_gap

exg = pd.DataFrame({"gene_id": E.index,
    "ex_early": (P[e[-1]] - P[e[0]]) / rsd,
    "ex_gap": (P[l[0]] - P[e[-1]]) / rsd,
    "ex_late": (P[l[-1]] - P[l[0]]) / rsd})
exg["ex_dev"] = exg.ex_early + exg.ex_gap
exg = exg.drop_duplicates("gene_id")
j = ev.merge(exg, on="gene_id", how="left", indicator=True)
jj = j[j._merge == "both"].copy()
jj.drop(columns=["_merge"]).to_csv(f"{OUT}/fig2a_event_coupling.tsv", sep="\t", index=False)
j[j._merge == "left_only"][["site_id", "cl", "gene_id"]].to_csv(f"{OUT}/fig2a_events_without_expr.tsv", sep="\t", index=False)
