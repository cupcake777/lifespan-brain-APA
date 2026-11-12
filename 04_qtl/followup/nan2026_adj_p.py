#!/usr/bin/env python3
# Reproduce the two adjusted Nan 2026 P values: Fig 3e (P = 0.004) and Fig 4e (P = 0.001).
# Inputs: nan2026_comparison.xlsx (lookup_all, from nan_compare.py and nan_summarize.py),
# braingvex_adult_bulk_comp.xlsx (prenatal_leads |z| and MAF; GxStage_leads postnatal |z|),
# gxstage_strict_classes_local_checks.tsv (193 HC3 core genes, nan_p_<celltype>, prenatal_maf).
import os
import pandas as pd, numpy as np, statsmodels.formula.api as smf, statsmodels.api as sm
from scipy.stats import fisher_exact

D = os.environ["APA_ROOT"] + "/data/04_qtl/"
rows = []
def add(test, model, n, OR, P, note=""): rows.append(dict(test=test, model=model, n=n, OR=OR, P=P, note=note))

# --- A: prenatal-specific vs prenatal-shared, gene-level replication in any cell type (Fig 3e) ---
L = pd.read_excel(D + "nan2026_comparison.xlsx", sheet_name="lookup_all")
L["o"] = (~L.same_transcript).astype(int)
L = L.sort_values(["o", "nan_p"]).drop_duplicates(["stage", "gene", "celltype"])
A = L[L.stage == "prenatal"].groupby(["gene", "group"]).agg(rep=("rep", "any")).reset_index()
A["ps"] = A.group.str.contains("absent|never").astype(int); A["rep"] = A.rep.astype(int)
pl = pd.read_excel(D + "braingvex_adult_bulk_comp.xlsx", sheet_name="prenatal_leads")
A = A.merge(pl[["group_id", "absz", "maf"]], left_on="gene", right_on="group_id")
t = A.groupby("ps").rep.agg(["sum", "size"])
cnt = f"prenatal-specific {t.loc[1,'sum']}/{t.loc[1,'size']} vs shared {t.loc[0,'sum']}/{t.loc[0,'size']}"
add("A Fig3e", "Fisher", len(A), np.nan, fisher_exact(pd.crosstab(A.ps, A.rep))[1], cnt)
m = smf.logit("rep ~ ps + absz + maf", A).fit(disp=0)
add("A Fig3e", "logit rep ~ group + |z|_prenatal + MAF  [Methods; text 0.004]", len(A), np.exp(m.params.ps), m.pvalues.ps, cnt)

# --- B: Reversed vs Amplified+Attenuated, gene x cell-type nominal P < 0.05 (Fig 4e) ---
gx = pd.read_csv(D + "gxstage_strict_classes_local_checks.tsv", sep="\t")
ct = ["Ast", "Exc", "Inh", "Mic", "Oli", "Opc"]
B = gx.melt(id_vars=["gene_id", "stage_dependent_class", "prenatal_maf"], value_vars=[f"nan_p_{c}" for c in ct],
            var_name="celltype", value_name="p").dropna(subset=["p"])
B["sig"] = (B.p < 0.05).astype(int); B["rev"] = (B.stage_dependent_class == "Reversed").astype(int)
gl = pd.read_excel(D + "braingvex_adult_bulk_comp.xlsx", sheet_name="GxStage_leads")
B = B.merge(gl[["gene", "absz_post"]], left_on="gene_id", right_on="gene")
s = B.groupby("stage_dependent_class").sig.agg(["sum", "size"])
cnt = "; ".join(f"{k} {r['sum']}/{r['size']} ({100*r['sum']/r['size']:.0f}%)" for k, r in s.iterrows()) + f"; genes {B.gene_id.nunique()}"
grp = pd.factorize(B.gene_id)[0]
add("B Fig4e", "Fisher, pairs (rev vs other)", len(B), np.nan, fisher_exact(pd.crosstab(B.rev, B.sig))[1], cnt)
for f, lab in [("sig ~ rev + prenatal_maf + C(celltype)", "MAF + cell type  [Methods; text 0.001]"),
               ("sig ~ rev + prenatal_maf + C(celltype) + absz_post", "MAF + cell type + |z|_postnatal")]:
    m = smf.logit(f, B).fit(disp=0)
    mc = smf.logit(f, B).fit(disp=0, cov_type="cluster", cov_kwds={"groups": grp})
    g = smf.gee(f, groups="gene_id", data=B, family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
    add("B Fig4e", "logit pairs, naive SE: " + lab, len(B), np.exp(m.params.rev), m.pvalues.rev)
    add("B Fig4e", "logit pairs, gene-clustered SE: " + lab, len(B), np.exp(mc.params.rev), mc.pvalues.rev)
    add("B Fig4e", "GEE exchangeable by gene: " + lab, len(B), np.exp(g.params.rev), g.pvalues.rev)
G = B.groupby(["gene_id", "rev"]).sig.max().reset_index()
add("B Fig4e", "Fisher, gene-level any cell type", len(G), np.nan, fisher_exact(pd.crosstab(G.rev, G.sig))[1],
    f"Reversed {G[G.rev==1].sig.sum()}/{(G.rev==1).sum()} vs other {G[G.rev==0].sig.sum()}/{(G.rev==0).sum()}")

out = pd.DataFrame(rows); out.to_csv(D + "r4_nan_adjusted_p.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(out.assign(OR=out.OR.round(3), P=out.P.map(lambda x: f"{x:.2g}")).to_string(index=False))
