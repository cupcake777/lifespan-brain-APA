# ED5d input: covariate-adjusted odds ratios of neurite-vs-soma gene sets per class (logistic, 11k-gene universe).
import gzip, zipfile, io, os, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm
warnings.filterwarnings("ignore")

ROOT = os.environ["APA_ROOT"] + "/"
D2 = ROOT + "data/02_traj"; D1 = ROOT + "data/01_quant"; REF = ROOT + "ref/"
OUT = ROOT + "data/figures/ed5"; os.makedirs(OUT, exist_ok=True)

ga = pd.read_csv(D2 + "/05_gene_aggregation.tsv", sep="\t")
ga["sym"] = ga.gene_symbol.astype(str).str.upper()
ga["cls"] = np.where(ga.gene_state == "C1-spec", "C1-specific",
             np.where(ga.gene_state == "C2-spec", "C2-specific", None))
ga.loc[ga.membership == "Shared", "cls"] = "Shared"
ga = ga[ga.cls.notna()].drop_duplicates("sym")

zt = zipfile.ZipFile(REF + "gky1270_supplemental_files.zip")
cio = pd.read_excel(io.BytesIO(zt.read("TableS2.xlsx")), sheet_name=0, skiprows=2)
ciolli = set(cio.Names.astype(str).str.split("|").str[-1].str.upper())
zs = zipfile.ZipFile(REF + "ScienceDirect_files_13Oct2026_13-32-29.121.zip")
tus = pd.read_excel(io.BytesIO(zs.read("1-s2.0-S0896627318302368-mmc5.xls")), sheet_name=2)
tushev = set(tus["gene.symbol"].astype(str).str.upper())

lengths = {}
with gzip.open(REF + "gencode.v40.annotation.gtf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "gene": continue
        name = f[8].split('gene_name "')[1].split('"')[0].upper() if 'gene_name "' in f[8] else None
        if name and name not in lengths:
            lengths[name] = int(f[4]) - int(f[3]) + 1
ga["gene_len"] = ga.sym.map(lengths)

E = pd.read_csv(D1 + "/expr_after_combat.tsv", sep="\t", index_col=0)
E.index = E.index.str.replace(r"\..*$", "", regex=True)
ga["mean_expr"] = ga.gene_id_base.map(E.mean(axis=1))
ga = ga.dropna(subset=["gene_len", "mean_expr"]).copy()
ga = ga[ga.n_sites > 0]
ga["log_n"] = np.log(ga.n_sites)
ga["log_len"] = np.log(ga.gene_len)

rows = []
for setname, S in [("Ciolli", ciolli), ("Tushev", tushev)]:
    ga["in_set"] = ga.sym.isin(S).astype(int)
    for cls in ["C1-specific", "C2-specific", "Shared"]:
        X = pd.DataFrame({"is_class": (ga.cls == cls).astype(float), "log_n": ga.log_n,
                          "log_len": ga.log_len, "mean_expr": ga.mean_expr})
        X = sm.add_constant(X)
        fit = sm.Logit(ga.in_set, X).fit(disp=0)
        b, se = fit.params["is_class"], fit.bse["is_class"]
        rows.append(dict(set=setname, cls=cls, n_class=int((ga.cls == cls).sum()),
                         n_set_in_class=int(((ga.cls == cls) & (ga.in_set == 1)).sum()),
                         OR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), P=fit.pvalues["is_class"]))
pd.DataFrame(rows).to_csv(OUT + "/ed5d_adjusted_or.tsv", sep="\t", index=False)
