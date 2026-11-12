#!/usr/bin/env python3
# S2 regulatory-element tables for Fig 2c: 3'UTR elements and CLIP peaks, quasi-Poisson rate ratios.
# Reconstructed from the Methods text; the original S2 producers are not in the release.
import os, re, gzip
import numpy as np
import pandas as pd
from scipy import stats
from patsy import cr
import statsmodels.api as sm

root = os.environ["APA_ROOT"]
OUT = os.environ.get("APA_S2_OUT", os.path.join(root, "data", "03_omics"))
P_DAPARS = os.environ.get("APA_S2_DAPARS", os.path.join(root, "raw", "all.apa.count"))
P_SITES = os.environ.get("APA_S2_SITES", os.path.join(root, "data", "02_traj", "02_sites_clustered.tsv"))
P_HG38 = os.environ.get("APA_S2_HG38", os.path.join(root, "ref", "hg38.fa"))
P_MATURE = os.environ.get("APA_S2_MATURE", os.path.join(root, "ref", "mature.fa"))
P_RPMM = os.environ.get("APA_S2_RPMM", os.path.join(root, "data", "03_omics", "miRNATissueAtlas2_human_miRNA_rpmm.parquet"))
P_POSTAR = os.environ.get("APA_S2_POSTAR", os.path.join(root, "ref", "human.txt.gz"))

RC = str.maketrans("ACGTN", "TGCAN")
MOTIFS = {"CPE_UUUUAU": "TTTTAT", "PUM_UGUAHAUA": "TGTA[ACT]ATA", "ARE_AUUUA": "ATTTA"}
PATS = {k: re.compile("(?=(%s))" % v) for k, v in MOTIFS.items()}
RBPS = {  # label -> (RBP name, tissue or None, method or None); AGO2 uses the Brain HITS-CLIP (CLIPper) subset
    "AGO2_brain": ("AGO2", "Brain", "HITS-CLIP,CLIPper"),
    "PUM1": ("PUM1", None, None),
    "PUM2": ("PUM2", None, "eCLIP"),
    "FMR1": ("FMR1", None, None),
}
COV_RBP = ("ELAVL1", None, "PAR-CLIP,PARalyzer")

# --- helpers ---
def fetch(fh, fai, chrom, s1, e1):
    row = fai.loc[chrom]
    def pos(i):
        return row.off + (i // row.lb) * row.lw + (i % row.lb)
    a = pos(s1 - 1); b = pos(e1 - 1) + 1
    fh.seek(a)
    return fh.read(b - a).replace(b"\n", b"").decode().upper()


def seed7(mseq):
    """7mer-m8 target site = reverse complement of miRNA nt 2-8 (0-based slice [1:8])."""
    if mseq is None or len(mseq) < 8:
        return None
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp[c] for c in reversed(mseq[1:8]))


def seed_code(s):
    v = 0
    for ch in s:
        v = v * 4 + "ACGT".index(ch)
    return v


def kmer_counts(seqs, lens, k):
    """Overlapping k-mer counts across all segments (gap 'N' separators prevent cross-segment matches)."""
    big = np.frombuffer(("N" * k).join(seqs).encode(), dtype=np.uint8)
    lut = np.full(256, 4, dtype=np.int64)
    for ch, v in zip(b"ACGT", range(4)):
        lut[ch] = v
    b = lut[big]
    n = len(b) - k + 1
    code = np.zeros(n, dtype=np.int64); valid = np.ones(n, dtype=bool)
    for j in range(k):
        bj = b[j:j + n]; valid &= bj < 4; code = code * 4 + np.where(bj < 4, bj, 0)
    return np.bincount(code[valid], minlength=4 ** k) / (lens.sum() / 1000.0)


def peaks_for(P, spec):
    rbp, tissue, method = spec
    d = P[P.rbp == rbp]
    if tissue:
        d = d[d.tissue == tissue]
    if method:
        d = d[d.method == method]
    return d


def overlap_counts(d, chrom, s1, e1):
    """Number of peaks overlapping each 1-based inclusive segment (0-based half-open peaks)."""
    out = np.zeros(len(s1), dtype=np.int64)
    pc = d.chrom.values.astype(str); ps_all = d.s.values; pe_all = d.e.values
    for c in np.unique(chrom):
        mp = pc == c
        if not mp.any():
            continue
        ps = np.sort(ps_all[mp]); pe = np.sort(pe_all[mp])
        idx = np.where(chrom == c)[0]
        a = s1[idx] - 1; b = e1[idx]
        out[idx] = np.searchsorted(ps, b, side="left") - np.searchsorted(pe, a, side="right")
    return out


def ns_design(x):
    """Natural cubic spline, 3 df, knots at tertiles of x (R ns(df = 3) space, intercept removed)."""
    q = np.quantile(x, [1 / 3, 2 / 3])
    B = np.asarray(cr(x, knots=q, lower_bound=x.min(), upper_bound=x.max()))
    return B[:, 1:]


def qp_fit(y, off, X):
    return sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(scale="X2")


def ci_p(b, se, df_):
    t = stats.t.ppf(0.975, df_)
    return np.exp(b), np.exp(b - t * se), np.exp(b + t * se), 2 * stats.t.sf(abs(b / se), df_)


def aUTR_model(m, mot, use_au):
    d = m[m.directional]
    y = np.concatenate([d["a_" + mot].values, d["c_" + mot].values]).astype(float)
    off = np.log(np.concatenate([d.alt_len.values, d.com_len.values]))
    isalt = np.concatenate([np.ones(len(d)), np.zeros(len(d))])
    cols = [np.ones(len(y)), isalt]
    if use_au:
        S_ = ns_design(np.concatenate([d.AU.values, d.AUc.values]))
        cols += [S_[:, i] for i in range(S_.shape[1])]
    r = qp_fit(y, off, np.column_stack(cols))
    return ci_p(r.params[1], r.bse[1], r.df_resid)


def group_model_alt(m, grp, ref, ycol, use_au=True, covariate=None):
    keep = grp | ref
    d = m[keep]
    y = d[ycol].values.astype(float) if isinstance(ycol, str) else ycol[keep].astype(float)
    off = np.log(d.alt_len.values)
    g = grp[keep].astype(float)
    cols = [np.ones(len(d)), g]
    if use_au:
        S_ = ns_design(d.AU.values)
        cols += [S_[:, i] for i in range(S_.shape[1])]
    if covariate is not None:
        cols.append(covariate[keep].astype(float))
    r = qp_fit(y, off, np.column_stack(cols))
    return ci_p(r.params[1], r.bse[1], r.df_resid)


def compare(out_df, ref, key_cols, val_cols):
    mm = out_df.merge(ref, on=key_cols, suffixes=("_rec", "_ref"))
    res = []
    for c in val_cols:
        a = mm[c + "_rec"].astype(float); b = mm[c + "_ref"].astype(float)
        rel = (a - b).abs() / b.abs().clip(lower=1e-300)
        res.append((c, float(rel.max()) if len(rel) else np.nan))
    return res

# --- reference S2 tables, read before the outputs below overwrite them ---
ref_files = {"S2_aUTR_vs_cUTR_elem": "S2_aUTR_vs_cUTR_elem.tsv",
             "S2_C1C2_vs_nondynamic_elem": "S2_C1C2_vs_nondynamic_elem.tsv",
             "S2_C1C2_vs_nondynamic_clip_adj": "S2_C1C2_vs_nondynamic_clip_adj.tsv"}
ref = {k: pd.read_csv(os.path.join(OUT, v), sep="\t") for k, v in ref_files.items() if os.path.exists(os.path.join(OUT, v))}

# --- events: DaPars2 proximal sites joined to the 63,097 clustered events ---
apa = pd.read_csv(P_DAPARS, sep="\t", usecols=[0, 1, 2, 3], dtype={"Loci": str})
apa["transcript_id"] = apa["Gene"].str.split("|").str[0]
apa["strand_d"] = apa["Gene"].str.split("|").str[3]
apa["chrom"] = apa["Loci"].str.split(":").str[0]
sp = apa["Loci"].str.split(":").str[1].str.split("-")
apa["S"] = sp.str[0].astype(np.int64)
apa["E"] = sp.str[1].astype(np.int64)
apa = apa.rename(columns={"Predicted_Proximal_APA": "P", "Loci": "locus"})
ev = pd.read_csv(P_SITES, sep="\t", low_memory=False)
m = ev.merge(apa[["transcript_id", "locus", "P", "strand_d", "chrom", "S", "E"]],
             on=["transcript_id", "locus"], how="left")
if m["P"].isna().any():
    raise RuntimeError("unmatched events between DaPars2 output and 02_sites_clustered")
if not (m["strand"] == m["strand_d"]).all():
    raise RuntimeError("strand disagreement")
plus = (m["strand"] == "+").values
P, S, E = m["P"].values, m["S"].values, m["E"].values
m["alt_s"] = np.where(plus, P, S);     m["alt_e"] = np.where(plus, E, P)
m["com_s"] = np.where(plus, S, P + 1); m["com_e"] = np.where(plus, P - 1, E)
m["alt_len"] = m["alt_e"] - m["alt_s"] + 1
m["com_len"] = m["com_e"] - m["com_s"] + 1
keep = (m["alt_len"] > 50) & (m["com_len"] > 50)   # 'length > 50' gives 62,918 events (text count)
m = m.loc[keep].reset_index(drop=True)
m["directional"] = m["directional_pass"].astype(str).str.lower().eq("true")
m["nondynamic"] = ~m["dynamic_any_pass"].astype(str).str.lower().eq("true")
print("events retained (both segments > 50 nt): %d; directional %d; C1 %d; C2 %d; non-dynamic %d" % (
    len(m), m.directional.sum(), (m.cluster_id == 1).sum(), (m.cluster_id == 2).sum(), m.nondynamic.sum()))

# --- sequences and motifs (sense strand, overlapping matches) ---
fai = pd.read_csv(P_HG38 + ".fai", sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                  names=["chr", "len", "off", "lb", "lw"]).set_index("chr")
alt, com = [], []
with open(P_HG38, "rb") as fh:
    for r in m.itertuples(index=False):
        a = fetch(fh, fai, r.chrom, r.alt_s, r.alt_e)
        c = fetch(fh, fai, r.chrom, r.com_s, r.com_e)
        if r.strand == "-":
            a = a.translate(RC)[::-1]; c = c.translate(RC)[::-1]
        alt.append(a); com.append(c)
m["alt_seq"] = alt; m["com_seq"] = com
m["AU"] = [(s.count("A") + s.count("T")) / len(s) for s in alt]
m["AUc"] = [(s.count("A") + s.count("T")) / len(s) for s in com]
for k, p in PATS.items():
    m["a_" + k] = np.fromiter((len(p.findall(s)) for s in m.alt_seq), np.int64, len(m))
    m["c_" + k] = np.fromiter((len(p.findall(s)) for s in m.com_seq), np.int64, len(m))

# --- groups ---
dirm = m.directional.values
C1 = (m.cluster_id == 1).values; C2 = (m.cluster_id == 2).values
ND = m.nondynamic.values

# --- Table 1: aUTR vs cUTR, directional events ---
rows = []
for mot in MOTIFS:
    crude = aUTR_model(m, mot, False); adj = aUTR_model(m, mot, True)
    d = m[dirm]
    ra = d["a_" + mot].sum() / d.alt_len.sum() * 1000; rc_ = d["c_" + mot].sum() / d.com_len.sum() * 1000
    rows.append(dict(motif=mot, **{"rate_aUTR_per_kb": ra, "rate_cUTR_per_kb": rc_,
                     "rate_ratio.Estimate": crude[0], "lo.Estimate": crude[1], "hi.Estimate": crude[2],
                     "p.Pr(>|t|)": crude[3], "rr_AUadj": adj[0], "lo_AUadj": adj[1], "hi_AUadj": adj[2], "p_AUadj": adj[3]}))
t1 = pd.DataFrame(rows)

# --- miRNA 7mer-m8 rows (set definitions unresolved; see notes) ---
# The mir7_brain_ctrl row (composition-matched control) is not reconstructed.
mat = {}; name = None
with open(P_MATURE) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
        elif name and name.startswith("hsa-"):
            mat[name] = line.upper().replace("U", "T")
pq = pd.read_parquet(P_RPMM)
mir_cols = [c for c in pq.columns if c.startswith("hsa-")]
row_brain = pq[pq["Tissue"] == "brain"][mir_cols].mean()
sets = {
    "brain_row_ge1000": row_brain[row_brain >= 1000].index.tolist(),
    "brain_row_top100": row_brain.dropna().sort_values(ascending=False).index[:100].tolist(),
}
mir_rows = []
mir_names = {"mir7_brain": sets["brain_row_ge1000"], "mir7_brain100": sets["brain_row_top100"]}
md = m[dirm]
for lab, names in mir_names.items():
    seeds = sorted({seed7(mat.get(n)) for n in names if seed7(mat.get(n))})
    ka = kmer_counts(md.alt_seq.tolist(), md.alt_len.values, 7)
    kc = kmer_counts(md.com_seq.tolist(), md.com_len.values, 7)
    a_rate = sum(ka[seed_code(s)] for s in seeds)
    c_rate = sum(kc[seed_code(s)] for s in seeds)
    mir_rows.append(dict(motif=lab, rate_aUTR_per_kb=a_rate, rate_cUTR_per_kb=c_rate, n_distinct_seeds=len(seeds)))
t1_mir = pd.DataFrame(mir_rows)

# --- Table 2: C1 / C2 vs non-dynamic (alternative segments, AU-adjusted) ---
rows2 = []
for mot in MOTIFS:
    for g, gm in [("C1", C1), ("C2", C2)]:
        est = group_model_alt(m, gm, ND, "a_" + mot, use_au=True)
        rows2.append(dict(motif=mot, group=g, **{"rate_ratio.Estimate": est[0], "lo.Estimate": est[1],
                                                 "hi.Estimate": est[2], "p.Pr(>|t|)": est[3]}))
t2 = pd.DataFrame(rows2)

# --- Table 3: CLIP peaks, accessibility adjusted with ELAVL1 peak count ---
want = {v[0] for v in RBPS.values()} | {COV_RBP[0]}
peak_rows = []
with gzip.open(P_POSTAR, "rt") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 8 and p[5] in want:
            peak_rows.append((p[0], int(p[1]), int(p[2]), p[5], p[6], p[7]))
P = pd.DataFrame(peak_rows, columns=["chrom", "s", "e", "rbp", "method", "tissue"])
ccov = overlap_counts(peaks_for(P, COV_RBP), m.chrom.values.astype(str), m.alt_s.values, m.alt_e.values)
rows3 = []
for lab, spec in RBPS.items():
    cnt = overlap_counts(peaks_for(P, spec), m.chrom.values.astype(str), m.alt_s.values, m.alt_e.values)
    for g, gm in [("C1", C1), ("C2", C2)]:
        est = group_model_alt(m, gm, ND, cnt, use_au=True, covariate=ccov)
        rows3.append(dict(rbp=lab, group=g, rate_ratio=est[0], lo=est[1], hi=est[2], p=est[3]))
t3 = pd.DataFrame(rows3)

# --- outputs ---
os.makedirs(OUT, exist_ok=True)
t1.to_csv(os.path.join(OUT, "S2_aUTR_vs_cUTR_elem.tsv"), sep="\t", index=False)
t2.to_csv(os.path.join(OUT, "S2_C1C2_vs_nondynamic_elem.tsv"), sep="\t", index=False)
t3.to_csv(os.path.join(OUT, "S2_C1C2_vs_nondynamic_clip_adj.tsv"), sep="\t", index=False)
t1_mir.to_csv(os.path.join(OUT, "S2_mirna_seed_rates_unresolved.tsv"), sep="\t", index=False)

# --- comparison with the reference tables ---
ver = []
for name, out_df, keys, vals in [
    ("S2_aUTR_vs_cUTR_elem", t1, ["motif"],
     ["rate_ratio.Estimate", "lo.Estimate", "hi.Estimate", "rr_AUadj", "lo_AUadj", "hi_AUadj"]),
    ("S2_C1C2_vs_nondynamic_elem", t2, ["motif", "group"],
     ["rate_ratio.Estimate", "lo.Estimate", "hi.Estimate"]),
    ("S2_C1C2_vs_nondynamic_clip_adj", t3, ["rbp", "group"],
     ["rate_ratio", "lo", "hi"]),
]:
    if name in ref:
        for col, mx in compare(out_df, ref[name], keys, vals):
            ver.append(dict(table=name, column=col, max_rel_diff=mx))
pd.DataFrame(ver).to_csv(os.path.join(OUT, "s2_verification.tsv"), sep="\t", index=False)
