# ATP2B1 saturation check: refit the HC3-core G x Stage leads after dropping bound PDUI values.
# Reads data/05_gxs and raw/genotype; writes data/05_gxs.
from __future__ import annotations

# --- load packages ---
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib
from scipy import stats

# --- paths and settings ---
APA_ROOT = Path(os.environ["APA_ROOT"])
D = APA_ROOT / "data" / "05_gxs"
PHENOTYPE_FILE = D / "pre.phenotype.combat.int.txt"
RAW_PDUI_FILE = D / "before.combat.txt"
COVARIATE_FILE = D / "covariates_for_qtl.txt"
STAGE_FILE = D / "LifeStage.txt"
CORE_FILE = D / "intx_hc3_core_groups.tsv"
GENOTYPE_PREFIX = APA_ROOT / "raw" / "genotype"
OUT_DIR = D
LOW, HIGH = 0.001, 0.999
MIN_STAGE_N, MIN_CARRIERS = 10, 3
FOCUS_GENE = "ATP2B1"

# --- helpers ---
def residualizer(cov: np.ndarray) -> np.ndarray:
    """Orthonormal basis of centered covariates; rank-revealing, drops collinear columns."""
    c = cov - cov.mean(axis=0, keepdims=True)
    keep = c.std(axis=0) > 0
    c = c[:, keep]
    u, s, _ = np.linalg.svd(c, full_matrices=False)
    rank = int((s > s.max() * max(c.shape) * np.finfo(float).eps).sum())
    return u[:, :rank]


def resid(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float) - np.mean(v)
    return v - q @ (q.T @ v)


def fit(y: np.ndarray, cols: list[np.ndarray], q: np.ndarray) -> dict:
    """OLS after covariate residualization (as in 04_hc3_refit.py) with HC3 covariance."""
    yr = resid(y, q)
    x = np.column_stack([resid(c, q) for c in cols])
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ yr
    e = yr - x @ beta
    df = int(y.size - q.shape[1] - 1 - x.shape[1])
    lev = np.einsum("ij,jk,ik->i", x, xtx_inv, x)
    w = (e / np.maximum(1 - lev, 1e-10)) ** 2
    hc3 = xtx_inv @ ((x.T * w) @ x) @ xtx_inv
    hc3 = (hc3 + hc3.T) / 2
    return {"beta": beta, "hc3": hc3, "df": df}


def stage_effects(f: dict) -> dict:
    """beta_pre, beta_post = beta_G + beta_GxS, HC3 CIs, P_GxStage (columns: G, stage, GxS)."""
    b, v, df = f["beta"], f["hc3"], f["df"]
    crit = stats.t.ppf(0.975, df)
    se_pre = np.sqrt(v[0, 0])
    b_post = b[0] + b[2]
    se_post = np.sqrt(max(v[0, 0] + v[2, 2] + 2 * v[0, 2], 0))
    p_gxs = 2 * stats.t.sf(abs(b[2] / np.sqrt(v[2, 2])), df)
    return {
        "beta_pre": b[0], "pre_ci_low": b[0] - crit * se_pre, "pre_ci_high": b[0] + crit * se_pre,
        "beta_post": b_post, "post_ci_low": b_post - crit * se_post, "post_ci_high": b_post + crit * se_post,
        "hc3_p_gxstage": p_gxs, "df": df,
    }


def carriers(g: np.ndarray) -> int:
    g = np.rint(g[np.isfinite(g)]).astype(int)
    if g.size == 0:
        return 0
    minor_is_alt = g.mean() / 2 <= 0.5
    return int(((g > 0) if minor_is_alt else (g < 2)).sum())


def read_rows(path: Path, ids: set[str], samples: list[str]) -> dict[str, np.ndarray]:
    out = {}
    for chunk in pd.read_csv(path, sep="\t", index_col=0, chunksize=4096):
        chunk.index = chunk.index.astype(str)
        hit = chunk.index.intersection(list(ids))
        for pid in hit:
            out[pid] = chunk.loc[pid, samples].to_numpy(float)
    missing = ids - set(out)
    if missing:
        raise RuntimeError(f"{path}: missing phenotypes {sorted(missing)[:5]}")
    return out


# --- load core table, samples, covariates and stage ---
core = pd.read_csv(CORE_FILE, sep="\t")
print(f"[input] core table {CORE_FILE} ({len(core)} groups)", flush=True)

psam = pd.read_csv(f"{GENOTYPE_PREFIX}.psam", sep="\t")
samples = psam.iloc[:, 0].astype(str).tolist()
cov = pd.read_csv(COVARIATE_FILE, sep="\t", index_col=0)
cov.columns = cov.columns.astype(str)
cov = cov.loc[:, samples].apply(pd.to_numeric).to_numpy(float).T
st = pd.read_csv(STAGE_FILE, sep="\t", usecols=[0, 1])
stage = st.set_index(st.columns[0]).iloc[:, 0]
stage.index = stage.index.astype(str)
stage = stage.loc[samples].to_numpy(float)
q_all = residualizer(cov)

# --- phenotypes, raw PDUI and genotypes of the leads ---
pids = set(core["lead_phenotype_id"].astype(str))
pheno = read_rows(PHENOTYPE_FILE, pids, samples)
raw = read_rows(RAW_PDUI_FILE, pids, samples)

variants = core["lead_variant_id"].astype(str).tolist()
index, i = {}, 0
wanted = set(variants)
with open(f"{GENOTYPE_PREFIX}.pvar") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        vid = line.split("\t", 3)[2]
        if vid in wanted:
            index[vid] = i
        i += 1
reader = pgenlib.PgenReader(f"{GENOTYPE_PREFIX}.pgen".encode())
geno_raw = np.empty((len(variants), len(samples)), dtype=np.int8)
reader.read_list(np.asarray([index[v] for v in variants], dtype=np.uint32), geno_raw)
reader.close()
geno = geno_raw.astype(float)
geno[geno < 0] = np.nan
print("[input] phenotypes, raw PDUI and genotypes loaded", flush=True)

# --- refit each lead on all samples, then after removing bound PDUI values ---
records = []
for k, row in core.reset_index(drop=True).iterrows():
    pid = str(row["lead_phenotype_id"])
    g = geno[k].copy()
    g[np.isnan(g)] = np.nanmean(g)
    y, r = pheno[pid], raw[pid]
    full = stage_effects(fit(y, [g, stage, g * stage], q_all))
    if abs(full["beta_pre"] - row["beta_pre"]) > 1e-4 or abs(full["beta_post"] - row["beta_post"]) > 1e-4:
        raise RuntimeError(f"refit does not reproduce core table for {row['group_id']}: wrong inputs")
    at_bound = (r <= LOW + 1e-12) | (r >= HIGH - 1e-12)
    keep = np.isfinite(r) & ~at_bound & np.isfinite(y)
    rec = {"group_id": row["group_id"], "gene_id": row["gene_id"], "lead_phenotype_id": pid,
           "lead_variant_id": row["lead_variant_id"], "point_estimate_class": row.get("point_estimate_class"),
           "strict_ci_direction_reversal": row.get("strict_ci_direction_reversal")}
    rec.update({f"full_{a}": b for a, b in full.items()})
    for s, lab in ((0, "pre"), (1, "post")):
        m = stage == s
        rec[f"{lab}_n"] = int(m.sum())
        rec[f"{lab}_n_at_bound"] = int((at_bound & m).sum())
        rec[f"{lab}_n_kept"] = int((keep & m).sum())
        rec[f"{lab}_carriers_kept"] = carriers(g[keep & m])
    ok = all(rec[f"{l}_n_kept"] >= MIN_STAGE_N and rec[f"{l}_carriers_kept"] >= MIN_CARRIERS for l in ("pre", "post"))
    rec["boundary_status"] = "complete" if ok else "not_evaluable"
    if ok:
        try:
            bf = stage_effects(fit(y[keep], [g[keep], stage[keep], g[keep] * stage[keep]], residualizer(cov[keep])))
            rec.update({f"bound_{a}": b for a, b in bf.items()})
            rec["bound_strict_reversal"] = bool(
                (bf["pre_ci_low"] > 0 and bf["post_ci_high"] < 0) or (bf["pre_ci_high"] < 0 and bf["post_ci_low"] > 0))
        except np.linalg.LinAlgError as exc:
            rec["boundary_status"] = f"not_estimable: {exc}"
    records.append(rec)

out = pd.DataFrame(records)
for col in ("bound_hc3_p_gxstage", "bound_strict_reversal"):
    if col not in out:
        out[col] = np.nan
out["bound_strict_reversal"] = out["bound_strict_reversal"].fillna(False).astype(bool)
out.to_csv(OUT_DIR / "boundary_refit_hc3_core.tsv", sep="\t", index=False)
done = out["boundary_status"] == "complete"
summary = pd.DataFrame([{
    "n_groups": len(out), "n_evaluable": int(done.sum()),
    "n_evaluable_p_gxstage_lt_0.05": int((out.loc[done, "bound_hc3_p_gxstage"] < 0.05).sum()),
    "n_strict_reversal_full": int(out["strict_ci_direction_reversal"].astype(str).eq("True").sum()),
    "n_strict_reversal_after_bound_removal": int(out.loc[done, "bound_strict_reversal"].sum()),
}])
summary.to_csv(OUT_DIR / "boundary_refit_sum.tsv", sep="\t", index=False)
print(summary.to_string(index=False), flush=True)

# --- ATP2B1 diagnostics ---
a = out.index[out["gene_id"] == FOCUS_GENE]
if len(a):
    k = int(a[0])
    pid = out.loc[k, "lead_phenotype_id"]
    g = geno[k].copy()
    g[np.isnan(g)] = np.nanmean(g)
    r, y = raw[pid], pheno[pid]
    diag = []
    for s, lab in ((0, "prenatal"), (1, "postnatal")):
        m = (stage == s) & np.isfinite(r)
        sat = (r[m] >= HIGH - 1e-12).astype(float)
        q_s = residualizer(cov[m])
        row = {"stage": lab, "n": int(m.sum()), "n_saturated": int(sat.sum()),
               "n_distinct_model_phenotype_among_saturated": int(np.unique(np.round(y[m][sat == 1], 8)).size)}
        for gg in (0, 1, 2):
            gm = np.rint(g[m]) == gg
            row[f"geno{gg}_n"] = int(gm.sum())
            row[f"geno{gg}_frac_saturated"] = float(sat[gm].mean()) if gm.any() else np.nan
        for name, yy in (("rawPDUI", r[m]), ("saturated", sat)):
            if np.var(yy) == 0:
                row[f"{name}_beta"] = np.nan
                row[f"{name}_note"] = "no variance"
                continue
            f1 = fit(yy, [g[m]], q_s)
            se = np.sqrt(f1["hc3"][0, 0])
            crit = stats.t.ppf(0.975, f1["df"])
            row.update({f"{name}_beta": f1["beta"][0], f"{name}_ci_low": f1["beta"][0] - crit * se,
                        f"{name}_ci_high": f1["beta"][0] + crit * se,
                        f"{name}_p": 2 * stats.t.sf(abs(f1["beta"][0] / se), f1["df"])})
        diag.append(row)
    d = pd.DataFrame(diag)
    d.to_csv(OUT_DIR / "atp2b1_stage_diagnostics.tsv", sep="\t", index=False)
    print(out.loc[[k]].T.to_string(), flush=True)
    print(d.T.to_string(), flush=True)
print(f"[done] outputs in {OUT_DIR}", flush=True)
