#!/usr/bin/env python3
# MAGMA .gsa.out helpers: parse VARIABLE rows, condition-run summary, and BH FDR checks.
# Imported by tests/disease/test_magma.py; no command-line entry point.
import os

import numpy as np
import pandas as pd


# --- helpers ---
def parse_gsa_out(path):
    """Return list of dict rows (string tokens kept as read) from a MAGMA .gsa.out file."""
    header = None
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()
            if tok[0] == "VARIABLE":
                header = tok
                continue
            if header is None:
                raise ValueError(f"no VARIABLE header before data in {path}")
            if len(tok) != len(header):
                raise ValueError(f"column count {len(tok)} != header {len(header)} in {path}: {s}")
            rows.append(dict(zip(header, tok)))
    if not rows:
        raise ValueError(f"no result rows in {path}")
    return rows


def bh(p):
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    cummin = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(cummin, 1.0)
    return q


def condition_summary(runs):
    """runs: list of (label, path, keep_variables); keep_variables=[] keeps all rows."""
    out = []
    for label, path, keep in runs:
        for r in parse_gsa_out(path):
            if keep and r["VARIABLE"] not in keep:
                continue
            out.append({
                "analysis": label,
                "variable": r["VARIABLE"],
                "type": r["TYPE"],
                "ngenes": r["NGENES"],
                "beta": r["BETA"],
                "beta_std": r["BETA_STD"],
                "se": r["SE"],
                "p": r["P"],
                "source": os.path.basename(path),
            })
    return pd.DataFrame(out)


def check_fdr(table, group_cols, p_col, within_col, global_col):
    """Recompute within-group and global BH from p and report the largest difference to the table."""
    df = table.copy()
    df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
    df["_within_recalc"] = np.nan
    if group_cols:
        for _, idx in df.groupby(group_cols).groups.items():
            df.loc[idx, "_within_recalc"] = bh(df.loc[idx, p_col].to_numpy())
    else:
        df["_within_recalc"] = bh(df[p_col].to_numpy())
    df["_global_recalc"] = bh(df[p_col].to_numpy())
    dw = (df["_within_recalc"] - pd.to_numeric(df[within_col])).abs().max()
    res = {"n_rows": len(df), "max_abs_diff_within": float(dw)}
    if global_col and global_col in df:
        res["max_abs_diff_global"] = float((df["_global_recalc"] - pd.to_numeric(df[global_col])).abs().max())
    return res
