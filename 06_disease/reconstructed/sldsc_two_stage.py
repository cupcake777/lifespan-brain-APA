#!/usr/bin/env python3
# Two-stage BH FDR for the binary apaQTL S-LDSC annotations (prenatal, postnatal) within each
# Population x Trait. Imported by fig5a_plot_data.py and tests/disease/test_sldsc_2stage.py.
# Run as a script to write data/06_disease/sldsc_binary_fdr_two_stage.tsv.
import os

import numpy as np
import pandas as pd

BINARY = {"prenatal_apa_qtl": "Prenatal", "postnatal_apa_qtl": "Postnatal"}


# --- helpers ---
def bh(p):
    """Benjamini-Hochberg q-values, same convention as statsmodels fdr_bh (monotone, capped at 1)."""
    p = np.asarray(p, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    cummin = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(cummin, 1.0)
    return q


def load_processed(path):
    df = pd.read_csv(path, sep="\t")
    for c in ["Enrichment", "Enrichment_p", "FDR_p"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def two_stage_fdr(df):
    """Two-stage BH within each Population x Trait over prenatal and postnatal binary annotations."""
    sub = df[df["Annotation"].isin(BINARY)].copy()
    sub["Stage"] = sub["Annotation"].map(BINARY)
    sub["FDR_two_stage"] = np.nan
    for _, idx in sub.groupby(["Population", "Trait"]).groups.items():
        sub.loc[idx, "FDR_two_stage"] = bh(sub.loc[idx, "Enrichment_p"].to_numpy())
    sub["n_family"] = sub.groupby(["Population", "Trait"])["Stage"].transform("size")
    return sub[["Population", "Trait", "Stage", "Annotation", "Enrichment", "Enrichment_p",
                "FDR_two_stage", "n_family"]].reset_index(drop=True)


def check_global76(df):
    """Recompute the one-BH-over-all-rows FDR and compare with the processed FDR_p column."""
    ok = df["Enrichment_p"].notna() & df["FDR_p"].notna()
    q = np.full(len(df), np.nan)
    q[ok.to_numpy()] = bh(df.loc[ok, "Enrichment_p"].to_numpy())
    obs = df["FDR_p"].to_numpy()
    rel = np.abs(q - obs) / np.maximum(np.abs(obs), 1e-300)
    return {
        "n_rows": int(len(df)),
        "n_compared": int(ok.sum()),
        "max_abs_diff": float(np.nanmax(np.abs(q - obs))),
        "max_rel_diff": float(np.nanmax(rel)),
    }


if __name__ == "__main__":
    APA = os.environ["APA_ROOT"]
    DATA = os.path.join(APA, "data/06_disease")
    df = load_processed(os.path.join(DATA, "sldsc_enrich_processed.tsv"))
    two = two_stage_fdr(df)
    two.to_csv(os.path.join(DATA, "sldsc_binary_fdr_two_stage.tsv"), sep="\t", index=False,
               float_format="%.6g")
