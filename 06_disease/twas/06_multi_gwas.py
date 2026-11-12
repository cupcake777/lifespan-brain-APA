# Build the TWAS master summary: one row per GWAS x stage.
# Uses the by-run association table when present, else the canonical ranked table.
import os
from pathlib import Path

import numpy as np
import pandas as pd

# --- config ---
DATA_DIR = Path(os.environ["APA_ROOT"]) / "data" / "06_disease"
SUMMARY_DIR = DATA_DIR / "by_run"
RANKED_ALL_PATH = DATA_DIR / "twas_ranked_all.tsv.gz"
MANIFEST_PATH = DATA_DIR / "gwas_manifest.tsv"
STAGES = ("prenatal", "postnatal")
EXPECTED_GWAS_COUNT = 19

# --- manifest ---
manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
if len(manifest) != EXPECTED_GWAS_COUNT:
    raise RuntimeError(f"Expected {EXPECTED_GWAS_COUNT} manifest rows, found {len(manifest)}")

# --- per GWAS x stage summary ---
consolidated = pd.read_csv(RANKED_ALL_PATH, sep="\t")
consolidated["gwas_id"] = consolidated["gwas_id"].astype(str)
consolidated["stage"] = consolidated["stage"].astype(str)

rows = []
for gwas_id, meta in manifest.set_index("gwas_id").to_dict("index").items():
    for stage in STAGES:
        path = SUMMARY_DIR / stage / f"{gwas_id}.association_ranked.tsv"
        if path.exists():
            df = pd.read_csv(path, sep="\t")
        else:
            df = consolidated.loc[
                consolidated["gwas_id"].eq(gwas_id) & consolidated["stage"].eq(stage)
            ].copy()
        if df.empty:
            continue

        p = pd.to_numeric(df["pvalue"], errors="coerce")
        z = pd.to_numeric(df["zscore"], errors="coerce")
        valid_p = p.notna() & np.isfinite(p) & (p >= 0) & (p <= 1)
        n_tested = int(valid_p.sum())
        threshold = (
            float(df["bonferroni_threshold"].dropna().astype(float).iloc[0])
            if "bonferroni_threshold" in df and df["bonferroni_threshold"].notna().any()
            else (0.05 / n_tested if n_tested else np.nan)
        )
        ranked = df.loc[valid_p].sort_values(["pvalue", "model_id"], na_position="last")
        top = ranked.iloc[0] if not ranked.empty else pd.Series(dtype=object)
        abs_ranked = df.loc[z.notna()].assign(abs_z=z.abs()).sort_values(
            ["abs_z", "model_id"], ascending=[False, True]
        )
        top_abs_z = float(abs_ranked.iloc[0]["abs_z"]) if not abs_ranked.empty else np.nan
        rows.append(
            {
                "gwas_id": gwas_id,
                "trait_label": meta.get("trait_label", ""),
                "ancestry": meta.get("ancestry", ""),
                "broad_category": meta.get("broad_category", ""),
                "stage": stage,
                "n_tested_features": n_tested,
                "bonferroni_threshold": threshold,
                "n_bonferroni_hits": int((valid_p & (p < threshold)).sum()) if np.isfinite(threshold) else 0,
                "n_positive_z_hits": int((valid_p & (p < threshold) & (z > 0)).sum()) if np.isfinite(threshold) else 0,
                "n_negative_z_hits": int((valid_p & (p < threshold) & (z < 0)).sum()) if np.isfinite(threshold) else 0,
                "top_feature": top.get("stable_feature_id", ""),
                "top_gene": top.get("gene_name", ""),
                "top_z": top.get("zscore", np.nan),
                "top_p": top.get("pvalue", np.nan),
                "top_abs_z": top_abs_z,
            }
        )

out = pd.DataFrame(rows)
expected_rows = len(manifest) * len(STAGES)
if len(out) != expected_rows:
    raise RuntimeError(f"Expected {expected_rows} GWAS-stage rows, found {len(out)}")
out.to_csv(DATA_DIR / "twas_master_summary.tsv", sep="\t", index=False)
print(f"master_rows={len(out)}")
