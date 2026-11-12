#!/usr/bin/env python3
# Cross-stage direction of interaction apaQTL effects (prenatal vs postnatal).
# Writes data/07_coloc/cross_stage_direction.tsv.
import os
from pathlib import Path

import numpy as np
import pandas as pd

# --- config ---

APA_ROOT = Path(os.environ["APA_ROOT"])
INTERACTION_SIG_QTL = APA_ROOT / "data/07_coloc/intx.signif_pairs.txt.gz"
CANDIDATES = APA_ROOT / "data/07_coloc/candidates_prioritized.tsv"
OUTPUT = APA_ROOT / "data/07_coloc/cross_stage_direction.tsv"

# --- inputs ---

inter = pd.read_csv(INTERACTION_SIG_QTL, sep="\t")
candidates = pd.read_csv(CANDIDATES, sep="\t")
# prenatal slope falls back to slope when b_g is absent
if "b_g" not in inter.columns:
    inter["b_g"] = np.nan
    b_g = pd.to_numeric(inter["slope"], errors="coerce")
else:
    b_g = pd.to_numeric(inter["b_g"], errors="coerce")
if "b_gi" not in inter.columns:
    inter["b_gi"] = np.nan
b_gi = pd.to_numeric(inter["b_gi"], errors="coerce")
inter["slope_prenatal"] = b_g
inter["slope_postnatal"] = b_g + b_gi

# --- direction ---

direction = []
for pre, post in zip(inter["slope_prenatal"], inter["slope_postnatal"], strict=False):
    if pd.isna(pre) or pd.isna(post):
        direction.append("not_classified")
    elif pre * post < 0:
        direction.append("Flip")
    elif abs(post) > abs(pre):
        direction.append("Amplifying")
    else:
        direction.append("Attenuating")
inter["cross_stage_direction"] = direction

substatus = []
for pre, post, dir_label in zip(
    inter["slope_prenatal"], inter["slope_postnatal"], inter["cross_stage_direction"], strict=False
):
    if pd.isna(pre) or pd.isna(post):
        substatus.append("none")
    elif dir_label == "Attenuating" and abs(post) < 0.1 * abs(pre):
        substatus.append("silenced_in_postnatal")
    elif dir_label == "Amplifying" and abs(pre) < 0.1 * abs(post):
        substatus.append("silenced_in_prenatal")
    else:
        substatus.append("none")
inter["cross_stage_silenced_substatus"] = substatus

# --- join ---

out = candidates[["candidate_id", "phenotype_id", "variant_id", "gene_name"]].merge(
    inter[
        [
            "phenotype_id",
            "variant_id",
            "b_g",
            "b_gi",
            "pval_gi" if "pval_gi" in inter.columns else "pval_nominal",
            "slope_prenatal",
            "slope_postnatal",
            "cross_stage_direction",
            "cross_stage_silenced_substatus",
        ]
    ].rename(columns={"pval_nominal": "pval_gi"}),
    on=["phenotype_id", "variant_id"],
    how="left",
)
out["interaction_apa_qtl_status"] = out["cross_stage_direction"].notna().map(
    {True: "significant", False: "not_significant"}
)

# --- write ---

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUTPUT, sep="\t", index=False)
