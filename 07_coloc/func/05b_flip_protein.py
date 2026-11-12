#!/usr/bin/env python3
# Summarize RNA flip candidates with available protein evidence.

import os
from pathlib import Path

import numpy as np
import pandas as pd


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
DATA = ROOT / "data" / "07_coloc"
CROSS_STAGE = DATA / "cross_stage_direction.tsv"
PROTEIN = DATA / "protein_evidence.tsv"
OUTPUT = DATA / "flip_with_protein_evidence.tsv"


# --- merge flip candidates with protein evidence ---

cross = pd.read_csv(CROSS_STAGE, sep="\t")
protein = pd.read_csv(PROTEIN, sep="\t") if PROTEIN.exists() else pd.DataFrame()
flip = cross.loc[cross["cross_stage_direction"].eq("Flip")].copy()
if protein.empty:
    merged = flip
    merged["protein_evidence"] = "MS_not_detected"
    merged["ms_log2_fold_per_allele"] = pd.NA
    merged["ms_pvalue"] = pd.NA
else:
    merged = flip.merge(protein, on=["candidate_id", "gene_name"], how="left")
    merged["protein_evidence"] = merged["protein_evidence"].fillna("MS_not_detected")

# protein evidence class for each flip candidate
evidence = merged["protein_evidence"].astype(str)
merged["protein_flip_concordance"] = np.select(
    [
        evidence.eq("MS_underpowered"),
        evidence.eq("MS_not_detected"),
        evidence.isin(["MS_higher_in_alt", "MS_lower_in_alt"]),
    ],
    ["protein_underpowered", "protein_not_detected", "protein_flip_partial"],
    default="protein_no_signal_either_stage",
)
out = merged.rename(
    columns={
        "slope_prenatal": "rna_flip_prenatal_slope",
        "slope_postnatal": "rna_flip_postnatal_slope",
        "ms_log2_fold_per_allele": "protein_prenatal_effect",
        "ms_pvalue": "protein_prenatal_pvalue",
    }
)
out["protein_postnatal_effect"] = pd.NA
out["protein_postnatal_pvalue"] = pd.NA
keep = [
    "candidate_id",
    "gene_name",
    "rna_flip_prenatal_slope",
    "rna_flip_postnatal_slope",
    "protein_prenatal_effect",
    "protein_prenatal_pvalue",
    "protein_postnatal_effect",
    "protein_postnatal_pvalue",
    "protein_flip_concordance",
]
out[[col for col in keep if col in out.columns]].to_csv(OUTPUT, sep="\t", index=False)
