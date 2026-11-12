#!/usr/bin/env python3
# Write the ranked candidate table from the functional support output.

import os
from pathlib import Path

import pandas as pd


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
DATA = ROOT / "data" / "07_coloc"
INPUT = DATA / "functional_support.tsv.gz"
OUTPUT = DATA / "ranked_candidates.tsv"


# --- rank candidates ---

df = pd.read_csv(INPUT, sep="\t")
keep = [
    "candidate_id",
    "gene_name",
    "stage",
    "tier",
    "framework_tier",
    "pip",
    "mechanism_class",
    "functional_support_category",
    "pas_disrupted",
    "density_supported_modules",
    "mirna_sites_per_100nt",
    "dur_rbp_hits_per_100nt",
    "variant_phastcons",
    "TotalScore",
    "TotalScore_status",
    "top_cell_types",
]
out = (
    df[keep]
    .sort_values(["TotalScore", "pip", "gene_name"], ascending=[False, False, True], kind="stable")
    .reset_index(drop=True)
)
out.index += 1
out.index.name = "rank"
out.to_csv(OUTPUT, sep="\t")
print(f"Written {len(out)} ranked candidates -> {OUTPUT}")
