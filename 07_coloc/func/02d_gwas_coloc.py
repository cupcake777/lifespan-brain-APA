#!/usr/bin/env python3
# Join apaQTL-GWAS coloc PP4 onto prioritized candidates; GWAS_coloc = PP4 >= 0.75.
# Writes data/07_coloc/gwas_coloc.tsv.
import os
from pathlib import Path

import pandas as pd

# --- config ---

APA_ROOT = Path(os.environ["APA_ROOT"])
CANDIDATES = APA_ROOT / "data/07_coloc/candidates_prioritized.tsv"
GWAS_COLOC = APA_ROOT / "data/07_coloc/coloc_run/all_colocalization.csv"
OUTPUT = APA_ROOT / "data/07_coloc/gwas_coloc.tsv"
PP4_THRESHOLD = 0.75

# --- inputs ---

candidates = pd.read_csv(CANDIDATES, sep="\t")
coloc = pd.read_csv(GWAS_COLOC)
coloc = coloc.rename(
    columns={
        "QTL_ID": "stage",
        "Phenotype": "phenotype_id",
        "PP4": "gwas_coloc_pp4",
        "GWAS_ID": "gwas_trait",
    }
)
coloc["gwas_coloc_pp4"] = pd.to_numeric(coloc["gwas_coloc_pp4"], errors="coerce")

# best PP4 per stage and phenotype
best = coloc.sort_values("gwas_coloc_pp4", ascending=False, kind="stable").drop_duplicates(
    ["stage", "phenotype_id"], keep="first"
)

# --- join ---

out = candidates[["candidate_id", "phenotype_id", "stage"]].merge(
    best[["stage", "phenotype_id", "gwas_coloc_pp4", "gwas_trait"]],
    on=["stage", "phenotype_id"],
    how="left",
)
out["GWAS_coloc"] = out["gwas_coloc_pp4"].ge(PP4_THRESHOLD).fillna(False)
out = out[["candidate_id", "GWAS_coloc", "gwas_coloc_pp4", "gwas_trait"]]

# --- write ---

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUTPUT, sep="\t", index=False)
