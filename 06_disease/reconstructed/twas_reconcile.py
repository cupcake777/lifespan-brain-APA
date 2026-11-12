#!/usr/bin/env python3
# Reconciles twas_master_summary.tsv Bonferroni hit counts with twas_ranked_all.tsv.gz per gwas_id x stage.
# Writes data/06_disease/twas_master_vs_ranked_reconcile.tsv.
import os

import pandas as pd

APA = os.environ["APA_ROOT"]
DATA = os.path.join(APA, "data/06_disease")

# --- master summary and ranked table ---
tm = pd.read_csv(os.path.join(DATA, "twas_master_summary.tsv"), sep="\t")
tm["n_bonferroni_hits"] = pd.to_numeric(tm["n_bonferroni_hits"]).astype(int)
tr = pd.read_csv(os.path.join(DATA, "twas_ranked_all.tsv.gz"), sep="\t",
                 usecols=["gwas_id", "stage", "bonferroni_significant"])
flag_written = tr["bonferroni_significant"].astype(str).str.lower().isin(["true", "1"])

# --- per-group counts and differences ---
ranked_grp = tr[flag_written].groupby(["gwas_id", "stage"]).size().rename("ranked_bonf").reset_index()
rg = tr.groupby(["gwas_id", "stage"]).size().rename("ranked_rows").reset_index()
rec = tm[["gwas_id", "stage", "n_tested_features", "bonferroni_threshold", "n_bonferroni_hits"]].merge(
    ranked_grp, on=["gwas_id", "stage"], how="outer").merge(rg, on=["gwas_id", "stage"], how="left")
rec["ranked_bonf"] = rec["ranked_bonf"].fillna(0).astype(int)
rec["ranked_rows"] = rec["ranked_rows"].fillna(0).astype(int)
rec["diff_master_minus_ranked"] = rec["n_bonferroni_hits"].fillna(0).astype(int) - rec["ranked_bonf"]
diff = rec[rec["diff_master_minus_ranked"] != 0]

print(f"[master] sum n_bonferroni_hits = {int(tm.n_bonferroni_hits.sum())}")
print(f"[ranked] sum Bonferroni rows = {int(flag_written.sum())}")
missing = sorted(set(map(tuple, tm[["gwas_id", "stage"]].values)) -
                 set(map(tuple, tr[["gwas_id", "stage"]].drop_duplicates().values)))
print(f"[groups in master but absent from ranked] {missing}")
print("[group count differences]")
print(diff[["gwas_id", "stage", "n_bonferroni_hits", "ranked_bonf", "diff_master_minus_ranked"]].to_string(index=False))

rec.to_csv(os.path.join(DATA, "twas_master_vs_ranked_reconcile.tsv"), sep="\t", index=False,
           float_format="%.6g")
