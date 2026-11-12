import os
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
DATA = ROOT / "data" / "06_disease"
TWAS_DIR = DATA / "twas"
WORKBOOK = ROOT / "pub/finalization/publication_v28/source_v27/Supplementary_Tables_v27.xlsx"
TARGETS = ["EAS_AD", "EAS_SCZ", "EUR_ANX", "EUR_MDD"]
STAGES = ["prenatal", "postnatal"]
OUTPUT_COLUMNS = [
    "gwas_id",
    "stage",
    "gene_name",
    "transcript_id",
    "original_phenotype_id",
    "zscore",
    "pvalue",
    "bonferroni_threshold",
    "bh_fdr_qvalue",
    "pred_perf_r2",
    "pred_perf_pval",
    "n_snps_in_model",
]

# --- loading v27 ST30 ---

workbook = load_workbook(WORKBOOK, read_only=True, data_only=True)
rows = list(workbook["ST30 TWAS"].values)
old = pd.DataFrame(rows[1:], columns=rows[0])
old = old[OUTPUT_COLUMNS]
old_affected = old[old["gwas_id"].isin(TARGETS)].copy()

# --- frozen-gwas TWAS replacement ---

new_tables = []
for stage in STAGES:
    model_map = pd.read_csv(
        TWAS_DIR / "inputs" / f"{stage}.model_id_map.tsv",
        sep="\t",
        usecols=["model_id", "original_phenotype_id", "transcript_id"],
    )
    for gwas_id in TARGETS:
        table = pd.read_csv(TWAS_DIR / "spredixcan" / stage / f"{gwas_id}.csv.gz")
        table["pvalue"] = pd.to_numeric(table["pvalue"], errors="coerce")
        n_valid = int(table["pvalue"].notna().sum())
        table["bonferroni_threshold"] = 0.05 / n_valid
        valid = table["pvalue"].dropna().sort_values()
        n = len(valid)
        table["bh_fdr_qvalue"] = float("nan")
        if n:
            ranked = valid * n / pd.Series(range(1, n + 1), index=valid.index)
            table.loc[ranked.index, "bh_fdr_qvalue"] = ranked.iloc[::-1].cummin().iloc[::-1].clip(upper=1)
        table = table[table["pvalue"] < table["bonferroni_threshold"]].copy()
        table = table.rename(columns={"gene": "model_id"})
        table = table.merge(model_map, on="model_id", how="left", validate="one_to_one")
        table["gwas_id"] = gwas_id
        table["stage"] = stage
        if table[["original_phenotype_id", "transcript_id"]].isna().any().any():
            raise RuntimeError(f"Missing TWAS annotation for {gwas_id} {stage}")
        new_tables.append(table[OUTPUT_COLUMNS])

new = pd.concat(new_tables, ignore_index=True)
replacement = pd.concat([old[~old["gwas_id"].isin(TARGETS)], new], ignore_index=True)
replacement = replacement.sort_values(["gwas_id", "stage", "pvalue"], kind="stable")

# --- manuscript-facing changes ---

key = ["gwas_id", "stage", "original_phenotype_id"]
old_delta = old_affected.rename(
    columns={
        "gene_name": "old_gene_name",
        "zscore": "old_zscore",
        "pvalue": "old_pvalue",
    }
)[key + ["old_gene_name", "old_zscore", "old_pvalue"]]
new_delta = new.rename(
    columns={
        "gene_name": "new_gene_name",
        "zscore": "new_zscore",
        "pvalue": "new_pvalue",
    }
)[key + ["new_gene_name", "new_zscore", "new_pvalue"]]
delta = old_delta.merge(new_delta, on=key, how="outer", validate="one_to_one", indicator=True)
delta["status"] = delta["_merge"].map(
    {"left_only": "removed", "right_only": "added", "both": "retained"}
)
delta = delta.drop(columns="_merge").sort_values(["gwas_id", "stage", "status", "new_pvalue"])

summary_rows = []
for gwas_id in TARGETS:
    for stage in STAGES:
        subset = delta[(delta["gwas_id"] == gwas_id) & (delta["stage"] == stage)]
        summary_rows.append(
            {
                "gwas_id": gwas_id,
                "stage": stage,
                "old_hits": int((subset["status"] != "added").sum()),
                "new_hits": int((subset["status"] != "removed").sum()),
                "retained": int((subset["status"] == "retained").sum()),
                "added": int((subset["status"] == "added").sum()),
                "removed": int((subset["status"] == "removed").sum()),
            }
        )
summary = pd.DataFrame(summary_rows)

totals = pd.DataFrame(
    [
        {
            "metric": "bonferroni_significant_associations",
            "v27": len(old),
            "frozen_replacement": len(replacement),
        },
        {
            "metric": "unique_prenatal_genes",
            "v27": old.loc[old["stage"] == "prenatal", "gene_name"].nunique(),
            "frozen_replacement": replacement.loc[
                replacement["stage"] == "prenatal", "gene_name"
            ].nunique(),
        },
        {
            "metric": "unique_postnatal_genes",
            "v27": old.loc[old["stage"] == "postnatal", "gene_name"].nunique(),
            "frozen_replacement": replacement.loc[
                replacement["stage"] == "postnatal", "gene_name"
            ].nunique(),
        },
    ]
)

# --- save ---

replacement.to_csv(DATA / "ST30_TWAS_replacement.tsv", sep="\t", index=False)
delta.to_csv(DATA / "TWAS_affected_delta.tsv", sep="\t", index=False)
summary.to_csv(DATA / "TWAS_affected_counts.tsv", sep="\t", index=False)
totals.to_csv(DATA / "TWAS_manuscript_totals.tsv", sep="\t", index=False)
