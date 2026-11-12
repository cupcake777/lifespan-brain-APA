#!/usr/bin/env python3
# Summarize prenatal and postnatal SPrediXcan TWAS results for one GWAS.
# Ranked and top-hit tables per stage, cross-stage comparison, QC metrics and top-10 reports.
import argparse
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

# --- config ---
APA = Path(os.environ["APA_ROOT"])
DATA_DIR = APA / "data" / "06_disease"
SPREDIXCAN_DIR = DATA_DIR / "twas"
INPUT_DIR = DATA_DIR / "twas" / "inputs"
QC_BY_RUN_DIR = APA / "data" / "notes" / "frozen_gwas_1029" / "qc" / "twas"
SUMMARY_DIR = DATA_DIR / "by_run"
REPORT_DIR = DATA_DIR / "reports"
QC_METRIC_DIR = DATA_DIR / "qc"
STAGES = ("prenatal", "postnatal")
parser = argparse.ArgumentParser()
parser.add_argument("--gwas-id", default="EAS_SCZ")
GWAS_ID = parser.parse_args().gwas_id
TOP_HIT_MIN_RANK = 50
N_TOP_REPORT = 10

ASSOC_COLUMNS = [
    "stage",
    "rank",
    "model_id",
    "gene_name",
    "gene_id_or_name",
    "transcript_id",
    "original_phenotype_id",
    "3utr_region",
    "strand",
    "chr",
    "group_left",
    "group_right",
    "transcript_left",
    "transcript_right",
    "zscore",
    "effect_size",
    "pvalue",
    "bh_fdr_qvalue",
    "bonferroni_threshold",
    "bonferroni_significant",
    "fdr_lt_0_05",
    "fdr_lt_0_10",
    "nominal_p_lt_0_05",
    "var_g",
    "pred_perf_r2",
    "pred_perf_pval",
    "pred_perf_qval",
    "n_snps_used",
    "n_snps_in_cov",
    "n_snps_in_model",
]

for path in (SUMMARY_DIR / "prenatal", SUMMARY_DIR / "postnatal", REPORT_DIR, QC_METRIC_DIR):
    path.mkdir(parents=True, exist_ok=True)

# --- per stage: ranking, q-values, hit tables ---
stage_data = {}
metrics = []
for stage in STAGES:
    assoc = pd.read_csv(SPREDIXCAN_DIR / f"{stage}_{GWAS_ID}.csv.gz")
    assoc = assoc.rename(columns={"gene": "model_id"})
    assoc["gwas_id"] = GWAS_ID
    assoc["stage"] = stage
    assoc["pvalue"] = pd.to_numeric(assoc["pvalue"], errors="coerce")

    # one annotation row per model_id
    ann = pd.read_csv(
        INPUT_DIR / f"{stage}.annotation.tsv",
        sep="\t",
        usecols=[
            "model_id",
            "original_phenotype_id",
            "transcript_id",
            "gene_id_or_name",
            "3utr_region",
            "strand",
            "chr",
            "group_left",
            "group_right",
            "transcript_left",
            "transcript_right",
        ],
        dtype=str,
    )
    ann = ann.drop_duplicates(subset=["model_id"], keep="first")
    merged = assoc.merge(ann, on="model_id", how="left", validate="m:1")

    # stable feature id: original_phenotype_id, else transcript|gene|3utr|strand
    primary = merged.get("original_phenotype_id", pd.Series("", index=merged.index)).fillna("")
    fallback = merged[["transcript_id", "gene_id_or_name", "3utr_region", "strand"]].fillna("").agg("|".join, axis=1)
    merged["stable_feature_id"] = primary.where(primary.ne(""), fallback)

    valid_p = merged["pvalue"].notna() & np.isfinite(merged["pvalue"]) & (merged["pvalue"] >= 0) & (merged["pvalue"] <= 1)
    n_valid = int(valid_p.sum())
    threshold = 0.05 / n_valid if n_valid else math.nan

    # BH q-values over valid p-values
    p = merged["pvalue"]
    q = pd.Series(np.nan, index=p.index, dtype=float)
    if n_valid:
        p_valid = p.loc[valid_p]
        order = np.argsort(p_valid.to_numpy(), kind="mergesort")
        sorted_p = p_valid.to_numpy()[order]
        n = len(sorted_p)
        raw = sorted_p * n / np.arange(1, n + 1)
        monotone = np.minimum.accumulate(raw[::-1])[::-1]
        monotone = np.clip(monotone, 0, 1)
        q_valid = pd.Series(index=p_valid.index[order], data=monotone)
        q.loc[q_valid.index] = q_valid
    merged["bh_fdr_qvalue"] = q

    merged["bonferroni_threshold"] = threshold
    merged["bonferroni_significant"] = valid_p & (merged["pvalue"] < threshold)
    merged["fdr_lt_0_05"] = merged["bh_fdr_qvalue"] < 0.05
    merged["fdr_lt_0_10"] = merged["bh_fdr_qvalue"] < 0.10
    merged["nominal_p_lt_0_05"] = valid_p & (merged["pvalue"] < 0.05)
    merged = merged.sort_values(["pvalue", "model_id"], na_position="last").reset_index(drop=True)
    merged["rank"] = np.arange(1, len(merged) + 1)

    metrics.append(
        {
            "gwas_id": GWAS_ID,
            "stage": stage,
            "association_rows": int(len(merged)),
            "valid_pvalues": n_valid,
            "bonferroni_threshold": threshold,
            "bonferroni_hits": int(merged["bonferroni_significant"].sum()),
            "fdr_lt_0_05_hits": int(merged["fdr_lt_0_05"].sum()),
            "fdr_lt_0_10_hits": int(merged["fdr_lt_0_10"].sum()),
            "nominal_p_lt_0_05_hits": int(merged["nominal_p_lt_0_05"].sum()),
            "missing_annotation_rows": int(merged["original_phenotype_id"].isna().sum()),
            "min_pvalue": float(p.loc[valid_p].min()) if n_valid else math.nan,
            "top_model_id": str(merged.iloc[0]["model_id"]) if len(merged) else "",
            "top_gene_name": str(merged.iloc[0]["gene_name"]) if len(merged) else "",
        }
    )

    cols = [c for c in ["gwas_id", *ASSOC_COLUMNS] if c in merged.columns] + ["stable_feature_id"]
    merged[cols].to_csv(SUMMARY_DIR / stage / f"{GWAS_ID}.association_ranked.tsv", sep="\t", index=False)
    top_mask = (
        (merged["rank"] <= TOP_HIT_MIN_RANK)
        | merged["bonferroni_significant"]
        | merged["fdr_lt_0_10"]
        | merged["nominal_p_lt_0_05"]
    )
    merged.loc[top_mask, cols].to_csv(SUMMARY_DIR / stage / f"{GWAS_ID}.top_hits.tsv", sep="\t", index=False)
    stage_data[stage] = merged

# --- cross-stage comparison ---
cols = [
    "stable_feature_id",
    "model_id",
    "gene_name",
    "gene_id_or_name",
    "transcript_id",
    "3utr_region",
    "strand",
    "chr",
    "group_left",
    "group_right",
    "rank",
    "zscore",
    "effect_size",
    "pvalue",
    "bh_fdr_qvalue",
    "bonferroni_significant",
    "fdr_lt_0_05",
    "fdr_lt_0_10",
    "nominal_p_lt_0_05",
    "pred_perf_r2",
    "n_snps_used",
]
prepared = {}
for stage, df in stage_data.items():
    sub = df[[c for c in cols if c in df.columns]]
    sub = sub.sort_values(["stable_feature_id", "pvalue", "model_id"], na_position="last")
    sub = sub.drop_duplicates(subset=["stable_feature_id"], keep="first")
    sub = sub.add_prefix(f"{stage}_")
    sub = sub.rename(columns={f"{stage}_stable_feature_id": "stable_feature_id"})
    prepared[stage] = sub

comp = prepared["prenatal"].merge(prepared["postnatal"], on="stable_feature_id", how="outer", indicator=True)
comp["presence"] = comp["_merge"].map(
    {
        "both": "both",
        "left_only": "prenatal_only",
        "right_only": "postnatal_only",
    }
)
comp = comp.drop(columns=["_merge"])
comp["effect_size_delta_post_minus_pre"] = pd.to_numeric(comp["postnatal_effect_size"], errors="coerce") - pd.to_numeric(
    comp["prenatal_effect_size"], errors="coerce"
)
pre_z = pd.to_numeric(comp["prenatal_zscore"], errors="coerce")
post_z = pd.to_numeric(comp["postnatal_zscore"], errors="coerce")
same_sign = pd.Series(np.sign(pre_z) == np.sign(post_z), index=comp.index, dtype="boolean")
same_sign.loc[pre_z.isna() | post_z.isna()] = pd.NA
comp["zscore_same_sign"] = same_sign
pre_p = pd.to_numeric(comp["prenatal_pvalue"], errors="coerce")
post_p = pd.to_numeric(comp["postnatal_pvalue"], errors="coerce")
comp["best_stage_by_pvalue"] = np.where(pre_p <= post_p, "prenatal", "postnatal")
comp.loc[pre_p.isna() & post_p.notna(), "best_stage_by_pvalue"] = "postnatal"
comp.loc[post_p.isna() & pre_p.notna(), "best_stage_by_pvalue"] = "prenatal"
comp.loc[pre_p.isna() & post_p.isna(), "best_stage_by_pvalue"] = ""
comp["best_pvalue"] = np.fmin(pre_p.fillna(np.inf), post_p.fillna(np.inf))
comp.loc[np.isinf(comp["best_pvalue"]), "best_pvalue"] = np.nan

comp = comp.sort_values(["best_pvalue", "stable_feature_id"], na_position="last").reset_index(drop=True)
comp.insert(0, "gwas_id", GWAS_ID)
comp.to_csv(SUMMARY_DIR / f"{GWAS_ID}.stage_comparison.tsv", sep="\t", index=False)

# --- QC metrics and reports ---
pd.DataFrame(metrics).to_csv(QC_METRIC_DIR / f"{GWAS_ID}.qc_metrics.tsv", sep="\t", index=False)

audit_tables = {
    stage: pd.read_csv(QC_BY_RUN_DIR / stage / f"{GWAS_ID}.gwas_alignment_audit.tsv", sep="\t") for stage in STAGES
}
rows = []
for m in metrics:
    audit = audit_tables[m["stage"]].iloc[0].to_dict()
    rows.append(
        {
            "gwas_id": GWAS_ID,
            "stage": m["stage"],
            "association_rows": m["association_rows"],
            "valid_pvalues": m["valid_pvalues"],
            "bonferroni_threshold": m["bonferroni_threshold"],
            "bonferroni_hits": m["bonferroni_hits"],
            "fdr_lt_0_05_hits": m["fdr_lt_0_05_hits"],
            "fdr_lt_0_10_hits": m["fdr_lt_0_10_hits"],
            "nominal_p_lt_0_05_hits": m["nominal_p_lt_0_05_hits"],
            "missing_annotation_rows": m["missing_annotation_rows"],
            "model_variant_count": audit.get("model_variant_count"),
            "matched_model_variant_count": audit.get("matched_model_variant_count"),
            "formatted_gwas_rows": audit.get("formatted_gwas_rows"),
            "exact_effect_allele_matches": audit.get("exact_effect_allele_matches"),
            "flipped_effect_allele_matches": audit.get("flipped_effect_allele_matches"),
            "allele_mismatch_events": audit.get("allele_mismatch_events"),
            "min_N": audit.get("min_N"),
            "max_N": audit.get("max_N"),
        }
    )
presence_counts = comp["presence"].value_counts(dropna=False).to_dict()
summary_rows = pd.DataFrame(rows)
summary_rows["stable_id_rule"] = "original_phenotype_id_else_transcript_gene_region_strand"
summary_rows["stage_comparison_both"] = presence_counts.get("both", 0)
summary_rows["stage_comparison_prenatal_only"] = presence_counts.get("prenatal_only", 0)
summary_rows["stage_comparison_postnatal_only"] = presence_counts.get("postnatal_only", 0)
both = comp[comp["presence"] == "both"]
summary_rows["shared_feature_same_sign_zscore"] = int(both["zscore_same_sign"].fillna(False).sum())
summary_rows["shared_feature_count"] = len(both)
summary_rows.to_csv(REPORT_DIR / f"{GWAS_ID}.summary_report.tsv", sep="\t", index=False)

top_rows = []
for stage in STAGES:
    top = stage_data[stage].head(N_TOP_REPORT)
    for _, row in top.iterrows():
        out = row.to_dict()
        out["gwas_id"] = GWAS_ID
        out["stage"] = stage
        top_rows.append(out)
pd.DataFrame(top_rows).to_csv(REPORT_DIR / f"{GWAS_ID}.top10_report.tsv", sep="\t", index=False)
