#!/usr/bin/env python3
# Bridge stage-specific APA-TWAS hits to HC3-core interaction-QTL evidence.
# Annotation and enrichment only; no interaction-component TWAS weights are built.
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# --- config ---
APA_ROOT = Path(os.environ["APA_ROOT"])
DATA = APA_ROOT / "data/06_disease"
TWAS_RANKED_PATH = DATA / "twas_ranked_all.tsv.gz"
TWAS_STAGE_COMPARISON_PATH = DATA / "twas_stage_comparison_all.tsv.gz"
SUMMARY_DIR = DATA / "by_run"
GWAS_MANIFEST_PATH = DATA / "gwas_manifest.tsv"
INTERACTION_CORE_PATH = APA_ROOT / "data/05_gxs/intx_hc3_core_groups.tsv"

ANNOTATED_TWAS_PATH = DATA / "twas_interaction_qtl_annotated.tsv.gz"
ENRICHMENT_PATH = DATA / "interaction_qtl_twas_enrichment_by_stage.tsv"
OVERLAP_HITS_PATH = DATA / "interaction_qtl_twas_overlap_hits.tsv"
STAGE_COMPARISON_PATH = DATA / "interaction_qtl_twas_stage_comparison.tsv"
STAGE_SPECIFICITY_PATH = DATA / "interaction_qtl_twas_stage_specificity.tsv"
REPORT_PATH = DATA / "interaction_qtl_twas_bridge_summary.md"
INPUT_AUDIT_PATH = DATA / "interaction_qtl_bridge_input_audit.tsv"

STAGES = ("prenatal", "postnatal")
MATCH_TYPES = ("interaction_core_gene", "exact_interaction_lead_feature")
HIT_DEFINITIONS = ("bonferroni_significant", "fdr_lt_0_10", "nominal_p_lt_0_05")
TOP_REPORT_ROWS = 20

TWAS_USECOLS = [
    "gwas_id",
    "stage",
    "rank",
    "model_id",
    "gene_name",
    "gene_id_or_name",
    "transcript_id",
    "stable_feature_id",
    "zscore",
    "effect_size",
    "pvalue",
    "bh_fdr_qvalue",
    "bonferroni_threshold",
    "bonferroni_significant",
    "fdr_lt_0_05",
    "fdr_lt_0_10",
    "nominal_p_lt_0_05",
    "pred_perf_r2",
    "n_snps_used",
]

INTERACTION_USECOLS = [
    "group_id",
    "gene_id",
    "stage_group_status",
    "lead_phenotype_id",
    "lead_variant_id",
    "group_q_eigenmt",
    "lead_pval_gi",
    "lead_b_gi",
    "beta_pre",
    "beta_gxstage",
    "beta_post",
    "hc3_pval_gxstage",
    "point_estimate_class",
    "strict_ci_direction_reversal",
    "hc3_robust_core",
    "genotype_support_status",
    "scaled_interaction_sign_consistent",
]

STAGE_COMPARISON_USECOLS = [
    "gwas_id",
    "stable_feature_id",
    "prenatal_gene_name",
    "prenatal_gene_id_or_name",
    "prenatal_transcript_id",
    "prenatal_rank",
    "prenatal_zscore",
    "prenatal_effect_size",
    "prenatal_pvalue",
    "prenatal_bonferroni_significant",
    "prenatal_fdr_lt_0_10",
    "prenatal_nominal_p_lt_0_05",
    "postnatal_gene_name",
    "postnatal_gene_id_or_name",
    "postnatal_transcript_id",
    "postnatal_rank",
    "postnatal_zscore",
    "postnatal_effect_size",
    "postnatal_pvalue",
    "postnatal_bonferroni_significant",
    "postnatal_fdr_lt_0_10",
    "postnatal_nominal_p_lt_0_05",
    "presence",
    "effect_size_delta_post_minus_pre",
    "zscore_same_sign",
    "best_stage_by_pvalue",
    "best_pvalue",
]

# --- helpers ---

def read_bool(series):
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def clean_key(series):
    return series.fillna("").astype(str).str.strip()


def neg_log10(value):
    p = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if not np.isfinite(p):
        return np.nan
    if p <= 0:
        return 300.0
    return float(-math.log10(p))


def pick_best_interaction_rows(core, key_col):
    ranked = core.sort_values(
        ["hc3_pval_gxstage", "group_q_eigenmt", "group_id"],
        na_position="last",
    )
    return ranked.drop_duplicates(key_col, keep="first")


def rename_interaction_columns(df, prefix):
    rename = {
        col: f"{prefix}_{col}"
        for col in df.columns
        if col not in {"stable_feature_id", "gene_key"}
    }
    return df.rename(columns=rename)


def read_twas_ranked(path):
    ranked = pd.read_csv(path, sep="\t", usecols=lambda col: col in TWAS_USECOLS)
    return ranked.reindex(columns=TWAS_USECOLS)


def read_stage_comparison(path):
    stage = pd.read_csv(path, sep="\t", usecols=lambda col: col in STAGE_COMPARISON_USECOLS)
    return stage.reindex(columns=STAGE_COMPARISON_USECOLS)


def coalesce_columns(df, exact_col, gene_col):
    exact = df.get(exact_col)
    gene = df.get(gene_col)
    if exact is None:
        return gene
    if gene is None:
        return exact
    return exact.where(exact.notna() & exact.astype(str).ne(""), gene)


def md_row(values):
    cells = [str(v).replace("|", "\\|").replace("\n", " ") for v in values]
    return "| " + " | ".join(cells) + " |"


def top_markdown_rows(df, columns, n):
    if df.empty:
        return "none"
    show = df[columns].head(n).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: f"{x:.3g}" if np.isfinite(x) else "")
    show = show.fillna("").astype(str)
    header = md_row(show.columns)
    separator = md_row(["---"] * len(show.columns))
    body = [md_row(row) for row in show.to_numpy(dtype=str)]
    return "\n".join([header, separator, *body])

# --- load inputs ---
manifest = pd.read_csv(GWAS_MANIFEST_PATH, sep="\t")

core = pd.read_csv(INTERACTION_CORE_PATH, sep="\t", usecols=INTERACTION_USECOLS)
core["hc3_robust_core"] = read_bool(core["hc3_robust_core"])
core["strict_ci_direction_reversal"] = read_bool(core["strict_ci_direction_reversal"])
gene_id = clean_key(core["gene_id"])
core["gene_key"] = gene_id.where(gene_id.ne(""), clean_key(core["group_id"]))
core["lead_phenotype_id"] = clean_key(core["lead_phenotype_id"])
core = core.loc[core["hc3_robust_core"]].copy()

feature_best = pick_best_interaction_rows(core, "lead_phenotype_id")
feature_best = feature_best.drop(columns=["gene_key"])
feature_best["stable_feature_id"] = feature_best["lead_phenotype_id"]
feature_best = rename_interaction_columns(feature_best, "exact")

gene_best = pick_best_interaction_rows(core, "gene_key")
gene_best = rename_interaction_columns(gene_best, "gene")

# --- load APA-TWAS ranked results ---
twas = read_twas_ranked(TWAS_RANKED_PATH)
supplemental_twas = []
present_pairs = set(zip(twas["gwas_id"].astype(str), twas["stage"].astype(str)))
for gwas_id in manifest["gwas_id"].dropna().astype(str).tolist():
    for stage in STAGES:
        if (gwas_id, stage) in present_pairs:
            continue
        path = SUMMARY_DIR / stage / f"{gwas_id}.association_ranked.tsv"
        if not path.exists():
            continue
        extra = read_twas_ranked(path)
        extra["gwas_id"] = gwas_id
        extra["stage"] = stage
        twas = pd.concat([twas, extra], ignore_index=True)
        supplemental_twas.append(str(path))
        present_pairs.add((gwas_id, stage))

gene_name = clean_key(twas["gene_name"])
twas["gene_key"] = gene_name.where(gene_name.ne(""), clean_key(twas["gene_id_or_name"]))
for col in HIT_DEFINITIONS + ("fdr_lt_0_05",):
    twas[col] = read_bool(twas[col])
for col in ["rank", "zscore", "effect_size", "pvalue", "bh_fdr_qvalue", "bonferroni_threshold", "pred_perf_r2", "n_snps_used"]:
    twas[col] = pd.to_numeric(twas[col], errors="coerce")

# --- annotate TWAS with interaction-QTL evidence ---
annotated = twas.merge(feature_best, on="stable_feature_id", how="left", validate="m:1")
annotated = annotated.merge(gene_best, on="gene_key", how="left", validate="m:1")
annotated["exact_interaction_lead_feature"] = annotated["exact_group_id"].notna()
annotated["interaction_core_gene"] = annotated["gene_group_id"].notna()
annotated["any_interaction_core"] = (
    annotated["exact_interaction_lead_feature"] | annotated["interaction_core_gene"]
)
annotated["interaction_match_level"] = np.select(
    [annotated["exact_interaction_lead_feature"], annotated["interaction_core_gene"]],
    ["exact_lead_feature", "core_gene"],
    default="none",
)
bridge_cols = [
    "group_id",
    "gene_id",
    "stage_group_status",
    "lead_phenotype_id",
    "lead_variant_id",
    "group_q_eigenmt",
    "lead_pval_gi",
    "lead_b_gi",
    "beta_pre",
    "beta_gxstage",
    "beta_post",
    "hc3_pval_gxstage",
    "point_estimate_class",
    "strict_ci_direction_reversal",
    "genotype_support_status",
    "scaled_interaction_sign_consistent",
]
for col in bridge_cols:
    annotated[f"bridge_{col}"] = coalesce_columns(annotated, f"exact_{col}", f"gene_{col}")
annotated["minus_log10_twas_p"] = annotated["pvalue"].map(neg_log10)
annotated["minus_log10_interaction_hc3_p"] = annotated["bridge_hc3_pval_gxstage"].map(neg_log10)
annotated["bridge_priority_score"] = (
    annotated["minus_log10_twas_p"] + annotated["minus_log10_interaction_hc3_p"]
)

# --- enrichment of interaction-QTL markers among TWAS hits ---
meta_dict = manifest.set_index("gwas_id").to_dict("index")
rows = []
for (gwas_id, stage), sub in annotated.groupby(["gwas_id", "stage"], sort=True):
    for marker_col in MATCH_TYPES:
        marker = read_bool(sub[marker_col])
        for hit_col in HIT_DEFINITIONS:
            hit = read_bool(sub[hit_col])
            a = int((marker & hit).sum())
            b = int((marker & ~hit).sum())
            c = int((~marker & hit).sum())
            d = int((~marker & ~hit).sum())
            try:
                oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
            except ValueError:
                oddsratio, pvalue = np.nan, np.nan
            marker_fraction = a / (a + b) if (a + b) else np.nan
            background_fraction = c / (c + d) if (c + d) else np.nan
            rows.append({
                "n_total": int(len(sub)),
                "n_marker": int(marker.sum()),
                "n_non_marker": int((~marker).sum()),
                "n_hits": int(hit.sum()),
                "n_marker_hits": a,
                "n_background_hits": c,
                "marker_hit_fraction": marker_fraction,
                "background_hit_fraction": background_fraction,
                "oddsratio": float(oddsratio) if np.isfinite(oddsratio) else oddsratio,
                "fisher_greater_pvalue": float(pvalue) if np.isfinite(pvalue) else np.nan,
                "gwas_id": gwas_id,
                "stage": stage,
                "marker_type": marker_col,
                "hit_definition": hit_col,
                "trait_label": meta_dict.get(gwas_id, {}).get("trait_label", ""),
                "ancestry": meta_dict.get(gwas_id, {}).get("ancestry", ""),
                "broad_category": meta_dict.get(gwas_id, {}).get("broad_category", ""),
            })
enrichment = pd.DataFrame(rows)
if not enrichment.empty:
    enrichment["fisher_bh_qvalue"] = np.nan
    for _, grp in enrichment.groupby(["marker_type", "hit_definition"]):
        p = pd.to_numeric(grp["fisher_greater_pvalue"], errors="coerce")
        valid_p = p.loc[p.notna() & np.isfinite(p) & (p >= 0) & (p <= 1)]
        if valid_p.empty:
            continue
        order = np.argsort(valid_p.to_numpy(), kind="mergesort")
        sorted_p = valid_p.to_numpy()[order]
        n = len(sorted_p)
        raw = sorted_p * n / np.arange(1, n + 1)
        adjusted = np.minimum.accumulate(raw[::-1])[::-1]
        enrichment.loc[valid_p.index[order], "fisher_bh_qvalue"] = np.clip(adjusted, 0, 1)
    enrichment = enrichment.sort_values(
        ["hit_definition", "marker_type", "fisher_greater_pvalue", "gwas_id", "stage"],
        na_position="last",
    ).reset_index(drop=True)

# --- overlap hits with interaction-QTL support ---
meta_cols = ["gwas_id", "trait_label", "ancestry", "broad_category"]
meta = manifest[meta_cols].drop_duplicates("gwas_id")
overlap = annotated.loc[
    annotated["any_interaction_core"]
    & (
        annotated["bonferroni_significant"]
        | annotated["fdr_lt_0_10"]
        | annotated["nominal_p_lt_0_05"]
    )
].copy()
overlap["twas_hit_tier"] = np.select(
    [
        overlap["bonferroni_significant"],
        overlap["fdr_lt_0_05"],
        overlap["fdr_lt_0_10"],
        overlap["nominal_p_lt_0_05"],
    ],
    ["bonferroni", "fdr_0.05", "fdr_0.10", "nominal_0.05"],
    default="not_nominal",
)
overlap = overlap.merge(meta, on="gwas_id", how="left", validate="m:1")
overlap_keep = [
    "gwas_id",
    "trait_label",
    "ancestry",
    "broad_category",
    "stage",
    "twas_hit_tier",
    "interaction_match_level",
    "rank",
    "model_id",
    "gene_name",
    "gene_id_or_name",
    "transcript_id",
    "stable_feature_id",
    "zscore",
    "effect_size",
    "pvalue",
    "bh_fdr_qvalue",
    "bonferroni_significant",
    "fdr_lt_0_10",
    "nominal_p_lt_0_05",
    "pred_perf_r2",
    "n_snps_used",
    "bridge_group_id",
    "bridge_gene_id",
    "bridge_stage_group_status",
    "bridge_point_estimate_class",
    "bridge_lead_phenotype_id",
    "bridge_lead_variant_id",
    "bridge_beta_pre",
    "bridge_beta_gxstage",
    "bridge_beta_post",
    "bridge_hc3_pval_gxstage",
    "bridge_strict_ci_direction_reversal",
    "bridge_genotype_support_status",
    "bridge_scaled_interaction_sign_consistent",
    "minus_log10_twas_p",
    "minus_log10_interaction_hc3_p",
    "bridge_priority_score",
]
overlap_hits = overlap[overlap_keep].sort_values(
    ["twas_hit_tier", "pvalue", "bridge_priority_score", "gwas_id", "stage"],
    ascending=[True, True, False, True, True],
    na_position="last",
).reset_index(drop=True)

# --- stage comparison of APA-TWAS effects ---
stage = read_stage_comparison(TWAS_STAGE_COMPARISON_PATH)
supplemental_stage = []
present_gwas = set(stage["gwas_id"].dropna().astype(str))
for gwas_id in manifest["gwas_id"].dropna().astype(str).tolist():
    if gwas_id in present_gwas:
        continue
    path = SUMMARY_DIR / f"{gwas_id}.stage_comparison.tsv"
    if not path.exists():
        continue
    extra = read_stage_comparison(path)
    extra["gwas_id"] = gwas_id
    stage = pd.concat([stage, extra], ignore_index=True)
    supplemental_stage.append(str(path))
    present_gwas.add(gwas_id)

for col in [
    "prenatal_bonferroni_significant",
    "prenatal_fdr_lt_0_10",
    "prenatal_nominal_p_lt_0_05",
    "postnatal_bonferroni_significant",
    "postnatal_fdr_lt_0_10",
    "postnatal_nominal_p_lt_0_05",
]:
    stage[col] = read_bool(stage[col])
for col in [
    "prenatal_rank",
    "prenatal_zscore",
    "prenatal_effect_size",
    "prenatal_pvalue",
    "postnatal_rank",
    "postnatal_zscore",
    "postnatal_effect_size",
    "postnatal_pvalue",
    "effect_size_delta_post_minus_pre",
    "best_pvalue",
]:
    stage[col] = pd.to_numeric(stage[col], errors="coerce")
prenatal_name = clean_key(stage["prenatal_gene_name"])
stage["gene_key"] = prenatal_name.where(prenatal_name.ne(""), clean_key(stage["postnatal_gene_name"]))
stage["gene_key"] = stage["gene_key"].where(
    stage["gene_key"].ne(""),
    clean_key(stage["prenatal_gene_id_or_name"]),
)
stage["gene_key"] = stage["gene_key"].where(
    stage["gene_key"].ne(""),
    clean_key(stage["postnatal_gene_id_or_name"]),
)
stage["delta_abs_z_post_minus_pre"] = stage["postnatal_zscore"].abs() - stage["prenatal_zscore"].abs()

stage_annotated = stage.merge(feature_best, on="stable_feature_id", how="left", validate="m:1")
stage_annotated = stage_annotated.merge(gene_best, on="gene_key", how="left", validate="m:1")
stage_annotated["exact_interaction_lead_feature"] = stage_annotated["exact_group_id"].notna()
stage_annotated["interaction_core_gene"] = stage_annotated["gene_group_id"].notna()
stage_annotated["interaction_match_level"] = np.select(
    [stage_annotated["exact_interaction_lead_feature"], stage_annotated["interaction_core_gene"]],
    ["exact_lead_feature", "core_gene"],
    default="none",
)
for col in [
    "group_id",
    "gene_id",
    "stage_group_status",
    "lead_phenotype_id",
    "beta_pre",
    "beta_gxstage",
    "beta_post",
    "hc3_pval_gxstage",
    "point_estimate_class",
    "strict_ci_direction_reversal",
]:
    stage_annotated[f"bridge_{col}"] = coalesce_columns(stage_annotated, f"exact_{col}", f"gene_{col}")
stage_keep = [
    "gwas_id",
    "interaction_match_level",
    "stable_feature_id",
    "gene_key",
    "prenatal_gene_name",
    "postnatal_gene_name",
    "presence",
    "prenatal_rank",
    "prenatal_zscore",
    "prenatal_effect_size",
    "prenatal_pvalue",
    "prenatal_bonferroni_significant",
    "prenatal_fdr_lt_0_10",
    "prenatal_nominal_p_lt_0_05",
    "postnatal_rank",
    "postnatal_zscore",
    "postnatal_effect_size",
    "postnatal_pvalue",
    "postnatal_bonferroni_significant",
    "postnatal_fdr_lt_0_10",
    "postnatal_nominal_p_lt_0_05",
    "effect_size_delta_post_minus_pre",
    "delta_abs_z_post_minus_pre",
    "zscore_same_sign",
    "best_stage_by_pvalue",
    "best_pvalue",
    "bridge_group_id",
    "bridge_gene_id",
    "bridge_stage_group_status",
    "bridge_point_estimate_class",
    "bridge_lead_phenotype_id",
    "bridge_beta_pre",
    "bridge_beta_gxstage",
    "bridge_beta_post",
    "bridge_hc3_pval_gxstage",
    "bridge_strict_ci_direction_reversal",
]
stage_overlap = stage_annotated.loc[
    stage_annotated["interaction_match_level"].ne("none"), stage_keep
].copy()
stage_overlap = stage_overlap.sort_values(
    ["best_pvalue", "gwas_id", "interaction_match_level"],
    na_position="last",
).reset_index(drop=True)

# --- stage specificity among interaction-overlap features ---
stage_specificity = pd.DataFrame()
if not stage_overlap.empty:
    rows = []
    for (gwas_id, match_level), sub in stage_overlap.groupby(["gwas_id", "interaction_match_level"], sort=True):
        nominal = pd.to_numeric(sub["best_pvalue"], errors="coerce") < 0.05
        bonf = read_bool(sub["prenatal_bonferroni_significant"]) | read_bool(sub["postnatal_bonferroni_significant"])
        rows.append({
            "gwas_id": gwas_id,
            "interaction_match_level": match_level,
            "n_overlap_features": int(len(sub)),
            "n_best_nominal_p_lt_0_05": int(nominal.sum()),
            "n_any_stage_bonferroni": int(bonf.sum()),
            "n_best_prenatal": int((sub["best_stage_by_pvalue"].eq("prenatal")).sum()),
            "n_best_postnatal": int((sub["best_stage_by_pvalue"].eq("postnatal")).sum()),
            "median_delta_abs_z_post_minus_pre": float(sub["delta_abs_z_post_minus_pre"].median(skipna=True)),
            "n_postnatal_larger_abs_z": int((sub["delta_abs_z_post_minus_pre"] > 0).sum()),
            "n_prenatal_larger_abs_z": int((sub["delta_abs_z_post_minus_pre"] < 0).sum()),
            "n_same_z_sign": int(read_bool(sub["zscore_same_sign"]).sum()),
        })
    stage_specificity = pd.DataFrame(rows).merge(meta, on="gwas_id", how="left", validate="m:1")
    stage_specificity = stage_specificity.sort_values(
        ["n_any_stage_bonferroni", "n_best_nominal_p_lt_0_05", "gwas_id", "interaction_match_level"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)

# --- input audit ---
ranked_gwas = sorted(twas["gwas_id"].dropna().astype(str).unique())
manifest_gwas = sorted(manifest["gwas_id"].dropna().astype(str).unique())
stage_gwas = sorted(stage["gwas_id"].dropna().astype(str).unique())
audit_rows = [
    ("twas_ranked_path", str(TWAS_RANKED_PATH)),
    ("twas_supplemental_ranked_files", ",".join(supplemental_twas)),
    ("twas_ranked_rows", len(twas)),
    ("twas_ranked_gwas_stage_pairs", twas[["gwas_id", "stage"]].drop_duplicates().shape[0]),
    ("twas_ranked_gwas_ids", ",".join(ranked_gwas)),
    ("gwas_manifest_path", str(GWAS_MANIFEST_PATH)),
    ("gwas_manifest_rows", len(manifest)),
    ("gwas_manifest_ids_not_in_twas_ranked", ",".join(sorted(set(manifest_gwas) - set(ranked_gwas)))),
    ("interaction_core_path", str(INTERACTION_CORE_PATH)),
    ("interaction_core_rows", len(core)),
    ("interaction_core_genes", core["gene_key"].nunique()),
    ("interaction_core_lead_features", core["lead_phenotype_id"].nunique()),
    ("annotated_exact_lead_feature_rows", int(annotated["exact_interaction_lead_feature"].sum())),
    ("annotated_core_gene_rows", int(annotated["interaction_core_gene"].sum())),
    ("stage_comparison_path", str(TWAS_STAGE_COMPARISON_PATH)),
    ("stage_comparison_supplemental_files", ",".join(supplemental_stage)),
    ("stage_comparison_gwas_ids", ",".join(stage_gwas)),
    ("stage_comparison_manifest_ids_not_present", ",".join(sorted(set(manifest_gwas) - set(stage_gwas)))),
    ("stage_comparison_overlap_rows", len(stage_overlap)),
]
audit = pd.DataFrame(audit_rows, columns=["metric", "value"])
audit.to_csv(INPUT_AUDIT_PATH, sep="\t", index=False)

# --- write tables ---
annotated_out_cols = [
    *TWAS_USECOLS,
    "interaction_match_level",
    "exact_interaction_lead_feature",
    "interaction_core_gene",
    "bridge_group_id",
    "bridge_gene_id",
    "bridge_stage_group_status",
    "bridge_point_estimate_class",
    "bridge_lead_phenotype_id",
    "bridge_lead_variant_id",
    "bridge_beta_pre",
    "bridge_beta_gxstage",
    "bridge_beta_post",
    "bridge_hc3_pval_gxstage",
    "bridge_strict_ci_direction_reversal",
    "bridge_genotype_support_status",
    "bridge_scaled_interaction_sign_consistent",
    "minus_log10_twas_p",
    "minus_log10_interaction_hc3_p",
    "bridge_priority_score",
]
annotated[annotated_out_cols].to_csv(ANNOTATED_TWAS_PATH, sep="\t", index=False, compression="gzip")
enrichment.to_csv(ENRICHMENT_PATH, sep="\t", index=False)
overlap_hits.to_csv(OVERLAP_HITS_PATH, sep="\t", index=False)
stage_overlap.to_csv(STAGE_COMPARISON_PATH, sep="\t", index=False)
stage_specificity.to_csv(STAGE_SPECIFICITY_PATH, sep="\t", index=False)

# --- report ---
total_pairs = annotated[["gwas_id", "stage"]].drop_duplicates().shape[0]
total_bonf = int(read_bool(annotated["bonferroni_significant"]).sum())
core_bonf = int(
    (read_bool(annotated["bonferroni_significant"]) & read_bool(annotated["interaction_core_gene"])).sum()
)
exact_bonf = int(
    (read_bool(annotated["bonferroni_significant"]) & read_bool(annotated["exact_interaction_lead_feature"])).sum()
)
enrich_top = enrichment.loc[
    enrichment["hit_definition"].eq("bonferroni_significant")
    & (enrichment["n_marker_hits"] > 0)
].sort_values(["fisher_greater_pvalue", "marker_type"], na_position="last")

report = [
    "# Interaction-QTL bridge for APA-TWAS",
    "",
    "This analysis annotates current stage-specific APA-TWAS results with HC3-core developmental interaction-QTL evidence from `data/05_gxs/intx_hc3_core_groups.tsv`.",
    "It is not an interaction-component TWAS model and should not be interpreted as new interaction TWAS weights.",
    "",
    "## Input scope",
    "",
    f"- TWAS rows: {len(annotated):,}",
    f"- GWAS x stage pairs in TWAS table: {total_pairs}",
    f"- Bonferroni-significant TWAS rows: {total_bonf}",
    f"- Bonferroni rows in HC3-core interaction genes: {core_bonf}",
    f"- Bonferroni rows matching exact interaction lead features: {exact_bonf}",
    "",
    "## Strongest Bonferroni enrichment rows",
    "",
    top_markdown_rows(
        enrich_top,
        [
            "gwas_id",
            "stage",
            "marker_type",
            "n_marker",
            "n_hits",
            "n_marker_hits",
            "oddsratio",
            "fisher_greater_pvalue",
            "fisher_bh_qvalue",
        ],
        TOP_REPORT_ROWS,
    ),
    "",
    "## Overlapping TWAS signals with interaction-QTL support",
    "",
    top_markdown_rows(
        overlap_hits,
        [
            "gwas_id",
            "stage",
            "twas_hit_tier",
            "interaction_match_level",
            "gene_name",
            "transcript_id",
            "zscore",
            "pvalue",
            "bridge_stage_group_status",
            "bridge_point_estimate_class",
            "bridge_hc3_pval_gxstage",
            "bridge_priority_score",
        ],
        TOP_REPORT_ROWS,
    ),
    "",
    "## Stage specificity among interaction-overlap features",
    "",
    top_markdown_rows(
        stage_specificity,
        [
            "gwas_id",
            "interaction_match_level",
            "n_overlap_features",
            "n_best_nominal_p_lt_0_05",
            "n_any_stage_bonferroni",
            "n_best_prenatal",
            "n_best_postnatal",
            "median_delta_abs_z_post_minus_pre",
        ],
        TOP_REPORT_ROWS,
    ),
    "",
    "## Input audit",
    "",
    top_markdown_rows(audit, ["metric", "value"], len(audit)),
    "",
]
REPORT_PATH.write_text("\n".join(report), encoding="utf-8")
