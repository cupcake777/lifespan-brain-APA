# Audit of the eigenMT group result: Meff vs tested pairs, multi-site aggregation,
# and groups significant only at group level. Writes interaction_eigenmt_* to data/05_gxs.
from pathlib import Path
import os

import numpy as np
import pandas as pd

# --- paths ---
APA_ROOT = Path(os.environ["APA_ROOT"])
D = APA_ROOT / "data" / "05_gxs"
GROUP_FILE = D / "interaction_group_summary.tsv"
EIGENMT_TOP_FILE = D / "interaction_group_eigenmt_top.tsv"
SIGNIFICANT_FILE = D / "interaction_significant_groups.tsv"
PREVIOUS_SIGNIFICANT_FILE = D / "previous115_interaction_significant_groups.tsv"
SITE_HIT_FILE = D / "sig.interaction.result"

# --- load ---
groups = pd.read_csv(GROUP_FILE, sep="\t")
eigenmt_top = pd.read_csv(EIGENMT_TOP_FILE, sep="\t")
significant = pd.read_csv(SIGNIFICANT_FILE, sep="\t")
previous = pd.read_csv(PREVIOUS_SIGNIFICANT_FILE, sep="\t", usecols=["group_id"])
site_hits = pd.read_csv(SITE_HIT_FILE, sep="\t", usecols=["phenotype_id"])

groups["primary_significant"] = groups["primary_significant"].astype(str).str.lower().eq("true")
site_hits["group_id"] = site_hits["phenotype_id"].str.split("|", regex=False).str[1]
eigenmt_ids = set(significant["group_id"].astype(str))
previous_ids = set(previous["group_id"].astype(str))
site_hit_ids = set(site_hits["group_id"].astype(str))

groups["previous_115"] = groups["group_id"].astype(str).isin(previous_ids)
groups["new_83"] = groups["primary_significant"] & ~groups["previous_115"]
groups["site_level_hit_group"] = groups["group_id"].astype(str).isin(site_hit_ids)
groups["group_only_eigenmt"] = groups["primary_significant"] & ~groups["site_level_hit_group"]

# --- 1. Meff relative to tested pairs ---
tested = groups.loc[groups["n_tested_pairs"] > 0].copy()
tested["meff_per_pair"] = tested["tests_emt"] / tested["n_tested_pairs"]
tested["group_factor"] = tested["tests_emt"] * tested["num_phenotypes"]
tested["group_factor_per_pair"] = tested["group_factor"] / tested["n_tested_pairs"]

ratio_columns = [
    "group_id", "gene_id", "n_group_phenotypes", "n_tested_pairs",
    "tests_emt", "num_phenotypes", "meff_per_pair", "group_factor",
    "group_factor_per_pair", "group_q_eigenmt", "primary_significant",
    "previous_115", "new_83", "site_level_hit_group", "group_only_eigenmt",
]
tested[ratio_columns].to_csv(D / "interaction_eigenmt_meff_pair_ratio.tsv", sep="\t", index=False)

scope_masks = {
    "tested_groups_with_pairs": np.ones(len(tested), dtype=bool),
    "eigenmt_198": tested["primary_significant"].to_numpy(),
    "previous_115": tested["previous_115"].to_numpy(),
    "new_83": tested["new_83"].to_numpy(),
}
quantiles = [("min", 0.00), ("p01", 0.01), ("p05", 0.05), ("p25", 0.25), ("median", 0.50),
             ("p75", 0.75), ("p95", 0.95), ("p99", 0.99), ("max", 1.00)]
ratio_summary_records = []
for scope, mask in scope_masks.items():
    subset = tested.loc[mask]
    for metric in ["meff_per_pair", "group_factor_per_pair"]:
        values = subset[metric].to_numpy(dtype=float)
        for label, probability in quantiles:
            ratio_summary_records.append({
                "scope": scope, "n_groups": len(subset), "metric": metric,
                "quantile": label, "value": float(np.quantile(values, probability)),
            })
pd.DataFrame(ratio_summary_records).to_csv(D / "interaction_eigenmt_meff_pair_ratio_summary.tsv",
                                           sep="\t", index=False)

# --- 2. multi-site aggregation check ---
top_check = tested.merge(
    eigenmt_top.rename(columns={"num_phenotypes": "top_num_phenotypes", "tests_emt": "top_tests_emt"}),
    on="group_id", how="left", validate="one_to_one",
)
top_check["expected_group_p"] = np.minimum(
    1.0, top_check["min_pval_gi"] * top_check["num_phenotypes"] * top_check["tests_emt"])
top_check["group_p_abs_diff"] = np.abs(top_check["expected_group_p"] - top_check["group_p_eigenmt"])
top_check["declared_site_count_match"] = top_check["n_group_phenotypes"].eq(top_check["num_phenotypes"])
top_check["tested_site_count_match"] = top_check["n_tested_phenotypes"].eq(top_check["num_phenotypes"])
top_check["top_site_count_match"] = top_check["num_phenotypes"].eq(top_check["top_num_phenotypes"])
top_check["top_meff_match"] = top_check["tests_emt"].eq(top_check["top_tests_emt"])

multisite = top_check.loc[top_check["n_group_phenotypes"] > 1, [
    "group_id", "gene_id", "n_group_phenotypes", "n_tested_phenotypes",
    "n_tested_pairs", "tests_emt", "num_phenotypes", "group_factor",
    "meff_per_pair", "group_factor_per_pair", "min_pval_gi",
    "expected_group_p", "group_p_eigenmt", "group_p_abs_diff",
    "group_q_eigenmt", "primary_significant", "new_83",
]].copy()
multisite["meff_aggregation"] = "one Meff from the shared group cis-genotype matrix"
multisite["phenotype_aggregation"] = "minimum p multiplied by num_phenotypes"
multisite.to_csv(D / "interaction_eigenmt_multisite_aggregation.tsv", sep="\t", index=False)

site_count_mismatches = top_check.loc[~top_check["declared_site_count_match"], [
    "group_id", "gene_id", "n_group_phenotypes", "n_tested_phenotypes",
    "num_phenotypes", "n_tested_pairs", "tests_emt", "group_q_eigenmt",
]].copy()
site_count_mismatches["interpretation"] = (
    "one declared phenotype was dropped as constant before eigenMT; "
    "num_phenotypes matches the nominally tested phenotype count")
site_count_mismatches.to_csv(D / "interaction_eigenmt_site_count_mismatches.tsv", sep="\t", index=False)

aggregation_audit = pd.DataFrame([
    ("tensorqtl_group_range_contract", "all phenotype cis ranges within a group must be identical"),
    ("meff_aggregation", "compute_tests once on the shared group genotype matrix after interaction MAF filtering"),
    ("phenotype_aggregation", "multiply the group minimum p by num_phenotypes"),
    ("group_p_formula", "min(1, min_p_GxStage * num_phenotypes * tests_emt)"),
    ("tested_groups_total", len(groups)),
    ("tested_groups_with_pairs", len(tested)),
    ("zero_pair_groups", int((groups["n_tested_pairs"] == 0).sum())),
    ("multi_site_groups", len(multisite)),
    ("multi_site_eigenmt_significant", int(multisite["primary_significant"].sum())),
    ("summary_top_row_mismatches", int(top_check["top_num_phenotypes"].isna().sum())),
    ("declared_site_count_mismatches", int((~top_check["declared_site_count_match"]).sum())),
    ("tested_site_count_mismatches", int((~top_check["tested_site_count_match"]).sum())),
    ("top_site_count_mismatches", int((~top_check["top_site_count_match"]).sum())),
    ("top_meff_mismatches", int((~top_check["top_meff_match"]).sum())),
    ("group_p_formula_mismatches_gt_1e-12", int((top_check["group_p_abs_diff"] > 1e-12).sum())),
], columns=["field", "value"])
aggregation_audit.to_csv(D / "interaction_eigenmt_aggregation_audit.tsv", sep="\t", index=False)

# --- 3. groups significant only at group level ---
group_only_columns = [
    "group_id", "gene_id", "group_q_eigenmt", "group_p_eigenmt",
    "n_group_phenotypes", "n_tested_pairs", "tests_emt", "num_phenotypes",
    "meff_per_pair", "group_factor_per_pair", "stage_group_status",
    "lead_phenotype_id", "lead_variant_id", "lead_b_g", "lead_b_gi",
]
group_only = tested.loc[tested["group_only_eigenmt"], group_only_columns].sort_values(
    ["group_q_eigenmt", "group_id"], kind="mergesort")
group_only.to_csv(D / "interaction_eigenmt_group_only_genes.tsv", sep="\t", index=False)

# --- gate ---
checks = {
    "site_level_hits_381": len(site_hits) == 381,
    "site_level_groups_268": len(site_hit_ids) == 268,
    "eigenmt_groups_198": len(eigenmt_ids) == 198,
    "previous_groups_115": len(previous_ids) == 115,
    "previous_115_retained": previous_ids.issubset(eigenmt_ids),
    "new_groups_83": int(groups["new_83"].sum()) == 83,
    "group_only_genes_10": len(group_only) == 10,
    "top_table_complete_for_tested_groups": top_check["top_num_phenotypes"].notna().all(),
    "tested_group_site_counts_match": top_check["tested_site_count_match"].all(),
    "top_site_counts_match": top_check["top_site_count_match"].all(),
    "top_meff_matches_summary": top_check["top_meff_match"].all(),
    "group_p_formula_reproduced": bool((top_check["group_p_abs_diff"] <= 1e-12).all()),
}
gate = pd.DataFrame({"check": list(checks), "pass": list(checks.values())})
gate.to_csv(D / "interaction_eigenmt_diagnostic_gate.tsv", sep="\t", index=False)
if not all(checks.values()):
    failed = [name for name, passed in checks.items() if not passed]
    raise RuntimeError(f"eigenMT diagnostic gate failed: {', '.join(failed)}")
print(f"[done] gate=PASS; site hits={len(site_hits)}; site groups={len(site_hit_ids)}; "
      f"eigenMT groups={len(eigenmt_ids)}; group-only={len(group_only)}")
