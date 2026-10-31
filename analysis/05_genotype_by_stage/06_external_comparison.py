# BrainGVEX (adult, postnatal) validation of the 193 HC3-core G x Stage leads:
# exact phenotype_id + variant_id matching, direction concordance, BH q over matched and
# over all 193 prespecified leads. Writes interaction_external_* to data/05_gxs.
from pathlib import Path
import csv
import gzip
import math
import os

import numpy as np
from scipy import stats

# --- paths and constants ---
APA_ROOT = Path(os.environ["APA_ROOT"])
D = APA_ROOT / "data" / "05_gxs"
Q4 = APA_ROOT / "data" / "04_qtl"
LEAD_FILE = D / "interaction_hc3_core_groups.tsv"
EXTERNAL_NOMINAL_FILE = Q4 / "eur_adult.cis_qtl_pairs.txt.gz"
EXTERNAL_SIGNIFICANT_FILE = Q4 / "eur_adult.significant_pairs.txt.gz"
EXTERNAL_SAMPLE_FILE = Q4 / "sample.Adult.final"
OUT_TABLE = D / "interaction_external_component.tsv"
OUT_SUMMARY = D / "interaction_external_component_summary.tsv"

SITE_HITS = 381
SITE_HIT_GROUPS = 268
EIGENMT_GROUPS = 198
HC3_CORE_GROUPS = 193


# --- helpers ---
def read_tsv(path):
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add_bh(rows, q_field, total_tests):
    finite = sorted((row for row in rows if row["external_p"] is not None), key=lambda row: float(row["external_p"]))
    running = 1.0
    for rank, row in reversed(list(enumerate(finite, start=1))):
        running = min(running, float(row["external_p"]) * total_tests / rank)
        row[q_field] = min(running, 1.0)


def format_value(value):
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.12g}"
    return str(value)


def association_metrics(rows):
    if not rows:
        return {"n": 0, "direction_concordant": 0, "direction_discordant": 0,
                "direction_concordant_fraction": None, "pearson_r": None, "pearson_p": None,
                "spearman_rho": None, "spearman_p": None, "direction_binomial_p_one_sided": None}
    local = np.asarray([float(row["local_post_beta"]) for row in rows])
    external = np.asarray([float(row["external_beta"]) for row in rows])
    concordant = int(np.sum(np.sign(local) == np.sign(external)))
    if len(rows) >= 2:
        pearson = stats.pearsonr(local, external)
        spearman = stats.spearmanr(local, external)
        pearson_r, pearson_p = float(pearson.statistic), float(pearson.pvalue)
        spearman_r, spearman_p = float(spearman.statistic), float(spearman.pvalue)
    else:
        pearson_r = pearson_p = spearman_r = spearman_p = None
    return {
        "n": len(rows),
        "direction_concordant": concordant,
        "direction_discordant": len(rows) - concordant,
        "direction_concordant_fraction": concordant / len(rows),
        "pearson_r": pearson_r, "pearson_p": pearson_p,
        "spearman_rho": spearman_r, "spearman_p": spearman_p,
        "direction_binomial_p_one_sided": float(stats.binomtest(concordant, len(rows), 0.5, alternative="greater").pvalue),
    }


# --- local HC3-core leads ---
leads = read_tsv(LEAD_FILE)
if len(leads) != HC3_CORE_GROUPS:
    raise ValueError(f"Expected {HC3_CORE_GROUPS} HC3-core leads, found {len(leads)}")
external_significant = {(row["phenotype_id"], row["variant_id"]) for row in read_tsv(EXTERNAL_SIGNIFICANT_FILE)}
external_sample_n = sum(1 for _ in EXTERNAL_SAMPLE_FILE.open())

lead_rows = {}
for row in leads:
    key = row["lead_phenotype_id"], row["lead_variant_id"]
    if key in lead_rows:
        raise ValueError(f"Duplicate lead key: {key}")
    lead_rows[key] = {
        "group_id": row["group_id"],
        "gene_id": row["gene_id"],
        "lead_phenotype_id": row["lead_phenotype_id"],
        "lead_variant_id": row["lead_variant_id"],
        "group_q_eigenmt": float(row["group_q_eigenmt"]),
        "local_pre_beta": float(row["beta_pre"]),
        "local_post_beta": float(row["beta_post"]),
        "local_gxstage_beta": float(row["beta_gxstage"]),
        "local_hc3_pre_ci95_low": float(row["hc3_pre_ci95_low"]),
        "local_hc3_pre_ci95_high": float(row["hc3_pre_ci95_high"]),
        "local_hc3_post_ci95_low": float(row["hc3_post_ci95_low"]),
        "local_hc3_post_ci95_high": float(row["hc3_post_ci95_high"]),
        "point_estimate_class": row["point_estimate_class"],
        "point_estimate_class_role": row["point_estimate_class_role"],
        "strict_ci_direction_reversal": row["strict_ci_direction_reversal"],
        "external_beta": None, "external_se": None, "external_p": None,
        "external_q_bh_matched": None, "external_q_bh_prespecified193": None,
        "external_af": None, "external_ma_count": None, "external_in_significant_pairs": False,
    }

# --- match to the external nominal pair table ---
with gzip.open(EXTERNAL_NOMINAL_FILE, "rt", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        key = row["phenotype_id"], row["variant_id"]
        if key not in lead_rows:
            continue
        target = lead_rows[key]
        target["external_beta"] = float(row["slope"])
        target["external_se"] = float(row["slope_se"])
        target["external_p"] = float(row["pval_nominal"])
        target["external_af"] = float(row["af"])
        target["external_ma_count"] = int(row["ma_count"])
        target["external_in_significant_pairs"] = key in external_significant

rows = list(lead_rows.values())
matched = [row for row in rows if row["external_p"] is not None]
add_bh(matched, "external_q_bh_matched", len(matched))
add_bh(matched, "external_q_bh_prespecified193", len(rows))

for row in rows:
    if row["external_beta"] is None:
        row.update({"external_direction_concordant": None, "external_nominal_support": None,
                    "external_q_support_matched": None, "external_q_support_prespecified193": None})
        continue
    local_sign = (float(row["local_post_beta"]) > 0) - (float(row["local_post_beta"]) < 0)
    external_sign = (float(row["external_beta"]) > 0) - (float(row["external_beta"]) < 0)
    row["external_direction_concordant"] = local_sign == external_sign
    row["external_nominal_support"] = float(row["external_p"]) < 0.05
    row["external_q_support_matched"] = float(row["external_q_bh_matched"]) < 0.05
    row["external_q_support_prespecified193"] = float(row["external_q_bh_prespecified193"]) < 0.05

# --- per-lead table ---
table_fields = [
    "group_id", "gene_id", "lead_phenotype_id", "lead_variant_id",
    "group_q_eigenmt", "local_pre_beta", "local_gxstage_beta",
    "local_post_beta", "local_hc3_pre_ci95_low", "local_hc3_pre_ci95_high",
    "local_hc3_post_ci95_low", "local_hc3_post_ci95_high",
    "point_estimate_class", "point_estimate_class_role",
    "strict_ci_direction_reversal", "external_beta", "external_se",
    "external_p", "external_q_bh_matched", "external_q_bh_prespecified193",
    "external_af", "external_ma_count", "external_in_significant_pairs",
    "external_direction_concordant", "external_nominal_support",
    "external_q_support_matched", "external_q_support_prespecified193",
    "external_source", "external_sample_n", "matching_key",
    "local_effect_source", "validation_scope", "result_hierarchy",
]
with OUT_TABLE.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=table_fields, delimiter="\t")
    writer.writeheader()
    for row in sorted(rows, key=lambda value: float(value["group_q_eigenmt"])):
        output = {field: row.get(field, "") for field in table_fields}
        output.update({
            "external_source": "BrainGVEX eur_adult apaQTL",
            "external_sample_n": external_sample_n,
            "matching_key": "phenotype_id + variant_id; exact REF/ALT string",
            "local_effect_source": "193-group HC3 core; direct beta_post",
            "validation_scope": "adult/postnatal component; not full GxStage replication",
            "result_hierarchy": "381 site-level / 268 groups -> 198 eigenMT -> 193 HC3 core",
        })
        writer.writerow({key: format_value(value) for key, value in output.items()})

# --- summary ---
matched_metrics = association_metrics(matched)
q_matched = [row for row in matched if row["external_q_support_matched"]]
q_prespecified = [row for row in matched if row["external_q_support_prespecified193"]]
q_matched_metrics = association_metrics(q_matched)
summary = [
    ("site_level_hits", SITE_HITS),
    ("represented_site_hit_groups", SITE_HIT_GROUPS),
    ("eigenmt_calibrated_groups", EIGENMT_GROUPS),
    ("hc3_robust_core_groups", HC3_CORE_GROUPS),
    ("lead_tests", len(rows)),
    ("external_sample_n", external_sample_n),
    ("exact_key_matched_leads", len(matched)),
    ("exact_key_unmatched_leads", len(rows) - len(matched)),
    ("external_nominal_p_lt_0.05", sum(bool(row["external_nominal_support"]) for row in matched)),
    ("external_bh_q_lt_0.05_matched", len(q_matched)),
    ("external_bh_q_lt_0.05_prespecified193", len(q_prespecified)),
]
summary.extend((f"all_matched_{key}", value) for key, value in matched_metrics.items())
summary.extend((f"q_matched_{key}", value) for key, value in q_matched_metrics.items())
summary.extend([
    ("q_matched_direction_discordant_genes", ",".join(
        str(row["gene_id"]) for row in q_matched if not row["external_direction_concordant"])),
    ("primary_multiple_testing", "BH over observed exact-key matches"),
    ("sensitivity_multiple_testing", "BH over 193 prespecified HC3-core leads with unmatched leads retained in denominator"),
    ("validation_scope", "adult/postnatal component; not full GxStage replication"),
])
with OUT_SUMMARY.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    for metric, value in summary:
        writer.writerow([metric, format_value(value)])
print(f"[done] matched={len(matched)}/{len(rows)}; direction={matched_metrics['direction_concordant']}/{len(matched)}; "
      f"q_matched={len(q_matched)}; q_prespecified193={len(q_prespecified)}")
