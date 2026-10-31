# HC3 refit of all eigenMT leads: strict direction reversal, boundary (PDUI 0.001/0.999)
# sensitivity, stage-specific genotype support and within-stage phenotype scaling.
# Writes interaction_hc3_* / interaction_refit_* / interaction_robustness_* to data/05_gxs.
from pathlib import Path
import os

import numpy as np
import pandas as pd
import pgenlib
from scipy import stats

# --- paths and settings ---
APA_ROOT = Path(os.environ["APA_ROOT"])
D = APA_ROOT / "data" / "05_gxs"
LEAD_FILE = D / "interaction_significant_groups.tsv"
PREVIOUS_LEAD_FILE = D / "previous115_interaction_significant_groups.tsv"
PREVIOUS_HC3_FILE = D / "previous115_interaction_robustness_summary.tsv"
PHENOTYPE_FILE = D / "pre.phenotype.combat.int.txt"
RAW_PDUI_FILE = D / "before.combat.txt"
COVARIATE_FILE = D / "covariates_for_qtl.txt"
STAGE_FILE = D / "LifeStage.txt"
GENOTYPE_PREFIX = APA_ROOT / "raw" / "genotype"
PHENOTYPE_CHUNK_SIZE = 4096
HC3_P_CUTOFF = 0.05
CI_LEVEL = 0.95
MIN_GENOTYPE_CELL_N = 20
VALIDATION_BETA_TOL = 2e-4
VALIDATION_SE_TOL = 2e-4
ORIGINAL_HC3_UNSTABLE = {"HSPA8", "MUC20-OT1", "EIF4G2"}
BOUNDARY_LOW = 0.001
BOUNDARY_HIGH = 0.999
BOUNDARY_ATOL = 1e-12
BOUNDARY_MIN_STAGE_N = 10
BOUNDARY_MIN_CARRIER_N = 3


# --- helpers ---
def make_residualizer(covariate_values: np.ndarray) -> np.ndarray:
    centered = covariate_values - covariate_values.mean(axis=0, keepdims=True)
    q, _ = np.linalg.qr(centered, mode="reduced")
    if np.linalg.matrix_rank(centered) != centered.shape[1]:
        raise ValueError("Covariate matrix is rank deficient")
    return q


def residualize_centered(values: np.ndarray, q: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    one_dimensional = values.ndim == 1
    if one_dimensional:
        values = values[None, :]
    residualized = values - (values @ q) @ q.T
    return residualized[0] if one_dimensional else residualized


def fit_ols(phenotype, raw_columns, q, covariate_count):
    y = np.asarray(phenotype, dtype=np.float64)
    y_res = residualize_centered(y - y.mean(), q)
    x_mat = np.column_stack([
        residualize_centered(np.asarray(column, dtype=np.float64) - np.mean(column), q)
        for column in raw_columns
    ])
    xtx_inv = np.linalg.inv(x_mat.T @ x_mat)
    beta = xtx_inv @ x_mat.T @ y_res
    residuals = y_res - x_mat @ beta
    degrees_of_freedom = int(y.size - covariate_count - 1 - x_mat.shape[1])

    sigma2 = float(residuals @ residuals / degrees_of_freedom)
    ordinary_cov = (xtx_inv * sigma2 + (xtx_inv * sigma2).T) / 2.0
    ordinary_se = np.sqrt(np.maximum(np.diag(ordinary_cov), 0.0))
    ordinary_t = beta / ordinary_se
    ordinary_p = 2.0 * stats.t.sf(np.abs(ordinary_t), degrees_of_freedom)

    leverage = np.einsum("ij,jk,ik->i", x_mat, xtx_inv, x_mat)
    hc3_weights = (residuals / np.maximum(1.0 - leverage, 1e-10)) ** 2
    meat = (x_mat.T * hc3_weights) @ x_mat
    hc3_cov = xtx_inv @ meat @ xtx_inv
    hc3_cov = (hc3_cov + hc3_cov.T) / 2.0
    hc3_se = np.sqrt(np.maximum(np.diag(hc3_cov), 0.0))
    hc3_t = beta / hc3_se
    hc3_p = 2.0 * stats.t.sf(np.abs(hc3_t), degrees_of_freedom)

    return {
        "beta": beta, "ordinary_cov": ordinary_cov, "ordinary_se": ordinary_se, "ordinary_p": ordinary_p,
        "hc3_cov": hc3_cov, "hc3_se": hc3_se, "hc3_p": hc3_p, "degrees_of_freedom": degrees_of_freedom,
    }


def point_estimate_class(beta_pre: float, beta_post: float) -> str:
    if beta_pre * beta_post < 0:
        return "Opposing"
    if abs(beta_post) > abs(beta_pre):
        return "Amplified"
    return "Attenuated"


def within_stage_standardize(phenotype, stage):
    result = np.asarray(phenotype, dtype=np.float64).copy()
    for stage_value in (0, 1):
        mask = stage == stage_value
        stage_sd = float(result[mask].std(ddof=1))
        if not np.isfinite(stage_sd) or stage_sd <= 0:
            raise ValueError(f"Stage {stage_value} phenotype variance is not positive")
        result[mask] = (result[mask] - result[mask].mean()) / stage_sd
    return result


def minor_allele_carrier_count(genotype: np.ndarray) -> int:
    observed = genotype[(genotype >= 0) & (genotype <= 2)].astype(int)
    if observed.size == 0:
        return 0
    allele_frequency = float(observed.sum() / (2 * observed.size))
    if allele_frequency <= 0.5:
        return int(np.sum(observed > 0))
    return int(np.sum(observed < 2))


# --- covariates, stage and sample order ---
leads = pd.read_csv(LEAD_FILE, sep="\t")
previous_ids = set(pd.read_csv(PREVIOUS_LEAD_FILE, sep="\t", usecols=["group_id"])["group_id"].astype(str))
previous_hc3 = pd.read_csv(PREVIOUS_HC3_FILE, sep="\t",
                           usecols=["group_id", "hc3_beta_gxstage", "hc3_se_gxstage", "hc3_pval_gxstage"])
leads["previous_115"] = leads["group_id"].astype(str).isin(previous_ids)
leads["new_83"] = ~leads["previous_115"]

covariates = pd.read_csv(COVARIATE_FILE, sep="\t", index_col=0)
covariates.index = covariates.index.astype(str)
covariates.columns = covariates.columns.astype(str)
covariates = covariates.apply(pd.to_numeric, errors="raise")

stage_table = pd.read_csv(STAGE_FILE, sep="\t", usecols=[0, 1])
stage_table.columns = ["sample", "stage"]
stage_series = stage_table.set_index(stage_table["sample"].astype(str))["stage"]

psam = pd.read_csv(f"{GENOTYPE_PREFIX}.psam", sep="\t")
sample_order = psam.iloc[:, 0].astype(str).tolist()
covariates = covariates.loc[:, sample_order]
stage = stage_series.loc[sample_order].to_numpy(dtype=np.float64)
covariate_values = covariates.to_numpy(dtype=np.float64).T
q = make_residualizer(covariate_values)

# --- lead genotypes (pgen) ---
variant_ids = leads["lead_variant_id"].astype(str).tolist()
target_variants = set(variant_ids)
variant_indices = {}
variant_index = 0
with open(f"{GENOTYPE_PREFIX}.pvar", "rt", encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        variant_id = line.rstrip("\n").split("\t")[2]
        if variant_id in target_variants:
            variant_indices[variant_id] = variant_index
        variant_index += 1

reader = pgenlib.PgenReader(f"{GENOTYPE_PREFIX}.pgen".encode())
raw_genotypes = np.empty((len(leads), len(sample_order)), dtype=np.int8)
reader.read_list(np.asarray([variant_indices[v] for v in variant_ids], dtype=np.uint32), raw_genotypes)
reader.close()
genotypes = raw_genotypes.astype(np.float64)
for row_index in range(len(leads)):
    missing = genotypes[row_index] < 0
    if missing.any():
        genotypes[row_index, missing] = genotypes[row_index, ~missing].mean()

# --- lead phenotypes (scaled PDUI) and raw PDUI for the boundary sensitivity ---
phenotype_ids = leads["lead_phenotype_id"].astype(str).tolist()
wanted_phenotypes = set(phenotype_ids)
phenotype_index_column = str(pd.read_csv(PHENOTYPE_FILE, sep="\t", nrows=0).columns[0])
phenotypes = {}
for chunk in pd.read_csv(PHENOTYPE_FILE, sep="\t", usecols=[phenotype_index_column] + sample_order,
                         index_col=0, dtype={s: np.float32 for s in sample_order},
                         chunksize=PHENOTYPE_CHUNK_SIZE):
    chunk.index = chunk.index.astype(str)
    for phenotype_id in chunk.index.intersection(wanted_phenotypes):
        phenotypes[phenotype_id] = chunk.loc[phenotype_id, sample_order].to_numpy(dtype=np.float64)

raw_index_column = str(pd.read_csv(RAW_PDUI_FILE, sep="\t", nrows=0).columns[0])
raw_header = pd.read_csv(RAW_PDUI_FILE, sep="\t", nrows=0)
missing_samples = sorted(set(sample_order) - set(map(str, raw_header.columns[1:])))
if missing_samples:
    raise RuntimeError(f"raw PDUI matrix is missing genotype samples: {missing_samples[:5]}")
raw_pdui = {}
for chunk in pd.read_csv(RAW_PDUI_FILE, sep="\t", usecols=[raw_index_column] + sample_order, index_col=0,
                         dtype={s: np.float64 for s in sample_order}, chunksize=PHENOTYPE_CHUNK_SIZE):
    chunk.index = chunk.index.astype(str)
    for phenotype_id in chunk.index.intersection(wanted_phenotypes):
        raw_pdui[phenotype_id] = chunk.loc[phenotype_id, sample_order].to_numpy(dtype=np.float64)
missing_targets = wanted_phenotypes.difference(raw_pdui)
if missing_targets:
    raise RuntimeError(f"raw PDUI matrix is missing target phenotypes: {sorted(missing_targets)[:5]}")

# --- HC3 refit per lead, with genotype support and within-stage scaling ---
records = []
validation_records = []
for row_number, lead in leads.iterrows():
    genotype = genotypes[row_number]
    phenotype = phenotypes[str(lead["lead_phenotype_id"])]
    fit = fit_ols(phenotype, [genotype, stage, genotype * stage], q, covariate_values.shape[1])

    beta = np.asarray(fit["beta"], dtype=float)
    ordinary_se = np.asarray(fit["ordinary_se"], dtype=float)
    hc3_cov = np.asarray(fit["hc3_cov"], dtype=float)
    hc3_se = np.asarray(fit["hc3_se"], dtype=float)
    hc3_p = np.asarray(fit["hc3_p"], dtype=float)

    # ordinary fit must reproduce the TensorQTL lead values
    canonical_beta = np.array([lead["lead_b_g"], lead["lead_b_i"], lead["lead_b_gi"]], dtype=float)
    canonical_se = np.array([lead["lead_b_g_se"], lead["lead_b_i_se"], lead["lead_b_gi_se"]], dtype=float)
    beta_diff = np.abs(beta - canonical_beta)
    se_diff = np.abs(ordinary_se - canonical_se)
    if not (np.all(beta_diff <= VALIDATION_BETA_TOL) and np.all(se_diff <= VALIDATION_SE_TOL)):
        raise RuntimeError(f"ordinary TensorQTL reconstruction failed for {lead['group_id']}")

    beta_pre = float(beta[0])
    beta_gxstage = float(beta[2])
    beta_post = beta_pre + beta_gxstage
    hc3_var_pre = float(hc3_cov[0, 0])
    hc3_var_gxstage = float(hc3_cov[2, 2])
    hc3_cov_pre_gxstage = float(hc3_cov[0, 2])
    hc3_var_post = hc3_var_pre + hc3_var_gxstage + 2.0 * hc3_cov_pre_gxstage
    if not np.isfinite(hc3_var_post) or hc3_var_post <= 0:
        raise RuntimeError(f"invalid HC3 postnatal variance for {lead['group_id']}: {hc3_var_post}")
    hc3_se_post = float(np.sqrt(hc3_var_post))

    critical = float(stats.t.ppf(0.5 + CI_LEVEL / 2.0, int(fit["degrees_of_freedom"])))
    pre_low = beta_pre - critical * float(hc3_se[0])
    pre_high = beta_pre + critical * float(hc3_se[0])
    post_low = beta_post - critical * hc3_se_post
    post_high = beta_post + critical * hc3_se_post
    descriptive_class = point_estimate_class(beta_pre, beta_post)
    strict_reversal = bool((pre_low > 0 and post_high < 0) or (pre_high < 0 and post_low > 0))
    hc3_robust = bool(np.isfinite(hc3_p[2]) and hc3_p[2] < HC3_P_CUTOFF)

    # genotype support (stage-specific counts, MAF, MAC)
    raw_genotype = raw_genotypes[row_number]
    support = {}
    all_counts = []
    for stage_value, stage_name in ((0, "prenatal"), (1, "postnatal")):
        observed = raw_genotype[stage == stage_value]
        observed = observed[(observed >= 0) & (observed <= 2)].astype(int)
        counts = np.bincount(observed, minlength=3)
        all_counts.extend(counts.tolist())
        allele_total = 2 * len(observed)
        alt_count = int(observed.sum())
        af = alt_count / allele_total
        support.update({
            f"{stage_name}_g0_n": int(counts[0]), f"{stage_name}_g1_n": int(counts[1]),
            f"{stage_name}_g2_n": int(counts[2]), f"{stage_name}_maf": min(af, 1.0 - af),
            f"{stage_name}_mac": min(alt_count, allele_total - alt_count),
        })
    support["min_genotype_cell_n"] = int(min(all_counts))
    support["genotype_support_status"] = (
        f"sufficient_at_{MIN_GENOTYPE_CELL_N}_development_threshold"
        if support["min_genotype_cell_n"] >= MIN_GENOTYPE_CELL_N
        else f"insufficient_at_{MIN_GENOTYPE_CELL_N}_development_threshold")

    # within-stage phenotype standardisation sensitivity
    try:
        scaled_phenotype = within_stage_standardize(phenotype, stage)
        scaled_fit = fit_ols(scaled_phenotype, [genotype, stage, genotype * stage], q, covariate_values.shape[1])
    except (np.linalg.LinAlgError, ValueError) as exc:
        scaling = {
            "scaled_beta_gxstage": np.nan, "scaled_se_gxstage": np.nan, "scaled_pval_gxstage": np.nan,
            "scaled_beta_pre": np.nan, "scaled_beta_post": np.nan,
            "scaled_interaction_sign_consistent": "not_estimable",
            "phenotype_scaling_status": "not_estimable", "phenotype_scaling_reason": str(exc),
        }
    else:
        scaled_beta_pre = float(scaled_fit["beta"][0])
        scaled_beta_gxstage = float(scaled_fit["beta"][2])
        scaling = {
            "scaled_beta_gxstage": scaled_beta_gxstage,
            "scaled_se_gxstage": float(scaled_fit["ordinary_se"][2]),
            "scaled_pval_gxstage": float(scaled_fit["ordinary_p"][2]),
            "scaled_beta_pre": scaled_beta_pre,
            "scaled_beta_post": scaled_beta_pre + scaled_beta_gxstage,
            "scaled_interaction_sign_consistent": "yes" if np.sign(scaled_beta_gxstage) == np.sign(beta_gxstage) else "no",
            "phenotype_scaling_status": "complete", "phenotype_scaling_reason": "",
        }

    record = {
        "group_id": lead["group_id"],
        "gene_id": lead["gene_id"],
        "previous_115": bool(lead["previous_115"]),
        "new_83": bool(lead["new_83"]),
        "group_q_eigenmt": float(lead["group_q_eigenmt"]),
        "lead_phenotype_id": lead["lead_phenotype_id"],
        "lead_variant_id": lead["lead_variant_id"],
        "degrees_of_freedom": int(fit["degrees_of_freedom"]),
        "beta_pre": beta_pre,
        "beta_gxstage": beta_gxstage,
        "beta_post": beta_post,
        "ordinary_se_gxstage": float(ordinary_se[2]),
        "hc3_se_pre": float(hc3_se[0]),
        "hc3_se_gxstage": float(hc3_se[2]),
        "hc3_var_pre": hc3_var_pre,
        "hc3_var_gxstage": hc3_var_gxstage,
        "hc3_cov_pre_gxstage": hc3_cov_pre_gxstage,
        "hc3_twice_cov_pre_gxstage": 2.0 * hc3_cov_pre_gxstage,
        "hc3_var_post": hc3_var_post,
        "hc3_se_post": hc3_se_post,
        "hc3_pval_gxstage": float(hc3_p[2]),
        "hc3_se_ratio": float(hc3_se[2] / ordinary_se[2]),
        "hc3_pre_ci95_low": pre_low,
        "hc3_pre_ci95_high": pre_high,
        "hc3_post_ci95_low": post_low,
        "hc3_post_ci95_high": post_high,
        "point_estimate_class": descriptive_class,
        "point_estimate_class_role": "descriptive only",
        "strict_ci_direction_reversal": strict_reversal,
        "strict_reversal_definition": "prenatal and postnatal HC3 95% CIs exclude zero in opposite directions",
        "hc3_robust_core": hc3_robust,
        "hc3_core_rule": "finite HC3 p_GxStage < 0.05",
    }
    record.update(support)
    record.update(scaling)
    records.append(record)
    validation_records.append({
        "group_id": lead["group_id"],
        "beta_g_abs_diff": beta_diff[0], "beta_stage_abs_diff": beta_diff[1], "beta_gxstage_abs_diff": beta_diff[2],
        "se_g_abs_diff": se_diff[0], "se_stage_abs_diff": se_diff[1], "se_gxstage_abs_diff": se_diff[2],
        "validation_status": "pass",
    })
    if row_number == 0 or (row_number + 1) % 25 == 0:
        print(f"[HC3] refit {row_number + 1}/{len(leads)}", flush=True)

hc3 = pd.DataFrame(records)
validation = pd.DataFrame(validation_records)

# --- previous-115 comparison and core set ---
old_comparison = hc3.loc[hc3["previous_115"], ["group_id", "beta_gxstage", "hc3_se_gxstage", "hc3_pval_gxstage"]].merge(
    previous_hc3, on="group_id", how="left", validate="one_to_one")
old_comparison["hc3_beta_abs_diff"] = np.abs(old_comparison["beta_gxstage"] - old_comparison["hc3_beta_gxstage"])
old_comparison["hc3_se_abs_diff"] = np.abs(old_comparison["hc3_se_gxstage_x"] - old_comparison["hc3_se_gxstage_y"])
old_comparison["hc3_p_abs_diff"] = np.abs(old_comparison["hc3_pval_gxstage_x"] - old_comparison["hc3_pval_gxstage_y"])
old_comparison.to_csv(D / "interaction_hc3_previous115_validation.tsv", sep="\t", index=False)

annotation_columns = ["group_id"] + [c for c in hc3.columns if c not in leads.columns]
annotated = leads.merge(hc3[annotation_columns], on="group_id", how="left", validate="one_to_one")
core = annotated.loc[annotated["hc3_robust_core"]].copy()
excluded = annotated.loc[~annotated["hc3_robust_core"]].copy()
unexpected_original_core = ORIGINAL_HC3_UNSTABLE.intersection(set(core["group_id"].astype(str)))
if unexpected_original_core:
    raise RuntimeError("previous HC3-unstable groups incorrectly entered the core: "
                       + ", ".join(sorted(unexpected_original_core)))

# --- PDUI-boundary sensitivity (drop raw PDUI == 0.001 or 0.999) for the HC3 core ---
if len(core) != 193:
    raise RuntimeError(f"expected 193 HC3-core groups for PDUI-boundary sensitivity, found {len(core)}")
lead_row_by_group = leads.reset_index().set_index("group_id")["index"].to_dict()
boundary_records = []
for _, group in core.iterrows():
    group_id = str(group["group_id"])
    lead_row_number = int(lead_row_by_group[group_id])
    phenotype_id = str(group["lead_phenotype_id"])
    raw_profile = raw_pdui[phenotype_id]
    phenotype = phenotypes[phenotype_id]
    genotype = genotypes[lead_row_number]
    raw_genotype = raw_genotypes[lead_row_number]
    boundary_mask = (np.isclose(raw_profile, BOUNDARY_LOW, atol=BOUNDARY_ATOL, rtol=0.0)
                     | np.isclose(raw_profile, BOUNDARY_HIGH, atol=BOUNDARY_ATOL, rtol=0.0))
    retained = np.isfinite(raw_profile) & ~boundary_mask & np.isfinite(phenotype)
    boundary_stage = stage[retained]
    pre_mask = retained & (stage == 0)
    post_mask = retained & (stage == 1)
    pre_n = int(pre_mask.sum())
    post_n = int(post_mask.sum())
    pre_carriers = minor_allele_carrier_count(raw_genotype[pre_mask])
    post_carriers = minor_allele_carrier_count(raw_genotype[post_mask])
    evaluable = bool(pre_n >= BOUNDARY_MIN_STAGE_N and post_n >= BOUNDARY_MIN_STAGE_N
                     and pre_carriers >= BOUNDARY_MIN_CARRIER_N and post_carriers >= BOUNDARY_MIN_CARRIER_N)
    record = {
        "group_id": group_id,
        "gene_id": group["gene_id"],
        "lead_phenotype_id": phenotype_id,
        "lead_variant_id": group["lead_variant_id"],
        "n_boundary_profiles": int(boundary_mask.sum()),
        "n_nonfinite_profiles": int((~np.isfinite(raw_profile)).sum()),
        "n_retained_profiles": int(retained.sum()),
        "prenatal_n_after_boundary": pre_n,
        "postnatal_n_after_boundary": post_n,
        "prenatal_minor_allele_carrier_n": pre_carriers,
        "postnatal_minor_allele_carrier_n": post_carriers,
        "evaluable_definition": "each stage n>=10 and minor-allele carrier n>=3 after removing raw PDUI exactly 0.001 or 0.999",
        "evaluable": evaluable,
        "boundary_hc3_pval_gxstage": np.nan,
        "interaction_p_lt_0_05": False,
        "fit_status": "not_evaluable" if not evaluable else "pending",
    }
    if evaluable:
        try:
            filtered_covariates = covariate_values[retained, :]
            filtered_q = make_residualizer(filtered_covariates)
            filtered_fit = fit_ols(phenotype[retained],
                                   [genotype[retained], boundary_stage, genotype[retained] * boundary_stage],
                                   filtered_q, filtered_covariates.shape[1])
            boundary_p = float(filtered_fit["hc3_p"][2])
            record["boundary_hc3_pval_gxstage"] = boundary_p
            record["interaction_p_lt_0_05"] = bool(np.isfinite(boundary_p) and boundary_p < HC3_P_CUTOFF)
            record["fit_status"] = "complete"
        except (np.linalg.LinAlgError, ValueError, FloatingPointError) as exc:
            record["fit_status"] = "not_estimable"
            record["fit_reason"] = str(exc)
    boundary_records.append(record)

boundary_sensitivity = pd.DataFrame(boundary_records)
boundary_summary = pd.DataFrame([{
    "n_input_hc3_core_groups": len(boundary_sensitivity),
    "n_evaluable_groups": int(boundary_sensitivity["evaluable"].sum()),
    "n_evaluable_interaction_p_lt_0_05": int(
        boundary_sensitivity.loc[boundary_sensitivity["evaluable"], "interaction_p_lt_0_05"].sum()),
    "n_not_evaluable_groups": int((~boundary_sensitivity["evaluable"]).sum()),
    "boundary_values": "raw PDUI == 0.001 or 0.999",
    "evaluable_definition": "each stage retains >=10 people and >=3 minor-allele carriers",
    "interaction_definition": "finite HC3 p_GxStage < 0.05 after profile removal",
    "source_raw_pdui": str(RAW_PDUI_FILE),
}])
boundary_sensitivity.to_csv(D / "interaction_hc3_boundary_sensitivity.tsv", sep="\t", index=False)
boundary_summary.to_csv(D / "interaction_hc3_boundary_sensitivity_summary.tsv", sep="\t", index=False)

# --- reversal counts by scope ---
strata = {
    "all_198": annotated,
    "previous_115": annotated.loc[annotated["previous_115"]],
    "new_83": annotated.loc[annotated["new_83"]],
    "hc3_robust_core": core,
    "hc3_robust_previous": core.loc[core["previous_115"]],
    "hc3_robust_new": core.loc[core["new_83"]],
}
reversal_summary = pd.DataFrame([{
    "scope": scope,
    "n_groups": len(table),
    "n_hc3_robust": int(table["hc3_robust_core"].sum()),
    "n_point_estimate_opposing": int(table["point_estimate_class"].eq("Opposing").sum()),
    "n_strict_ci_direction_reversal": int(table["strict_ci_direction_reversal"].sum()),
    "strict_reversal_definition": "prenatal and postnatal HC3 95% CIs exclude zero in opposite directions",
} for scope, table in strata.items()])

# --- write outputs ---
annotated.to_csv(D / "interaction_hc3_summary.tsv", sep="\t", index=False)
annotated.to_csv(D / "interaction_robustness_summary.tsv", sep="\t", index=False)
core.to_csv(D / "interaction_hc3_core_groups.tsv", sep="\t", index=False)
excluded.to_csv(D / "interaction_hc3_excluded_groups.tsv", sep="\t", index=False)
validation.to_csv(D / "interaction_refit_validation.tsv", sep="\t", index=False)
reversal_summary.to_csv(D / "interaction_hc3_reversal_summary.tsv", sep="\t", index=False)
print(f"[done] all={len(annotated)}; HC3 core={len(core)}; excluded={len(excluded)}; "
      f"strict reversals={int(core['strict_ci_direction_reversal'].sum())}; "
      f"scaling complete={int(annotated['phenotype_scaling_status'].eq('complete').sum())}")
