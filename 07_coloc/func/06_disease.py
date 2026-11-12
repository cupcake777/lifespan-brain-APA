#!/usr/bin/env python3
# Internal eQTL contract, then GWAS disease chains and association evidence tiers.
# Reads data/07_coloc tables and func/config.yml; writes data/07_coloc tables 04, 05 and 07.

from __future__ import annotations

import gzip
import json
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


# --- helpers ---

def _expand_env(obj):
    # expand ${VAR} in config values (paths set via config/paths.sh)
    if isinstance(obj, dict):
        return {k: _expand_env(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env(v) for v in obj]
    return os.path.expandvars(obj) if isinstance(obj, str) else obj


def normalize_gene_id(value: object) -> str:
    return str(value).split(".", 1)[0]


def read_samples(path: Path) -> set[str]:
    with path.open(encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def as_bool(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
RESULTS = ROOT / "data/07_coloc"
RESULTS.mkdir(parents=True, exist_ok=True)
CONFIG_PATH = Path(__file__).resolve().parent / "config.yml"
EQTL_PATH = RESULTS / "internal_eqtl.tsv"
MECHANISM_PATH = RESULTS / "02_signal_mechanisms.tsv.gz"
OUTPUT_PATH = RESULTS / "07_internal_eqtl_integration_contract.tsv.gz"
SUMMARY_PATH = RESULTS / "internal_summary.tsv"
CONTRACT_PATH = RESULTS / "internal_contract.json"
CONSEQUENCES = RESULTS / "03_event_consequences.tsv.gz"
ASSOCIATIONS = RESULTS / "01_significant_qtl_universe.tsv.gz"
EXTERNAL_QTL = RESULTS / "external_all.tsv"
OUT_DISEASE = RESULTS / "04_disease_chains_results.tsv.gz"
OUT_TIERS = RESULTS / "05_evidence_tiers_results.tsv.gz"
ANALYSIS_VERSION = "framework_2026-07-25_v5"
PP4_THRESHOLD = 0.75

MECHANISM_ORDER = {
    "direct_high_confidence": 0,
    "local_regulatory_supported": 1,
    "direct_exploratory": 2,
    "predicted_only": 3,
    "association_mechanism_unresolved": 4,
}
CONSEQUENCE_ORDER = {
    "mediation_supported": 0,
    "conditional_association_only": 1,
    "predicted_consequence_only": 2,
    "consequence_unresolved": 3,
}
DISEASE_ORDER = {
    "complete_high_confidence_chain": 0,
    "high_confidence_locus_incomplete_mechanism": 1,
    "exploratory_disease_chain": 2,
    "regional_disease_overlap_only": 3,
    "no_disease_evidence": 4,
}

with CONFIG_PATH.open(encoding="utf-8") as handle:
    config = _expand_env(yaml.safe_load(handle))

# --- internal eQTL contract ---

gtf_path = Path(config["reference_resources"]["gene_annotation_gtf"]["path"])
symbol_to_ids = defaultdict(set)
id_to_symbol = {}
transcript_to_gene = {}
opener = gzip.open if gtf_path.suffix == ".gz" else open
with opener(gtf_path, "rt", encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] not in {"gene", "transcript"}:
            continue
        gene_id_match = re.search(r'gene_id "([^"]+)"', fields[8])
        gene_name_match = re.search(r'gene_name "([^"]+)"', fields[8])
        if not gene_id_match:
            continue
        gene_id = normalize_gene_id(gene_id_match.group(1))
        gene_name = gene_name_match.group(1) if gene_name_match else gene_id
        if fields[2] == "gene":
            symbol_to_ids[gene_name].add(gene_id)
            symbol_to_ids[gene_id].add(gene_id)
            id_to_symbol[gene_id] = gene_name
        else:
            transcript_match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if transcript_match:
                transcript_to_gene[normalize_gene_id(transcript_match.group(1))] = gene_id

mechanisms = pd.read_csv(
    MECHANISM_PATH,
    sep="\t",
    usecols=["candidate_id", "signal_id", "stage", "phenotype_id", "gene_name"],
).drop_duplicates("candidate_id")

candidate_ids = []
mapping_status = []
for row in mechanisms.itertuples(index=False):
    transcript_id = normalize_gene_id(str(row.phenotype_id).split("|", 1)[0])
    transcript_gene_id = transcript_to_gene.get(transcript_id)
    ids = [transcript_gene_id] if transcript_gene_id else sorted(symbol_to_ids.get(str(row.gene_name), set()))
    candidate_ids.append(";".join(ids) if ids else pd.NA)
    mapping_status.append(
        "exact_transcript_to_gene"
        if transcript_gene_id
        else "exact_unique_symbol"
        if len(ids) == 1
        else "ambiguous_multiple_gencode_ids"
        if ids
        else "unmapped"
    )
mechanisms["candidate_gene_id"] = candidate_ids
mechanisms["candidate_gene_mapping_status"] = mapping_status

eqtl = pd.read_csv(EQTL_PATH, sep="\t").rename(columns={"phenotype_id": "egene_id"})
if set(eqtl["analysis_version"].dropna().astype(str)) != {ANALYSIS_VERSION}:
    raise ValueError("Internal eQTL coloc table is not a current-v5 raw calculation")
if eqtl["pipeline_run_id"].nunique() != 1:
    raise ValueError("Internal eQTL coloc table contains multiple pipeline runs")
if set(eqtl["candidate_id"].dropna()) != set(mechanisms["candidate_id"]):
    raise ValueError("Internal eQTL coloc candidate coverage does not match current mechanism signals")
for stage in ("prenatal", "postnatal"):
    stage_table = eqtl.loc[eqtl["stage"].eq(stage)]
    expected_apa = str(config["inputs"][f"{stage}_nominal_apaqtl"]["path"])
    expected_eqtl = str(config["inputs"][f"{stage}_internal_eqtl_nominal"]["path"])
    if set(stage_table["apa_nominal_path"].dropna().astype(str)) != {expected_apa}:
        raise ValueError(f"Internal eQTL coloc APA provenance is incorrect for {stage}")
    if set(stage_table["eqtl_nominal_path"].dropna().astype(str)) != {expected_eqtl}:
        raise ValueError(f"Internal eQTL coloc eQTL provenance is incorrect for {stage}")

table = eqtl.merge(
    mechanisms,
    on=["candidate_id", "stage"],
    how="left",
    validate="many_to_one",
)
if table["gene_name"].isna().any():
    raise ValueError("At least one internal eQTL candidate lacks an exact mechanism-table mapping")

table["egene_id_normalized"] = table["egene_id"].map(normalize_gene_id)
table["egene_symbol"] = table["egene_id_normalized"].map(id_to_symbol)
candidate_sets = table["candidate_gene_id"].fillna("").map(
    lambda value: set(str(value).split(";")) if value else set()
)
table["gene_mapping_scope"] = [
    "same_gene" if egene in candidate_set else "nearby_egene_context"
    for egene, candidate_set in zip(table["egene_id_normalized"], candidate_sets)
]

sample_metadata = {}
for stage in ("prenatal", "postnatal"):
    apa_path = Path(config["inputs"][f"{stage}_qtl_samples"]["path"])
    eqtl_path = Path(config["inputs"][f"{stage}_internal_eqtl_samples"]["path"])
    apa_samples = read_samples(apa_path)
    eqtl_samples = read_samples(eqtl_path)
    overlap = apa_samples.intersection(eqtl_samples)
    sample_metadata[stage] = {
        "apa_sample_n": len(apa_samples),
        "eqtl_sample_n": len(eqtl_samples),
        "overlap_sample_n": len(overlap),
        "overlap_fraction_apa": len(overlap) / len(apa_samples),
        "overlap_fraction_eqtl": len(overlap) / len(eqtl_samples),
        "apa_sample_path": str(apa_path),
        "eqtl_sample_path": str(eqtl_path),
    }

for column in (
    "apa_sample_n", "eqtl_sample_n", "overlap_sample_n",
    "overlap_fraction_apa", "overlap_fraction_eqtl",
):
    table[column] = table["stage"].map({stage: values[column] for stage, values in sample_metadata.items()})

input_support = table["eQTL_coloc"].astype(str).str.lower().isin({"true", "1"})
table["same_gene_coloc_supported"] = (
    table["gene_mapping_scope"].eq("same_gene")
    & input_support
    & pd.to_numeric(table["eqtl_coloc_pp4"], errors="coerce").ge(PP4_THRESHOLD)
)
table["sample_overlap_status"] = "substantial_same_cohort_overlap"
table["evidence_role"] = table["gene_mapping_scope"].map({
    "same_gene": "primary_same_cohort_same_gene_molecular_sharing",
    "nearby_egene_context": "nearby_egene_context_only",
})
table["orthogonal_evidence"] = False
table["primary_functional_support"] = table["same_gene_coloc_supported"]
table["standalone_consequence_eligible"] = False
table["tier_eligible"] = table["primary_functional_support"]
table["integration_status"] = pd.Series(
    "same_gene_not_colocalized", index=table.index, dtype="object"
)
table.loc[
    table["gene_mapping_scope"].eq("nearby_egene_context"), "integration_status"
] = "context_only_not_transferable_to_candidate_gene"
table.loc[table["same_gene_coloc_supported"], "integration_status"] = "primary_same_cohort_support"
table["interpretation_limit"] = (
    "primary_same_cohort_genetic_sharing_not_independent_replication_or_mediation"
)
table["recommended_sensitivity"] = "coloc_prior_sensitivity"
table["gencode_release"] = config["reference_resources"]["gene_annotation_gtf"]["release"]
table["analysis_version"] = ANALYSIS_VERSION
table["eqtl_pipeline_run_id"] = table["pipeline_run_id"]
table["eqtl_apa_nominal_path"] = table["apa_nominal_path"]
table["eqtl_nominal_path"] = table["eqtl_nominal_path"]
table["config_path"] = str(CONFIG_PATH)

output_columns = [
    "candidate_id", "signal_id", "stage", "phenotype_id", "gene_name",
    "candidate_gene_id", "candidate_gene_mapping_status", "egene_id",
    "egene_id_normalized", "egene_symbol", "gene_mapping_scope", "eQTL_coloc",
    "eqtl_coloc_pp4", "eqtl_coloc_status", "same_gene_coloc_supported",
    "apa_sample_n", "eqtl_sample_n", "overlap_sample_n", "overlap_fraction_apa",
    "overlap_fraction_eqtl", "sample_overlap_status", "evidence_role",
    "orthogonal_evidence", "primary_functional_support",
    "standalone_consequence_eligible", "tier_eligible", "integration_status",
    "interpretation_limit", "recommended_sensitivity", "gencode_release",
    "analysis_version", "eqtl_pipeline_run_id", "eqtl_apa_nominal_path", "eqtl_nominal_path", "config_path",
]
table[output_columns].to_csv(OUTPUT_PATH, sep="\t", index=False, compression="gzip")

summary_rows = []
for stage, frame in [("all", table), *table.groupby("stage", sort=True)]:
    summary_rows.extend([
        {"stage": stage, "metric": "rows", "value": len(frame)},
        {"stage": stage, "metric": "candidates", "value": frame["candidate_id"].nunique()},
        {"stage": stage, "metric": "same_gene_rows", "value": frame["gene_mapping_scope"].eq("same_gene").sum()},
        {"stage": stage, "metric": "same_gene_supported_rows", "value": frame["same_gene_coloc_supported"].sum()},
        {
            "stage": stage,
            "metric": "same_gene_supported_candidates",
            "value": frame.loc[frame["same_gene_coloc_supported"], "candidate_id"].nunique(),
        },
        {
            "stage": stage,
            "metric": "candidates_without_same_gene_row",
            "value": frame["candidate_id"].nunique()
            - frame.loc[frame["gene_mapping_scope"].eq("same_gene"), "candidate_id"].nunique(),
        },
        {"stage": stage, "metric": "primary_functional_support_rows", "value": frame["primary_functional_support"].sum()},
        {"stage": stage, "metric": "standalone_consequence_eligible_rows", "value": frame["standalone_consequence_eligible"].sum()},
        {"stage": stage, "metric": "tier_eligible_rows", "value": frame["tier_eligible"].sum()},
        {"stage": stage, "metric": "orthogonal_evidence_rows", "value": frame["orthogonal_evidence"].sum()},
    ])
summary = pd.DataFrame(summary_rows)
summary["analysis_version"] = ANALYSIS_VERSION
summary["pipeline_run_id"] = eqtl["pipeline_run_id"].iloc[0]
summary.to_csv(SUMMARY_PATH, sep="\t", index=False)

contract = {
    "analysis_version": ANALYSIS_VERSION,
    "evidence_role": "primary_same_cohort_molecular_sharing",
    "same_gene_mapping": "exact candidate-gene to eGene Ensembl ID via GENCODE v40",
    "nearby_gene_policy": "retain as context only; never transfer to the APA candidate gene",
    "orthogonal_evidence": False,
    "tier_eligible_rule": "supported exact same-gene rows support the functional tier but are not standalone consequences",
    "standalone_consequence_eligible": False,
    "pp4_reporting_threshold": PP4_THRESHOLD,
    "sample_overlap": sample_metadata,
    "recommended_robustness": ["coloc prior sensitivity"],
    "inputs": {
        "eqtl_coloc": str(EQTL_PATH),
        "eqtl_pipeline_run_id": eqtl["pipeline_run_id"].iloc[0],
        "raw_provenance": {
            stage: {
                "apa_nominal_path": str(eqtl.loc[eqtl["stage"].eq(stage), "apa_nominal_path"].iloc[0]),
                "eqtl_nominal_path": str(eqtl.loc[eqtl["stage"].eq(stage), "eqtl_nominal_path"].iloc[0]),
            }
            for stage in ("prenatal", "postnatal")
        },
        "mechanisms": str(MECHANISM_PATH),
        "gencode_gtf": str(gtf_path),
    },
    "outputs": {
        "mapped_table": str(OUTPUT_PATH),
        "summary": str(SUMMARY_PATH),
    },
}
CONTRACT_PATH.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

# --- disease integration ---

thresholds = config["candidate_sets"]["evidence_thresholds"]
pp4_exploratory = float(thresholds["pp4_exploratory"]["value"])
pp4_high = float(thresholds["pp4_high_confidence"]["value"])
pip_high = float(thresholds["pip_high_confidence"]["value"])

mechanisms = pd.read_csv(MECHANISM_PATH, sep="\t")
consequences = pd.DataFrame(columns=["signal_id", "candidate_id", "downstream_consequence_evidence_label"])
if CONSEQUENCES.exists():
    consequence_table = pd.read_csv(CONSEQUENCES, sep="\t")
    if not consequence_table.empty:
        consequence_table["_order"] = consequence_table["downstream_consequence_evidence_label"].map(CONSEQUENCE_ORDER).fillna(99)
        consequences = (
            consequence_table.sort_values(["_order", "indirect_fdr"], na_position="last", kind="stable")
            .drop_duplicates("candidate_id")
            .drop(columns="_order")
        )
if consequences["downstream_consequence_evidence_label"].eq("orthogonal_support_only").any():
    raise ValueError("Stale orthogonal-support rows are not allowed in the v5 consequence input")
stale_columns = [
    column for column in consequences.columns
    if column.startswith("orthogonal_") or column.startswith("same_gene_external_")
]
if stale_columns:
    raise ValueError(f"Stale external evidence columns are not allowed in v5: {stale_columns}")

internal_eqtl = pd.read_csv(OUTPUT_PATH, sep="\t", low_memory=False)
if set(internal_eqtl["analysis_version"].dropna().astype(str)) != {ANALYSIS_VERSION}:
    raise ValueError("Internal eQTL contract does not match current-v5")
internal_support = (
    internal_eqtl.loc[internal_eqtl["gene_mapping_scope"].eq("same_gene")]
    .assign(eqtl_coloc_pp4=lambda frame: pd.to_numeric(frame["eqtl_coloc_pp4"], errors="coerce"))
    .sort_values(["primary_functional_support", "eqtl_coloc_pp4"], ascending=[False, False], kind="stable")
    .drop_duplicates("candidate_id")
    [["candidate_id", "primary_functional_support", "eqtl_coloc_pp4"]]
    .rename(columns={
        "primary_functional_support": "internal_same_cohort_eqtl_support",
        "eqtl_coloc_pp4": "internal_same_cohort_eqtl_pp4",
    })
)

external_qtl = pd.read_csv(EXTERNAL_QTL, sep="\t", low_memory=False)
if set(external_qtl["analysis_version"].dropna().astype(str)) != {ANALYSIS_VERSION}:
    raise ValueError("External-QTL contract does not match current-v5")
external_qtl["external_validation_eligible"] = external_qtl["external_validation_eligible"].map(as_bool)
external_qtl["PP.H4"] = pd.to_numeric(external_qtl["PP.H4"], errors="coerce")
external_support = (
    external_qtl.sort_values(
        ["external_validation_eligible", "PP.H4"], ascending=[False, False], kind="stable"
    )
    .drop_duplicates("candidate_id")
    [["candidate_id", "external_validation_eligible", "source", "qtl_type", "PP.H4", "evidence_role"]]
    .rename(columns={
        "external_validation_eligible": "eligible_external_qtl_support",
        "source": "strongest_external_qtl_source",
        "qtl_type": "strongest_external_qtl_type",
        "PP.H4": "strongest_external_qtl_pp4",
        "evidence_role": "strongest_external_qtl_evidence_role",
    })
)

# --- disease comparisons ---

coloc = pd.read_csv(config["inputs"]["gwas_coloc_results"]["path"])
disease = coloc.rename(columns={
    "GWAS_ID": "gwas_trait",
    "QTL_ID": "stage",
    "Locus": "locus",
    "Phenotype": "phenotype_id",
    "PP0": "gwas_coloc_pp0",
    "PP1": "gwas_coloc_pp1",
    "PP2": "gwas_coloc_pp2",
    "PP3": "gwas_coloc_pp3",
    "PP4": "gwas_coloc_pp4",
})
posterior_columns = [f"gwas_coloc_pp{index}" for index in range(5)]
for column in posterior_columns:
    disease[column] = pd.to_numeric(disease[column], errors="coerce")
disease["disease_comparison_id"] = (
    disease["gwas_trait"].astype(str) + "::" + disease["stage"].astype(str) + "::"
    + disease["phenotype_id"].astype(str) + "::" + disease["locus"].astype(str)
)
if disease["disease_comparison_id"].duplicated().any():
    raise ValueError("Source coloc results are not unique by trait, stage, phenotype, and locus")

frames = []
for path_text in config["inputs"]["coloc_prior_sensitivity_results"]["paths"]:
    path = Path(path_text)
    if not path.exists():
        continue
    frame = pd.read_csv(path, sep="\t")
    if "phenotype_id" not in frame.columns and "candidate_id" in frame.columns:
        frame["phenotype_id"] = frame["candidate_id"].astype(str).str.split("::", n=2).str[1]
    if "stage" not in frame.columns and "candidate_id" in frame.columns:
        frame["stage"] = frame["candidate_id"].astype(str).str.split("::", n=2).str[0]
    if "gwas_trait" not in frame.columns:
        frame["gwas_trait"] = "EUR_SCZ"
    locus_column = next(
        (
            column
            for column in ["historical_gwas_locus_label", "historical_gwas_hg19_locus_label", "gwas_locus"]
            if column in frame.columns
        ),
        None,
    )
    frame["locus"] = frame[locus_column] if locus_column else pd.NA
    robust_column = next(
        (column for column in ["prior_robust_local_coloc", "default_pairwise_support_prior_robust"] if column in frame.columns),
        None,
    )
    shared_column = next(
        (column for column in ["candidate_in_conditional_95pct_H4_set", "candidate_in_shared_credible_set"] if column in frame.columns),
        None,
    )
    frame["coloc_prior_sensitivity_acceptable"] = frame[robust_column].map(as_bool) if robust_column else False
    frame["proposed_variant_in_shared_set"] = frame[shared_column].map(as_bool) if shared_column else pd.NA
    keep = [
        "gwas_trait", "stage", "phenotype_id", "locus", "candidate_id",
        "coloc_prior_sensitivity_acceptable", "proposed_variant_in_shared_set",
    ]
    frames.append(frame[[column for column in keep if column in frame.columns]])
if not frames:
    prior = pd.DataFrame(columns=[
        "gwas_trait", "stage", "phenotype_id", "locus", "candidate_id",
        "coloc_prior_sensitivity_acceptable", "proposed_variant_in_shared_set",
    ])
else:
    combined = pd.concat(frames, ignore_index=True, sort=False)
    combined["coloc_prior_sensitivity_acceptable"] = combined["coloc_prior_sensitivity_acceptable"].map(as_bool)
    combined["proposed_variant_in_shared_set"] = combined["proposed_variant_in_shared_set"].map(as_bool)
    prior = combined.groupby(
        ["gwas_trait", "stage", "phenotype_id", "locus", "candidate_id"], as_index=False, dropna=False
    ).agg(
        coloc_prior_sensitivity_acceptable=("coloc_prior_sensitivity_acceptable", "max"),
        proposed_variant_in_shared_set=("proposed_variant_in_shared_set", "max"),
    )

mechanism_columns = [
    "signal_id", "candidate_id", "variant_id", "pip", "finemap_tier",
    "snp_to_apa_evidence_label", "mechanism_star_eligible", "direct_pas_evidence",
    "interaction_group_id", "interaction_gene_id",
    "interaction_representative_phenotype_id", "interaction_representative_variant_id",
    "interaction_previous_115", "interaction_new_83",
    "gxstage_class", "gxstage_classification_status", "gxstage_class_source",
    "pattern_tau", "pattern_delta", "pattern_method",
    "interaction_pair_signal_assignment_method", "interaction_pair_signal_lead_variant",
    "interaction_pair_ld_r2_to_signal_lead", "interaction_pair_cs_id",
    "gxstage_mechanism_variant_link", "gxstage_mechanism_variant_link_r2",
]
candidate_evidence = prior.merge(
    mechanisms[mechanism_columns],
    on="candidate_id",
    how="left",
    validate="many_to_one",
)
candidate_evidence = candidate_evidence.merge(
    consequences[[column for column in consequences.columns if column in {
        "candidate_id", "downstream_consequence_evidence_label", "indirect_effect", "indirect_ci_low",
        "indirect_ci_high", "indirect_fdr", "outcome_type", "outcome_id",
    }]],
    on="candidate_id",
    how="left",
    validate="many_to_one",
)
candidate_evidence = candidate_evidence.merge(
    internal_support, on="candidate_id", how="left", validate="many_to_one"
)
candidate_evidence = candidate_evidence.merge(
    external_support, on="candidate_id", how="left", validate="many_to_one"
)
candidate_evidence["_mechanism_order"] = candidate_evidence["snp_to_apa_evidence_label"].map(MECHANISM_ORDER).fillna(99)
candidate_evidence["_consequence_order"] = candidate_evidence["downstream_consequence_evidence_label"].map(CONSEQUENCE_ORDER).fillna(99)
candidate_evidence["_prior_order"] = ~candidate_evidence["coloc_prior_sensitivity_acceptable"].map(as_bool)
candidate_evidence["_shared_order"] = ~candidate_evidence["proposed_variant_in_shared_set"].map(as_bool)
candidate_evidence["_pip_order"] = -pd.to_numeric(candidate_evidence["pip"], errors="coerce").fillna(-1)
candidate_keys = ["gwas_trait", "stage", "phenotype_id", "locus"]
candidate_lists = candidate_evidence.groupby(candidate_keys, as_index=False, dropna=False).agg(
    locus_candidate_count=("candidate_id", "nunique"),
    locus_candidate_ids=("candidate_id", lambda values: ";".join(sorted(set(values.dropna().astype(str))))),
)
candidate_best = (
    candidate_evidence.sort_values(
        candidate_keys + ["_prior_order", "_shared_order", "_mechanism_order", "_consequence_order", "_pip_order", "candidate_id"],
        kind="stable",
    )
    .drop_duplicates(candidate_keys)
    .drop(columns=["_prior_order", "_shared_order", "_mechanism_order", "_consequence_order", "_pip_order"])
    .merge(candidate_lists, on=candidate_keys, how="left", validate="one_to_one")
)
disease = disease.merge(candidate_best, on=candidate_keys, how="left", validate="one_to_one")
disease = disease.rename(columns={"variant_id": "proposed_functional_variant", "pip": "proposed_variant_pip"})
disease["locus_candidate_link_status"] = np.where(
    disease["candidate_id"].notna(),
    "exact_candidate_mapping",
    "no_exact_candidate_mapping",
)
disease["coloc_prior_sensitivity_status"] = np.where(
    disease["coloc_prior_sensitivity_acceptable"].map(as_bool), "acceptable", "not_available_or_not_robust"
)

disease[["population", "trait"]] = disease["gwas_trait"].astype(str).str.split("_", n=1, expand=True)
disease["sldsc_annotation"] = disease["stage"].astype(str) + "_apa_qtl"
sldsc = pd.read_csv(config["inputs"]["sldsc_enrichment_results"]["path"], sep="\t").rename(columns={
    "Trait": "trait", "Population": "population", "Annotation": "sldsc_annotation",
    "Enrichment": "sldsc_enrichment", "Enrichment_std_error": "sldsc_enrichment_se",
    "Enrichment_p": "sldsc_enrichment_p", "FDR_p": "sldsc_enrichment_fdr",
})
disease = disease.merge(
    sldsc[[
        "trait", "population", "sldsc_annotation", "sldsc_enrichment", "sldsc_enrichment_se",
        "sldsc_enrichment_p", "sldsc_enrichment_fdr",
    ]],
    on=["trait", "population", "sldsc_annotation"],
    how="left",
    validate="many_to_one",
)

disease["coloc_posterior_sum"] = disease[posterior_columns].sum(axis=1, min_count=5)
disease["coloc_posterior_sum_status"] = np.where(
    disease["coloc_posterior_sum"].sub(1).abs().le(1e-6),
    "pass",
    "fail",
)
if disease["coloc_posterior_sum_status"].eq("fail").any():
    raise ValueError("At least one coloc posterior vector does not sum to one")

disease["coloc_evidence_tier"] = "not_supported"
disease.loc[disease["gwas_coloc_pp4"].gt(pp4_exploratory), "coloc_evidence_tier"] = "exploratory"
disease.loc[disease["gwas_coloc_pp4"].ge(pp4_high), "coloc_evidence_tier"] = "high_pp4_prior_unresolved"
disease.loc[
    disease["gwas_coloc_pp4"].ge(pp4_high) & disease["coloc_prior_sensitivity_acceptable"].map(as_bool),
    "coloc_evidence_tier",
] = "high_confidence_prior_robust"

labels = []
for _, row in disease.iterrows():
    high_coloc = row["coloc_evidence_tier"] == "high_confidence_prior_robust"
    high_variant = pd.notna(row.get("proposed_variant_pip")) and row["proposed_variant_pip"] >= pip_high
    mechanism = row.get("snp_to_apa_evidence_label")
    consequence = row.get("downstream_consequence_evidence_label")
    strong_mechanism = mechanism in {"direct_high_confidence", "local_regulatory_supported"}
    strong_consequence = consequence == "mediation_supported"
    shared = as_bool(row.get("proposed_variant_in_shared_set"))
    if high_coloc and high_variant and strong_mechanism and strong_consequence and shared:
        labels.append("complete_high_confidence_chain")
        continue
    if high_coloc:
        labels.append("high_confidence_locus_incomplete_mechanism")
        continue
    partial = (
        mechanism in {"direct_exploratory", "predicted_only", "local_regulatory_supported"}
        or consequence in {"conditional_association_only", "predicted_consequence_only", "mediation_supported"}
        or as_bool(row.get("internal_same_cohort_eqtl_support"))
        or as_bool(row.get("eligible_external_qtl_support"))
    )
    if row["gwas_coloc_pp4"] > 0.5 and partial:
        labels.append("exploratory_disease_chain")
    elif row["gwas_coloc_pp4"] > 0.5:
        labels.append("regional_disease_overlap_only")
    else:
        labels.append("no_disease_evidence")
disease["disease_chain_evidence_label"] = labels
disease["analysis_version"] = ANALYSIS_VERSION
disease["config_path"] = str(CONFIG_PATH)
disease["pp4_exploratory_threshold"] = pp4_exploratory
disease["pp4_high_confidence_threshold"] = pp4_high
disease["disease_formal_unit"] = "trait_x_stage_x_phenotype_x_locus"

# --- association evidence tiers ---

associations = pd.read_csv(ASSOCIATIONS, sep="\t", low_memory=False)
tier_mechanisms = mechanisms[[
    "signal_id", "candidate_id", "variant_id", "snp_to_apa_evidence_label", "mechanism_star_eligible"
]].rename(columns={
    "candidate_id": "mechanism_candidate_id",
    "variant_id": "mechanism_variant_id",
})
tier_table = associations.merge(
    tier_mechanisms,
    on="signal_id",
    how="left",
)
tier_table = tier_table.merge(
    consequences[[column for column in consequences.columns if column in {
        "signal_id", "downstream_consequence_evidence_label", "indirect_fdr",
    }]],
    on="signal_id",
    how="left",
)
disease_best = disease.copy()
disease_best["_disease_order"] = disease_best["disease_chain_evidence_label"].map(DISEASE_ORDER).fillna(99)
disease_best = disease_best.loc[disease_best["signal_id"].notna()].sort_values(
    ["_disease_order", "gwas_coloc_pp4"], ascending=[True, False], kind="stable"
).drop_duplicates("signal_id")
tier_table = tier_table.merge(
    disease_best[[
        "signal_id", "gwas_trait", "gwas_coloc_pp4", "coloc_evidence_tier",
        "disease_chain_evidence_label", "proposed_variant_in_shared_set",
        "internal_same_cohort_eqtl_support", "internal_same_cohort_eqtl_pp4",
        "eligible_external_qtl_support", "strongest_external_qtl_source",
        "strongest_external_qtl_type", "strongest_external_qtl_pp4",
        "strongest_external_qtl_evidence_role",
    ]],
    on="signal_id",
    how="left",
)
tier_table["evidence_tier"] = "D"
exploratory = (
    tier_table["finemap_tier"].isin(["exploratory", "high_confidence"])
    | tier_table["snp_to_apa_evidence_label"].isin(["direct_exploratory", "predicted_only", "local_regulatory_supported"])
    | tier_table["disease_chain_evidence_label"].isin(["exploratory_disease_chain", "regional_disease_overlap_only"])
)
tier_table.loc[exploratory, "evidence_tier"] = "C"
tier_table.loc[
    tier_table["disease_chain_evidence_label"].eq("high_confidence_locus_incomplete_mechanism"),
    "evidence_tier",
] = "B"
tier_table.loc[
    tier_table["disease_chain_evidence_label"].eq("complete_high_confidence_chain"),
    "evidence_tier",
] = "A"
tier_table["disease_chain_star_eligible"] = tier_table["disease_chain_evidence_label"].eq("complete_high_confidence_chain")
tier_table["composite_total_score"] = "disabled"
tier_table["analysis_version"] = ANALYSIS_VERSION

# --- output ---

disease.sort_values(["gwas_trait", "stage", "phenotype_id", "locus"], kind="stable").to_csv(
    OUT_DISEASE, sep="\t", index=False, compression="gzip"
)
tier_table.sort_values(["evidence_tier", "stage", "phenotype_id", "variant_pos"], kind="stable").to_csv(
    OUT_TIERS, sep="\t", index=False, compression="gzip"
)
