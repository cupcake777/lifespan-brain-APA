#!/usr/bin/env python3
# Validate the functional-QTL package and write its validation table and manifest.
# Inputs are the module 07 tables in data/07_coloc and the cluster 05_func work directory.

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pandas as pd
import yaml


# --- helpers ---

def read(path: Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def propagated_signal_metadata_pass(
    frame: pd.DataFrame,
    source: pd.DataFrame,
    fields: list[str],
) -> bool:
    required = {"signal_id", *fields}
    if not required.issubset(frame.columns) or not required.issubset(source.columns):
        return False
    if source["signal_id"].duplicated().any():
        return False
    unlinked = frame["signal_id"].isna()
    if unlinked.any() and frame.loc[unlinked, fields].notna().any().any():
        return False
    linked = frame.loc[~unlinked, ["signal_id", *fields]].merge(
        source[["signal_id", *fields]],
        on="signal_id",
        how="left",
        validate="many_to_one",
        suffixes=("_observed", "_expected"),
        indicator=True,
    )
    if not linked["_merge"].eq("both").all():
        return False
    return all(
        (
            linked[f"{column}_observed"].eq(linked[f"{column}_expected"])
            | (
                linked[f"{column}_observed"].isna()
                & linked[f"{column}_expected"].isna()
            )
        ).all()
        for column in fields
    )


def normalized_row_counts(frame: pd.DataFrame, columns: list[str]) -> dict[tuple[str, ...], int]:
    if not set(columns).issubset(frame.columns):
        return {}
    normalized = frame[columns].fillna("").astype(str)
    return normalized.value_counts(sort=False).to_dict()


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
FUNC_DIR = Path(__file__).resolve().parent
RESULTS = ROOT / "data/07_coloc"
FIGURES = ROOT / "figures/07_coloc"
CLUSTER_WORK = Path(os.environ["APA_HPC_HOME"]) / "project/work/apa/05_func/work"
MECHANISM_WORK = CLUSTER_WORK / "mechanism"
FIGURE_WORK = CLUSTER_WORK / "figures"
VERSION = "framework_2026-07-25_v5"
FINALIZER_VERSION = "10_finalize"
GXSTAGE_CLASSES = {"Opposing", "Amplified", "Attenuated", "Uncertain"}
GXSTAGE_METADATA_FIELDS = [
    "interaction_group_id",
    "interaction_gene_id",
    "interaction_representative_phenotype_id",
    "interaction_representative_variant_id",
    "interaction_previous_115",
    "interaction_new_83",
    "gxstage_class",
    "gxstage_classification_status",
    "gxstage_class_source",
    "pattern_tau",
    "pattern_delta",
    "pattern_method",
]
GXSTAGE_LINK_FIELDS = [
    "interaction_pair_signal_assignment_method",
    "interaction_pair_signal_lead_variant",
    "interaction_pair_ld_r2_to_signal_lead",
    "interaction_pair_cs_id",
    "gxstage_mechanism_variant_link",
    "gxstage_mechanism_variant_link_r2",
]
GXSTAGE_FIELDS = [*GXSTAGE_METADATA_FIELDS, *GXSTAGE_LINK_FIELDS]


# --- reads ---

universe = read(RESULTS / "01_significant_qtl_universe.tsv.gz")
mechanisms = read(RESULTS / "02_signal_mechanisms.tsv.gz")
consequences = read(RESULTS / "03_event_consequences.tsv.gz")
excluded = read(RESULTS / "05_excluded_or_blocked.tsv.gz")
enrichment = read(RESULTS / "06_functional_enrichment.tsv.gz")
gxstage_tests = read(CLUSTER_WORK / "enrichment" / "gxstage.tsv")
enrichment_summary = read(CLUSTER_WORK / "enrichment" / "summary.tsv")
eqtl_raw = read(CLUSTER_WORK / "coloc" / "internal_eqtl.tsv")
internal = read(RESULTS / "07_internal_eqtl_integration_contract.tsv.gz")
pas_details = read(MECHANISM_WORK / "pas.tsv.gz")
rbp_details = read(MECHANISM_WORK / "rbp.tsv.gz")
external = read(CLUSTER_WORK / "coloc" / "external_all.tsv")
external_summary = read(CLUSTER_WORK / "coloc" / "external_summary.tsv")
source_audit = read(CLUSTER_WORK / "coloc" / "external_sources.tsv")
published = read(CLUSTER_WORK / "coloc" / "published.tsv")
disease = read(RESULTS / "04_disease_chains_results.tsv.gz")
tiers = read(RESULTS / "05_evidence_tiers_results.tsv.gz")

coloc_local = {"case_summary.tsv", "case_pairs.tsv", "novel_summary.tsv"}
coloc_names = [
    "case_summary.tsv",
    "case_pairs.tsv",
    "case_sensitivity.tsv",
    "case_posteriors.tsv",
    "case_competitors.tsv",
    "novel_summary.tsv",
    "novel_sensitivity.tsv",
    "novel_posteriors.tsv",
    "novel_competitors.tsv",
]
coloc_tables = {
    name: read((RESULTS if name in coloc_local else CLUSTER_WORK / "coloc") / name)
    for name in coloc_names
}
case_layout = read(FIGURE_WORK / "case_layout.tsv")
case_source = read(RESULTS / "case_source.tsv")
support_layout = read(FIGURE_WORK / "support_layout.tsv")
support_source = read(FIGURE_WORK / "support_source.tsv")
with (FUNC_DIR / "config.yml").open(encoding="utf-8") as handle:
    config = yaml.safe_load(handle)
interaction_group_source = read(Path(os.path.expandvars(config["inputs"]["interaction_group_results"]["path"])))
gxstage_contract_source = read(Path(os.path.expandvars(config["inputs"]["gxstage_classes"]["path"])))
interaction_posterior_source = read(Path(os.path.expandvars(config["inputs"]["interaction_posterior_results"]["path"])))

# --- checks ---

internal_primary = as_bool(internal["primary_functional_support"])
internal_same_gene = internal["gene_mapping_scope"].eq("same_gene")
internal_standalone = as_bool(internal["standalone_consequence_eligible"])
internal_tier = as_bool(internal["tier_eligible"])

external_supported = as_bool(external["coloc_supported"])
external_eligible = as_bool(external["external_validation_eligible"])
eligible_apa = external_eligible & external["qtl_type"].eq("apaQTL")
eligible_pqtl = external_eligible & external["qtl_type"].eq("pQTL")
candidate_eligible = external.groupby("candidate_id")["external_validation_eligible"].apply(
    lambda values: as_bool(values).any()
)
summary_eligible = external_summary.groupby("candidate_id")["any_external_qtl_coloc"].apply(
    lambda values: as_bool(values).any()
)

posterior_columns = [f"gwas_coloc_pp{index}" for index in range(5)]
posterior_sum = disease[posterior_columns].apply(pd.to_numeric, errors="coerce").sum(axis=1, min_count=5)
case_support_actual = {
    row.gene_name: (
        as_bool(pd.Series([row.internal_support_panel])).iloc[0],
        as_bool(pd.Series([row.external_overlap_event_panel])).iloc[0],
    )
    for row in case_layout.itertuples(index=False)
}
case_support_expected = {
    # Same-gene APA interval overlap: CEP41's local
    # chr7:130393771-130399039 event overlaps GTEx Brain_BA9
    # chr7:130393771-130400821 (PP.H4 ~ 1.0, 2603 shared SNPs) -> external TRUE.
    "CEP41": (False, True),
    "GTDC1": (True, False),
    "STMN3": (False, True),
    "ZSWIM7": (True, True),
}
support_figure_expected = {
    "cell_reference_reconstruction",
    "cell_sample_composition",
    "cell_sample_fit",
    "cell_covariate_capture",
    "cell_qtl_residual_signal",
    "internal_eqtl_sharing",
    "external_qtl_support",
    "sldsc_disease_enrichment",
}
case_figure_paths = [
    FIGURES / str(name)
    for column in ["png_file", "pdf_file", "svg_file"]
    for name in case_layout[column]
]
support_figure_paths = [
    FIGURES / str(name)
    for column in ["png_file", "pdf_file", "svg_file"]
    for name in support_layout[column]
]
internal_config = config["parameters"]["internal_eqtl_integration"]
mechanism_candidate_ids = set(mechanisms["candidate_id"])
mechanism_run_ids = set(mechanisms["pipeline_run_id"].dropna().astype(str))
local_candidate_ids = set(
    mechanisms.loc[
        as_bool(mechanisms["local_regulatory_eligible"])
        & as_bool(mechanisms["allele_harmonized"]),
        "candidate_id",
    ]
)

group_source_columns = [
    "group_id", "gene_id", "lead_phenotype_id", "lead_variant_id", "previous_115", "new_83",
]
group_source_key_unique = bool(
    not interaction_group_source["group_id"].duplicated().any()
    and not interaction_group_source.duplicated(["lead_phenotype_id", "lead_variant_id"]).any()
)
interaction_group_map = interaction_group_source[group_source_columns].rename(columns={
    "group_id": "interaction_group_id",
    "gene_id": "interaction_gene_id",
    "lead_phenotype_id": "phenotype_id",
    "lead_variant_id": "variant_id",
    "previous_115": "interaction_previous_115",
    "new_83": "interaction_new_83",
})
interaction_group_map["interaction_representative_phenotype_id"] = interaction_group_map["phenotype_id"]
interaction_group_map["interaction_representative_variant_id"] = interaction_group_map["variant_id"]

gxstage_contract_columns = [
    "lead_phenotype_id", "lead_variant_id", "final_pattern", "classification_status",
    "pattern_tau", "pattern_delta", "pattern_method",
]
gxstage_contract = gxstage_contract_source.rename(columns={
    "lead_phenotype_id": "phenotype_id",
    "lead_variant_id": "variant_id",
    "final_pattern": "gxstage_class",
    "classification_status": "gxstage_classification_status",
})
contract_key_unique = bool(
    not gxstage_contract.duplicated(["phenotype_id", "variant_id"]).any()
)
contract_formal_class_pass = bool(
    gxstage_contract["gxstage_classification_status"].eq("formal").all()
    and gxstage_contract["gxstage_class"].isin(GXSTAGE_CLASSES).all()
)

posterior_source_columns = [
    "group_id", "lead_phenotype_id", "lead_variant_id", "posterior_pattern",
    "hc3_pattern_tau", "hc3_pattern_delta",
]
posterior_source = interaction_posterior_source[posterior_source_columns].rename(columns={
    "group_id": "posterior_interaction_group_id",
    "lead_phenotype_id": "phenotype_id",
    "lead_variant_id": "variant_id",
    "posterior_pattern": "source_posterior_pattern",
    "hc3_pattern_tau": "source_pattern_tau",
    "hc3_pattern_delta": "source_pattern_delta",
})
if posterior_source["posterior_interaction_group_id"].duplicated().any():
    raise ValueError("Interaction posterior source is not unique by group_id")
if posterior_source.duplicated(["phenotype_id", "variant_id"]).any():
    raise ValueError("Interaction posterior source is not unique by representative pair")
gxstage_source_check = (
    gxstage_contract[[
        "phenotype_id", "variant_id", "gxstage_class", "pattern_tau", "pattern_delta",
    ]]
    .merge(
        interaction_group_map[["phenotype_id", "variant_id", "interaction_group_id"]],
        on=["phenotype_id", "variant_id"],
        how="left",
        validate="one_to_one",
    )
    .merge(
        posterior_source,
        on=["phenotype_id", "variant_id"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )
)
gxstage_contract_source_match_pass = bool(
    gxstage_source_check["interaction_group_id"].notna().all()
    and gxstage_source_check["_merge"].eq("both").all()
    and gxstage_source_check["interaction_group_id"].astype(str).eq(
        gxstage_source_check["posterior_interaction_group_id"].astype(str)
    ).all()
    and gxstage_source_check["gxstage_class"].astype(str).eq(
        gxstage_source_check["source_posterior_pattern"].astype(str)
    ).all()
    and pd.to_numeric(gxstage_source_check["pattern_tau"], errors="coerce").notna().all()
    and pd.to_numeric(gxstage_source_check["source_pattern_tau"], errors="coerce").notna().all()
    and pd.to_numeric(gxstage_source_check["pattern_delta"], errors="coerce").notna().all()
    and pd.to_numeric(gxstage_source_check["source_pattern_delta"], errors="coerce").notna().all()
    and pd.Series(
        pd.to_numeric(gxstage_source_check["pattern_tau"], errors="coerce").to_numpy()
        - pd.to_numeric(gxstage_source_check["source_pattern_tau"], errors="coerce").to_numpy()
    ).abs().le(1e-12).all()
    and pd.Series(
        pd.to_numeric(gxstage_source_check["pattern_delta"], errors="coerce").to_numpy()
        - pd.to_numeric(gxstage_source_check["source_pattern_delta"], errors="coerce").to_numpy()
    ).abs().le(1e-12).all()
)

expected_interaction_rows = (
    universe[["stage", "phenotype_id", "variant_id", "signal_id"]]
    .drop_duplicates()
    .merge(
        interaction_group_map,
        on=["phenotype_id", "variant_id"],
        how="inner",
        validate="many_to_one",
    )
    .merge(
        gxstage_contract[[
            "phenotype_id", "variant_id", "gxstage_class", "gxstage_classification_status",
            "pattern_tau", "pattern_delta", "pattern_method",
        ]],
        on=["phenotype_id", "variant_id"],
        how="left",
        validate="many_to_one",
    )
)
current_group_without_pattern = expected_interaction_rows["gxstage_class"].isna()
expected_interaction_rows.loc[
    current_group_without_pattern, "gxstage_classification_status"
] = "posterior_class_not_available"
expected_interaction_rows["gxstage_class_source"] = pd.Series(
    "03_interaction_previous115_posterior_pattern_representative_pair",
    index=expected_interaction_rows.index,
)
expected_interaction_rows.loc[
    current_group_without_pattern, "gxstage_class_source"
] = "03_interaction_current_hc3_core_without_formal_posterior_pattern"

expected_interaction_signal_links = expected_interaction_rows[
    ["signal_id", *GXSTAGE_METADATA_FIELDS]
].drop_duplicates()
expected_interaction_signal_unique = not expected_interaction_signal_links["signal_id"].duplicated().any()
expected_formal_signal_links = expected_interaction_signal_links.loc[
    expected_interaction_signal_links["gxstage_classification_status"].eq("formal")
    & expected_interaction_signal_links["gxstage_class"].isin(GXSTAGE_CLASSES)
].copy()
expected_formal_signal_ids = set(expected_formal_signal_links["signal_id"].dropna().astype(str))
expected_gxstage_groups = expected_formal_signal_links[GXSTAGE_METADATA_FIELDS].drop_duplicates()
expected_gxstage_group_unique = not expected_gxstage_groups["interaction_group_id"].duplicated().any()
expected_flip_n = int(expected_gxstage_groups["gxstage_class"].eq("Opposing").sum())
expected_nonflip_n = len(expected_gxstage_groups) - expected_flip_n
expected_current_group_n = expected_interaction_rows["interaction_group_id"].nunique()
expected_unclassified_group_n = expected_interaction_rows.loc[
    expected_interaction_rows["gxstage_classification_status"].eq("posterior_class_not_available"),
    "interaction_group_id",
].nunique()

mechanism_interaction = mechanisms.loc[
    mechanisms.get("interaction_group_id", pd.Series(index=mechanisms.index, dtype=object)).notna()
].copy()
mechanism_interaction_ids = set(
    mechanism_interaction.get("signal_id", pd.Series(dtype=object)).dropna().astype(str)
)
mechanism_formal = mechanism_interaction.loc[
    mechanism_interaction.get("gxstage_classification_status", pd.Series(index=mechanism_interaction.index, dtype=object))
    .astype(str).str.strip().eq("formal")
    & mechanism_interaction.get("gxstage_class", pd.Series(index=mechanism_interaction.index, dtype=object))
    .astype(str).str.strip().isin(GXSTAGE_CLASSES)
].copy()
mechanism_formal_ids = set(mechanism_formal.get("signal_id", pd.Series(dtype=object)).dropna().astype(str))

summary_values = dict(zip(enrichment_summary.get("metric", []), enrichment_summary.get("value", [])))
expected_blocked = consequences.loc[
    ~consequences["model_status"].eq("ok"),
    ["candidate_id", "stage", "phenotype_id", "variant_id", "model_status"],
].rename(columns={"model_status": "reason"})
observed_blocked = excluded.loc[
    excluded["step"].eq("mediation_model"),
    ["candidate_id", "stage", "phenotype_id", "variant_id", "reason"],
]
blocked_match_columns = ["candidate_id", "stage", "phenotype_id", "variant_id", "reason"]
eqtl_raw_run_ids = set(eqtl_raw["pipeline_run_id"].dropna().astype(str))
eqtl_raw_version_pass = set(eqtl_raw["analysis_version"].dropna().astype(str)) == {VERSION}
eqtl_raw_candidate_coverage_pass = set(eqtl_raw["candidate_id"]) == mechanism_candidate_ids
eqtl_raw_provenance_pass = bool(
    eqtl_raw[["apa_nominal_path", "eqtl_nominal_path", "apa_nominal_bytes", "eqtl_nominal_bytes"]]
    .notna()
    .all()
    .all()
)
eqtl_raw_mechanism_run_link_pass = eqtl_raw_run_ids == mechanism_run_ids or (
    eqtl_raw_version_pass
    and eqtl_raw_candidate_coverage_pass
    and eqtl_raw_provenance_pass
)
support_external = support_source.loc[support_source["figure"].eq("external QTL support")].copy()
support_external_values = pd.to_numeric(
    support_external.get("value", pd.Series(dtype=float)),
    errors="coerce",
).fillna(0)
support_external_overlap_eligible = support_external.loc[
    support_external["metric"].astype(str).eq("eligible support rows")
    & support_external["feature"].astype(str).str.endswith(" / Same-gene APA overlap"),
    "value",
]
support_external_pqtl_eligible = support_external.loc[
    support_external["metric"].astype(str).eq("eligible support rows")
    & support_external["feature"].astype(str).eq("ROSMAP DLPFC pQTL / Same-gene protein"),
    "value",
]

checks = {
    "pipeline_version": VERSION,
    "universe_rows": len(universe),
    "universe_version_pass": set(universe["analysis_version"].dropna().astype(str)) == {VERSION},
    "universe_key_unique_pass": bool(
        not universe.duplicated(["stage", "phenotype_id", "variant_id"]).any()
    ),
    "universe_signal_assignment_pass": bool(as_bool(universe["independent_signal_assigned"]).all()),
    "interaction_group_source_key_unique_pass": group_source_key_unique,
    "gxstage_contract_key_unique_pass": contract_key_unique,
    "gxstage_contract_formal_class_pass": contract_formal_class_pass,
    "gxstage_contract_source_match_pass": gxstage_contract_source_match_pass,
    "universe_interaction_group_columns_pass": set(GXSTAGE_METADATA_FIELDS).issubset(universe.columns),
    "universe_interaction_group_accounting_pass": (
        normalized_row_counts(
            universe.loc[
                universe.get("interaction_group_id", pd.Series(index=universe.index, dtype=object)).notna()
            ],
            ["stage", "phenotype_id", "variant_id", "signal_id", *GXSTAGE_METADATA_FIELDS],
        )
        == normalized_row_counts(
            expected_interaction_rows,
            ["stage", "phenotype_id", "variant_id", "signal_id", *GXSTAGE_METADATA_FIELDS],
        )
    ),
    "current_hc3_overlap_group_count": expected_current_group_n,
    "formal_gxstage_group_count": len(expected_gxstage_groups),
    "unclassified_current_hc3_overlap_group_count": expected_unclassified_group_n,
    "mechanism_rows": len(mechanisms),
    "mechanism_version_pass": set(mechanisms["analysis_version"].dropna().astype(str)) == {VERSION},
    "gxstage_formal_group_unique_pass": expected_gxstage_group_unique,
    "mechanism_gxstage_columns_pass": set(GXSTAGE_FIELDS).issubset(mechanisms.columns),
    "mechanism_interaction_group_coverage_pass": (
        expected_interaction_signal_unique
        and mechanism_interaction_ids
        == set(expected_interaction_signal_links["signal_id"].dropna().astype(str))
        and propagated_signal_metadata_pass(
            mechanism_interaction,
            expected_interaction_signal_links,
            GXSTAGE_METADATA_FIELDS,
        )
    ),
    "mechanism_formal_gxstage_coverage_pass": (
        expected_gxstage_group_unique
        and mechanism_formal_ids == expected_formal_signal_ids
        and propagated_signal_metadata_pass(
            mechanism_formal,
            expected_formal_signal_links,
            GXSTAGE_METADATA_FIELDS,
        )
    ),
    "mechanism_variant_link_status_pass": (
        not mechanism_interaction.empty
        and mechanism_interaction["gxstage_mechanism_variant_link"].isin({
            "same_variant",
            "same_susie_credible_set",
            "interaction_variant_ld_to_selected_signal_representative",
            "same_stage_qtl_signal_pairwise_link_unresolved",
        }).all()
    ),
    "consequence_rows": len(consequences),
    "consequence_gxstage_columns_pass": set(GXSTAGE_FIELDS).issubset(consequences.columns),
    "consequence_gxstage_propagation_pass": propagated_signal_metadata_pass(
        consequences, mechanisms, GXSTAGE_FIELDS
    ),
    "mechanism_run_id_pass": mechanisms["pipeline_run_id"].nunique() == 1,
    "pas_detail_version_pass": set(pas_details["analysis_version"].dropna().astype(str)) == {VERSION},
    "pas_detail_run_link_pass": set(pas_details["pipeline_run_id"].dropna().astype(str)) == mechanism_run_ids,
    "pas_detail_candidate_coverage_pass": set(pas_details["candidate_id"]) == mechanism_candidate_ids,
    "rbp_detail_version_pass": set(rbp_details["analysis_version"].dropna().astype(str)) == {VERSION},
    "rbp_detail_run_link_pass": set(rbp_details["pipeline_run_id"].dropna().astype(str)) == mechanism_run_ids,
    "rbp_detail_candidate_coverage_pass": set(rbp_details["candidate_id"]) == local_candidate_ids,
    "eqtl_raw_version_pass": eqtl_raw_version_pass,
    "eqtl_raw_run_id_pass": eqtl_raw["pipeline_run_id"].nunique() == 1,
    "eqtl_raw_mechanism_run_link_pass": eqtl_raw_mechanism_run_link_pass,
    "eqtl_raw_candidate_coverage_pass": eqtl_raw_candidate_coverage_pass,
    "eqtl_raw_provenance_pass": eqtl_raw_provenance_pass,
    "internal_rows": len(internal),
    "internal_candidate_count": internal["candidate_id"].nunique(),
    "internal_primary_support_rows": int(internal_primary.sum()),
    "internal_primary_only_same_gene_pass": bool((~internal_primary | internal_same_gene).all()),
    "internal_standalone_consequence_zero_pass": bool(not internal_standalone.any()),
    "internal_tier_matches_primary_support_pass": bool(internal_tier.equals(internal_primary)),
    "internal_candidate_coverage_pass": internal["candidate_id"].nunique() == mechanisms["candidate_id"].nunique(),
    "internal_version_pass": set(internal["analysis_version"].dropna().astype(str)) == {VERSION},
    "internal_raw_run_link_pass": set(internal["eqtl_pipeline_run_id"].dropna().astype(str)) == set(eqtl_raw["pipeline_run_id"].dropna().astype(str)),
    "external_rows": len(external),
    "external_candidate_count": external["candidate_id"].nunique(),
    "external_source_count": source_audit["source"].nunique(),
    "external_source_contract_pass": set(source_audit["source"]) == {
        "GTEx_v11_Brain_Cortex_apaQTL",
        "GTEx_v11_Brain_BA9_apaQTL",
        "ROSMAP_DLPFC_pQTL",
    },
    "external_cache_force_rebuild_pass": set(source_audit["external_cache_mode"].dropna().astype(str)) == {"force_rebuild"},
    "external_candidate_coverage_pass": external["candidate_id"].nunique() == mechanisms["candidate_id"].nunique(),
    "external_eligible_implies_numerical_support_pass": bool((~external_eligible | external_supported).all()),
    "external_apa_eligible_overlap_pass": bool(
        external.loc[eligible_apa, "event_match_scope"].eq("same_gene_apa_overlap").all()
    ),
    "external_pqtl_eligible_same_gene_pass": bool(
        external.loc[eligible_pqtl, "event_match_scope"].eq("same_gene_protein").all()
    ),
    "external_other_apa_context_not_eligible_pass": bool(
        not external.loc[
            external["event_match_scope"].eq("same_gene_other_event"),
            "external_validation_eligible",
        ].pipe(as_bool).any()
    ),
    "external_summary_candidate_source_unique_pass": bool(
        not external_summary.duplicated(["candidate_id", "source"]).any()
    ),
    "external_summary_eligibility_consistency_pass": bool(
        candidate_eligible.sort_index().equals(summary_eligible.sort_index())
    ),
    "external_version_pass": set(external["analysis_version"].dropna().astype(str)) == {VERSION},
    "same_gene_apa_overlap_test_rows": int(external["event_match_scope"].eq("same_gene_apa_overlap").sum()),
    "same_gene_apa_overlap_eligible_rows": int(eligible_apa.sum()),
    "same_gene_other_event_context_rows": int(external["event_match_scope"].eq("same_gene_other_event").sum()),
    "pqtl_eligible_rows": int(eligible_pqtl.sum()),
    "published_rows": len(published),
    "published_version_pass": set(published["analysis_version"].dropna().astype(str)) == {VERSION},
    "published_total_score_absent_pass": "TotalScore" not in published.columns,
    "coloc_file_count": len(coloc_tables),
    "coloc_nonempty_pass": all(not frame.empty for frame in coloc_tables.values()),
    "coloc_version_pass": all(
        "analysis_version" in frame.columns
        and set(frame["analysis_version"].dropna().astype(str)) == {VERSION}
        for frame in coloc_tables.values()
    ),
    "disease_rows": len(disease),
    "disease_ids_unique_pass": bool(not disease["disease_comparison_id"].duplicated().any()),
    "disease_posterior_sum_pass": bool(posterior_sum.sub(1).abs().le(1e-6).all()),
    "disease_candidate_link_rows": int(disease["candidate_id"].notna().sum()),
    "disease_prior_robust_rows": int(as_bool(disease["coloc_prior_sensitivity_acceptable"]).sum()),
    "disease_version_pass": set(disease["analysis_version"].dropna().astype(str)) == {VERSION},
    "disease_gxstage_columns_pass": set(GXSTAGE_FIELDS).issubset(disease.columns),
    "disease_gxstage_propagation_pass": propagated_signal_metadata_pass(
        disease, mechanisms, GXSTAGE_FIELDS
    ),
    "tier_rows": len(tiers),
    "tier_version_pass": set(tiers["analysis_version"].dropna().astype(str)) == {VERSION},
    "tier_row_count_pass": len(tiers) == 321097,
    "tier_universe_row_match_pass": len(tiers) == len(universe),
    "excluded_version_pass": set(excluded["analysis_version"].dropna().astype(str)) <= {VERSION},
    "excluded_schema_pass": {"stage", "phenotype_id", "variant_id", "step", "reason"}.issubset(excluded.columns),
    "excluded_mediation_model_accounting_pass": (
        normalized_row_counts(expected_blocked, blocked_match_columns)
        == normalized_row_counts(observed_blocked, blocked_match_columns)
        and excluded.loc[excluded["step"].eq("mediation_model"), "severity"].eq("blocked").all()
        and excluded.loc[excluded["step"].eq("mediation_model"), "source_row_id"].fillna("").eq("").all()
    ),
    "enrichment_rows": len(enrichment),
    "enrichment_version_pass": set(enrichment["analysis_version"].dropna().astype(str)) == {VERSION},
    "gxstage_test_interaction_group_unit_pass": (
        not gxstage_tests.empty
        and {"analysis_unit", "flip_n", "nonflip_n"}.issubset(gxstage_tests.columns)
        and gxstage_tests["analysis_unit"].eq("interaction_group").all()
    ),
    "gxstage_test_group_counts_pass": (
        not gxstage_tests.empty
        and {"flip_n", "nonflip_n"}.issubset(gxstage_tests.columns)
        and pd.to_numeric(gxstage_tests["flip_n"], errors="coerce").eq(expected_flip_n).all()
        and pd.to_numeric(gxstage_tests["nonflip_n"], errors="coerce").eq(expected_nonflip_n).all()
    ),
    "enrichment_summary_gxstage_counts_pass": (
        pd.to_numeric(pd.Series([summary_values.get("gxstage_group_n")]), errors="coerce").eq(len(expected_gxstage_groups)).all()
        and pd.to_numeric(pd.Series([summary_values.get("flip_group_n")]), errors="coerce").eq(expected_flip_n).all()
        and pd.to_numeric(pd.Series([summary_values.get("nonflip_group_n")]), errors="coerce").eq(expected_nonflip_n).all()
    ),
    "internal_config_primary_support_pass": (
        internal_config.get("primary_functional_support") == "exact_same_gene_PP4_gte_0.75"
        and internal_config.get("tier_eligible") == "exact_same_gene_primary_support_only"
        and internal_config.get("standalone_consequence_eligible") is False
        and internal_config.get("independent_replication") is False
    ),
    "case_figure_rows": len(case_layout),
    "case_figure_set_pass": set(case_layout["gene_name"].astype(str)) == set(case_support_expected),
    "case_figure_support_contract_pass": case_support_actual == case_support_expected,
    "case_layout_pass": bool(
        case_layout["automated_layout_status"].eq("pass").all()
        and as_bool(case_layout["files_exist"]).all()
        and case_layout["text_overlap_count"].eq(0).all()
        and case_layout["tick_label_overlap_count"].eq(0).all()
        and case_layout["clipped_text_count"].eq(0).all()
    ),
    "case_figure_files_pass": all(path.exists() and path.stat().st_size > 0 for path in case_figure_paths),
    "case_source_contract_pass": (
        set(case_source["association_analysis_version"].dropna().astype(str)) == {VERSION}
        and set(case_source["consequence_analysis_version"].dropna().astype(str)) == {VERSION}
        and set(case_source["internal_eqtl_analysis_version"].dropna().astype(str)) == {VERSION}
        and set(case_source["external_analysis_version"].dropna().astype(str)) == {VERSION}
        and set(case_source["pas_detail_analysis_version"].dropna().astype(str)) == {VERSION}
        and set(case_source["rbp_detail_analysis_version"].dropna().astype(str)) == {VERSION}
        and set(case_source["mirna_analysis_version"].dropna().astype(str)) == {VERSION}
    ),
    "support_figure_rows": len(support_layout),
    "support_figure_set_pass": set(support_layout["figure"].astype(str)) == support_figure_expected,
    "support_layout_pass": bool(
        support_layout["automated_status"].eq("pass").all()
        and support_layout["missing_plot_values"].eq(0).all()
        and support_layout["label_overlap_count"].eq(0).all()
        and support_layout["clipping_status"].eq("pass").all()
    ),
    "support_figure_files_pass": all(path.exists() and path.stat().st_size > 0 for path in support_figure_paths),
    "support_source_contract_pass": (
        not support_source.empty
        and set(support_source["figure_contract_version"].dropna().astype(str)) == {VERSION}
        and set(support_source["figure"].astype(str)) == {
            "cell reference reconstruction",
            "cell sample composition",
            "cell sample fit",
            "cell covariate capture",
            "cell QTL residual signal",
            "internal eQTL sharing",
            "external QTL support",
            "sLDSC disease enrichment",
        }
    ),
    "support_external_qtl_support_pass": (
        not support_external.empty
        and support_external_values.notna().all()
        and int(pd.to_numeric(support_external_overlap_eligible, errors="coerce").fillna(0).sum())
        == int(eligible_apa.sum())
        and int(pd.to_numeric(support_external_pqtl_eligible, errors="coerce").fillna(0).sum())
        == int(eligible_pqtl.sum())
    ),
    "support_sldsc_76_unique_tests_pass": (
        len(support_source.loc[support_source["figure"].eq("sLDSC disease enrichment")]) == 76
        and support_source.loc[
            support_source["figure"].eq("sLDSC disease enrichment"), "feature"
        ].nunique() == 76
    ),
}

failed = []
for name, value in checks.items():
    if name.endswith("_pass"):
        try:
            passed = bool(value)
        except (TypeError, ValueError):
            passed = False
        if not passed:
            failed.append(name)
if failed:
    raise SystemExit("Current-v5 result validation failed: " + ", ".join(failed))

validation = pd.DataFrame([checks])
validation.insert(0, "analysis_version", VERSION)
validation.insert(1, "finalizer_version", FINALIZER_VERSION)
validation_path = RESULTS / "00_validation.tsv"
temporary = validation_path.with_suffix(".tsv.tmp")
validation.to_csv(temporary, sep="\t", index=False)
temporary.replace(validation_path)

# --- manifest ---

evidence_paths = sorted(
    (
        path for path in CLUSTER_WORK.rglob("*")
        if path.is_file()
        and "cache" not in path.relative_to(CLUSTER_WORK).parts
        and "ld" not in path.relative_to(CLUSTER_WORK).parts
        and not path.name.endswith(".tmp")
        and not path.name.endswith((".err", ".out"))
    ),
    key=str,
)

manifest_groups = {
    "contract": [
        ROOT / "data/notes/framework.md",
        FUNC_DIR / "config.yml",
    ],
    "workflow": [
        FUNC_DIR / "run.sh",
        FUNC_DIR / "01_prepare.py",
        FUNC_DIR / "02_mechanism.py",
        FUNC_DIR / "03_deconvolution.R",
        FUNC_DIR / "04_consequence.py",
        FUNC_DIR / "05_coloc.R",
        FUNC_DIR / "06_disease.py",
        FUNC_DIR / "07_enrichment.R",
        FUNC_DIR / "08_plot_cases.py",
        FUNC_DIR / "09_plot.R",
        FUNC_DIR / "10_finalize.py",
    ],
    "result": [
        RESULTS / "00_validation.tsv",
        RESULTS / "01_significant_qtl_universe.tsv.gz",
        RESULTS / "02_signal_mechanisms.tsv.gz",
        RESULTS / "03_event_consequences.tsv.gz",
        RESULTS / "04_disease_chains_results.tsv.gz",
        RESULTS / "05_evidence_tiers_results.tsv.gz",
        RESULTS / "05_excluded_or_blocked.tsv.gz",
        RESULTS / "06_functional_enrichment.tsv.gz",
        RESULTS / "07_internal_eqtl_integration_contract.tsv.gz",
    ],
    "final_figure": [*case_figure_paths, *support_figure_paths],
    "evidence": evidence_paths,
}

expected_results = {
    Path(path).name
    for role in ("result", "final_figure")
    for path in manifest_groups[role]
}
present_results = {path.name for path in RESULTS.iterdir() if path.is_file()}
unexpected_results = sorted(present_results - expected_results - {"00_manifest.tsv"})
missing_results = sorted(expected_results - present_results)
if unexpected_results or missing_results:
    raise SystemExit(
        f"Result directory contract failed; missing={missing_results}; unexpected={unexpected_results}"
    )
manifest_rows = []
for role, paths in manifest_groups.items():
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            raise SystemExit(f"Final artifact missing or empty: {path}")
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        manifest_rows.append({
            "role": role,
            "path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
            "bytes": path.stat().st_size,
            "sha256": digest.hexdigest(),
            "analysis_version": VERSION,
            "status": "pass",
        })
manifest = pd.DataFrame(manifest_rows).drop_duplicates("path")
manifest_path = RESULTS / "00_manifest.tsv"
temporary_manifest = manifest_path.with_suffix(".tsv.tmp")
manifest.to_csv(temporary_manifest, sep="\t", index=False)
temporary_manifest.replace(manifest_path)

print("Validated canonical final package")
print(validation.to_string(index=False))
print(f"Written {len(manifest)} manifest rows -> {manifest_path}")
