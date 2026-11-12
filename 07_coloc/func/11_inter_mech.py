#!/usr/bin/env python3
# Score interaction-group variant events for PAS, RBP, RNA structure, APARENT and CLIP.
# Writes the variant, group and coverage tables and the coverage figure.

# --- packages ---
from __future__ import annotations
import ast
from bisect import bisect_right
from collections import Counter
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import numpy as np
import pandas as pd
import pyBigWig
import pysam
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --- helpers ---

def _expand_env(obj):
    """Expand ${VAR} in config values (paths set via config/paths.sh)."""
    if isinstance(obj, dict):
        return {k: _expand_env(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env(v) for v in obj]
    return os.path.expandvars(obj) if isinstance(obj, str) else obj


# --- config ---
ROOT = Path(os.environ["APA_ROOT"])
RESULTS = ROOT / "data/07_coloc"
FIGURES = ROOT / "figures/07_coloc"
CLUSTER = Path(os.environ["APA_HPC_HOME"]) / "project/work/apa/05_func"
SOURCE = Path(__file__).resolve().parent / "02_mechanism.py"
SOURCE_SHA256 = "d95e12ed09b19175434d2485d9de725c75cd8ee48f31efaabdc13701e08d1a3d"
PARENT_SOURCE_SHA256 = "13f5f382ebc99c44607e87d2bf8ae7069b0b47d0de9458aa7b47fa5bb400446c"
WITHIN_GROUP_P = 0.05
RBP_DELTA_THRESHOLD = 0.20
POSTAR_RBP_LIMIT = 20
POSTAR_CHUNK_SIZE = 250000
APARENT_MODEL_NAME = "aparent_large_lessdropout_all_libs_no_sampleweights.h5"
APARENT_WINDOW_HALF = 102
RBPDB_HUMAN_PROTEINS = ROOT / "ref/RBPDB_v1.3.1_proteins_human_2012-11-21.csv"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42})

# --- reuse mechanism functions from 02_mechanism.py ---
source_tree = ast.parse(SOURCE.read_text())
functions = [node for node in source_tree.body if isinstance(node, ast.FunctionDef)]
function_module = ast.Module(body=[ast.ImportFrom(module="__future__",
    names=[ast.alias(name="annotations")], level=0), *functions], type_ignores=[])
exec(compile(ast.fix_missing_locations(function_module), str(SOURCE), "exec"), globals())

with (Path(__file__).resolve().parent / "config.yml").open() as handle:
    config = _expand_env(yaml.safe_load(handle))
core = pd.read_csv(ROOT / "data/05_gxs/intx_hc3_core_groups.tsv", sep="\t")
assert core.group_id.is_unique and len(core) == 193
assert core.stage_group_status.eq("Interaction-only group").sum() == 77
anchors = pd.read_csv(config["inputs"]["apa_phenotype_annotation"]["path"], sep="\t",
                      usecols=["phenotype_id", "Predicted_Proximal_APA"]).drop_duplicates("phenotype_id")
group_columns = ["group_id", "lead_phenotype_id", "lead_variant_id", "stage_group_status",
                 "tests_emt", "num_phenotypes", "beta_pre", "beta_post", "strict_ci_direction_reversal"]

# --- all HC3 leads and within-group supported variant-event candidates ---
parts = []
pair_file = CLUSTER.parent / "03_interaction/results/02_group_eigenmt/interaction_pair_summary.tsv.gz"
for chunk in pd.read_csv(pair_file, sep="\t", chunksize=250000,
                         usecols=["group_id", "phenotype_id", "variant_id", "pval_gi", "b_g", "b_gi", "af"]):
    block = chunk.merge(core[group_columns], on="group_id", how="inner", validate="many_to_one")
    block["within_group_eigenmt_p"] = np.minimum(1, block.pval_gi*block.tests_emt*block.num_phenotypes)
    block["is_hc3_lead_pair"] = block.phenotype_id.eq(block.lead_phenotype_id) & block.variant_id.eq(block.lead_variant_id)
    parts.append(block.loc[block.within_group_eigenmt_p.le(WITHIN_GROUP_P) | block.is_hc3_lead_pair])
table = pd.concat(parts, ignore_index=True).drop_duplicates(["group_id", "phenotype_id", "variant_id"])
assert table.is_hc3_lead_pair.sum() == 193
table["matches_lead_event"] = table.phenotype_id.eq(table.lead_phenotype_id)
phenotype_parts = table.phenotype_id.str.split("|", expand=True, regex=False)
variant_parts = table.variant_id.str.split(":", expand=True)
coordinates = phenotype_parts[2].str.extract(r"(?P<phenotype_chr>[^:]+):(?P<phenotype_start>\d+)-(?P<phenotype_end>\d+)")
table = pd.concat([table, coordinates], axis=1)
table[["phenotype_start", "phenotype_end"]] = table[["phenotype_start", "phenotype_end"]].astype(int)
table["phenotype_strand"] = phenotype_parts[3]
table["variant_chr"] = variant_parts[0]
table["variant_pos"] = variant_parts[1].astype(int)
table["variant_ref"] = variant_parts[2]
table["variant_alt"] = variant_parts[3]
table["candidate_id"] = "HC3::"+table.group_id+"::"+table.phenotype_id+"::"+table.variant_id
table = table.merge(anchors, on="phenotype_id", how="left", validate="many_to_one")
table["proximal_pas_pos"] = pd.to_numeric(table.Predicted_Proximal_APA, errors="coerce")
table["distal_pas_pos"] = np.where(table.phenotype_strand.eq("+"), table.phenotype_end, table.phenotype_start)
table["pas_anchor_resolution"] = "proximal_pipeline_distal_boundary_approximation"
table["pip"] = np.nan
table["candidate_role"] = np.where(table.is_hc3_lead_pair, "HC3_lead",
    "exploratory_within_group_eigenMT_P_le_0.05_not_finemapped")
fasta = pysam.FastaFile(config["reference_resources"]["genome_fasta"]["path"])
table["allele_harmonized"] = [fasta.fetch(row.variant_chr, row.variant_pos-1,
    row.variant_pos-1+len(row.variant_ref)).upper() == row.variant_ref
    for row in table.itertuples(index=False)]

# --- direct PAS and eligibility under the existing anchor definitions ---
pas_rows = []
motifs = set(config["parameters"]["pas_hexamers"])
for row in table.to_dict("records"):
    series = pd.Series(row)
    direct = []
    for label in ["proximal", "distal"]:
        position = series[f"{label}_pas_pos"]
        if pd.notna(position) and series.allele_harmonized:
            direct.append(direct_pas_change(series, int(position), label, fasta, motifs))
        else:
            direct.append({"pas_core_change": "not_scored", "variant_relative_to_pas_nt": np.nan,
                "pas_name": label, "predicted_pdui_direction": "not_available"})
    pas_rows.append({"candidate_id": series.candidate_id, **collapse_pas_results(direct),
                     "pas_status": "scored_available_anchors" if series.allele_harmonized else "reference_conflict"})
table = table.merge(pd.DataFrame(pas_rows), on="candidate_id", validate="one_to_one")
table["local_regulatory_eligible"] = False
for column in ["proximal_variant_distance_nt", "distal_variant_distance_nt"]:
    distance = pd.to_numeric(table[column], errors="coerce")
    table["local_regulatory_eligible"] |= distance.between(-100, -1) | distance.between(10, 50)
local = table.loc[table.local_regulatory_eligible & table.allele_harmonized].copy()

# --- motif scores with separate prenatal and postnatal expression support ---
rbp_motifs = parse_cisbp_meme(Path(config["reference_resources"]["cisbp_rna_motifs"]["path"]))
expression = load_rbp_expression(Path(config["reference_resources"]["pfc_rbp_expression"]["path"]))
assert expression is not None
for stage in ["prenatal", "postnatal"]:
    rows = [score_local_rbp(pd.Series(row), fasta, rbp_motifs, expression,
                            min_mean_tpm=1.0, stage=stage) for row in local.to_dict("records")]
    scored = pd.DataFrame(rows)
    if len(scored):
        scored = scored.rename(columns={c: f"{stage}_{c}" for c in scored if c != "candidate_id"})
        table = table.merge(scored, on="candidate_id", how="left", validate="one_to_one")
table["rbp_expression_supported"] = False
for stage in ["prenatal", "postnatal"]:
    column = f"{stage}_v2_rbp_variant_delta_support"
    if column in table:
        table["rbp_expression_supported"] |= table[column].fillna(False).astype(bool)

# --- RNAfold, APARENT, conservation and CLIP ---
structure = score_local_structure(local, fasta,
    float(config["parameters"]["rna_structure"]["local_delta_mfe_threshold_kcal_mol"]))
table = table.merge(structure, on="candidate_id", how="left", validate="one_to_one")
table["local_structure_change"] = table.local_structure_change.fillna(False).astype(bool)
near_aparent = table.allele_harmonized & table.variant_ref.str.len().eq(1) & table.variant_alt.str.len().eq(1)
distance_pre = (table.variant_pos-table.proximal_pas_pos).abs()
distance_post = (table.variant_pos-table.distal_pas_pos).abs()
near_aparent &= distance_pre.le(APARENT_WINDOW_HALF) | distance_post.le(APARENT_WINDOW_HALF)
aparent = score_aparent_variants(table.loc[near_aparent], fasta,
                               Path(config["reference_resources"]["aparent_model"]["path"]))
if len(aparent):
    table = table.merge(aparent, on="candidate_id", how="left", validate="one_to_one")
conservation = annotate_conservation(table,
    Path(config["reference_resources"]["phastcons_bigwig"]["path"]))
table = table.merge(conservation, on="candidate_id", how="left", validate="one_to_one")
# --- recheck the original stage candidates in the same corrected CLIP scan ---
stage_columns = ["candidate_id", "signal_id", "independent_signal_assigned", "allele_harmonized",
    "phenotype_id", "phenotype_chr", "phenotype_start", "phenotype_end", "phenotype_strand", "variant_pos"]
stage = pd.read_csv(RESULTS / "01_significant_qtl_universe.tsv.gz", sep="\t", usecols=stage_columns)
stage = stage.loc[stage.independent_signal_assigned.fillna(False) & stage.allele_harmonized.fillna(False)]
stage = stage.merge(anchors, on="phenotype_id", how="left", validate="many_to_one")
stage["proximal_pas_pos"] = pd.to_numeric(stage.Predicted_Proximal_APA, errors="coerce")
stage["distal_pas_pos"] = np.where(stage.phenotype_strand.eq("+"), stage.phenotype_end, stage.phenotype_start)
sign = np.where(stage.phenotype_strand.eq("+"), 1, -1)
stage_local = pd.Series(False, index=stage.index)
for label in ["proximal", "distal"]:
    distance = (stage.variant_pos-stage[f"{label}_pas_pos"])*sign
    stage_local |= distance.between(-100, -1) | distance.between(10, 50)
stage = stage.loc[stage_local & stage.proximal_pas_pos.notna()].copy()
clip = stream_postar_overlaps(pd.concat([local, stage], ignore_index=True),
    Path(config["reference_resources"]["postar3_clip_peaks"]["path"])/"ALL_RBP_binding_sites_merged.csv")
stage_clip = stage[["candidate_id", "signal_id", "phenotype_id", "variant_pos"]].merge(
    clip, on="candidate_id", validate="one_to_one")
stage_clip.to_csv(RESULTS/"11_stage_clip_recheck.tsv", sep="\t", index=False)
table = table.merge(clip, on="candidate_id", how="left", validate="one_to_one")
table["clip_overlap"] = table.postar_variant_overlap_count.fillna(0).gt(0)
table["mechanism_source_sha256"] = SOURCE_SHA256
table["mechanism_parent_sha256"] = PARENT_SOURCE_SHA256
table["interpretation"] = "sequence_prediction_or_CLIP_overlap; no interaction credible set or disease mediation"
# --- group-level coverage while retaining exact-event distinctions ---
table["aparent_scored"] = table.APARENT_predicted_PDUI_delta.notna()
table.to_csv(RESULTS/"11_interaction_variant_mechanisms.tsv.gz", sep="\t", index=False)
evidence = ["local_regulatory_eligible", "aparent_scored", "direct_pas_evidence",
            "rbp_expression_supported", "local_structure_change", "clip_overlap"]
summary = core[group_columns].copy()
for scope, subset in [("all_group_events", table), ("lead_event", table.loc[table.matches_lead_event]),
                       ("lead_variant", table.loc[table.is_hc3_lead_pair])]:
    grouped = subset.groupby("group_id")[evidence].max().add_prefix(scope+"_")
    counts = subset.groupby("group_id").size().rename(scope+"_n_candidates")
    summary = summary.merge(grouped, on="group_id", how="left", validate="one_to_one")
    summary = summary.merge(counts, on="group_id", how="left", validate="one_to_one")
summary.to_csv(RESULTS/"11_interaction_group_mechanisms.tsv", sep="\t", index=False)
coverage = []
for status, block in summary.groupby("stage_group_status"):
    for scope in ["all_group_events", "lead_event", "lead_variant"]:
        for feature in evidence:
            hits = int(block[scope+"_"+feature].fillna(False).sum())
            coverage.append({"detection_status": status, "scope": scope, "evidence": feature,
                             "measure_role": "model_coverage" if feature in evidence[:2] else "sequence_evidence",
                             "n_groups": len(block), "n_supported": hits, "fraction": hits/len(block)})
coverage = pd.DataFrame(coverage)
coverage.to_csv(RESULTS/"11_interaction_mechanism_coverage.tsv", sep="\t", index=False)

# --- figure: exact lead-event annotation coverage by detection status ---
order = ["Interaction-only group", "Prenatal detected group", "Postnatal detected group", "Shared group"]
labels = ["Interaction only (77)", "Prenatal (15)", "Postnatal (59)", "Shared (42)"]
fig, axes = plt.subplots(1, 2, figsize=(7.8, 2.8), width_ratios=[2, 4], layout="constrained", sharey=True)
lead_coverage = coverage.loc[coverage.scope.eq("lead_event")]
for ax, features, titles, title in [
    (axes[0], evidence[:2], ["PAS local", "APARENT"], "Model coverage"),
    (axes[1], evidence[2:], ["PAS", "RBP motif", "RNAfold", "CLIP"], "Sequence predictions / overlap")]:
    fractions = lead_coverage.pivot(index="detection_status", columns="evidence", values="fraction").loc[order, features]
    counts = lead_coverage.pivot(index="detection_status", columns="evidence", values="n_supported").loc[order, features]
    rendered = ax.imshow(fractions.to_numpy()*100, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(np.arange(len(features)), titles)
    ax.set_yticks(np.arange(4), labels)
    ax.tick_params(length=0)
    for row_index in range(4):
        for column_index in range(len(features)):
            ax.text(column_index, row_index, str(counts.iloc[row_index, column_index]),
                    ha="center", va="center", fontsize=8)
    ax.set_title(title)
    for spine in ax.spines.values():
        spine.set_visible(False)
colorbar = fig.colorbar(rendered, ax=axes, shrink=.85, pad=.035)
colorbar.set_label("HC3 groups (%)")
for extension in ["png", "pdf"]:
    fig.savefig(FIGURES/f"11_interaction_mechanism_coverage.{extension}", dpi=300)
