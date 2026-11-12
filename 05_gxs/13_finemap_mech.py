# Mechanism annotation (PAS, RBP motifs, RNAfold, APARENT, conservation, CLIP) of the usable
# cohort-LD SuSiE credible-set variants of the G x Stage groups.
# Reuses the mechanism functions of the cluster 05_func 02_mechanism_v2.py (executed from its AST).
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
import tempfile

import numpy as np
import pandas as pd
import pyBigWig
import pysam
import yaml

# --- paths and settings ---
APA_ROOT = Path(os.environ["APA_ROOT"])
D = APA_ROOT / "data" / "05_gxs"
ROOT = Path(f"{os.environ['APA_HPC_HOME']}/project/work/apa/05_func")
SOURCE = ROOT / "scripts/02_mechanism_v2.py"
SOURCE_SHA256 = "d95e12ed09b19175434d2485d9de725c75cd8ee48f31efaabdc13701e08d1a3d"
POSTAR_CHUNK_SIZE = 250000
POSTAR_RBP_LIMIT = 20
APARENT_WINDOW_HALF = 102
APARENT_MODEL_NAME = "aparent_large_lessdropout_all_libs_no_sampleweights.h5"
RBP_DELTA_THRESHOLD = 0.20
RBPDB_HUMAN_PROTEINS = ROOT / "ref/RBPDB_v1.3.1_proteins_human_2012-11-21.csv"

# --- reuse the mechanism functions from 02_mechanism_v2.py (function definitions only) ---
source_tree = ast.parse(SOURCE.read_text())
functions = [node for node in source_tree.body if isinstance(node, ast.FunctionDef)]
function_module = ast.Module(body=[ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0),
                                   *functions], type_ignores=[])
exec(compile(ast.fix_missing_locations(function_module), str(SOURCE), "exec"), globals())

with (ROOT / "config.yml").open() as handle:
    config = yaml.safe_load(handle)
anchors = pd.read_csv(config["inputs"]["apa_phenotype_annotation"]["path"], sep="\t",
                      usecols=["phenotype_id", "Predicted_Proximal_APA"]).drop_duplicates("phenotype_id")

# --- POSTAR chunk-boundary check (overlaps must not be lost or doubled across chunks) ---
with tempfile.TemporaryDirectory() as check_dir:
    check_file = Path(check_dir) / "clip.csv"
    check_columns = [f"column_{i}" for i in range(22)]
    for position, name in [(0, "RBP"), (4, "l1"), (20, "RBP_Query"), (21, "CLIP_Technology")]:
        check_columns[position] = name
    check_rows = []
    for coordinate in [None, "chr1:10-20", "chr1:12-18"]:
        row = [""] * 22
        row[0], row[4], row[20], row[21] = "PUM2", coordinate, "PUM2", "eCLIP"
        check_rows.append(row)
    pd.DataFrame(check_rows, columns=check_columns).to_csv(check_file, index=False)
    POSTAR_CHUNK_SIZE = 1
    checked = stream_postar_overlaps(pd.DataFrame([dict(candidate_id="check", phenotype_chr="chr1",
                                                        phenotype_start=1, phenotype_end=30, variant_pos=15)]),
                                     check_file)
    assert checked.postar_variant_overlap_count.iloc[0] == 2
    POSTAR_CHUNK_SIZE = 250000
print("POSTAR_CHUNK_CHECK PASS", flush=True)

# --- usable cohort-LD credible-set candidates ---
runs = pd.read_csv(D / "susie_run_summary.tsv", sep="\t")
usable = runs.loc[runs.usable, ["group_id", "stage"]]
pips = pd.read_csv(D / "susie_variant_pips.tsv", sep="\t")
table = pips.loc[pips.cs_id.gt(0)].merge(usable, on=["group_id", "stage"], how="inner", validate="many_to_one")
assert table.fit_usable.all()
table["candidate_id"] = ("GXS_SUSIE_V2::" + table.group_id + "::" + table.stage
                         + "::" + table.phenotype_id + "::" + table.variant_id)
assert table.candidate_id.is_unique
table["candidate_role"] = "usable_cohort_ld_credible_set_member"
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
table = table.merge(anchors, on="phenotype_id", how="left", validate="many_to_one")
table["proximal_pas_pos"] = pd.to_numeric(table.Predicted_Proximal_APA, errors="coerce")
table["distal_pas_pos"] = np.where(table.phenotype_strand.eq("+"), table.phenotype_end, table.phenotype_start)
table["pas_anchor_resolution"] = "proximal_pipeline_distal_boundary_approximation"
fasta = pysam.FastaFile(config["reference_resources"]["genome_fasta"]["path"])
table["allele_harmonized"] = [fasta.fetch(row.variant_chr, row.variant_pos - 1,
                                          row.variant_pos - 1 + len(row.variant_ref)).upper() == row.variant_ref
                              for row in table.itertuples(index=False)]
print(f"CANDIDATES {len(table)} credible-set variants; {table.group_id.nunique()} groups; "
      f"reference conflicts {(~table.allele_harmonized).sum()}", flush=True)

# --- direct PAS change and local eligibility ---
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

# --- RBP motif scores with separate prenatal and postnatal expression support ---
rbp_motifs = parse_cisbp_meme(Path(config["reference_resources"]["cisbp_rna_motifs"]["path"]))
expression = load_rbp_expression(Path(config["reference_resources"]["pfc_rbp_expression"]["path"]))
for stage in ["prenatal", "postnatal"]:
    rows = [score_local_rbp(pd.Series(row), fasta, rbp_motifs, expression, min_mean_tpm=1.0, stage=stage)
            for row in local.to_dict("records")]
    scored = pd.DataFrame(rows)
    if len(scored):
        scored = scored.rename(columns={c: f"{stage}_{c}" for c in scored if c != "candidate_id"})
        table = table.merge(scored, on="candidate_id", how="left", validate="one_to_one")
    print(f"RBP {stage}: {len(rows)} local candidates", flush=True)
table["rbp_expression_supported"] = False
for stage in ["prenatal", "postnatal"]:
    column = f"{stage}_v2_rbp_variant_delta_support"
    if column in table:
        table["rbp_expression_supported"] |= table[column].fillna(False).astype(bool)

# --- RNAfold, APARENT, conservation and CLIP ---
structure = score_local_structure(local, fasta, float(config["parameters"]["rna_structure"]["local_delta_mfe_threshold_kcal_mol"]))
table = table.merge(structure, on="candidate_id", how="left", validate="one_to_one")
table["local_structure_change"] = table.local_structure_change.fillna(False).astype(bool)
near_aparent = table.allele_harmonized & table.variant_ref.str.len().eq(1) & table.variant_alt.str.len().eq(1)
distance_pre = (table.variant_pos - table.proximal_pas_pos).abs()
distance_post = (table.variant_pos - table.distal_pas_pos).abs()
near_aparent &= distance_pre.le(APARENT_WINDOW_HALF) | distance_post.le(APARENT_WINDOW_HALF)
aparent = score_aparent_variants(table.loc[near_aparent], fasta, Path(config["reference_resources"]["aparent_model"]["path"]))
if len(aparent):
    table = table.merge(aparent, on="candidate_id", how="left", validate="one_to_one")
conservation = annotate_conservation(table, Path(config["reference_resources"]["phastcons_bigwig"]["path"]))
table = table.merge(conservation, on="candidate_id", how="left", validate="one_to_one")
clip = stream_postar_overlaps(local, Path(config["reference_resources"]["postar3_clip_peaks"]["path"]) / "ALL_RBP_binding_sites_merged.csv")
table = table.merge(clip, on="candidate_id", how="left", validate="one_to_one")
table["clip_overlap"] = table.postar_variant_overlap_count.fillna(0).gt(0)
table["mechanism_source_sha256"] = SOURCE_SHA256
table["aparent_scored"] = table.APARENT_predicted_PDUI_delta.notna()

# --- save ---
table.to_csv(D / "susie_group_variant_mechanisms.tsv.gz", sep="\t", index=False)
evidence = ["local_regulatory_eligible", "aparent_scored", "direct_pas_evidence",
            "rbp_expression_supported", "local_structure_change", "clip_overlap"]
print(table.groupby("stage")[evidence].sum(), flush=True)
print("RNAfold:", subprocess.check_output(["RNAfold", "--version"], text=True).strip(), flush=True)
print("COMPLETE: usable cohort-LD credible-set mechanism annotation", flush=True)
