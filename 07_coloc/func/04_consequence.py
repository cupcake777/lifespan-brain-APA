#!/usr/bin/env python3
# Consequence opportunities per APA event (miRNA seed, RBP, ARE, GRE, URE).
# Signal-level mediation of genotype effects through APA, per stage and outcome.

from __future__ import annotations

import csv
import gzip
import math
import multiprocessing
import os
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import statsmodels.api as sm
import yaml
from scipy.linalg import LinAlgError, lstsq
from scipy.stats import norm, pearsonr
from statsmodels.stats.multitest import multipletests


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
CONFIG_PATH = Path(__file__).resolve().parent / "config.yml"
RESULTS = ROOT / "data" / "07_coloc"
UNIVERSE = RESULTS / "01_significant_qtl_universe.tsv.gz"
SIGNALS = RESULTS / "02_signal_mechanisms.tsv.gz"
EVENT_ANNOTATIONS = RESULTS / "annotations.tsv.gz"
SUMMARY = RESULTS / "summary.tsv"
OUTPUT = RESULTS / "03_event_consequences.tsv.gz"
EXCLUDED_OUTPUT = RESULTS / "05_excluded_or_blocked.tsv.gz"
ANALYSIS_VERSION = "framework_2026-07-25_v5"
GRE_MOTIFS = ("GTGTG", "TGTGT", "GTGTT", "TGTTT")
URE_MOTIFS = ("TTTTT", "TTTGT", "TTTCT", "TGTGT")
GXSTAGE_FIELDS = [
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
    "interaction_pair_signal_assignment_method",
    "interaction_pair_signal_lead_variant",
    "interaction_pair_ld_r2_to_signal_lead",
    "interaction_pair_cs_id",
    "gxstage_mechanism_variant_link",
    "gxstage_mechanism_variant_link_r2",
]
EXCLUDED_FIELDS = [
    "candidate_id",
    "source_row_id",
    "stage",
    "phenotype_id",
    "variant_id",
    "step",
    "reason",
    "severity",
    "analysis_version",
    "config_path",
]


# --- helpers ---

def _expand_env(obj):
    # expand ${VAR} in config values
    if isinstance(obj, dict):
        return {k: _expand_env(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env(v) for v in obj]
    return os.path.expandvars(obj) if isinstance(obj, str) else obj


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def overlapping_count(sequence: str, motif: str) -> int:
    count = 0
    start = 0
    while True:
        index = sequence.find(motif, start)
        if index < 0:
            return count
        count += 1
        start = index + 1


def normalize_sample(sample: str) -> str:
    return str(sample).removeprefix("Intensity ").replace("_", "-")


def load_selected_matrix(
    path: Path,
    wanted: set[str],
    key_column: int,
    value_start: int,
    strip_key_version: bool = False,
) -> tuple[list[str], dict[str, pd.Series]]:
    rows = {}
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open(encoding="utf-8")) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = [normalize_sample(sample) for sample in header[value_start:]]
        for fields in reader:
            if len(fields) <= key_column:
                continue
            key = fields[key_column].split(".")[0] if strip_key_version else fields[key_column]
            if key not in wanted:
                continue
            values = pd.to_numeric(pd.Series(fields[value_start:], index=samples), errors="coerce")
            rows[key] = values.groupby(level=0).mean()
    return samples, rows


def motif_consensus(name: str, width: int, rows: list[list[float]]) -> str | None:
    if name and width and len(rows) >= width:
        return "".join("ACGT"[int(np.argmax(row))] for row in rows[:width])
    return None


def fit_mediation(
    g: pd.Series,
    mediator: pd.Series,
    outcome: pd.Series,
    covariates: pd.DataFrame,
    other_signals: pd.DataFrame,
    bootstrap_n: int,
    seed: int,
    min_n: int,
) -> dict:
    model_parts = [g.rename("G"), mediator.rename("M"), outcome.rename("Y"), covariates]
    if not other_signals.empty:
        model_parts.append(other_signals)
    frame = pd.concat(model_parts, axis=1, join="inner").replace([np.inf, -np.inf], np.nan).dropna()
    covariate_columns = [column for column in frame.columns if column not in {"G", "M", "Y"}]
    frame = frame.loc[:, ["G", "M", "Y", *covariate_columns]]
    counts = frame["G"].round().astype(int).value_counts().sort_index()
    base = {
        "complete_case_n": len(frame),
        "genotype_counts": ";".join(f"{genotype}:{count}" for genotype, count in counts.items()),
        "covariates": ";".join(covariate_columns),
        "model_status": "ok",
    }
    if len(frame) < min_n or counts.size < 2:
        base["model_status"] = "insufficient_complete_cases_or_genotype_variation"
        return base

    c = frame[covariate_columns].to_numpy(dtype=float) if covariate_columns else np.empty((len(frame), 0))
    g_values = frame["G"].to_numpy(dtype=float)
    m_values = frame["M"].to_numpy(dtype=float)
    y_values = frame["Y"].to_numpy(dtype=float)
    x_m = sm.add_constant(np.column_stack([g_values, c]), has_constant="add")
    x_y = sm.add_constant(np.column_stack([g_values, m_values, c]), has_constant="add")
    model_m = sm.OLS(m_values, x_m).fit()
    model_y = sm.OLS(y_values, x_y).fit()
    model_total = sm.OLS(y_values, x_m).fit()
    a = model_m.params[1]
    b = model_y.params[2]
    indirect = a * b
    sobel_se = np.sqrt((b * model_m.bse[1]) ** 2 + (a * model_y.bse[2]) ** 2)
    sobel_z = indirect / sobel_se if sobel_se > 0 else np.nan

    _, _, full_rank_m, _ = lstsq(x_m, m_values, lapack_driver="gelsy", check_finite=False)
    _, _, full_rank_y, _ = lstsq(x_y, y_values, lapack_driver="gelsy", check_finite=False)
    rng = np.random.default_rng(seed)
    bootstrap = np.empty(bootstrap_n, dtype=float)
    valid_draw_n = 0
    attempt_n = 0
    rank_reject_n = 0
    solver_error_n = 0
    max_attempt_n = bootstrap_n * 10
    while valid_draw_n < bootstrap_n and attempt_n < max_attempt_n:
        attempt_n += 1
        sampled = rng.integers(0, len(frame), len(frame))
        try:
            coef_m, _, rank_m, _ = lstsq(
                x_m[sampled], m_values[sampled], lapack_driver="gelsy", check_finite=False
            )
            coef_y, _, rank_y, _ = lstsq(
                x_y[sampled], y_values[sampled], lapack_driver="gelsy", check_finite=False
            )
        except LinAlgError:
            solver_error_n += 1
            continue
        if rank_m < full_rank_m or rank_y < full_rank_y:
            rank_reject_n += 1
            continue
        indirect_boot = coef_m[1] * coef_y[2]
        if not np.isfinite(indirect_boot):
            solver_error_n += 1
            continue
        bootstrap[valid_draw_n] = indirect_boot
        valid_draw_n += 1

    if valid_draw_n == bootstrap_n:
        alpha = 0.025
        ci_low, ci_high = np.quantile(bootstrap, [alpha, 1 - alpha])
        lower_tail = (int(np.sum(bootstrap <= 0)) + 1) / (bootstrap_n + 1)
        upper_tail = (int(np.sum(bootstrap >= 0)) + 1) / (bootstrap_n + 1)
        bootstrap_p = min(1.0, 2 * min(lower_tail, upper_tail))
    else:
        base["model_status"] = "bootstrap_insufficient_estimable_resamples"
        ci_low = np.nan
        ci_high = np.nan
        bootstrap_p = np.nan

    mediator_residual = sm.OLS(m_values, x_m).fit().resid
    outcome_residual = sm.OLS(y_values, x_m).fit().resid
    residual_r, residual_p = pearsonr(mediator_residual, outcome_residual)
    return {
        **base,
        "a_effect": a,
        "a_se": model_m.bse[1],
        "b_effect": b,
        "b_se": model_y.bse[2],
        "b_pvalue": model_y.pvalues[2],
        "direct_effect": model_y.params[1],
        "direct_effect_pvalue": model_y.pvalues[1],
        "total_effect": model_total.params[1],
        "total_effect_pvalue": model_total.pvalues[1],
        "indirect_effect": indirect,
        "indirect_ci_low": ci_low,
        "indirect_ci_high": ci_high,
        "indirect_bootstrap_pvalue": bootstrap_p,
        "bootstrap_valid_draw_n": valid_draw_n,
        "bootstrap_attempt_n": attempt_n,
        "bootstrap_rank_reject_n": rank_reject_n,
        "bootstrap_solver_error_n": solver_error_n,
        "sobel_z": sobel_z,
        "sobel_pvalue": 2 * norm.sf(abs(sobel_z)) if pd.notna(sobel_z) else np.nan,
        "residual_association_r": residual_r,
        "residual_association_pvalue": residual_p,
    }


# --- loading tested events and references ---

associations = pd.read_csv(
    UNIVERSE,
    sep="\t",
    usecols=[
        "phenotype_id", "gene_name", "phenotype_chr", "phenotype_start",
        "phenotype_end", "phenotype_strand",
    ],
).drop_duplicates("phenotype_id")
with CONFIG_PATH.open(encoding="utf-8") as handle:
    config = _expand_env(yaml.safe_load(handle))
event_info = pd.read_csv(
    config["inputs"]["apa_phenotype_annotation"]["path"],
    sep="\t",
    usecols=["phenotype_id", "Predicted_Proximal_APA"],
).drop_duplicates("phenotype_id")
events = associations.merge(event_info, on="phenotype_id", how="left", validate="one_to_one")

# miRNA seed sites: reverse complement of mature nt 2-8
mirna_path = Path(config["reference_resources"]["mirbase_mature_fasta"]["path"])
records: list[tuple[str, str]] = []
header = ""
sequence: list[str] = []
with mirna_path.open(encoding="utf-8", errors="ignore") as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header and sequence:
                records.append((header, "".join(sequence)))
            header = line[1:]
            sequence = []
        else:
            sequence.append(line.upper().replace("U", "T"))
    if header and sequence:
        records.append((header, "".join(sequence)))

mirna_sites: dict[str, str] = {}
for header_text, mature in records:
    name = header_text.split()[0]
    if not name.startswith("hsa-") or len(mature) < 8:
        continue
    seed = mature[1:8]
    mirna_sites.setdefault(revcomp(seed), name)

# CIS-BP RBP consensus motifs from letter-probability matrices
cisbp_path = Path(config["reference_resources"]["cisbp_rna_motifs"]["path"])
rbp_motifs: list[tuple[str, str]] = []
current_name = ""
width = 0
rows: list[list[float]] = []
in_matrix = False
with cisbp_path.open(encoding="utf-8", errors="ignore") as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if line.startswith("MOTIF "):
            consensus = motif_consensus(current_name, width, rows)
            if consensus is not None:
                rbp_motifs.append((current_name, consensus))
            fields = line.split(maxsplit=2)
            current_name = fields[2] if len(fields) > 2 else fields[1]
            width = 0
            rows = []
            in_matrix = False
        elif line.startswith("letter-probability matrix"):
            width_text = line.split("w=")[-1].split()[0]
            width = int(width_text)
            rows = []
            in_matrix = True
        elif in_matrix and line:
            values = line.split()
            if len(values) == 4:
                rows.append([float(value) for value in values])
                if len(rows) == width:
                    consensus = motif_consensus(current_name, width, rows)
                    if consensus is not None:
                        rbp_motifs.append((current_name, consensus))
                    current_name = ""
                    width = 0
                    rows = []
                    in_matrix = False
consensus = motif_consensus(current_name, width, rows)
if consensus is not None:
    rbp_motifs.append((current_name, consensus))

# ARE patterns: RNAfold table plus fixed AUUUA-type motifs
are_table = pd.read_csv(config["reference_resources"]["viennarna_rnafold"]["path"], usecols=["Pattern"])
are_patterns = {
    value
    for value in are_table["Pattern"].dropna().astype(str).str.upper()
    if len(value) <= 13 and not (set(value) - {"A", "C", "G", "T"})
}
are_patterns.update({"ATTTA", "TATTTAT", "AATTTA", "TTATTTATT"})

fasta = pysam.FastaFile(config["reference_resources"]["genome_fasta"]["path"])


# --- event-wise DUR annotation ---

rows = []
for event in events.itertuples(index=False):
    proximal = pd.to_numeric(event.Predicted_Proximal_APA, errors="coerce")
    distal = event.phenotype_end if event.phenotype_strand == "+" else event.phenotype_start
    status = "ok"
    sequence = ""
    dur_start = pd.NA
    dur_end = pd.NA
    if pd.isna(proximal) or pd.isna(distal):
        status = "missing_pas_anchor"
    else:
        dur_start = int(min(proximal, distal))
        dur_end = int(max(proximal, distal))
        if dur_end <= dur_start:
            status = "nonpositive_dur_length"
        else:
            chrom = str(event.phenotype_chr)
            sequence = fasta.fetch(chrom, dur_start - 1, dur_end).upper()
            if event.phenotype_strand == "-":
                sequence = revcomp(sequence)

    mirna_counter: Counter[str] = Counter()
    rbp_counter: Counter[str] = Counter()
    if status == "ok":
        for motif, name in mirna_sites.items():
            count = overlapping_count(sequence, motif)
            if count:
                mirna_counter[name] += count
        for name, motif in rbp_motifs:
            count = overlapping_count(sequence, motif)
            if count:
                rbp_counter[name] += count

    are_total = sum(overlapping_count(sequence, motif) for motif in are_patterns) if status == "ok" else 0
    gre_total = sum(overlapping_count(sequence, motif) for motif in GRE_MOTIFS) if status == "ok" else 0
    ure_total = sum(overlapping_count(sequence, motif) for motif in URE_MOTIFS) if status == "ok" else 0
    length = len(sequence)
    rows.append({
        "event_annotation_id": event.phenotype_id,
        "phenotype_id": event.phenotype_id,
        "gene_name": event.gene_name,
        "phenotype_chr": event.phenotype_chr,
        "phenotype_strand": event.phenotype_strand,
        "proximal_pas_pos": proximal,
        "distal_pas_pos": distal,
        "dur_start": dur_start,
        "dur_end": dur_end,
        "dur_length": length,
        "event_annotation_status": status,
        "predicted_consequence_available": status == "ok",
        "mirna_total_sites": int(sum(mirna_counter.values())),
        "mirna_unique_families": int(len(mirna_counter)),
        "top_mirna_families": ";".join(f"{name}:{count}" for name, count in mirna_counter.most_common(10)),
        "dur_rbp_consensus_total_hits": int(sum(rbp_counter.values())),
        "dur_rbp_consensus_unique_rbps": int(len(rbp_counter)),
        "top_dur_rbp_consensus_hits": ";".join(f"{name}:{count}" for name, count in rbp_counter.most_common(10)),
        "are_total_hits": int(are_total),
        "gre_total_hits": int(gre_total),
        "ure_total_hits": int(ure_total),
    })
fasta.close()


# --- tested-universe density thresholds ---

annotations = pd.DataFrame(rows)
valid_length = annotations["dur_length"].replace(0, np.nan)
for prefix, count_column in {
    "mirna": "mirna_total_sites",
    "rbp": "dur_rbp_consensus_total_hits",
    "are": "are_total_hits",
    "gre": "gre_total_hits",
    "ure": "ure_total_hits",
}.items():
    density_column = f"{prefix}_hits_per_100nt"
    threshold_column = f"{prefix}_density_threshold_per_100nt"
    support_column = f"{prefix}_density_supported"
    annotations[density_column] = annotations[count_column] / valid_length * 100
    positive = pd.to_numeric(annotations[density_column], errors="coerce")
    positive = positive.loc[positive.gt(0)]
    threshold = float(positive.quantile(0.75)) if not positive.empty else math.nan
    annotations[threshold_column] = threshold
    annotations[support_column] = annotations[density_column].ge(threshold) & annotations[count_column].gt(0)

annotations["pas_anchor_resolution"] = "proximal_pipeline_distal_boundary_approximation"
annotations["genome_fasta_release"] = config["reference_resources"]["genome_fasta"]["release"]
annotations["mirna_resource_release"] = config["reference_resources"]["mirbase_mature_fasta"]["release"]
annotations["rbp_resource_release"] = config["reference_resources"]["cisbp_rna_motifs"]["release"]
annotations["analysis_version"] = ANALYSIS_VERSION
annotations["config_path"] = str(CONFIG_PATH)
annotations.to_csv(EVENT_ANNOTATIONS, sep="\t", index=False, compression="gzip")

summary = pd.DataFrame({
    "metric": ["tested_events", "annotated_events", "missing_or_invalid_events", "mirna_seed_patterns", "rbp_consensus_motifs"],
    "value": [len(annotations), annotations["predicted_consequence_available"].sum(), (~annotations["predicted_consequence_available"]).sum(), len(mirna_sites), len(rbp_motifs)],
})
summary.to_csv(SUMMARY, sep="\t", index=False)
print(f"Written {len(annotations)} event-wise consequence annotations -> {EVENT_ANNOTATIONS}")


# --- loading framework inputs ---

inputs = config["inputs"]
mediation_config = config["parameters"]["mediation"]
signals = pd.read_csv(SIGNALS, sep="\t")
signals = signals.drop_duplicates("signal_id")
missing_gxstage_fields = [column for column in GXSTAGE_FIELDS if column not in signals.columns]
if missing_gxstage_fields:
    raise ValueError(f"Signal mechanisms lack GxStage fields: {missing_gxstage_fields}")
signal_gxstage = signals[["signal_id", *GXSTAGE_FIELDS]].copy()

# GTF gene_id (version stripped) and gene_name -> gene_id
gene_map = {}
with gzip.open(Path(config["reference_resources"]["gene_annotation_gtf"]["path"]), "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        gene_id_match = re.search(r'gene_id "([^"]+)"', fields[8])
        gene_name_match = re.search(r'gene_name "([^"]+)"', fields[8])
        if gene_id_match:
            gene_id = gene_id_match.group(1).split(".")[0]
            gene_map[gene_id] = gene_id
            if gene_name_match:
                gene_map[gene_name_match.group(1)] = gene_id

signals["expression_gene_id"] = signals["gene_name"].astype(str).map(gene_map)
signals.loc[signals["gene_name"].astype(str).str.startswith("ENSG"), "expression_gene_id"] = (
    signals.loc[signals["gene_name"].astype(str).str.startswith("ENSG"), "gene_name"].astype(str).str.split(".").str[0]
)

variant_ids = set(signals["representative_variant_id"].dropna().astype(str))
phenotype_ids = set(signals["phenotype_id"].dropna().astype(str))
expression_ids = set(signals["expression_gene_id"].dropna().astype(str))
_, genotypes = load_selected_matrix(Path(inputs["genotype_dosage"]["path"]), variant_ids, 0, 1)
_, expression = load_selected_matrix(
    Path(inputs["bulk_rna_tpm"]["path"]), expression_ids, 0, 1, strip_key_version=True
)

# protein: best-supported row per gene, log2(intensity + 1)
table = pd.read_csv(Path(inputs["protein_abundance"]["path"]), sep="\t")
table = table.loc[table["Gene name"].astype(str).isin(set(signals["gene_name"].astype(str)))].copy()
table["Unique peptides"] = pd.to_numeric(table["Unique peptides"], errors="coerce")
table = table.sort_values("Unique peptides", ascending=False, kind="stable").drop_duplicates("Gene name")
intensity_columns = [column for column in table.columns if column.startswith("Intensity ")]
protein = {}
for _, row in table.iterrows():
    gene = str(row["Gene name"])
    values = pd.to_numeric(row[intensity_columns], errors="coerce")
    values.index = [normalize_sample(column) for column in intensity_columns]
    protein[gene] = np.log2(values.groupby(level=0).mean().where(values.groupby(level=0).mean().gt(0)) + 1)

stage_data = {}
for stage in ["prenatal", "postnatal"]:
    _, mediator = load_selected_matrix(Path(inputs[f"{stage}_qtl_phenotypes"]["path"]), phenotype_ids, 3, 4)
    covariates = pd.read_csv(Path(inputs[f"{stage}_mediation_covariates"]["path"]), sep="\t", index_col=0).T
    covariates.index = [normalize_sample(sample) for sample in covariates.index]
    covariates = covariates.apply(pd.to_numeric, errors="coerce")
    covariates = covariates.loc[:, covariates.nunique(dropna=True).gt(1)]
    stage_data[stage] = {"mediator": mediator, "covariates": covariates}


# --- event-level predicted consequence annotations ---

predicted = pd.read_csv(EVENT_ANNOTATIONS, sep="\t")
if predicted.duplicated("phenotype_id").any():
    raise ValueError("Event consequence annotations are not unique by phenotype_id")


# --- genotype-aware mediation models ---

event_signal_variants = (
    signals.groupby(["stage", "phenotype_id"], sort=False)["representative_variant_id"]
    .apply(list)
    .to_dict()
)


def fit_signal(signal_index: int) -> list[dict]:
    signal = signal_records[signal_index]
    stage = signal.stage
    mediator = stage_data[stage]["mediator"].get(signal.phenotype_id)
    genotype = genotypes.get(signal.representative_variant_id)
    if mediator is None or genotype is None:
        return []
    event_variants = event_signal_variants[(stage, signal.phenotype_id)]
    other = {
        f"other_signal_{variant}": genotypes[variant]
        for variant in event_variants
        if variant != signal.representative_variant_id and variant in genotypes
    }
    other_signals = pd.DataFrame(other)
    outcomes = []
    if pd.notna(signal.expression_gene_id) and signal.expression_gene_id in expression:
        outcomes.append(("rna_expression", signal.expression_gene_id, np.log2(expression[signal.expression_gene_id] + 1)))
    if signal.gene_name in protein:
        outcomes.append(("protein_abundance", signal.gene_name, protein[signal.gene_name]))
    signal_rows = []
    for outcome_type, outcome_id, outcome in outcomes:
        model = fit_mediation(
            genotype,
            mediator,
            outcome,
            stage_data[stage]["covariates"],
            other_signals,
            int(mediation_config["bootstrap_n"]),
            int(mediation_config["bootstrap_seed"]),
            int(mediation_config["minimum_complete_case_n"]),
        )
        signal_rows.append({
            "signal_id": signal.signal_id,
            "candidate_id": signal.candidate_id,
            "stage": stage,
            "phenotype_id": signal.phenotype_id,
            "gene_name": signal.gene_name,
            "variant_id": signal.variant_id,
            "genotype_variant_id": signal.representative_variant_id,
            "outcome_type": outcome_type,
            "outcome_id": outcome_id,
            "outcome_scope": "same_gene",
            "mediator_scale": inputs[f"{stage}_qtl_phenotypes"]["scale"],
            **model,
        })
    return signal_rows


rows = []
signal_records = list(signals.itertuples(index=False))
worker_n = max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
with ProcessPoolExecutor(max_workers=worker_n, mp_context=multiprocessing.get_context("fork")) as executor:
    for signal_rows in executor.map(fit_signal, range(len(signal_records)), chunksize=1):
        rows.extend(signal_rows)

results = pd.DataFrame(rows)
if results.empty:
    results = pd.DataFrame(columns=[
        "signal_id", "candidate_id", "stage", "phenotype_id", "gene_name", "variant_id",
        "outcome_type", "outcome_id", "outcome_scope", "model_status", "indirect_bootstrap_pvalue",
    ])
else:
    # BH within stage x outcome_type
    indirect_fdr = pd.Series(np.nan, index=results.index, dtype=float)
    for _, index in results.groupby(["stage", "outcome_type"], dropna=False).groups.items():
        values = pd.to_numeric(results.loc[index, "indirect_bootstrap_pvalue"], errors="coerce")
        valid = values.notna()
        if valid.any():
            indirect_fdr.loc[values.index[valid]] = multipletests(values.loc[valid], method="fdr_bh")[1]
    results["indirect_fdr"] = indirect_fdr
    results = results.merge(predicted, on=["phenotype_id", "gene_name"], how="left", validate="many_to_one")
    results["predicted_consequence_available"] = results["predicted_consequence_available"].fillna(False).astype(bool)
    fdr = float(mediation_config["fdr"])
    labels = []
    for _, row in results.iterrows():
        ci_supported = pd.notna(row.get("indirect_ci_low")) and row["indirect_ci_low"] * row["indirect_ci_high"] > 0
        if ci_supported and pd.notna(row.get("indirect_fdr")) and row["indirect_fdr"] <= fdr:
            labels.append("mediation_supported")
        elif pd.notna(row.get("b_pvalue")) and row["b_pvalue"] < 0.05:
            labels.append("conditional_association_only")
        elif bool(row.get("predicted_consequence_available", False)):
            labels.append("predicted_consequence_only")
        else:
            labels.append("consequence_unresolved")
    results["downstream_consequence_evidence_label"] = labels

results = results.merge(signal_gxstage, on="signal_id", how="left", validate="many_to_one")

results["bootstrap_n"] = int(mediation_config["bootstrap_n"])
results["bootstrap_pvalue_estimator"] = mediation_config["bootstrap_pvalue_estimator"]
results["multiple_testing_family"] = results.get("stage", pd.Series(dtype=str)).astype(str) + "::" + results.get("outcome_type", pd.Series(dtype=str)).astype(str)
results["causal_limitation"] = mediation_config["causal_interpretation"]
results["analysis_version"] = ANALYSIS_VERSION
results["config_path"] = str(CONFIG_PATH)

if EXCLUDED_OUTPUT.exists():
    excluded = pd.read_csv(EXCLUDED_OUTPUT, sep="\t")
    missing_excluded_fields = [column for column in EXCLUDED_FIELDS if column not in excluded.columns]
    if missing_excluded_fields:
        raise ValueError(f"Excluded/blocked table lacks fields: {missing_excluded_fields}")
    excluded = excluded.loc[~excluded["step"].eq("mediation_model"), EXCLUDED_FIELDS].copy()
else:
    excluded = pd.DataFrame(columns=EXCLUDED_FIELDS)
blocked_models = results.loc[
    ~results["model_status"].eq("ok"),
    ["candidate_id", "stage", "phenotype_id", "variant_id", "model_status"],
].copy()
blocked_models["source_row_id"] = pd.NA
blocked_models["step"] = "mediation_model"
blocked_models["reason"] = blocked_models.pop("model_status")
blocked_models["severity"] = "blocked"
blocked_models["analysis_version"] = ANALYSIS_VERSION
blocked_models["config_path"] = str(CONFIG_PATH)
blocked_models = blocked_models[EXCLUDED_FIELDS]
excluded = pd.concat([excluded, blocked_models], ignore_index=True)

results.to_csv(OUTPUT, sep="\t", index=False, compression="gzip")
excluded.to_csv(EXCLUDED_OUTPUT, sep="\t", index=False, compression="gzip")
print(f"Written {len(results)} stage-event-outcome models -> {OUTPUT}")
print(f"Recorded {len(blocked_models)} blocked mediation models -> {EXCLUDED_OUTPUT}")
