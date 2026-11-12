#!/usr/bin/env python3
# Draw the numbered APA mechanism case figures.
# Writes the case source/layout tables and locus caches to data/07_coloc.
from __future__ import annotations
import os
from datetime import datetime, timezone
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator

try:
    import pgenlib
except ImportError:
    pgenlib = None

import pysam


# --- helpers ---

def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def parse_variant(variant_id: object) -> dict[str, object]:
    chrom, pos, ref, alt = str(variant_id).split(":")
    return {"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt}


def source_signature(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "source_path": str(path),
        "source_bytes": stat.st_size,
        "source_mtime_unix": stat.st_mtime,
    }


def load_filtered_nominal(
    source_path: Path,
    cache_path: Path,
    meta_path: Path,
    phenotype_ids: set[str],
    usecols: list[str],
    chunksize: int,
) -> pd.DataFrame:
    if not phenotype_ids:
        return pd.DataFrame(columns=usecols)
    cache_current = False
    if cache_path.exists() and meta_path.exists() and cache_path.stat().st_size > 0:
        meta = pd.read_csv(meta_path, sep="\t")
        if not meta.empty:
            observed = meta.iloc[0].to_dict()
            signature = source_signature(source_path)
            cached_ids = set(str(observed.get("phenotype_ids", "")).split(";")) - {""}
            cache_current = (
                str(observed.get("source_path")) == signature["source_path"]
                and int(observed.get("source_bytes", -1)) == int(signature["source_bytes"])
                and abs(float(observed.get("source_mtime_unix", -1)) - float(signature["source_mtime_unix"])) <= 1e-3
                and phenotype_ids.issubset(cached_ids)
            )
    if cache_current:
        cached = pd.read_csv(cache_path, sep="\t", usecols=usecols)
        return cached.loc[cached["phenotype_id"].isin(phenotype_ids)].copy()

    rows: list[pd.DataFrame] = []
    for chunk in pd.read_csv(source_path, sep="\t", usecols=usecols, chunksize=chunksize):
        subset = chunk.loc[chunk["phenotype_id"].isin(phenotype_ids)].copy()
        if not subset.empty:
            rows.append(subset)
    output = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=usecols)
    output.to_csv(cache_path, sep="\t", index=False)
    meta = pd.DataFrame([
        {
            **source_signature(source_path),
            "phenotype_ids": ";".join(sorted(phenotype_ids)),
            "rows": len(output),
            "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
        }
    ])
    meta.to_csv(meta_path, sep="\t", index=False)
    return output


def ld_bin(value: object) -> str:
    if pd.isna(value):
        return "Missing"
    value = float(value)
    if value >= 0.8:
        return "0.8-1.0"
    if value >= 0.6:
        return "0.6-0.8"
    if value >= 0.4:
        return "0.4-0.6"
    if value >= 0.2:
        return "0.2-0.4"
    return "0.0-0.2"


def read_dosage(reader: object, variant_index: int, sample_n: int) -> np.ndarray:
    dosage = np.empty(sample_n, dtype=np.float32)
    reader.read_dosages(variant_index, dosage)
    dosage = dosage.astype(float)
    dosage[dosage < 0] = np.nan
    return dosage


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def short_mirna_name(name: object) -> str:
    return str(name).removeprefix("hsa-")


def set_clean_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(TEXT_MID)
    axis.spines["bottom"].set_color(TEXT_MID)
    axis.tick_params(axis="both", colors=TEXT_DARK, labelsize=7, width=0.7, length=3)


def add_panel_label(axis: plt.Axes, label: str) -> None:
    # Label inside the top-left corner, clear of tick labels and titles.
    axis.text(
        0.015,
        0.985,
        label,
        transform=axis.transAxes,
        fontsize=12,
        fontweight="normal",
        ha="left",
        va="top",
        color="black",
        gid="check_text",
    )


def annotate_lead_label(
    axis: plt.Axes,
    data: pd.DataFrame,
    lead_variant: str,
    x_col: str,
    y_col: str,
    label: str,
    va: str = "bottom",
    ha: str = "center",
    offset_frac: float = 0.06,
    margin_frac: float = 0.18,
) -> None:
    # Put the lead label inside the axes, expanding the y-limits when needed.
    lead = data.loc[data["variant_id"].eq(lead_variant)]
    if lead.empty:
        return
    y_lo, y_hi = axis.get_ylim()
    y_range = y_hi - y_lo
    y_value = float(lead[y_col].iloc[0])
    if va == "bottom":
        label_y = y_value + offset_frac * y_range
        required = label_y + margin_frac * y_range
        if required > y_hi:
            axis.set_ylim(top=required)
    else:
        label_y = y_value - offset_frac * y_range
        required = label_y - margin_frac * y_range
        if required < y_lo:
            axis.set_ylim(bottom=required)
    axis.text(
        float(lead[x_col].iloc[0]),
        label_y,
        label,
        ha=ha,
        va=va,
        fontsize=9,
        fontweight="bold",
        color="black",
        gid="check_text",
        clip_on=True,
    )


def pas_motif_label(values: pd.Series) -> str | None:
    if values.empty or pd.isna(values.iloc[0]):
        return None
    motifs = [motif.strip() for motif in str(values.iloc[0]).split(";") if motif.strip()]
    return "/".join(motifs[:2]) if motifs else None


def plot_locus(axis: plt.Axes, data: pd.DataFrame, y_col: str, lead_variant: str, title: str, ylabel: str) -> None:
    set_clean_axis(axis)
    for label in LD_ORDER:
        subset = data.loc[data["ld_bin"].eq(label)]
        if subset.empty:
            continue
        axis.scatter(
            subset["position"],
            subset[y_col],
            s=28,
            color=LD_COLORS[label],
            edgecolor="#555555",
            linewidth=0.25,
            alpha=0.95 if label != "Missing" else 0.75,
            zorder=2 if label != "Missing" else 1,
        )
    lead = data.loc[data["variant_id"].eq(lead_variant)]
    if not lead.empty:
        axis.scatter(
            lead["position"],
            lead[y_col],
            s=120,
            marker="D",
            color=LEAD_PURPLE,
            edgecolor="black",
            linewidth=0.8,
            zorder=4,
        )
    axis.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=8, color="black", gid="check_text")
    axis.set_ylabel(ylabel, fontsize=11, color=TEXT_DARK)
    axis.set_xlabel("")
    axis.ticklabel_format(style="plain", axis="x", useOffset=False)
    axis.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))


# --- config ---

DATA_DIR = Path(os.environ["APA_ROOT"]) / "data/07_coloc"
FIGURE_DIR = Path(os.environ["APA_ROOT"]) / "figures/07_coloc"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

SIGNAL_PATH = DATA_DIR / "02_signal_mechanisms.tsv.gz"
CONSEQUENCE_PATH = DATA_DIR / "03_event_consequences.tsv.gz"
PAS_PATH = DATA_DIR / "pas.tsv.gz"
RBP_PATH = DATA_DIR / "rbp.tsv.gz"
INTERNAL_PATH = DATA_DIR / "07_internal_eqtl_integration_contract.tsv.gz"
INTERNAL_RAW_PATH = DATA_DIR / "internal_eqtl.tsv"
EXTERNAL_PATH = DATA_DIR / "external_all.tsv"
LOCAL_APA_NOMINAL = DATA_DIR / "local_apa_postnatal.tsv.gz"
ANNOTATIONS_PATH = DATA_DIR / "annotations.tsv.gz"
MIRNA_FASTA = Path(os.environ["APA_ROOT"]) / "ref/mature.fa"
GENOME_FASTA = Path(os.environ["APA_ROOT"]) / "ref/hg38.fa"

LD_PGEN = Path(os.environ["APA_ROOT"]) / "data/05_gxs/pre_input/genotype.pgen"
LD_PVAR = Path(os.environ["APA_ROOT"]) / "data/05_gxs/pre_input/genotype.pvar"

OUT_SOURCE = DATA_DIR / "case_source.tsv"
OUT_LAYOUT = DATA_DIR / "case_layout.tsv"
OUT_APA_CACHE = DATA_DIR / "case_apa_locus.tsv.gz"
OUT_APA_META = DATA_DIR / "case_apa_locus.meta.tsv"
OUT_EQTL_CACHE = DATA_DIR / "case_eqtl_locus.tsv.gz"
OUT_EQTL_META = DATA_DIR / "case_eqtl_locus.meta.tsv"
OUT_LD_CACHE = DATA_DIR / "case_locus_ld.tsv.gz"
OUT_LOCUS_SOURCE = DATA_DIR / "case_locus_source.tsv.gz"

TARGETS = [
    {"figure_code": "01", "gene_name": "CEP41", "variant_id": "chr7:130393789:C:T"},
    {"figure_code": "02", "gene_name": "GTDC1", "variant_id": "chr2:143938128:C:G"},
    {"figure_code": "03", "gene_name": "STMN3", "variant_id": "chr20:63641230:G:T"},
    {
        "figure_code": "04",
        "gene_name": "ZSWIM7",
        "variant_id": "chr17:15976596:C:T",
        "display_variant_label": "rs1045599",
    },
]

TRACK_GREY = "#C9D1DA"
DUR_GREEN = "#79B876"
PROX_PURPLE = "#7B2D8E"
DIST_BLUE = "#3B83B5"
VARIANT_RED = "#B31921"
TEXT_DARK = "#38424E"
TEXT_MID = "#5F6975"
LOSS_RED = "#C51B24"
GAIN_GREEN = "#087F5B"
LEAD_PURPLE = "#6F3A8A"

LD_COLORS = {
    "0.8-1.0": "#FF1A1A",
    "0.6-0.8": "#FFB22E",
    "0.4-0.6": "#128A22",
    "0.2-0.4": "#9AD3E8",
    "0.0-0.2": "#111B8A",
    "Missing": "#D0D0D0",
}
LD_ORDER = ["Missing", "0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]

TOP_MIRNA_FAMILIES = 10
MAX_SITES_PER_FAMILY = 6
MAX_MIRNA_LABELS = 4
MIRNA_SITE_COLORS = {"8mer": "#C0392B", "7mer-m8": "#D9822B", "7mer-a1": "#2F80ED"}


# --- load current example rows ---

signal_columns = [
    "candidate_id",
    "signal_id",
    "stage",
    "gene_name",
    "phenotype_id",
    "variant_id",
    "phenotype_chr",
    "phenotype_start",
    "phenotype_end",
    "phenotype_strand",
    "proximal_pas_pos",
    "distal_pas_pos",
    "dur_length",
    "pval_nominal",
    "beta_alt_PDUI",
    "pip",
    "apa_direction",
    "snp_to_apa_evidence_label",
    "direct_pas_evidence",
    "local_rbp_motif_change",
    "local_structure_change",
    "analysis_version",
]
consequence_columns = [
    "candidate_id",
    "outcome_type",
    "outcome_scope",
    "complete_case_n",
    "direct_effect_pvalue",
    "indirect_fdr",
    "residual_association_pvalue",
    "downstream_consequence_evidence_label",
    "analysis_version",
]

target_frame = pd.DataFrame(TARGETS)
target_keys = set(target_frame[["gene_name", "variant_id"]].apply(tuple, axis=1))
signal = pd.read_csv(SIGNAL_PATH, sep="\t", usecols=signal_columns)
signal = signal.loc[
    signal["stage"].eq("postnatal")
    & signal[["gene_name", "variant_id"]].apply(tuple, axis=1).isin(target_keys)
].copy()
signal = signal.drop_duplicates(["gene_name", "variant_id"])

if len(signal) != len(TARGETS):
    raise RuntimeError("The configured case examples are not uniquely available in the signal table")

source = signal.merge(target_frame, on=["gene_name", "variant_id"], how="left", validate="one_to_one")
source = source.rename(columns={"analysis_version": "association_analysis_version"})

consequence = pd.read_csv(CONSEQUENCE_PATH, sep="\t", usecols=consequence_columns)
consequence = consequence.loc[
    consequence["outcome_type"].eq("rna_expression")
    & consequence["outcome_scope"].eq("same_gene")
    & consequence["candidate_id"].isin(source["candidate_id"])
].copy()
consequence = consequence.drop_duplicates("candidate_id").rename(
    columns={"analysis_version": "consequence_analysis_version"}
)
source = source.merge(consequence, on="candidate_id", how="left", validate="one_to_one")
if source["downstream_consequence_evidence_label"].isna().any():
    raise RuntimeError("A configured case example has no same-gene RNA consequence row")

pas = pd.read_csv(
    PAS_PATH,
    sep="\t",
    usecols=["candidate_id", "pas_label", "pas_pos", "ref_pas_motifs", "analysis_version"],
)
pas = pas.loc[pas["candidate_id"].isin(source["candidate_id"])].copy()
mechanism_version = source["association_analysis_version"].iloc[0]
if set(pas["analysis_version"].dropna().astype(str)) != {mechanism_version}:
    raise RuntimeError("PAS motif detail table is not from the current mechanism analysis")
if set(pas["candidate_id"]) != set(source["candidate_id"]):
    raise RuntimeError("A configured case example lacks regenerated PAS motif detail")

rbp = pd.read_csv(
    RBP_PATH,
    sep="\t",
    usecols=[
        "candidate_id",
        "rbp_name",
        "delta_norm_score",
        "delta_direction",
        "analysis_version",
        "rbp_expression_supported",
        "rbp_expression_gene_symbol",
    ],
)
if set(rbp["analysis_version"].dropna().astype(str)) != {mechanism_version}:
    raise RuntimeError("RBP motif detail table is not from the current mechanism analysis")
rbp = rbp.loc[
    rbp["candidate_id"].isin(source["candidate_id"])
    & rbp["delta_direction"].isin(["gain", "loss"])
    & rbp["rbp_expression_supported"].astype(str).eq("True")
].copy()

annotations = pd.read_csv(
    ANNOTATIONS_PATH,
    sep="\t",
    usecols=[
        "phenotype_id",
        "mirna_total_sites",
        "mirna_unique_families",
        "top_mirna_families",
        "analysis_version",
    ],
)
if set(annotations["analysis_version"].dropna().astype(str)) != {mechanism_version}:
    raise RuntimeError("miRNA annotation table is not from the current mechanism analysis")
annotations = annotations.loc[annotations["phenotype_id"].isin(source["phenotype_id"])].copy()
annotations = annotations.drop_duplicates("phenotype_id").rename(
    columns={"analysis_version": "mirna_analysis_version"}
)
source = source.merge(annotations, on="phenotype_id", how="left", validate="one_to_one")
if source["top_mirna_families"].isna().any():
    raise RuntimeError("A configured case example lacks miRNA annotation rows")

internal_contract = pd.read_csv(
    INTERNAL_PATH,
    sep="\t",
    usecols=[
        "candidate_id",
        "egene_id",
        "gene_mapping_scope",
        "eqtl_coloc_pp4",
        "primary_functional_support",
        "eqtl_nominal_path",
        "eqtl_sample_n",
        "analysis_version",
    ],
)
internal_contract["primary_functional_support"] = internal_contract["primary_functional_support"].map(as_bool)
internal_contract = internal_contract.loc[
    internal_contract["candidate_id"].isin(source["candidate_id"])
    & internal_contract["gene_mapping_scope"].eq("same_gene")
    & internal_contract["primary_functional_support"]
].copy()
internal_contract["eqtl_coloc_pp4"] = pd.to_numeric(internal_contract["eqtl_coloc_pp4"], errors="coerce")
internal_contract = (
    internal_contract.sort_values(["candidate_id", "eqtl_coloc_pp4", "egene_id"], ascending=[True, False, True], kind="stable")
    .drop_duplicates("candidate_id")
    .rename(
        columns={
            "egene_id": "eqtl_phenotype_id",
            "analysis_version": "internal_eqtl_analysis_version",
        }
    )
)

internal_raw = pd.read_csv(
    INTERNAL_RAW_PATH,
    sep="\t",
    usecols=[
        "candidate_id",
        "phenotype_id",
        "stage",
        "eQTL_coloc",
        "eqtl_coloc_pp4",
        "eqtl_coloc_status",
        "n_shared_snps",
        "top_apa_snp",
        "top_eqtl_snp",
        "apa_nominal_path",
        "eqtl_nominal_path",
        "eqtl_sample_n",
        "analysis_version",
    ],
)
internal_raw["eQTL_coloc"] = internal_raw["eQTL_coloc"].map(as_bool)
internal_raw = internal_raw.loc[
    internal_raw["candidate_id"].isin(internal_contract["candidate_id"])
    & internal_raw["eQTL_coloc"]
].copy()
internal_raw["eqtl_coloc_pp4"] = pd.to_numeric(internal_raw["eqtl_coloc_pp4"], errors="coerce")
internal_raw = internal_raw.rename(columns={"phenotype_id": "eqtl_phenotype_id"})
internal_selected = internal_contract.merge(
    internal_raw[
        [
            "candidate_id",
            "eqtl_phenotype_id",
            "eQTL_coloc",
            "eqtl_coloc_status",
            "n_shared_snps",
            "top_apa_snp",
            "top_eqtl_snp",
            "apa_nominal_path",
        ]
    ],
    on=["candidate_id", "eqtl_phenotype_id"],
    how="left",
    validate="one_to_one",
)
if internal_selected["n_shared_snps"].isna().any():
    raise RuntimeError("A primary same-gene eQTL support row lacks raw coloc row details")
source = source.merge(internal_selected, on="candidate_id", how="left")

external = pd.read_csv(
    EXTERNAL_PATH,
    sep="\t",
    usecols=[
        "candidate_id",
        "source",
        "event_match_scope",
        "PP.H4",
        "external_validation_eligible",
        "analysis_version",
    ],
)
external["external_validation_eligible"] = external["external_validation_eligible"].map(as_bool)
external = external.loc[
    external["candidate_id"].isin(source["candidate_id"])
    & external["external_validation_eligible"]
    & external["event_match_scope"].eq("same_gene_apa_overlap")
].copy()
external["source_short"] = external["source"].str.replace("GTEx_v11_Brain_", "", regex=False)
external = (
    external.groupby("candidate_id", as_index=False)
    .agg(
        external_sources=("source_short", lambda values: ", ".join(sorted(values))),
        external_apa_pp4=("PP.H4", "max"),
        external_analysis_version=("analysis_version", "first"),
    )
)
source = source.merge(external, on="candidate_id", how="left")

source["pas_detail_analysis_version"] = mechanism_version
source["rbp_detail_analysis_version"] = mechanism_version
source["gene_order"] = source["figure_code"].astype(int)
source = source.sort_values("gene_order").reset_index(drop=True)

# --- build locus-level plotting data ---

nominal_usecols = ["phenotype_id", "variant_id", "af", "pval_nominal", "slope", "slope_se"]
apa_loci = load_filtered_nominal(
    LOCAL_APA_NOMINAL,
    OUT_APA_CACHE,
    OUT_APA_META,
    set(source["phenotype_id"].astype(str)),
    nominal_usecols,
    chunksize=700_000,
)
apa_loci = apa_loci.merge(
    source[["candidate_id", "phenotype_id"]],
    on="phenotype_id",
    how="inner",
    validate="many_to_one",
)

eqtl_loci_parts: list[pd.DataFrame] = []
for eqtl_path_text, group in source.dropna(subset=["eqtl_phenotype_id"]).groupby("eqtl_nominal_path", dropna=True):
    eqtl_path = Path(str(eqtl_path_text))
    eqtl_loci_parts.append(
        load_filtered_nominal(
            eqtl_path,
            OUT_EQTL_CACHE,
            OUT_EQTL_META,
            set(group["eqtl_phenotype_id"].astype(str)),
            nominal_usecols,
            chunksize=1_000_000,
        )
    )
eqtl_loci = pd.concat(eqtl_loci_parts, ignore_index=True) if eqtl_loci_parts else pd.DataFrame(columns=nominal_usecols)
if not eqtl_loci.empty:
    eqtl_loci = eqtl_loci.merge(
        source.loc[
            source["eqtl_phenotype_id"].notna(),
            ["candidate_id", "eqtl_phenotype_id"],
        ].rename(columns={"eqtl_phenotype_id": "phenotype_id"}),
        on="phenotype_id",
        how="inner",
        validate="many_to_one",
    )

for frame in [apa_loci, eqtl_loci]:
    if not frame.empty:
        parsed = frame["variant_id"].map(parse_variant)
        frame["position"] = [item["pos"] for item in parsed]
        frame["neglog10p"] = -np.log10(pd.to_numeric(frame["pval_nominal"], errors="coerce").clip(lower=1e-300, upper=1.0))

case_variant_sets: dict[str, tuple[str, set[str]]] = {}
for row in source.itertuples(index=False):
    variants = set(apa_loci.loc[apa_loci["candidate_id"].eq(row.candidate_id), "variant_id"].astype(str))
    if pd.notna(row.eqtl_phenotype_id):
        variants.update(eqtl_loci.loc[eqtl_loci["candidate_id"].eq(row.candidate_id), "variant_id"].astype(str))
    case_variant_sets[row.candidate_id] = (str(row.variant_id), variants)

rows: list[dict[str, object]] = []
if pgenlib is None or not LD_PGEN.exists() or not LD_PVAR.exists():
    for candidate_id, (lead, variants) in case_variant_sets.items():
        for variant in sorted(variants):
            rows.append({
                "candidate_id": candidate_id,
                "variant_id": variant,
                "ld_r2_to_case_variant": np.nan,
                "ld_status": "pgenlib_or_reference_missing",
            })
    ld_table = pd.DataFrame(rows)
else:
    needed = {lead for lead, _ in case_variant_sets.values()}
    for _, variants in case_variant_sets.values():
        needed.update(variants)
    index_map: dict[str, int] = {}
    variant_index = 0
    with LD_PVAR.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] in needed:
                index_map[fields[2]] = variant_index
                if len(index_map) == len(needed):
                    break
            variant_index += 1
    reader = pgenlib.PgenReader(bytes(str(LD_PGEN), "utf8"))
    sample_n = reader.get_raw_sample_ct()
    dosage_cache: dict[str, np.ndarray] = {}

    try:
        for candidate_id, (lead, variants) in case_variant_sets.items():
            lead_index = index_map.get(lead)
            lead_dosage = None
            if lead_index is not None:
                lead_dosage = read_dosage(reader, lead_index, sample_n)
            for variant in sorted(variants):
                if lead_index is None or variant not in index_map:
                    rows.append({
                        "candidate_id": candidate_id,
                        "variant_id": variant,
                        "ld_r2_to_case_variant": np.nan,
                        "ld_status": "variant_absent_from_ld_reference",
                    })
                    continue
                if variant == lead:
                    r2 = 1.0
                else:
                    if variant not in dosage_cache:
                        dosage_cache[variant] = read_dosage(reader, index_map[variant], sample_n)
                    mask = np.isfinite(lead_dosage) & np.isfinite(dosage_cache[variant])
                    r2 = np.nan
                    if mask.sum() >= 3:
                        a = lead_dosage[mask]
                        b = dosage_cache[variant][mask]
                        if np.nanstd(a) != 0 and np.nanstd(b) != 0:
                            r2 = float(np.corrcoef(a, b)[0, 1] ** 2)
                rows.append({
                    "candidate_id": candidate_id,
                    "variant_id": variant,
                    "ld_r2_to_case_variant": r2,
                    "ld_status": "ok" if np.isfinite(r2) else "not_estimable",
                })
    finally:
        reader.close()

    ld_table = pd.DataFrame(rows)
    ld_table.to_csv(OUT_LD_CACHE, sep="\t", index=False)
for frame in [apa_loci, eqtl_loci]:
    if not frame.empty:
        frame.drop(columns=[column for column in ["ld_r2_to_case_variant", "ld_status", "ld_bin"] if column in frame.columns], inplace=True)
        frame = frame.merge(ld_table, on=["candidate_id", "variant_id"], how="left")

apa_loci = apa_loci.merge(ld_table, on=["candidate_id", "variant_id"], how="left")
eqtl_loci = eqtl_loci.merge(ld_table, on=["candidate_id", "variant_id"], how="left") if not eqtl_loci.empty else eqtl_loci
apa_loci["ld_bin"] = apa_loci["ld_r2_to_case_variant"].map(ld_bin)
if not eqtl_loci.empty:
    eqtl_loci["ld_bin"] = eqtl_loci["ld_r2_to_case_variant"].map(ld_bin)

source["apa_locus_rows"] = source["candidate_id"].map(apa_loci.groupby("candidate_id").size()).fillna(0).astype(int)
source["eqtl_locus_rows"] = source["candidate_id"].map(eqtl_loci.groupby("candidate_id").size()).fillna(0).astype(int)
source["shared_locus_rows"] = 0

default_labels = []
for variant_id in source["variant_id"]:
    parsed = parse_variant(variant_id)
    default_labels.append(f"{parsed['chrom'].replace('chr', '')}:{parsed['pos']} {parsed['ref']}>{parsed['alt']}")
source["display_variant_label"] = source["display_variant_label"].fillna(pd.Series(default_labels, index=source.index))
source["figure_base"] = source.apply(
    lambda row: f"{row.figure_code}_case_{str(row.gene_name).lower()}_apa_eqtl_mechanism"
    if pd.notna(row.eqtl_phenotype_id)
    else f"{row.figure_code}_case_{str(row.gene_name).lower()}_apa_mechanism",
    axis=1,
)
source["pas_rbp_plotting_scope"] = (
    "RBP motif deltas are expression-supported, variant-window ranked labels; "
    "motif genomic positions are not in rbp.tsv.gz"
)

# --- miRNA seed-site scanning for the structure lane ---

records: dict[str, str] = {}
header: str | None = None
chunks: list[str] = []
with MIRNA_FASTA.open("r", encoding="utf-8", errors="ignore") as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper().replace("U", "T")
            header = line[1:]
            chunks = []
        else:
            chunks.append(line)
if header is not None:
    records[header] = "".join(chunks).upper().replace("U", "T")

families: dict[str, dict] = {}
for header_text, seq in records.items():
    mirna_name = header_text.split()[0]
    if not mirna_name.startswith("hsa-") or len(seq) < 8:
        continue
    seed = seq[1:8]
    families.setdefault(seed, {
        "seed": seed,
        "family_id": f"seed_{seed}",
        "representative_mirna": mirna_name,
    })
for entry in families.values():
    seed_rc = revcomp(entry["seed"])
    sixmer_rc = revcomp(entry["seed"][:6])
    entry["site_8mer"] = "A" + seed_rc
    entry["site_7mer_m8"] = seed_rc
    entry["site_7mer_a1"] = "A" + sixmer_rc
mirna_families = sorted(families.values(), key=lambda item: item["representative_mirna"])

case_mirna_sites: dict[str, pd.DataFrame] = {}
with pysam.FastaFile(str(GENOME_FASTA)) as fasta:
    for row in source.itertuples(index=False):
        strand = str(row.phenotype_strand)
        dur_low = int(min(row.proximal_pas_pos, row.distal_pas_pos))
        dur_high = int(max(row.proximal_pas_pos, row.distal_pas_pos))
        sequence = fasta.fetch(str(row.phenotype_chr), dur_low - 1, dur_high).upper()
        if strand == "-":
            sequence = revcomp(sequence)
        selected: list[str] = []
        if pd.notna(row.top_mirna_families):
            for part in str(row.top_mirna_families).split(";"):
                part = part.strip()
                if not part:
                    continue
                selected.append(part.split(":", 1)[0])
                if len(selected) >= TOP_MIRNA_FAMILIES:
                    break
        wanted = set(selected)
        site_rows: list[dict] = []
        for family in mirna_families:
            name = family["representative_mirna"]
            if name not in wanted:
                continue
            fam_rows: list[dict] = []
            for site_type, motif_key in [
                ("8mer", "site_8mer"),
                ("7mer-m8", "site_7mer_m8"),
                ("7mer-a1", "site_7mer_a1"),
            ]:
                motif = family[motif_key]
                starts: list[int] = []
                start = sequence.find(motif)
                while start != -1:
                    starts.append(start)
                    start = sequence.find(motif, start + 1)
                for pos in starts:
                    fam_rows.append(
                        {
                            "representative_mirna": name,
                            "family_id": family["family_id"],
                            "site_type": site_type,
                            "motif": motif,
                            "dur_start": pos,
                            "dur_end": pos + len(motif),
                        }
                    )
            fam_rows = sorted(
                fam_rows,
                key=lambda r: ({"8mer": 0, "7mer-m8": 1, "7mer-a1": 2}[r["site_type"]], r["dur_start"]),
            )[:MAX_SITES_PER_FAMILY]
            site_rows.extend(fam_rows)
        sites = pd.DataFrame(site_rows)
        if not sites.empty:
            motif_len = sites["dur_end"] - sites["dur_start"]
            if strand == "+":
                sites["site_x"] = dur_low + sites["dur_start"] + motif_len / 2.0
            else:
                sites["site_x"] = dur_high - sites["dur_start"] - motif_len / 2.0
            present = set(sites["representative_mirna"])
            rank = {name: index for index, name in enumerate(selected) if name in present}
            sites["family_rank"] = sites["representative_mirna"].map(rank)
            sites["site_type_rank"] = sites["site_type"].map({"8mer": 0, "7mer-m8": 1, "7mer-a1": 2})
            sites = sites.sort_values(["family_rank", "site_type_rank", "dur_start"], kind="stable")
        case_mirna_sites[row.candidate_id] = sites

source["mirna_plotting_scope"] = (
    "Per-site seed matches (8mer/7mer-m8/7mer-a1) for the top "
    f"{TOP_MIRNA_FAMILIES} miRNA families from top_mirna_families, capped at "
    f"{MAX_SITES_PER_FAMILY} sites per family; scanned on the DUR sequence from hg38.fa"
)
source["n_plotted_mirna_sites"] = source["candidate_id"].map(
    lambda value: 0 if case_mirna_sites[value].empty else int(len(case_mirna_sites[value]))
)
source["n_plotted_mirna_families"] = source["candidate_id"].map(
    lambda value: 0 if case_mirna_sites[value].empty else int(case_mirna_sites[value]["representative_mirna"].nunique())
)

all_locus = pd.concat(
    [
        apa_loci.assign(association_layer="APA"),
        eqtl_loci.assign(association_layer="same_gene_eQTL") if not eqtl_loci.empty else pd.DataFrame(),
    ],
    ignore_index=True,
)
all_locus.to_csv(OUT_LOCUS_SOURCE, sep="\t", index=False)


# --- draw one figure for each case example ---

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)

layout_rows: list[dict[str, object]] = []
for row in source.itertuples(index=False):
    has_internal = pd.notna(row.eqtl_phenotype_id)
    has_external = pd.notna(row.external_apa_pp4)
    case_apa = apa_loci.loc[apa_loci["candidate_id"].eq(row.candidate_id)].copy()
    case_eqtl = eqtl_loci.loc[eqtl_loci["candidate_id"].eq(row.candidate_id)].copy() if has_internal else pd.DataFrame()
    case_apa = case_apa.sort_values(["ld_bin", "position"], kind="stable")
    case_eqtl = case_eqtl.sort_values(["ld_bin", "position"], kind="stable") if not case_eqtl.empty else case_eqtl
    case_mirna = case_mirna_sites.get(row.candidate_id, pd.DataFrame())

    if has_internal:
        shared = case_apa[["variant_id", "neglog10p", "ld_r2_to_case_variant", "ld_bin"]].rename(
            columns={"neglog10p": "apa_neglog10p"}
        ).merge(
            case_eqtl[["variant_id", "neglog10p"]].rename(columns={"neglog10p": "eqtl_neglog10p"}),
            on="variant_id",
            how="inner",
        )
        source.loc[source["candidate_id"].eq(row.candidate_id), "shared_locus_rows"] = len(shared)
        figure = plt.figure(figsize=(10.8, 7.35))
        grid = figure.add_gridspec(
            3,
            2,
            height_ratios=[1.0, 1.0, 1.42],
            width_ratios=[1.28, 0.95],
            hspace=0.58,
            wspace=0.26,
        )
        apa_axis = figure.add_subplot(grid[0, 0])
        eqtl_axis = figure.add_subplot(grid[1, 0], sharex=apa_axis)
        shared_axis = figure.add_subplot(grid[0:2, 1])
        structure_axis = figure.add_subplot(grid[2, :])
        add_panel_label(apa_axis, "a")
        add_panel_label(eqtl_axis, "b")
        add_panel_label(structure_axis, "c")
        add_panel_label(shared_axis, "d")
        plot_locus(apa_axis, case_apa, "neglog10p", row.variant_id, "apaQTL association", r"$-\log_{10}(P)$")
        plot_locus(eqtl_axis, case_eqtl, "neglog10p", row.variant_id, "eQTL association", r"$-\log_{10}(P)$")
        axis = shared_axis
        lead_variant = row.variant_id
        pp4 = float(row.eqtl_coloc_pp4)
        n_shared = int(row.n_shared_snps)
        set_clean_axis(axis)
        for label in LD_ORDER:
            subset = shared.loc[shared["ld_bin"].eq(label)]
            if subset.empty:
                continue
            axis.scatter(
                subset["apa_neglog10p"],
                subset["eqtl_neglog10p"],
                s=30,
                color=LD_COLORS[label],
                edgecolor="#555555",
                linewidth=0.25,
                alpha=0.95 if label != "Missing" else 0.75,
                zorder=2 if label != "Missing" else 1,
            )
        lead = shared.loc[shared["variant_id"].eq(lead_variant)]
        if not lead.empty:
            axis.scatter(
                lead["apa_neglog10p"],
                lead["eqtl_neglog10p"],
                s=130,
                marker="D",
                color=LEAD_PURPLE,
                edgecolor="black",
                linewidth=0.8,
                zorder=4,
            )
        axis.set_title(
            "Shared apaQTL\u2013eQTL association",
            loc="left",
            fontsize=13,
            fontweight="bold",
            pad=8,
            color="black",
            gid="check_text",
        )
        axis.set_xlabel(r"$-\log_{10}(\mathrm{APA}\ P)$", fontsize=11, color=TEXT_DARK)
        axis.set_ylabel(r"$-\log_{10}(\mathrm{eQTL}\ P)$", fontsize=11, color=TEXT_DARK)
        axis.text(
            0.98,
            0.05,
            f"PP4 = {pp4:.3f}\n{n_shared:,} shared variants",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            linespacing=1.4,
            fontsize=8.5,
            color=TEXT_DARK,
            gid="check_text",
        )
        handles = [
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=LD_COLORS[label],
                   markeredgecolor="#555555", markersize=6.5, label=label)
            for label in ["0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0.0-0.2", "Missing"]
        ]
        legend = axis.legend(
            handles=handles,
            title=r"LD r$^2$",
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.02, 0.92),
            borderaxespad=0,
            fontsize=8,
            title_fontsize=8.5,
            labelspacing=1.2,
        )
        legend.get_title().set_gid("check_text")
        for text in legend.get_texts():
            text.set_gid("check_text")
        plt.setp(apa_axis.get_xticklabels(), visible=False)
        plt.setp(eqtl_axis.get_xticklabels(), visible=False)
        eqtl_axis.set_xlabel("")
        variant_label = str(row.display_variant_label)
        annotate_lead_label(apa_axis, case_apa, row.variant_id, "position", "neglog10p", variant_label)
        annotate_lead_label(eqtl_axis, case_eqtl, row.variant_id, "position", "neglog10p", variant_label)
        annotate_lead_label(
            shared_axis,
            shared,
            row.variant_id,
            "apa_neglog10p",
            "eqtl_neglog10p",
            variant_label,
            va="top",
            ha="right",
        )
    else:
        figure = plt.figure(figsize=(10.8, 5.8))
        grid = figure.add_gridspec(2, 1, height_ratios=[1.0, 1.25], hspace=0.74)
        apa_axis = figure.add_subplot(grid[0, 0])
        structure_axis = figure.add_subplot(grid[1, 0])
        add_panel_label(apa_axis, "a")
        add_panel_label(structure_axis, "b")
        plot_locus(apa_axis, case_apa, "neglog10p", row.variant_id, "apaQTL association", r"$-\log_{10}(P)$")
        annotate_lead_label(
            apa_axis,
            case_apa,
            row.variant_id,
            "position",
            "neglog10p",
            str(row.display_variant_label),
        )

    figure.suptitle(
        f"{row.gene_name} apaQTL\u2013eQTL evidence" if has_internal else f"{row.gene_name} apaQTL evidence",
        x=0.05,
        y=0.985,
        ha="left",
        fontsize=17,
        fontweight="bold",
        gid="check_text",
    )
    # Final layout before the structure lane (miRNA labels measure text boxes).
    figure.subplots_adjust(left=0.08, right=0.86 if has_internal else 0.91, top=0.88 if has_internal else 0.84, bottom=0.27)
    axis = structure_axis
    start = float(row.phenotype_start)
    end = float(row.phenotype_end)
    proximal = float(row.proximal_pas_pos)
    distal = float(row.distal_pas_pos)
    variant = parse_variant(row.variant_id)
    variant_position = float(variant["pos"])
    span = end - start

    # Top raised for the miRNA lane; bottom kept at -1.22 for the RBP lane.
    axis.set_ylim(-1.22, 1.78)
    axis.set_yticks([])
    for spine in ["left", "right", "top"]:
        axis.spines[spine].set_visible(False)
    axis.spines["bottom"].set_color(TEXT_MID)
    axis.tick_params(axis="x", labelsize=8, colors=TEXT_DARK, length=3, width=0.7)
    axis.xaxis.set_major_locator(MaxNLocator(nbins=2, integer=True))
    axis.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))

    if row.phenotype_strand == "-":
        shared_start, shared_end = proximal, end
        dur_start, dur_end = distal, proximal
        axis.set_xlim(end + 0.02 * span, start - 0.16 * span)
    else:
        shared_start, shared_end = start, proximal
        dur_start, dur_end = proximal, distal
        axis.set_xlim(start - 0.16 * span, end + 0.02 * span)

    axis.add_patch(Rectangle((min(shared_start, shared_end), 0.25), abs(shared_end - shared_start), 0.20, color=TRACK_GREY, linewidth=0))
    axis.add_patch(Rectangle((min(dur_start, dur_end), 0.22), abs(dur_end - dur_start), 0.26, color=DUR_GREEN, linewidth=0))
    axis.axvline(proximal, ymin=0.54, ymax=0.83, color=PROX_PURPLE, linewidth=2.2)
    axis.axvline(distal, ymin=0.54, ymax=0.83, color=DIST_BLUE, linewidth=2.2)
    axis.axvline(variant_position, ymin=0.54, ymax=0.83, color=VARIANT_RED, linewidth=2.2)

    proximal_label = pas_motif_label(
        pas.loc[
            pas["candidate_id"].eq(row.candidate_id) & pas["pas_label"].eq("proximal"), "ref_pas_motifs"
        ]
    )
    distal_label = pas_motif_label(
        pas.loc[
            pas["candidate_id"].eq(row.candidate_id) & pas["pas_label"].eq("distal"), "ref_pas_motifs"
        ]
    )
    variant_near_proximal = abs(variant_position - proximal) <= max(80.0, 0.03 * span)
    # Near-proximal variant label drops below the RBP rows.
    variant_label_y = -0.80 if variant_near_proximal else -0.07

    label_offset = max(0.065 * span, 240.0)
    proximal_label_x = proximal - label_offset if row.phenotype_strand == "+" else proximal + label_offset
    distal_label_x = distal + label_offset if row.phenotype_strand == "+" else distal - label_offset
    proximal_text = "proximal PAS" if proximal_label is None else f"proximal PAS\n{proximal_label}"
    distal_text = "distal PAS" if distal_label is None else f"distal PAS\n{distal_label}"
    proximal_annotation = axis.annotate(
        proximal_text,
        xy=(proximal, 0.92),
        xytext=(proximal_label_x, 0.76),
        textcoords="data",
        ha="right",
        va="center",
        fontsize=7.6,
        linespacing=1.25,
        color=PROX_PURPLE,
        arrowprops={"arrowstyle": "-", "color": PROX_PURPLE, "linewidth": 0.8, "shrinkA": 2, "shrinkB": 4},
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5, "alpha": 0.92},
    )
    distal_annotation = axis.annotate(
        distal_text,
        xy=(distal, 0.92),
        xytext=(distal_label_x, 0.76),
        textcoords="data",
        ha="left",
        va="center",
        fontsize=7.6,
        linespacing=1.25,
        color=DIST_BLUE,
        arrowprops={"arrowstyle": "-", "color": DIST_BLUE, "linewidth": 0.8, "shrinkA": 2, "shrinkB": 4},
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5, "alpha": 0.92},
    )
    proximal_annotation.set_gid("check_text")
    distal_annotation.set_gid("check_text")
    axis.text(
        variant_position,
        variant_label_y,
        str(row.display_variant_label),
        ha="left" if row.phenotype_strand == "-" else "right",
        va="top",
        fontsize=8,
        color=VARIANT_RED,
        fontweight="bold",
        gid="check_text",
    )
    axis.text(
        (shared_start + shared_end) / 2,
        0.05,
        "shared 3'UTR",
        ha="center",
        va="top",
        fontsize=8,
        color=TEXT_MID,
        gid="check_text",
    )
    axis.text(
        (dur_start + dur_end) / 2,
        0.05,
        f"distal 3'UTR {int(abs(dur_end - dur_start)):,} nt",
        ha="center",
        va="top",
        fontsize=8,
        color="#3F6B4E",
        gid="check_text",
    )
    axis.text(
        0.0,
        1.02,
        "APA structure and regulatory features",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=13,
        fontweight="bold",
        color="black",
        gid="check_text",
    )
    axis.hlines(-0.17, min(dur_start, dur_end), max(dur_start, dur_end), color="#E0E4E8", linewidth=1.2)

    top_rbp = rbp.loc[rbp["candidate_id"].eq(row.candidate_id)].copy()
    losses = top_rbp.loc[top_rbp["delta_direction"].eq("loss")].nsmallest(3, "delta_norm_score")
    gains = top_rbp.loc[top_rbp["delta_direction"].eq("gain")].nlargest(3, "delta_norm_score")
    top_rbp = pd.concat([losses, gains], ignore_index=True)
    top_rbp["rbp_label"] = top_rbp["rbp_expression_gene_symbol"].where(
        top_rbp["rbp_expression_gene_symbol"].notna(),
        top_rbp["rbp_name"].map(
            lambda name: str(name)[1:].split(")", 1)[0]
            if str(name).startswith("(") and ")" in str(name)
            else str(name).split("_", 1)[0]
        ),
    )
    dur_low = min(dur_start, dur_end)
    dur_high = max(dur_start, dur_end)
    dur_width = dur_high - dur_low
    loss_rows = top_rbp.loc[top_rbp["delta_direction"].eq("loss")].reset_index(drop=True)
    gain_rows = top_rbp.loc[top_rbp["delta_direction"].eq("gain")].reset_index(drop=True)
    loss_positions = np.linspace(dur_high - 0.18 * dur_width, dur_high - 0.44 * dur_width, max(len(loss_rows), 1))
    gain_positions = np.linspace(dur_high - 0.64 * dur_width, dur_high - 0.92 * dur_width, max(len(gain_rows), 1))

    for index, (item, x_pos) in enumerate(zip(loss_rows.itertuples(index=False), loss_positions)):
        label_y = -0.43 - 0.11 * (index % 2)
        axis.vlines(x_pos, -0.34, -0.17, color=LOSS_RED, linewidth=1.0)
        axis.scatter([x_pos], [-0.17], marker="D", s=28, color=LOSS_RED, zorder=3)
        axis.text(x_pos, label_y, item.rbp_label, ha="center", va="top", fontsize=7.5, color=LOSS_RED, fontweight="bold", gid="check_text")
    for index, (item, x_pos) in enumerate(zip(gain_rows.itertuples(index=False), gain_positions)):
        label_y = -0.43 - 0.11 * (index % 2)
        axis.vlines(x_pos, -0.34, -0.17, color=GAIN_GREEN, linewidth=1.0)
        axis.scatter([x_pos], [-0.17], marker="D", s=28, color=GAIN_GREEN, zorder=3)
        axis.text(x_pos, label_y, item.rbp_label, ha="center", va="top", fontsize=7.5, color=GAIN_GREEN, fontweight="bold", gid="check_text")

    # miRNA seed-site lane above the gene track; unplaceable labels are skipped.
    if case_mirna is not None and not case_mirna.empty:
        axis.hlines(1.12, min(dur_start, dur_end), max(dur_start, dur_end), color="#E0E4E8", linewidth=1.2)
        for site in case_mirna.itertuples(index=False):
            axis.vlines(site.site_x, 1.12, 1.26, color=MIRNA_SITE_COLORS[site.site_type], linewidth=1.3, zorder=3)
        mirna_handles = [
            Line2D(
                [0],
                [0],
                marker="v",
                linestyle="none",
                markerfacecolor=color,
                markeredgecolor=color,
                markersize=5.5,
                label=site_type,
            )
            for site_type, color in MIRNA_SITE_COLORS.items()
        ]
        # Measure after a fresh draw; unrendered text reports stale extents.
        axis.figure.canvas.draw()
        renderer = axis.figure.canvas.get_renderer()
        axis_window = axis.get_window_extent(renderer=renderer)
        for name in case_mirna["representative_mirna"].drop_duplicates().head(MAX_MIRNA_LABELS).tolist():
            x = float(case_mirna.loc[case_mirna["representative_mirna"].eq(name), "site_x"].median())
            placed = False
            for label_y in (1.34, 1.45):
                candidate = axis.text(
                    x,
                    label_y,
                    short_mirna_name(name),
                    ha="center",
                    va="bottom",
                    fontsize=7.5,
                    color=TEXT_DARK,
                    gid="check_text",
                )
                axis.figure.canvas.draw()
                window_box = candidate.get_window_extent(renderer=renderer)
                collides = any(
                    window_box.overlaps(other.get_window_extent(renderer=renderer))
                    for other in axis.figure.findobj(plt.Text)
                    if other.get_gid() == "check_text"
                    and other.get_visible()
                    and other.get_text()
                    and other is not candidate
                )
                if not collides and axis_window.x0 - 1 <= window_box.x0 and window_box.x1 <= axis_window.x1 + 1:
                    axis.vlines(x, 1.26, label_y - 0.04, color="#94A3B8", linewidth=0.7)
                    placed = True
                    break
                candidate.remove()
            if not placed:
                print(
                    f"[layout-check] miRNA label skipped (no collision-free row): {short_mirna_name(name)}",
                    file=sys.stderr,
                )

    axis.set_xlabel("Genomic coordinate (GRCh38)", fontsize=10, color=TEXT_DARK)
    if case_mirna is not None and not case_mirna.empty:
        mirna_legend = axis.figure.legend(
            handles=mirna_handles,
            title="miRNA seed site",
            loc="lower center",
            frameon=False,
            ncol=3,
            bbox_to_anchor=(0.38, 0.02),
            borderaxespad=0.0,
            fontsize=8,
            title_fontsize=8.5,
            labelspacing=0.45,
            handletextpad=0.55,
            columnspacing=1.1,
        )
        mirna_legend.get_title().set_gid("check_text")
        for text in mirna_legend.get_texts():
            text.set_gid("check_text")
    legend_handles = [
        Line2D([0], [0], marker="D", linestyle="-", color=LOSS_RED, markerfacecolor=LOSS_RED, markersize=5, label="Loss"),
        Line2D([0], [0], marker="D", linestyle="-", color=GAIN_GREEN, markerfacecolor=GAIN_GREEN, markersize=5, label="Gain"),
    ]
    legend = axis.figure.legend(
        handles=legend_handles,
        title="RBP motif",
        loc="lower center",
        frameon=False,
        ncol=2,
        bbox_to_anchor=(0.64, 0.02),
        borderaxespad=0.0,
        fontsize=8,
        title_fontsize=8.5,
        labelspacing=0.45,
        handletextpad=0.55,
    )
    legend.get_title().set_gid("check_text")
    for text in legend.get_texts():
        text.set_gid("check_text")
    base = row.figure_base
    output_png = FIGURE_DIR / f"{base}.png"
    output_pdf = FIGURE_DIR / f"{base}.pdf"
    output_svg = FIGURE_DIR / f"{base}.svg"
    figure.savefig(output_png, dpi=360, facecolor="white")
    figure.savefig(output_pdf, facecolor="white")
    figure.savefig(output_svg, facecolor="white")

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    canvas = figure.bbox
    checked = [
        text
        for text in figure.findobj(plt.Text)
        if text.get_gid() == "check_text" and text.get_visible() and text.get_text()
    ]
    boxes = [(text, text.get_window_extent(renderer=renderer)) for text in checked]
    # Pairwise overlap over all checked texts, across axes.
    text_overlap = 0
    for index, (first_text, first_box) in enumerate(boxes):
        for second_text, second_box in boxes[index + 1 :]:
            if first_box.overlaps(second_box):
                text_overlap += 1
                print(
                    f"[layout-check] text-text overlap: {first_text.get_text()!r} {first_box} "
                    f"vs {second_text.get_text()!r} {second_box}",
                    file=sys.stderr,
                )
    tick_boxes = []
    for axis in figure.axes:
        for ticker, limits in ((axis.xaxis, axis.get_xlim()), (axis.yaxis, axis.get_ylim())):
            lo, hi = min(limits), max(limits)
            for tick in ticker.get_major_ticks():
                # Skip ticks outside the axis limits; they never render.
                if not (lo <= tick.get_loc() <= hi):
                    continue
                for label in (tick.label1, tick.label2):
                    if not (label.get_visible() and label.get_text()):
                        continue
                    tick_boxes.append((axis, label, label.get_window_extent(renderer=renderer)))
    tick_overlap = 0
    for index, (first_axis, first_tick, first_box) in enumerate(tick_boxes):
        for second_axis, second_tick, second_box in tick_boxes[index + 1 :]:
            if first_axis is second_axis and first_box.overlaps(second_box):
                tick_overlap += 1
                print(
                    f"[layout-check] tick-tick overlap: {first_tick.get_text()!r} {first_box} "
                    f"vs {second_tick.get_text()!r} {second_box}",
                    file=sys.stderr,
                )
    text_tick_overlap = 0
    for text, text_box in boxes:
        for axis, tick, tick_box in tick_boxes:
            if text.axes is axis and text_box.overlaps(tick_box):
                text_tick_overlap += 1
                print(
                    f"[layout-check] text-tick overlap: {text.get_text()!r} {text_box} "
                    f"vs {tick.get_text()!r} {tick_box}",
                    file=sys.stderr,
                )
    clipped = sum(
        box.x0 < canvas.x0 or box.y0 < canvas.y0 or box.x1 > canvas.x1 or box.y1 > canvas.y1
        for _, box in boxes
    )
    checked_count = len(checked)

    layout_rows.append(
        {
            "figure": base,
            "gene_name": row.gene_name,
            "png_file": output_png.name,
            "pdf_file": output_pdf.name,
            "svg_file": output_svg.name,
            "width_in": figure.get_figwidth(),
            "height_in": figure.get_figheight(),
            "internal_support_panel": has_internal,
            "external_overlap_event_panel": has_external,
            "apa_locus_rows": int(len(case_apa)),
            "eqtl_locus_rows": int(len(case_eqtl)) if has_internal else 0,
            "shared_locus_rows": int(len(shared)) if has_internal else 0,
            "mirna_site_count": int(len(case_mirna)),
            "mirna_family_count": int(case_mirna["representative_mirna"].nunique()) if not case_mirna.empty else 0,
            "checked_text_count": checked_count,
            "text_overlap_count": text_overlap,
            "tick_label_overlap_count": tick_overlap,
            "text_tick_overlap_count": text_tick_overlap,
            "clipped_text_count": clipped,
            "files_exist": all(path.exists() and path.stat().st_size > 0 for path in [output_png, output_pdf, output_svg]),
            "automated_layout_status": "pass" if text_overlap == 0 and tick_overlap == 0 and text_tick_overlap == 0 and clipped == 0 else "fail",
            "human_review_status": "not_required",
            "human_review_note": "Automated text, tick, clipping, locus-row, and file checks completed.",
        }
    )
    plt.close(figure)

layout = pd.DataFrame(layout_rows)
source.to_csv(OUT_SOURCE, sep="\t", index=False)
layout.to_csv(OUT_LAYOUT, sep="\t", index=False)
if not layout["automated_layout_status"].eq("pass").all():
    raise RuntimeError("Case-figure layout validation failed")

print(layout.to_string(index=False))
