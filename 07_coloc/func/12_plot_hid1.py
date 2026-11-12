#!/usr/bin/env python3
# Plot the HID1 PAS-loss apaQTL: distal PAS track and genotype-stratified PDUI per life stage.
# Writes the figure (PNG, PDF, SVG), the plotted sample table, a stage summary and a layout check.

# --- packages ---
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
import pandas as pd
import pgenlib
import pysam


# --- config ---
ROOT = Path(os.environ["APA_ROOT"])
DATA = ROOT / "data/07_coloc"
FIGURES = ROOT / "figures/07_coloc"
PRE_INPUT = ROOT / "data/05_gxs/pre_input"

PHENOTYPE_ID = "ENST00000528902.1|HID1|chr17:74961597-74961989|-"
VARIANT_ID = "chr17:74961620:T:C"
PROXIMAL_PAS = 74961838
DISTAL_PAS = 74961597
EVENT_END = 74961989
TRACK_RIGHT = 74961562

SIGNAL_FILE = DATA / "02_signal_mechanisms.tsv.gz"
PAS_FILE = ROOT / "data/04_qtl/pas.tsv.gz"
PHENOTYPE_FILE = ROOT / "data/05_gxs/phenotype.bed.gz"
PGEN_FILE = PRE_INPUT / "genotype.pgen"
PVAR_FILE = PRE_INPUT / "genotype.pvar"
PSAM_FILE = PRE_INPUT / "genotype.psam"
STAGE_FILE = PRE_INPUT / "LifeStage.txt"

FIG_STEM = FIGURES / "13_case_hid1_pas_loss_mechanism"
SOURCE_FILE = DATA / "hid1_mechanism_samples.tsv.gz"
SUMMARY_FILE = DATA / "hid1_mechanism_summary.tsv"
LAYOUT_FILE = DATA / "hid1_mechanism_layout_validation.tsv"

FIGURE_WIDTH = 7.2
FIGURE_HEIGHT = 6.2
DPI = 420

INK = "#20262C"
MUTED = "#59636D"
GRID = "#DDE3E8"
SHARED_UTR = "#C8D0D8"
DISTAL_UTR = "#6FB276"
PROXIMAL_COLOR = "#7B3294"
DISTAL_COLOR = "#2C7FB8"
VARIANT_COLOR = "#C62828"
GENOTYPE_COLORS = ["#B9BEC4", "#56B4E9", "#0072B2"]
STAGE_ORDER = ["Prenatal", "Postnatal"]
GENOTYPE_ORDER = ["T/T", "T/C", "C/C"]


# --- helpers ---
def add_panel_label(axis, label):
    axis.text(
        -0.055,
        1.04,
        label,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color=INK,
        gid="layout_text",
    )


def plot_boxplot(axis, samples, stage, signal, panel):
    subset = samples.loc[samples["stage"].eq(stage)].copy()
    positions = np.arange(3)
    arrays = [
        subset.loc[subset["genotype"].eq(genotype), "normalized_pdui"].to_numpy()
        for genotype in GENOTYPE_ORDER
    ]
    nonempty = [(pos, values) for pos, values in zip(positions, arrays) if len(values)]
    box = axis.boxplot(
        [values for _, values in nonempty],
        positions=[pos for pos, _ in nonempty],
        widths=0.56,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": INK, "linewidth": 1.2},
        whiskerprops={"color": MUTED, "linewidth": 0.8},
        capprops={"color": MUTED, "linewidth": 0.8},
        boxprops={"edgecolor": MUTED, "linewidth": 0.8},
    )
    for patch, (position, _) in zip(box["boxes"], nonempty):
        patch.set_facecolor(GENOTYPE_COLORS[position])
        patch.set_alpha(0.72)

    rng = np.random.default_rng(1709 if stage == "Prenatal" else 1710)
    for position, values in zip(positions, arrays):
        jitter = rng.uniform(-0.17, 0.17, size=len(values))
        axis.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=10,
            facecolor=GENOTYPE_COLORS[position],
            edgecolor="white",
            linewidth=0.25,
            alpha=0.60,
            zorder=3,
        )

    counts = [len(values) for values in arrays]
    axis.set_xticks(positions, [f"{label}\nn={count}" for label, count in zip(GENOTYPE_ORDER, counts)])
    axis.set_xlim(-0.55, 2.55)
    axis.set_title(stage, loc="left", fontsize=10, fontweight="bold", color=INK, pad=10,
                   gid="layout_text")
    p_value = float(signal["pval_nominal"])
    p_text = f"{p_value:.2e}" if p_value < 1e-3 else f"{p_value:.3g}"
    axis.text(
        1.0,
        1.02,
        f"P = {p_text}; beta_C = {float(signal['beta_alt_PDUI']):.2f}",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color=MUTED,
        gid="layout_text",
    )
    add_panel_label(axis, panel)
    axis.grid(axis="y", color=GRID, linewidth=0.55, zorder=0)
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["left", "bottom"]].set_color(MUTED)
    axis.tick_params(axis="both", colors=INK, labelsize=7.5)
    axis.set_axisbelow(True)


# --- load signal rows ---
signal_columns = [
    "stage",
    "phenotype_id",
    "gene_name",
    "variant_id",
    "pval_nominal",
    "beta_alt_PDUI",
    "pip",
    "direct_pas_evidence",
    "distal_pas_core_change",
    "direct_pas_predicted_pdui_direction",
    "APARENT_predicted_PDUI_delta",
    "snp_to_apa_evidence_label",
]
signal_table = pd.read_csv(SIGNAL_FILE, sep="\t", usecols=signal_columns)
selected = signal_table.loc[
    signal_table["phenotype_id"].eq(PHENOTYPE_ID)
    & signal_table["variant_id"].eq(VARIANT_ID)
    & signal_table["stage"].isin(["prenatal", "postnatal"])
].copy()
if len(selected) != 2 or selected["stage"].nunique() != 2:
    raise RuntimeError("Expected exactly one HID1 mechanism row per life stage")
signals = selected.set_index("stage").loc[["prenatal", "postnatal"]].reset_index()

# --- load distal PAS change ---
pas = pd.read_csv(PAS_FILE, sep="\t")
pas_selected = pas.loc[
    pas["candidate_id"].str.contains(PHENOTYPE_ID, regex=False)
    & pas["candidate_id"].str.endswith(VARIANT_ID)
    & pas["pas_label"].eq("distal")
].drop_duplicates(
    subset=[
        "pas_window_ref_seq",
        "pas_window_alt_seq",
        "ref_pas_motifs",
        "pas_core_change",
    ]
)
if len(pas_selected) != 1:
    raise RuntimeError("Expected one unique HID1 distal-PAS sequence change")
pas_row = pas_selected.iloc[0]

# --- load phenotype ---
with pysam.TabixFile(str(PHENOTYPE_FILE)) as tabix:
    header = list(tabix.header)[-1].lstrip("#").split("\t")
    records = [
        record.split("\t")
        for record in tabix.fetch("17", 74961596, 74961989)
        if PHENOTYPE_ID in record
    ]
if len(records) != 1:
    raise RuntimeError("Expected one HID1 phenotype row in the indexed BED")
row = records[0]
phenotype = pd.DataFrame(
    {
        "sample": header[4:],
        "normalized_pdui": pd.to_numeric(pd.Series(row[4:]), errors="coerce").to_numpy(),
    }
)

# --- load genotype ---
variant_index = None
index = 0
with PVAR_FILE.open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 3 and fields[2] == VARIANT_ID:
            variant_index = index
            break
        index += 1
if variant_index is None:
    raise RuntimeError(f"Variant absent from PVAR: {VARIANT_ID}")

psam = pd.read_csv(PSAM_FILE, sep="\t", dtype=str)
reader = pgenlib.PgenReader(bytes(str(PGEN_FILE), "utf8"))
try:
    if reader.get_raw_sample_ct() != len(psam):
        raise RuntimeError("PGEN and PSAM sample counts differ")
    dosage = np.empty(len(psam), dtype=np.float32)
    reader.read_dosages(variant_index, dosage)
finally:
    reader.close()
dosage = dosage.astype(float)
dosage[dosage < 0] = np.nan
hardcall = np.where(
    dosage < 0.5,
    0,
    np.where(dosage < 1.5, 1, 2),
)
labels = pd.Series(hardcall).map(dict(enumerate(GENOTYPE_ORDER)))
genotype = pd.DataFrame(
    {
        "sample": psam["#IID"].astype(str),
        "alt_c_dosage": dosage,
        "genotype": labels,
    }
)

# --- assemble samples ---
stage = pd.read_csv(STAGE_FILE, sep="\t", dtype={"sample": str})
stage["stage"] = stage["LifeStage"].map({0: "Prenatal", 1: "Postnatal"})
samples = phenotype.merge(genotype, on="sample", validate="one_to_one")
samples = samples.merge(stage[["sample", "stage"]], on="sample", validate="one_to_one")
samples = samples.dropna(subset=["normalized_pdui", "alt_c_dosage", "genotype", "stage"]).copy()
samples["stage"] = pd.Categorical(samples["stage"], STAGE_ORDER, ordered=True)
samples["genotype"] = pd.Categorical(samples["genotype"], GENOTYPE_ORDER, ordered=True)
samples = samples.sort_values(["stage", "genotype", "sample"], kind="stable")
samples.to_csv(SOURCE_FILE, sep="\t", index=False)

# --- stage summary ---
summary = (
    samples.groupby(["stage", "genotype"], observed=False)["normalized_pdui"]
    .agg(n="size", median="median", mean="mean", sd="std")
    .reset_index()
)
signal_summary = signals.rename(columns={"stage": "stage_key"}).copy()
signal_summary["stage"] = signal_summary["stage_key"].str.capitalize()
summary = summary.merge(
    signal_summary[["stage", "pval_nominal", "beta_alt_PDUI", "pip", "APARENT_predicted_PDUI_delta"]],
    on="stage",
    how="left",
    validate="many_to_one",
)
summary["phenotype_id"] = PHENOTYPE_ID
summary["variant_id"] = VARIANT_ID
summary["phenotype_source"] = str(PHENOTYPE_FILE)
summary["genotype_source"] = str(PGEN_FILE)
summary.to_csv(SUMMARY_FILE, sep="\t", index=False)

# --- figure ---
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)
figure = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), facecolor="white")
grid = figure.add_gridspec(2, 2, height_ratios=[0.88, 1.35], hspace=0.54, wspace=0.26)
track_axis = figure.add_subplot(grid[0, :])
prenatal_axis = figure.add_subplot(grid[1, 0])
postnatal_axis = figure.add_subplot(grid[1, 1], sharey=prenatal_axis)

# PAS track
aparent_delta = float(signals["APARENT_predicted_PDUI_delta"].dropna().iloc[0])
track_axis.set_xlim(EVENT_END + 10, TRACK_RIGHT)
track_axis.set_ylim(0.0, 1.35)
track_axis.set_yticks([])
track_axis.spines[["top", "left", "right"]].set_visible(False)
track_axis.spines["bottom"].set_color(MUTED)
track_axis.tick_params(axis="x", colors=MUTED, labelsize=7, length=3)
track_axis.ticklabel_format(axis="x", style="plain", useOffset=False)
track_axis.set_xticks([74961950, 74961850, 74961750, 74961650])
track_axis.set_xlabel("Genomic coordinate (GRCh38; minus strand)", fontsize=8, color=INK, labelpad=6)
track_axis.set_title(
    "HID1 APA structure and PAS variant",
    loc="left",
    fontsize=10,
    fontweight="bold",
    color=INK,
    pad=10,
    gid="layout_text",
)
add_panel_label(track_axis, "a")

track_y = 0.60
track_h = 0.14
track_axis.add_patch(Rectangle((PROXIMAL_PAS, track_y), EVENT_END - PROXIMAL_PAS, track_h,
                               facecolor=SHARED_UTR, edgecolor="none"))
track_axis.add_patch(Rectangle((DISTAL_PAS, track_y), PROXIMAL_PAS - DISTAL_PAS, track_h,
                               facecolor=DISTAL_UTR, edgecolor="none"))
track_axis.plot([EVENT_END, DISTAL_PAS], [track_y + track_h / 2] * 2, color=INK, linewidth=0.7, zorder=0)
for xpos in np.linspace(EVENT_END - 18, DISTAL_PAS + 14, 13):
    track_axis.plot([xpos, xpos - 8], [track_y + track_h / 2, track_y + track_h / 2],
                    color=INK, linewidth=0.8)

track_axis.vlines(PROXIMAL_PAS, track_y - 0.08, track_y + track_h + 0.22,
                  color=PROXIMAL_COLOR, linewidth=2.2)
track_axis.vlines(DISTAL_PAS, track_y - 0.08, track_y + track_h + 0.22,
                  color=DISTAL_COLOR, linewidth=2.2)
track_axis.scatter([74961620], [0.99], marker="D", s=58, facecolor=VARIANT_COLOR,
                   edgecolor="white", linewidth=0.7, zorder=4)
track_axis.vlines(74961620, track_y + track_h, 0.95, color=VARIANT_COLOR,
                  linewidth=1.0, linestyles="dotted")

track_axis.text((EVENT_END + PROXIMAL_PAS) / 2, track_y - 0.10, "shared 3' UTR",
                ha="center", va="top", fontsize=7.5, color=MUTED, gid="layout_text")
track_axis.text((PROXIMAL_PAS + DISTAL_PAS) / 2, track_y - 0.10, "distal 3' UTR 241 nt",
                ha="center", va="top", fontsize=7.5, color="#3E7145", gid="layout_text")
track_axis.text(PROXIMAL_PAS, 0.88, "proximal PAS", ha="center", va="bottom",
                fontsize=7.5, color=PROXIMAL_COLOR, gid="layout_text")
track_axis.text(DISTAL_PAS, 0.88, "distal PAS", ha="center", va="bottom",
                fontsize=7.5, color=DISTAL_COLOR, gid="layout_text")
track_axis.text(74961620, 1.10, "chr17:74961620 T>C", ha="center", va="bottom",
                fontsize=8, fontweight="bold", color=VARIANT_COLOR, gid="layout_text")

ref_motif = str(pas_row["ref_pas_motifs"])
alt_motif = "AATAGG"
track_axis.text(0.01, 0.08, f"T allele   {ref_motif}   PAS",
                transform=track_axis.transAxes, ha="left", va="bottom", fontsize=8,
                family="monospace", color=INK, gid="layout_text")
track_axis.text(0.01, 0.005, f"C allele   {alt_motif}   loss",
                transform=track_axis.transAxes, ha="left", va="bottom", fontsize=8,
                family="monospace", color=VARIANT_COLOR, gid="layout_text")
track_axis.text(0.99, 0.045, f"APARENT delta PDUI = {aparent_delta:.3f}",
                transform=track_axis.transAxes, ha="right", va="bottom", fontsize=7.5,
                color=MUTED, gid="layout_text")

plot_boxplot(prenatal_axis, samples, "Prenatal", signals.iloc[0], "b")
plot_boxplot(postnatal_axis, samples, "Postnatal", signals.iloc[1], "c")
prenatal_axis.set_ylabel("Normalized PDUI", fontsize=8.5, color=INK)
postnatal_axis.tick_params(axis="y", labelleft=False)
figure.supxlabel(f"Genotype at {VARIANT_ID}", x=0.52, y=0.035, fontsize=8.5, color=INK)
figure.suptitle(
    "HID1 PAS-loss apaQTL",
    x=0.09,
    y=0.958,
    ha="left",
    fontsize=15,
    fontweight="bold",
    color=INK,
    gid="layout_text",
)
figure.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.14)

DATA.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)
figure.savefig(FIG_STEM.with_suffix(".png"), dpi=DPI, facecolor="white")
figure.savefig(FIG_STEM.with_suffix(".pdf"), facecolor="white")
figure.savefig(FIG_STEM.with_suffix(".svg"), facecolor="white")

# --- layout check ---
figure.canvas.draw()
renderer = figure.canvas.get_renderer()
figure_box = figure.bbox
texts = [
    text for text in figure.findobj(match=Text)
    if text.get_visible() and text.get_text().strip()
]
boxes = [(text, text.get_window_extent(renderer=renderer)) for text in texts]
clipped_labels = [
    text.get_text().replace("\n", "/")
    for text, box in boxes
    if box.x0 < figure_box.x0 - 1
    or box.y0 < figure_box.y0 - 1
    or box.x1 > figure_box.x1 + 1
    or box.y1 > figure_box.y1 + 1
]
clipped = len(clipped_labels)
tagged = [(text, box) for text, box in boxes if text.get_gid() == "layout_text"]
overlaps = 0
for first_index, (_, first_box) in enumerate(tagged):
    for _, second_box in tagged[first_index + 1:]:
        if first_box.overlaps(second_box):
            overlaps += 1

genotype_counts = (
    samples.groupby(["stage", "genotype"], observed=False)
    .size()
    .reindex(pd.MultiIndex.from_product([STAGE_ORDER, GENOTYPE_ORDER]), fill_value=0)
)
files = [FIG_STEM.with_suffix(suffix) for suffix in [".png", ".pdf", ".svg"]]
status = "pass" if overlaps == 0 and clipped == 0 and all(path.exists() and path.stat().st_size > 0 for path in files) else "fail"
layout = pd.DataFrame(
    [
        {
            "figure": FIG_STEM.name,
            "width_in": figure.get_figwidth(),
            "height_in": figure.get_figheight(),
            "stage_panels": 2,
            "genotype_groups_expected": 6,
            "genotype_groups_present": int((genotype_counts > 0).sum()),
            "sample_n": len(samples),
            "visible_text_count": len(texts),
            "checked_annotation_count": len(tagged),
            "annotation_overlap_count": overlaps,
            "clipped_text_count": clipped,
            "clipped_labels": "|".join(clipped_labels),
            "minimum_text_pt": min(text.get_fontsize() for text in texts),
            "png_pdf_svg_exist": all(path.exists() and path.stat().st_size > 0 for path in files),
            "automated_layout_status": status,
        }
    ]
)
layout.to_csv(LAYOUT_FILE, sep="\t", index=False)
if not layout["automated_layout_status"].eq("pass").all():
    raise RuntimeError(f"HID1 figure layout validation failed:\n{layout.to_string(index=False)}")
plt.close(figure)
