#!/usr/bin/env python3
# Draw the ATP2B1 aUTR figure of miRNA 8mer sites and RBP motifs.
# Writes PNG, PDF and SVG figures and a layout validation table.

# --- packages ---
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


# --- config ---
DATA = Path(os.environ["APA_ROOT"]) / "data/07_coloc"
FIGURES = Path(os.environ["APA_ROOT"]) / "figures/07_coloc"
ANNOTATION = DATA / "12_autr_mirna_rbp_annotation.tsv"
OUT_STEM = FIGURES / "12_autr_atp2b1_sequence_candidates"
FIGSIZE = (10, 3.4)
ATP2B1_SITE = "ENST00000551310.1|ATP2B1|chr12:89655583-89656107|-"
AUTR_LENGTH = 525
MIRNA_COLOR = "#0F766E"
RBP_COLOR = "#7B3FA1"
BAR_COLOR = "#D8DEE6"


# --- loading data ---
annotations = pd.read_csv(ANNOTATION, sep="\t")
atp2b1 = annotations.loc[annotations["site_id"].eq(ATP2B1_SITE)]
mirna = atp2b1.loc[
    atp2b1["annotation_type"].eq("mirna_seed_match")
    & atp2b1["strict_mirna_candidate"].astype(str).str.lower().eq("true")
].copy()
rbp = atp2b1.loc[atp2b1["annotation_type"].eq("rbp_motif")].copy()
assert len(mirna) == 10 and len(rbp) == 4


mirna["x"] = mirna["positions"].map(
    lambda value: int([part for part in str(value).split(";") if part.startswith("8mer")][0].split(":")[1].split(",")[0])
)
mirna["host"] = mirna["host_genes"].map(lambda host_text: next(
    (host for host in str(host_text).split(";") if host and not host.startswith("ENSG")),
    str(host_text).split(";")[0],
))
rbp["x_list"] = rbp["positions"].map(lambda value: [int(v) for v in str(value).split(",")])


# --- figure ---
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})
fig, ax = plt.subplots(figsize=FIGSIZE)
ax.set_xlim(0, AUTR_LENGTH + 1)
ax.set_ylim(-2.4, 3.1)

ax.add_patch(plt.Rectangle((1, -0.14), AUTR_LENGTH - 1, 0.28, facecolor=BAR_COLOR, edgecolor="#AAB7C4", lw=0.8))
for boundary, label in ((1, "proximal cleavage site"), (AUTR_LENGTH, "distal cleavage site")):
    ax.vlines(boundary, -0.35, 0.35, color="#52606D", lw=1.4)
    ax.annotate(
        label,
        xy=(boundary, -0.42),
        xytext=(boundary, -0.75),
        ha="left" if boundary == 1 else "right",
        fontsize=8,
        color="#52606D",
        arrowprops=dict(arrowstyle="-", lw=0.6, color="#AAB7C4"),
    )

for index, (_, row) in enumerate(mirna.sort_values("x").iterrows()):
    level = (1.05, 1.7, 2.35)[index % 3]
    ax.vlines(row["x"], 0.16, level - 0.12, color=MIRNA_COLOR, lw=1.0, alpha=0.75)
    ax.plot(row["x"], level - 0.12, marker="o", ms=4.5, color=MIRNA_COLOR, clip_on=False)
    ax.text(row["x"], level + 0.05, f"{row['annotation_name'].replace('hsa-', '')}\nhost {row['host']}",
            ha="center", va="bottom", fontsize=7.6, color=MIRNA_COLOR, linespacing=1.25)

rbp_levels = {"RBM5": -1.05, "FUS": -1.75, "HNRNPR": -1.05, "EIF2S1": -1.75}
for _, row in rbp.iterrows():
    level = rbp_levels[row["annotation_name"]]
    for position in row["x_list"]:
        ax.vlines(position, -0.16, level + 0.20, color=RBP_COLOR, lw=1.0, alpha=0.75)
        ax.plot(position, level + 0.20, marker="^", ms=4.5, color=RBP_COLOR, clip_on=False)
    ax.text(
        row["x_list"][0],
        level - 0.05,
        f"{row['annotation_name']}\nTPM {row['prenatal_mean_tpm']:.0f}/{row['postnatal_mean_tpm']:.0f}",
        ha="center", va="top", fontsize=7.6, color=RBP_COLOR, linespacing=1.25,
    )

ax.text(
    -0.012, 1.65, "miRNA\n8mer sites", ha="right", va="bottom", fontsize=8.5,
    color=MIRNA_COLOR, fontweight="bold", transform=ax.get_yaxis_transform(),
)
ax.text(
    -0.012, -1.05, "RBP\nconsensus", ha="right", va="center", fontsize=8.5,
    color=RBP_COLOR, fontweight="bold", transform=ax.get_yaxis_transform(),
)
ax.set_yticks([])
ax.set_xticks([1, 100, 200, 300, 400, 500, 525])
ax.set_xlabel("aUTR position (nt, transcript orientation) · ATP2B1 ENST00000551310.1 · chr12:89655583-89656107 · minus strand")
ax.spines[["top", "right", "left"]].set_visible(False)

fig.subplots_adjust(left=0.13, right=0.985, top=0.97, bottom=0.20)
fig.savefig(f"{OUT_STEM}.png", dpi=300)
fig.savefig(f"{OUT_STEM}.pdf")
fig.savefig(f"{OUT_STEM}.svg")


# --- layout check ---
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
figure_bbox = fig.bbox
clipped = []
for text_object in ax.texts:
    bbox = text_object.get_window_extent(renderer=renderer)
    if not figure_bbox.contains(bbox.x0, bbox.y0) or not figure_bbox.contains(bbox.x1, bbox.y1):
        clipped.append(text_object.get_text())
assert not clipped, f"Text outside figure bounds: {clipped}"
text_bboxes = [text_object.get_window_extent(renderer=renderer) for text_object in ax.texts]
overlaps = sum(
    1
    for first in range(len(text_bboxes))
    for second in range(first + 1, len(text_bboxes))
    if text_bboxes[first].overlaps(text_bboxes[second])
)
assert overlaps == 0, f"{overlaps} overlapping label pairs"
validation = DATA / "12_autr_atp2b1_sequence_candidates_layout_validation.tsv"
validation.write_text(
    "check\tvalue\n"
    f"figure_width_in\t{FIGSIZE[0]}\n"
    f"figure_height_in\t{FIGSIZE[1]}\n"
    f"strict_mirna_features\t{len(mirna)}\n"
    f"rbp_features\t{len(rbp)}\n"
    f"clipped_text\t{len(clipped)}\n"
    f"overlapping_label_pairs\t{overlaps}\n"
    "status\tPASS\n",
    encoding="utf-8",
)
