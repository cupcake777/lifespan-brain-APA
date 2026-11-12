#!/usr/bin/env python3
# HID1 distal PAS mechanism figure: PAS core sequence (a) and APA evidence track (b).
# Figure goes to figures/07_coloc, layout check table to data/07_coloc.

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --- config ---
APA_ROOT = Path(os.environ["APA_ROOT"])
FIG_DIR = APA_ROOT / "figures/07_coloc"
DATA_DIR = APA_ROOT / "data/07_coloc"
FIG_STEM = FIG_DIR / "HID1_sequence_mechanism"
FIGSIZE = (10.8, 5.7)

NAVY = "#172B45"
SLATE = "#4B5563"
LIGHT = "#D7DEE7"
PURPLE = "#7B3FA1"
GREEN = "#0B7D5B"
RED = "#C7253E"

VARIANT_POS = 74961620
PAS_POS = 74961597
REF_MOTIF = "AATAGA"
ALT_MOTIF = "AATAGG"


# --- helpers ---
def draw_sequence(ax, y, label, motif, motif_color, detail):
    ax.text(0.02, y, label, transform=ax.transAxes, fontsize=10, fontweight="bold",
            color=motif_color, va="center")
    ax.text(0.12, y, "...TAGAA", transform=ax.transAxes, fontsize=12,
            family="DejaVu Sans Mono", color=NAVY, va="center")
    ax.text(0.31, y, f"[{motif}]", transform=ax.transAxes, fontsize=12,
            family="DejaVu Sans Mono", fontweight="bold", color=motif_color,
            va="center", bbox=dict(facecolor=motif_color, alpha=0.10,
                                   edgecolor="none", pad=2.5))
    ax.text(0.51, y, "ATGCAA...", transform=ax.transAxes, fontsize=12,
            family="DejaVu Sans Mono", color=NAVY, va="center")
    ax.text(0.98, y, detail, transform=ax.transAxes, fontsize=9.5,
            color=motif_color, va="center", ha="right")


# --- figure setup ---
FIG_DIR.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

fig = plt.figure(figsize=FIGSIZE, facecolor="white")
grid = fig.add_gridspec(2, 1, height_ratios=[0.9, 1.25], hspace=0.38,
                        left=0.08, right=0.97, top=0.79, bottom=0.16)
ax_seq = fig.add_subplot(grid[0])
ax_track = fig.add_subplot(grid[1])

fig.text(0.05, 0.955, "HID1 apaQTL sequence mechanism", fontsize=21,
         fontweight="bold", color=NAVY, va="top")
fig.text(0.05, 0.895,
         "chr17:74,961,620 T>C  |  minus strand  |  transcript allele A>G  |  ENST00000528902.1",
         fontsize=10.5, color=SLATE, va="top")

# --- panel a ---
# PAS core sequence, REF and ALT motifs
ax_seq.set_axis_off()
ax_seq.text(-0.045, 1.08, "a", transform=ax_seq.transAxes, fontsize=13,
            fontweight="bold", va="top")
ax_seq.text(0.00, 1.08, "Distal PAS core sequence", transform=ax_seq.transAxes,
            fontsize=14, fontweight="bold", va="top")

draw_sequence(ax_seq, 0.62, "REF", REF_MOTIF, PURPLE,
              "PAS-list hexamer")
draw_sequence(ax_seq, 0.25, "ALT", ALT_MOTIF, RED,
              "core motif loss")
ax_seq.annotate("", xy=(0.46, 0.35), xytext=(0.46, 0.53),
                xycoords=ax_seq.transAxes,
                arrowprops=dict(arrowstyle="-|>", color=SLATE, lw=1.2))
ax_seq.text(0.485, 0.44, "A>G", transform=ax_seq.transAxes, fontsize=9.5,
            color=RED, fontweight="bold", va="center")
ax_seq.plot([0, 1], [0.02, 0.02], transform=ax_seq.transAxes,
            color=LIGHT, lw=0.8, clip_on=False)

# --- panel b ---
# APA structure and regulatory evidence
ax_track.text(-0.045, 1.16, "b", transform=ax_track.transAxes, fontsize=13,
              fontweight="bold", va="top")
ax_track.text(0.00, 1.16, "APA structure and regulatory evidence",
              transform=ax_track.transAxes, fontsize=14, fontweight="bold",
              va="top")

ax_track.set_xlim(74961635, 74961585)
ax_track.set_ylim(-1.15, 1.05)
ax_track.set_yticks([])
ax_track.spines["left"].set_visible(False)
ax_track.spines["bottom"].set_color(SLATE)
ax_track.tick_params(axis="x", colors=SLATE, labelsize=8)
ax_track.set_xticks([74961630, 74961620, 74961610, 74961600, 74961590])
ax_track.ticklabel_format(style="plain", axis="x", useOffset=False)
ax_track.set_xlabel("Genomic coordinate (GRCh38; transcript runs rightward)",
                    fontsize=9.5, color=SLATE, labelpad=7)

ax_track.plot([74961633, 74961588], [0.36, 0.36], color=LIGHT, lw=8,
              solid_capstyle="butt", zorder=1)
ax_track.plot([VARIANT_POS, PAS_POS], [0.36, 0.36], color="#78B875", lw=8,
              solid_capstyle="butt", zorder=2)
ax_track.text((VARIANT_POS + PAS_POS) / 2, 0.13,
              "PAS core-to-cleavage interval", color=SLATE, fontsize=8.5,
              ha="center")

ax_track.vlines(VARIANT_POS, 0.22, 0.72, color=RED, lw=3, zorder=4)
ax_track.scatter([VARIANT_POS], [0.36], marker="D", s=42, color=RED,
                 edgecolor="white", linewidth=0.6, zorder=5)
ax_track.annotate("variant\n74,961,620 T>C",
                  xy=(VARIANT_POS, 0.70), xytext=(VARIANT_POS + 5.8, 0.94),
                  fontsize=9, fontweight="bold", color=RED, ha="center",
                  arrowprops=dict(arrowstyle="-", color=RED, lw=1.1))

ax_track.vlines(PAS_POS, 0.22, 0.72, color=PURPLE, lw=3, zorder=4)
ax_track.annotate("distal PAS\n74,961,597",
                  xy=(PAS_POS, 0.70), xytext=(PAS_POS - 5.0, 0.94),
                  fontsize=9, fontweight="bold", color=PURPLE, ha="center",
                  arrowprops=dict(arrowstyle="-", color=PURPLE, lw=1.1))
ax_track.text((VARIANT_POS + PAS_POS) / 2, 0.55, "23 nt", fontsize=8.5,
              color=SLATE, ha="center", va="bottom")

ax_track.plot([VARIANT_POS, VARIANT_POS], [0.18, -0.04], color=LIGHT, lw=1.2)
ax_track.text(0.09, 0.42, "Variant-local RBP motif score", transform=ax_track.transAxes,
              color=SLATE, fontsize=8.5, ha="left")
ax_track.scatter([0.10, 0.10], [0.26, 0.11], transform=ax_track.transAxes,
                 marker="D", s=34, color=[GREEN, RED], zorder=4)
ax_track.text(0.13, 0.26, "FUS  +0.250 gain", transform=ax_track.transAxes,
              color=GREEN, fontsize=8.7, fontweight="bold", ha="left", va="center")
ax_track.text(0.13, 0.11, "KHDRBS1  −0.173 loss", transform=ax_track.transAxes,
              color=RED, fontsize=8.7, fontweight="bold", ha="left", va="center")

ax_track.text(0.58, 0.39, "Distal PAS usage ↓  →  shorter 3′ UTR / PDUI",
              transform=ax_track.transAxes, color=NAVY, fontsize=10.2,
              fontweight="bold", ha="left")
ax_track.text(0.58, 0.23, "prenatal  PIP 0.401  ·  shortening",
              transform=ax_track.transAxes, color=RED, fontsize=9.2, ha="left")
ax_track.text(0.58, 0.09, "postnatal  PIP 0.466  ·  shortening",
              transform=ax_track.transAxes, color=RED, fontsize=9.2, ha="left")

fig.text(0.08, 0.035,
         "Direct-exploratory · RNAfold ΔMFE −1.3 kcal/mol · POSTAR3 overlap 0 · causal cleavage and binding effects remain untested",
         fontsize=8.2, color=SLATE, va="center")

# --- save ---
fig.savefig(f"{FIG_STEM}.png", dpi=300, facecolor="white")
fig.savefig(f"{FIG_STEM}.pdf", facecolor="white")
fig.savefig(f"{FIG_STEM}.svg", facecolor="white")

# layout check: no text outside the canvas
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
figure_bbox = fig.bbox
texts = [text for ax in [ax_seq, ax_track] for text in ax.texts if text.get_text()]
clipped = []
for text in texts:
    bbox = text.get_window_extent(renderer=renderer)
    if (bbox.x0 < figure_bbox.x0 or bbox.y0 < figure_bbox.y0 or
            bbox.x1 > figure_bbox.x1 or bbox.y1 > figure_bbox.y1):
        clipped.append(text.get_text())
assert not clipped, f"Text outside figure bounds: {clipped}"

validation = DATA_DIR / "HID1_sequence_mechanism_layout_validation.tsv"
validation.write_text(
    "check\tvalue\n"
    f"figure_width_in\t{FIGSIZE[0]}\n"
    f"figure_height_in\t{FIGSIZE[1]}\n"
    f"text_objects\t{len(texts)}\n"
    f"clipped_text\t{len(clipped)}\n"
    "status\tPASS\n",
    encoding="utf-8",
)
