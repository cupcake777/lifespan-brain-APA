# Fig 4 panel c: gating and saturation mechanism schematics (conceptual; no data values on the curves).
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon
import matplotlib.image as mpimg

for f in ["arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"]:
    fm.fontManager.addfont(os.path.join(os.environ["APA_FONT_DIR"], f))
plt.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                     "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "pdf.fonttype": 42})
FIG_DIR = os.path.join(os.environ["APA_ROOT"], "figures")
HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, GREY_CDS, GREY_UTR, GREY_LT = "#222222", "#5A5A5A", "#6E6E6E", "#C9C9C9", "#E6E6E6"
PRE, POST = "#E3A73F", "#6BAED6"
PRE_TXT, POST_TXT = "#B07D22", "#3C84B8"
W_MM, H_MM = 96, 52
mm = 1 / 25.4
fig = plt.figure(figsize=(W_MM * mm, H_MM * mm), dpi=600)

# switch icons; the saturated green of the ON icon is recoloured to the postnatal blue
ICON = os.path.join(HERE, "assets")
sw_on, sw_off = mpimg.imread(os.path.join(ICON, "switch_on.png")), mpimg.imread(os.path.join(ICON, "switch_off.png"))
_g = (sw_on[..., 1] > sw_on[..., 0] + 0.15) & (sw_on[..., 1] > sw_on[..., 2] + 0.05)
sw_on = sw_on.copy(); sw_on[_g, 0], sw_on[_g, 1], sw_on[_g, 2] = 0x6B / 255, 0xAE / 255, 0xD6 / 255


def gene_model(ax, x0, y, utr_len, h=0.9, dim=False, pas=None):
    """Promoter arrow, CDS and grey 3' UTR; pas lists (x, label) cleavage marks."""
    a = 0.35 if dim else 1.0
    ax.plot([x0 - 0.8, x0 + 3.2 + utr_len + 0.4], [y, y], color=MUTED, lw=0.6, alpha=a, solid_capstyle="butt")
    ax.plot([x0 - 0.6, x0 - 0.6, x0 - 0.1], [y, y + 0.9, y + 0.9], color=INK, lw=0.7, alpha=a)
    ax.add_patch(Polygon([[x0 - 0.15, y + 0.65], [x0 - 0.15, y + 1.15], [x0 + 0.2, y + 0.9]], closed=True, color=INK, alpha=a, lw=0))
    ax.add_patch(FancyBboxPatch((x0, y - h / 2), 3.0, h, boxstyle="round,pad=0,rounding_size=0.15", fc=GREY_CDS, ec="none", alpha=a))
    ax.add_patch(FancyBboxPatch((x0 + 3.0, y - h * 0.32), utr_len, h * 0.64, boxstyle="round,pad=0,rounding_size=0.12",
                                fc=GREY_UTR, ec="none", alpha=a))
    for xp, lab in (pas or []):
        ax.add_patch(Polygon([[xp - 0.22, y + h * 0.95], [xp + 0.22, y + h * 0.95], [xp, y + h * 0.42]], closed=True, color=INK, lw=0))
        ax.text(xp, y + h * 1.15, lab, ha="center", va="bottom", fontsize=5.5, color=INK)


def mini_axes(rect):
    ax = fig.add_axes(rect)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    ax.set_xlim(-0.3, 2.3); ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["0", "1", "2"])
    ax.set_xlabel("Allele dosage", labelpad=1); ax.tick_params(labelsize=6, pad=1)
    return ax


def eg(xr, gene):
    tg = fig.text(xr, 0.95, gene, fontsize=6, style="italic", color=MUTED, va="top", ha="right")
    fig.canvas.draw(); bb = tg.get_window_extent(); x_left = fig.transFigure.inverted().transform((bb.x0, bb.y0))[0]
    fig.text(x_left - 0.004, 0.95, "e.g.", fontsize=6, color=MUTED, va="top", ha="right")


# --- gating ---
axg = fig.add_axes([0.02, 0.60, 0.46, 0.30]); axg.set_xlim(0, 20); axg.set_ylim(0, 6); axg.axis("off")
fig.text(0.02, 0.95, "Gating", fontsize=7, fontweight="bold", color=INK, va="top")
eg(0.47, "TESPA1")
gene_model(axg, 2.0, 4.3, 4.5, dim=True)
axg.imshow(sw_off, extent=(13.0, 16.6, 3.55, 5.05), aspect="auto", zorder=3)
axg.text(17.1, 4.3, "Prenatal", va="center", fontsize=6, color=PRE_TXT)
gene_model(axg, 2.0, 1.4, 4.5)
axg.imshow(sw_on, extent=(13.0, 16.6, 0.65, 2.15), aspect="auto", zorder=3)
axg.text(17.1, 1.4, "Postnatal", va="center", fontsize=6, color=POST_TXT)
a1 = mini_axes([0.09, 0.14, 0.29, 0.40])
a1.plot([0, 1, 2], [0.45, 0.45, 0.45], color=PRE, lw=1.0, ls=(0, (2.5, 1.5)), marker="o", ms=2.6)
a1.plot([0, 1, 2], [0.25, 0.50, 0.75], color=POST, lw=1.0, marker="o", ms=2.6)
a1.set_ylim(0, 1); a1.set_yticks([]); a1.set_ylabel("PDUI", labelpad=2)
a1.text(2.25, 0.75, "Effect\nvisible", fontsize=5.5, color=POST_TXT, va="center", ha="left", linespacing=0.95)
a1.text(2.25, 0.42, "No\neffect", fontsize=5.5, color=PRE_TXT, va="center", ha="left", linespacing=0.95)

# --- saturation ---
axs = fig.add_axes([0.52, 0.60, 0.46, 0.30]); axs.set_xlim(0, 20); axs.set_ylim(0, 6); axs.axis("off")
fig.text(0.52, 0.95, "Saturation", fontsize=7, fontweight="bold", color=INK, va="top")
eg(0.97, "ATP2B1")
gene_model(axs, 1.5, 2.6, 9.5, pas=[(6.4, "pPAS"), (13.6, "dPAS")])
# usage bars under the model: prenatal mixed, postnatal all distal
for yy, frac, col, lab in [(1.0, 0.55, PRE, "Prenatal"), (0.15, 1.0, POST, "Postnatal")]:
    axs.add_patch(Rectangle((4.5, yy), 9.5, 0.45, fc=GREY_LT, ec="none"))
    axs.add_patch(Rectangle((4.5, yy), 9.5 * frac, 0.45, fc=col, ec="none", alpha=0.85))
    axs.text(14.4, yy + 0.22, lab, va="center", fontsize=5.5, color=POST_TXT if lab == "Postnatal" else PRE_TXT)
axs.text(4.2, 0.62, "Distal\nusage", ha="right", va="center", fontsize=5.5, color=MUTED, linespacing=0.9)
a2 = mini_axes([0.59, 0.14, 0.29, 0.40])
a2.axhspan(1.0, 1.12, color=GREY_LT, lw=0)
a2.axhline(1.0, color=MUTED, lw=0.5)
a2.plot([0, 1, 2], [0.35, 0.55, 0.75], color=PRE, lw=1.0, marker="o", ms=2.6)
a2.plot([0, 1, 2], [0.985, 0.985, 0.985], color=POST, lw=1.0, marker="o", ms=2.6, clip_on=False)
a2.set_ylim(0, 1.12); a2.set_yticks([0, 1]); a2.set_yticklabels(["0", "1"]); a2.set_ylabel("PDUI", labelpad=1)
a2.text(2.36, 0.985, "Ceiling", fontsize=5.5, color=POST_TXT, va="center", ha="left")
a2.text(2.25, 0.72, "Effect", fontsize=5.5, color=PRE_TXT, va="center", ha="left")

fig.savefig(os.path.join(FIG_DIR, "fig4_c.png"), dpi=600, transparent=False, facecolor="white")
fig.savefig(os.path.join(FIG_DIR, "fig4_c.pdf"), facecolor="white")
