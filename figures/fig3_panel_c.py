# Fig 3 panel c: CEP41 distal hexamer variant chr7:130393789 C>T; allele sequences and per-genotype coverage tracks.
import glob, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Rectangle
from dna_features_viewer import GraphicRecord, GraphicFeature

ROOT = os.environ["APA_ROOT"]
SEQ_TSV = os.path.join(ROOT, "data/figures/fig3/fig3_cep41_seq.tsv")
DEPTH_DIR = os.path.join(ROOT, "data/withheld/fig3_cep41_depth")
FIG_DIR = os.path.join(ROOT, "figures")
for f in ["arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"]:
    fm.fontManager.addfont(os.path.join(os.environ["APA_FONT_DIR"], f))
plt.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.linewidth": 0.35, "pdf.fonttype": 42})
INK, MUTED, EXON, CDS = "#222222", "#5A5A5A", "#8A9099", "#3E434A"
STAGE = {"Prenatal": "#E3A73F", "Postnatal": "#6BAED6"}
GT = ["C/C", "C/T", "T/T"]

W0, W1 = 130393500, 130399300
X = lambda p: W1 - np.asarray(p)          # minus strand: 5' on the left
UTR_S, CDS_S, CDS_E = 130393771, 130398894, 130399039
PPAS, DPAS, VAR = 130396874, 130393771, 130393789

seq = pd.read_csv(SEQ_TSV, sep="\t").gene_strand_seq[0].replace("T", "U")   # index 12 = variant
dep_files = {(os.path.basename(f).split("_")[1], int(os.path.basename(f).split("_")[2])): f
             for f in glob.glob(os.path.join(DEPTH_DIR, "CEP41_*_primary_*.depth.tsv"))}
assert len(dep_files) == 6, len(dep_files)
pos = np.arange(W0, W1 + 1)
cutr = (pos >= PPAS) & (pos <= CDS_S - 1)
autr = (pos >= DPAS) & (pos < PPAS)


def track(st, k):
    d = pd.read_csv(dep_files[(st, k)], sep="\t", header=None, names=["c", "p", "d"])
    v = np.zeros(len(pos)); idx = d.p.values - W0; ok = (idx >= 0) & (idx < len(pos)); v[idx[ok]] = d.d.values[ok]
    v = pd.Series(v).rolling(41, center=True, min_periods=1).mean().values
    return v / v[cutr].mean()


MM = 1 / 25.4
fig = plt.figure(figsize=(108 * MM, 60 * MM), dpi=600)
L, Wd = 0.155, 0.735
# --- allele sequences ---
ax_s = fig.add_axes([0.50, 0.875, 0.39, 0.115]); ax_s.set_xlim(0, 1); ax_s.set_ylim(0, 1); ax_s.axis("off")
lo, hi = 6, 19; step = 1 / (hi - lo + 1)
for r, (lab, base) in enumerate([("C allele", "G"), ("T allele", "A")]):
    yy = 0.72 - r * 0.44
    ax_s.text(-0.02, yy, lab, ha="right", va="center", fontsize=5.8, color=INK)
    for i, k in enumerate(range(lo, hi)):
        ch = base if k == 12 else seq[k]; inhex = 11 <= k <= 16
        ax_s.text((i + 0.5) * step, yy, ch, ha="center", va="center", family="Arial", fontsize=6.4,
                  color=INK if inhex else "#9A9A9A", weight="bold" if inhex else "normal")
    ax_s.add_patch(Rectangle(((12 - lo) * step, yy - 0.2), step, 0.4, fill=False, lw=0.6, ec=INK))
ax_s.add_patch(Rectangle(((11 - lo) * step, 0.06), 6 * step, 0.88, fc="#EBEBEB", ec="none", zorder=-1))
# --- gene model ---
ax_g = fig.add_axes([L, 0.70, Wd, 0.14]); ax_g.set_xlim(0, W1 - W0); ax_g.set_ylim(-1, 1); ax_g.axis("off")
feats = [GraphicFeature(start=int(X(CDS_S - 1)), end=int(X(UTR_S)), strand=0, color=EXON, linewidth=0, thickness=5),
         GraphicFeature(start=int(X(CDS_E)), end=int(X(CDS_S)), strand=0, color=CDS, linewidth=0, thickness=9)]
sub = ax_g.inset_axes([0, 0.15, 1, 0.7]); GraphicRecord(sequence_length=W1 - W0, features=feats).plot(ax=sub, with_ruler=False, draw_line=False)
sub.set_xlim(0, W1 - W0); sub.axis("off")
ax_g.plot([0, X(CDS_E)], [0, 0], color="#888888", lw=0.5)
for p_, lab in [(PPAS, "pPAS"), (DPAS, "dPAS")]:
    ax_g.plot([X(p_)], [0.62], marker="v", ms=2.6, color=INK, mec="none"); ax_g.text(X(p_), 0.78, lab, ha="center", va="bottom", fontsize=5.8)
ax_g.text(-60, 0, "CEP41", ha="right", va="center", fontsize=6, style="italic")
# --- coverage tracks ---
th, gap, ggap = 0.085, 0.012, 0.03; y = 0.655; boxes = []
for st in ["Prenatal", "Postnatal"]:
    for k in range(3):
        b = y - th; boxes.append((st, k, b)); y -= th + gap
    y -= ggap - gap
ymax = max(np.nanmax(track(st, k)) for st, k, _ in boxes) * 1.04   # shared scale, no clipping
for st, k, b in boxes:
    a = fig.add_axes([L, b, Wd, th]); v = track(st, k)
    a.fill_between(X(pos), 0, v, color=STAGE[st], lw=0)
    a.axvspan(X(PPAS), X(DPAS), color="#F2F2F2", zorder=-2, lw=0)
    a.axvline(X(VAR), color=INK, lw=0.45, ls=(0, (2, 1.5)))
    a.set_xlim(0, W1 - W0); a.set_ylim(0, ymax)
    for sp in ["top", "right", "bottom"]: a.spines[sp].set_visible(False)
    a.set_xticks([]); a.set_yticks([0, 1]); a.tick_params(axis="y", length=1.5, width=0.35, pad=1, labelsize=5.2)
    a.text(-0.075, 0.42, GT[k], transform=a.transAxes, ha="right", va="center", fontsize=6)
    a.text(1.012, 0.42, f"{v[autr].mean() / v[cutr].mean():.2f}", transform=a.transAxes, ha="left", va="center", fontsize=5.8)
for st in ["Prenatal", "Postnatal"]:
    bs = [b for s, _, b in boxes if s == st]; y0, y1 = min(bs), max(bs) + th
    fig.lines.append(plt.Line2D([0.036, 0.036], [y0 + 0.004, y1 - 0.004], transform=fig.transFigure, color=STAGE[st], lw=1.6))
    fig.text(0.018, (y0 + y1) / 2, st, rotation=90, ha="center", va="center", fontsize=6)
top = max(b for _, _, b in boxes) + th
fig.text(L, top + 0.006, "Relative coverage", ha="left", va="bottom", fontsize=5.8, color=MUTED)
fig.text(L + Wd + 0.008, top + 0.006, "aUTR/\ncUTR", ha="left", va="bottom", fontsize=5.4, color=MUTED)
yb = min(b for _, _, b in boxes) - 0.035; xf = lambda xx: L + Wd * xx / (W1 - W0)
fig.lines.append(plt.Line2D([xf(0), xf(1000)], [yb, yb], transform=fig.transFigure, color=INK, lw=0.6))
fig.text(xf(500), yb - 0.008, "1 kb", ha="center", va="top", fontsize=5.8)
fig.savefig(os.path.join(FIG_DIR, "fig3_c.png"), dpi=600, facecolor="white")
fig.savefig(os.path.join(FIG_DIR, "fig3_c.pdf"), facecolor="white")
