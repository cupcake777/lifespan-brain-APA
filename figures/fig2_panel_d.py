# Fig 2 panel d: PIK3CA transcripts with polyA-site, APA-site, miR-7 and PUM motif marks (chr3 locus).
import os, re, gzip, csv, io
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicRecord, GraphicFeature
from matplotlib import font_manager

ROOT = os.environ["APA_ROOT"]
FIG_DIR = os.path.join(ROOT, "figures")
DATA = os.path.join(ROOT, "data/figures/fig2")
MINUS = "\u2212"
W0, W1 = 179198000, 179241000
C1, C2 = "#4A7BA7", "#D17A4A"
TEAL = "#7A9E9F"; INK = "#222222"
EXON = "#5A626A"; CDS = "#24282E"
for _f in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
    font_manager.fontManager.addfont(os.path.join(os.environ["APA_FONT_DIR"], _f))
FONT = "Arial"

# --- transcripts (GENCODE v40 PIK3CA exons and CDS) ---
tx = {}
with gzip.open(os.path.join(ROOT, "ref/gencode.v40.annotation.gtf.gz"), "rt") as f:
    for line in f:
        if line.startswith("#") or not line.startswith("chr3\t"):
            continue
        p = line.split("\t")
        if 'gene_name "PIK3CA"' not in p[8] or p[2] not in ("exon", "CDS"):
            continue
        tid = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
        nm = re.search(r'transcript_name "([^"]+)"', p[8])
        d = tx.setdefault(tid, {"name": nm.group(1) if nm else tid, "exons": [], "cds": []})
        d["exons" if p[2] == "exon" else "cds"].append((int(p[3]), int(p[4])))
for d in tx.values():
    d["exons"].sort(); d["cds"].sort()
order = sorted(tx, key=lambda t: int(re.sub(r".*-", "", tx[t]["name"])))

# --- polyASite2 atlas PAS in the window ---
atlas = []
with gzip.open(os.path.join(ROOT, "ref/polyasite2/atlas.clusters.2.0.GRCh38.96.bed.gz"), "rt") as f:
    for line in f:
        p = line.split("\t")
        if p[0] == "3" and p[9] == "TE" and W0 < int(p[1]) < W1:
            atlas.append(int(p[1]))

# --- APA sites and the representative site (S1 table, read once) ---
S1_PATH = os.path.join(DATA, "S1_pair_results_371.tsv")
with open(S1_PATH) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
apa_ends = sorted(int(re.sub(r".*:(\d+)-(\d+)\|.*", r"\2", r["site_id"])) for r in rows
                  if r["gene_symbol"] == "PIK3CA" and float(r["q"]) < 0.05)
REP_END = 179235371
REP_SITE = "ENST00000675467.1|PIK3CA|chr3:179234094-179235371|+"
_rep = [r for r in rows if r["site_id"] == REP_SITE]
assert len(_rep) == 1, len(_rep)
T_LABEL = f"t = {float(_rep[0]['t_apa_M1']):.1f}".replace("-", MINUS)

# --- miR-7 and PUM motif scan on hg38 chr3 (+ strand; index 0 = OFF) ---
fa = os.path.join(ROOT, "ref/hg38.fa")
fai = {}
with open(fa + ".fai") as f:
    for line in f:
        p = line.split("\t")
        fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
length, offset, lb, lw = fai["chr3"]
OFF = 179234365
p0, p1 = OFF - 1, W1 - 1
b0 = offset + (p0 // lb) * lw + (p0 % lb)
b1 = offset + ((p1 - 1) // lb) * lw + ((p1 - 1) % lb) + 1
with open(fa, "rb") as f:
    f.seek(b0)
    seq = f.read(b1 - b0).decode().replace("\n", "").upper()
assert len(seq) == W1 - OFF, len(seq)
mir7_hits = [m.start() for m in re.finditer(r"GTCTTCCA", seq)]
pum_hits = [m.start() for m in re.finditer(r"TGTA[ACT]ATA", seq)]
mir7_site = OFF + mir7_hits[0]   # first hit
pum_site = OFF + pum_hits[0]     # first hit

# --- render ---
n = len(order)
fig, axs = plt.subplots(n, 1, figsize=(3.80, 3.10), dpi=600, sharex=True,
                        gridspec_kw={"hspace": 0.24, "left": 0.12, "right": 0.98, "top": 0.90, "bottom": 0.17})
for i, tid in enumerate(order):
    ax = axs[i]
    d = tx[tid]
    feats = []
    for s, e in d["exons"]:
        if e < W0 or s > W1: continue
        feats.append(GraphicFeature(start=max(s, W0) - W0, end=min(e, W1) - W0, strand=0,
                                    color=EXON, thickness=9.5, linecolor=EXON, linewidth=0))
    for s, e in d["cds"]:
        if e < W0 or s > W1: continue
        feats.append(GraphicFeature(start=max(s, W0) - W0, end=min(e, W1) - W0, strand=0,
                                    color=CDS, thickness=9.5, linecolor=CDS, linewidth=0))
    rec = GraphicRecord(sequence_length=W1 - W0, features=feats)
    rec.plot(ax=ax, with_ruler=False, draw_line=True)
    if ax.lines:
        ax.lines[0].set_color("#888888")
        ax.lines[0].set_linewidth(0.5)
    ax.set_xlim(0, W1 - W0)
    ax.set_ylim(-1.15, 1.15)
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(bottom=False, left=False, labelbottom=False)
    t_name = d["name"].replace("PIK3CA-", "")
    ax.text(-650, 0, t_name, fontsize=6.8, family=FONT, color=INK, va="center", ha="right", weight="bold")


def vline_all(x, color, lw, ls="-", zorder=5, alpha=0.9):
    for ax in axs:
        ax.axvline(x - W0, color=color, lw=lw, ls=ls, zorder=zorder, alpha=alpha)


for e in atlas:
    vline_all(e, TEAL, 0.60, alpha=0.85)
for e in apa_ends:
    vline_all(e, C1, 0.85, ls=(0, (3, 2)), zorder=6)
vline_all(mir7_site, C1, 1.6, zorder=7)
vline_all(pum_site, C2, 1.6, zorder=7)

idx_208 = [i for i, t in enumerate(order) if "208" in tx[t]["name"]][0]
axs[idx_208].annotate(T_LABEL, xy=(REP_END - W0, 0.5), xytext=(REP_END - W0 + 400, 1.0),
                      arrowprops=dict(arrowstyle="->", color=INK, lw=0.6),
                      fontsize=5.6, family=FONT, color=INK, va="bottom", weight="bold")

ax_leg = axs[0]
leg_y = 2.45
ax_leg.plot([600, 2200], [leg_y, leg_y], color=C1, lw=0.9, ls=(0, (3, 2)), clip_on=False)
ax_leg.text(2600, leg_y, "APA site", fontsize=6.4, family=FONT, color=C1, va="center", clip_on=False, weight="bold")
ax_leg.plot([11500, 13100], [leg_y, leg_y], color=TEAL, lw=0.75, clip_on=False)
ax_leg.text(13500, leg_y, "Atlas PAS", fontsize=6.4, family=FONT, color=TEAL, va="center", clip_on=False, weight="bold")
ax_leg.plot([21500, 22600], [leg_y, leg_y], color=C1, lw=1.6, clip_on=False)
ax_leg.plot([22800, 23900], [leg_y, leg_y], color=C2, lw=1.6, clip_on=False)
ax_leg.text(24300, leg_y, "miR-7 / PUM motif", fontsize=6.4, family=FONT, color=INK, va="center", clip_on=False, weight="bold")

axs[-1].set_axis_on()
ticks_genomic = [179200000, 179210000, 179220000, 179230000, 179240000]
axs[-1].set_xlim(0, W1 - W0)
axs[-1].set_xticks([t - W0 for t in ticks_genomic])
axs[-1].set_xticklabels(["179.20", "179.21", "179.22", "179.23", "179.24"], fontsize=6.4, family=FONT, color=INK)
axs[-1].tick_params(axis="x", length=3.0, width=0.5, pad=2.0, bottom=True, labelbottom=True)
axs[-1].spines["bottom"].set_visible(True)
axs[-1].spines["bottom"].set_linewidth(0.5)
axs[-1].spines["bottom"].set_color(INK)
for sp in ["top", "left", "right"]:
    axs[-1].spines[sp].set_visible(False)
axs[-1].set_yticks([])
axs[-1].set_xlabel("chr3 (Mb, GRCh38)", fontsize=6.8, family=FONT, color=INK, labelpad=2.5)

# --- crop outer white margin ---
buf = io.BytesIO()
fig.savefig(buf, format="png", dpi=600, facecolor="white")
plt.close(fig)
buf.seek(0)
im = Image.open(buf)
arr = np.array(im.convert("L"))
non_white_rows = np.where((arr < 252).sum(axis=1) > 0)[0]
non_white_cols = np.where((arr < 252).sum(axis=0) > 0)[0]
r0 = max(0, non_white_rows.min() - 20); r1 = min(arr.shape[0], non_white_rows.max() + 20)
c0 = max(0, non_white_cols.min() - 25); c1 = min(arr.shape[1], non_white_cols.max() + 25)
im.crop((c0, r0, c1, r1)).save(os.path.join(FIG_DIR, "fig2_d.png"))
