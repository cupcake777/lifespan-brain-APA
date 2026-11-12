#!/usr/bin/env python3
# Figure 5c disease-linked gene matrix: APA, CHB eQTL, GTEx brain and TWAS evidence per gene.
# Writes the source table, figure png/pdf and layout check.

import csv
import gzip
import os
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# --- config ---
APA_ROOT = Path(os.environ["APA_ROOT"])
DATA = APA_ROOT / "data/07_coloc"
TWAS = APA_ROOT / "data/06_disease/twas_ranked_all.tsv.gz"
FIG = APA_ROOT / "figures/07_coloc"
PREFERRED = [
    ("SCAF1", "EUR_SCZ"), ("NRBP2", "EUR_SCZ"), ("MTRFR", "EAS_SCZ"),
    ("ATP2B1", "EAS_SCZ"), ("OGFOD2", "EAS_SCZ"), ("ARL3", "EAS_SCZ"),
]


# --- helpers ---
def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# --- loading data ---
genes = read_rows(DATA / "gene_gwas_classification_1mb.tsv")
signals = read_rows(DATA / "signal_lookup_1mb.tsv")
chb_rows = read_rows(DATA / "chb_stage_eqtl_coloc_context.tsv")
gtex_rows = read_rows(DATA / "gtex_brain_eqtl_coloc_tissues.tsv")
# expected 87 gene-GWAS rows over 72 genes
assert len(genes) == 87 and len({row["gene"] for row in genes}) == 72

signals_by_pair = defaultdict(list)
for signal in signals:
    signals_by_pair[(signal["gene"], signal["GWAS"])].append(signal)
chb_by_signal = defaultdict(set)
for row in chb_rows:
    if row["status"] == "ok":
        chb_by_signal[row["signal_id"]].add(row["chb_stage"])
gtex_by_signal = defaultdict(set)
for row in gtex_rows:
    if row["status"] == "ok":
        gtex_by_signal[row["signal_id"]].add(row["tissue"])

gene_names = {row["gene"] for row in genes}
twas = defaultdict(lambda: {"tested": set(), "significant": set(), "min_p": 1.0})
with gzip.open(TWAS, "rt", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["gene_name"] not in gene_names:
            continue
        item = twas[(row["gene_name"], row["gwas_id"])]
        item["tested"].add(row["stage"])
        if row["bonferroni_significant"] == "True":
            item["significant"].add(row["stage"])
        if row["pvalue"]:
            item["min_p"] = min(item["min_p"], float(row["pvalue"]))


# --- gene by GWAS evidence ---

rows = []
for gene in genes:
    key = gene["gene"], gene["GWAS"]
    signal_group = signals_by_pair[key]
    twas_row = twas.get(key)
    chb_coverage = set().union(*(chb_by_signal[s["signal_id"]] for s in signal_group))
    gtex_coverage = set().union(*(gtex_by_signal[s["signal_id"]] for s in signal_group))
    rows.append({
        "gene": gene["gene"], "GWAS": gene["GWAS"],
        "apa_stages": ";".join(stage for stage in ("prenatal", "postnatal", "interaction") if stage in {s["stage"] for s in signal_group}),
        "apa_max_pp4": max(float(s["apa_pp4"]) for s in signal_group),
        "n_apa_signals": len(signal_group),
        "chb_max_pp4": gene["chb_max_pp4"],
        "chb_evaluable_stages": len(chb_coverage),
        "gtex_brain_max_pp4": gene["gtex_brain_max_pp4"],
        "gtex_evaluable_tissues": len(gtex_coverage),
        "twas_tested_stages": ";".join(sorted(twas_row["tested"])) if twas_row else "",
        "twas_bonferroni_stages": ";".join(sorted(twas_row["significant"])) if twas_row else "",
        "twas_min_p": twas_row["min_p"] if twas_row else "",
        "merged_1mb_locus_id": gene["merged_1mb_locus_id"],
        "candidate_status": "APA-only" if gene["apa_only_gene_gwas"] == "True" else "APA+eQTL",
        "all_gwas_gene_status": "APA-only" if gene["apa_only_gene_all_gwas"] == "True" else "APA+eQTL",
    })

by_key = {(row["gene"], row["GWAS"]): row for row in rows}
selected = [by_key[key] for key in PREFERRED]
for status in ("APA-only", "APA+eQTL"):
    choices = sorted(
        (row for row in rows if row["candidate_status"] == status and row["gene"] not in {item["gene"] for item in selected}),
        key=lambda row: (-row["apa_max_pp4"], row["gene"], row["GWAS"]),
    )
    selected.extend(choices[:2])
assert len(selected) == 10 and len({row["gene"] for row in selected}) == 10
for index, row in enumerate(selected, 1):
    row["panel_rank"] = index

fields = list(rows[0]) + ["panel_rank"]
with (DATA / "Figure5c_disease_gene_matrix.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


# --- figure 5c ---

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
fig, ax = plt.subplots(figsize=(14, 5.2))
fig.patch.set_facecolor("white")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

header = [(.025, "Gene"), (.145, "GWAS"), (.25, "APA stage"), (.355, "1 Mb locus"),
          (.49, "APA-GWAS PP4"), (.61, "CHB eQTL PP4"), (.735, "GTEx brain PP4"),
          (.855, "TWAS"), (.93, "Status")]
for x, label in header:
    ax.text(x, .93, label, fontsize=10.5, fontweight="bold", va="center")
ax.plot([.02, .985], [.885, .885], color="#91999B", linewidth=.8)

for index, row in enumerate(selected):
    y = .84 - index * .075
    if index % 2 == 1:
        ax.add_patch(Rectangle((.02, y - .035), .965, .07, facecolor="#F4F6F6", edgecolor="none"))
    ax.add_patch(Rectangle((.02, y - .034), .004, .068,
                           facecolor="#C47735" if row["candidate_status"] == "APA-only" else "#238078", edgecolor="none"))
    ax.text(.03, y, row["gene"], va="center", fontsize=12, fontweight="bold")
    ax.text(.145, y, row["GWAS"].replace("_", " "), va="center", fontsize=11)
    stage_text = "/".join({"prenatal": "Pr", "postnatal": "Po", "interaction": "Int"}[stage]
                          for stage in row["apa_stages"].split(";"))
    ax.text(.265, y, stage_text, va="center", fontsize=10.5)
    locus = row["merged_1mb_locus_id"].split(";")[0]
    locus = locus.split("_chr", 1)[-1].replace("_L", ":L")
    ax.text(.355, y, f"chr{locus}", va="center", fontsize=10.5)

    for x, field, coverage, denominator in (
        (.49, "apa_max_pp4", "n_apa_signals", None),
        (.61, "chb_max_pp4", "chb_evaluable_stages", 2),
        (.735, "gtex_brain_max_pp4", "gtex_evaluable_tissues", 13),
    ):
        value = row[field]
        if value == "":
            face, label, ink = "#E7E9E8", "NA", "#586366"
        else:
            numeric = float(value)
            face = plt.cm.cividis(numeric)
            label = f"{numeric:.2f}"
            ink = "white" if numeric < .52 else "#202625"
        ax.add_patch(Rectangle((x, y - .028), .083, .056, facecolor=face, edgecolor="white", linewidth=.5))
        ax.text(x + .0415, y + .006, label, ha="center", va="center", fontsize=10.5, color=ink, fontweight="bold")
        if denominator is not None:
            ax.text(x + .0415, y - .018, f"{row[coverage]}/{denominator}", ha="center", va="center", fontsize=7.5, color=ink)

    twas_sig = row["twas_bonferroni_stages"]
    if twas_sig:
        twas_text = "/".join("Pr" if stage == "prenatal" else "Po" for stage in twas_sig.split(";"))
        twas_color = "#146D66"
    elif row["twas_tested_stages"]:
        twas_text, twas_color = "none", "#596466"
    else:
        twas_text, twas_color = "NA", "#8A9496"
    ax.text(.855, y + .007, twas_text, va="center", fontsize=10, color=twas_color, fontweight="bold" if twas_sig else "normal")
    ax.text(.855, y - .02, f"{len(row['twas_tested_stages'].split(';')) if row['twas_tested_stages'] else 0}/2", va="center", fontsize=7.5, color="#596466")
    ax.text(.93, y, row["candidate_status"], va="center", fontsize=9.5,
            color="#9A531E" if row["candidate_status"] == "APA-only" else "#146D66", fontweight="bold")

ax.text(.025, .055, "Pr prenatal   Po postnatal   Int interaction   TWAS: Bonferroni hits by tested stage; none = tested without a hit", fontsize=9, color="#455052")
ax.text(.025, .018, "Status uses PP4 >= 0.7 in matched-stage CHB and GTEx brain; NA is untested. Cross-stage CHB is supplemental.", fontsize=9, color="#455052")

fig.canvas.draw()
canvas = fig.canvas.get_renderer()
frame = fig.bbox
outside = sum(not frame.contains(*text.get_window_extent(canvas).get_points()[0]) or
              not frame.contains(*text.get_window_extent(canvas).get_points()[1])
              for text in ax.texts)
assert outside == 0
FIG.mkdir(parents=True, exist_ok=True)
fig.savefig(FIG / "Figure5c_disease_gene_matrix.png", dpi=220, facecolor="white")
fig.savefig(FIG / "Figure5c_disease_gene_matrix.pdf", facecolor="white")
plt.close(fig)

with (DATA / "Figure5c_layout_validation.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["panel_rows", "source_rows", "text_outside_canvas", "canvas_inches", "status"])
    writer.writerow([len(selected), len(rows), outside, "14x5.2", "PASS"])
