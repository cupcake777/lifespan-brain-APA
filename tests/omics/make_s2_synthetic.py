#!/usr/bin/env python3
# Synthetic inputs for s2_element.py: 1,500 DaPars2-style events (C1 / C2 / non-dynamic), genome, miRNAs, CLIP peaks.
# Planted: PUM motif TGTA[ACT]ATA at high density in C1 alternative segments; AGO2 brain CLIP peaks over them.
# Usage: make_s2_synthetic.py <outdir>
import os
import sys
import gzip
import random

import numpy as np
import pandas as pd

out = sys.argv[1]
os.makedirs(out, exist_ok=True)

# --- simulate events ---
rng = random.Random(20261019)
nprng = np.random.default_rng(20261019)
chroms = {"chr1": 400000, "chr2": 400000}
seqs = {c: list("".join(rng.choice("ACGT") for _ in range(L))) for c, L in chroms.items()}
n_ev = 1500
rows = []
for i in range(n_ev):
    chrom = rng.choice(list(chroms))
    strand = rng.choice(["+", "-"])
    S = rng.randint(2000, chroms[chrom] - 3000)
    Ulen = rng.randint(400, 1500)
    E = S + Ulen - 1
    u = rng.random()
    cl = 1 if u < 0.40 else (2 if u < 0.70 else 0)
    A = rng.randint(120, 700)
    # alternative segment: [P, E] on + strand, [S, P] on - strand (1-based, inclusive)
    P = E - A + 1 if strand == "+" else S + A - 1
    rows.append(dict(i=i, chrom=chrom, strand=strand, S=S, E=E, P=P, cl=cl))
ev = pd.DataFrame(rows)

# --- plant PUM motifs ---
for r in ev.itertuples():
    a0, a1 = (r.P, r.E) if r.strand == "+" else (r.S, r.P)
    if r.cl == 1:
        for _ in range(max(1, (a1 - a0) // 120)):
            pos = rng.randint(a0, a1 - 9)
            motif = "TGTA" + rng.choice("ACT") + "ATA"
            seqs[r.chrom][pos - 1: pos - 1 + 8] = list(motif)

# --- write genome fasta and index ---
fa = os.path.join(out, "hg38_synthetic.fa")
with open(fa, "wb") as f, open(fa + ".fai", "w") as fai:
    off = 0
    for c, L in chroms.items():
        hdr = (">%s\n" % c).encode()
        f.write(hdr); off += len(hdr)
        s = "".join(seqs[c]).encode()
        for k in range(0, L, 50):
            f.write(s[k:k + 50] + b"\n")
        fai.write("%s\t%d\t%d\t50\t51\n" % (c, L, off))
        off += ((L + 49) // 50) * 51

# --- DaPars2 table ---
# columns used: Gene, Predicted_Proximal_APA, Loci
gene = ["ENST%08d.1|ENSG%08d|%s|%s" % (r.i, r.i, r.chrom, r.strand) for r in ev.itertuples()]
loci = ["%s:%d-%d" % (r.chrom, r.S, r.E) for r in ev.itertuples()]
apa = pd.DataFrame({"Gene": gene, "fit_value": nprng.uniform(10, 1000, len(ev)).round(1),
                    "Predicted_Proximal_APA": ev.P.values, "Loci": loci})
apa.to_csv(os.path.join(out, "all.apa.count"), sep="\t", index=False)

# --- event table with class flags ---
tx = [g.split("|")[0] for g in gene]
sites = pd.DataFrame({
    "site_id": [t + "|" + l + "|" + s for t, l, s in zip(tx, loci, ev.strand.values)],
    "transcript_id": tx,
    "gene_id": [g.split("|")[1] for g in gene],
    "locus": loci,
    "strand": ev.strand.values,
    "cluster_id": [{0: np.nan, 1: 1.0, 2: 2.0}[c] for c in ev.cl.values],
})
sites["directional_pass"] = ev.cl.isin([1, 2]).values
sites["dynamic_any_pass"] = ev.cl.isin([1, 2]).values | (nprng.random(len(ev)) < 0.3)
sites["association_pass"] = True
sites.to_csv(os.path.join(out, "02_sites_clustered.tsv"), sep="\t", index=False)

# --- miRNAs and brain expression ---
# one 'brain' row in rpmm
mir = {"hsa-miR-1000-5p": "UGAGGUAGUAGGUUGUAUAGUU", "hsa-miR-2000-3p": "UAGCAGCACAUCAUGAUUUGA",
       "hsa-miR-3000-5p": "UUGGCAUUCUGUCCACCUCCA", "hsa-miR-4000-5p": "AAAGUGCUUACAGUGCAGGUAG"}
with open(os.path.join(out, "mature.fa"), "w") as f:
    for k, v in mir.items():
        f.write(">%s\n%s\n" % (k, v))
pq = {"Tissue": ["brain"], "Organ_system": ["brain"], "Metric": ["mean"]}
for k in mir:
    pq[k] = [2500.0 if k in ("hsa-miR-1000-5p", "hsa-miR-4000-5p") else 300.0]
pd.DataFrame(pq).to_parquet(os.path.join(out, "rpmm.parquet"))

# --- CLIP peaks ---
# POSTAR-style; columns 1-8 used: chrom, start, end, -, strand, RBP, method, tissue
peaks = []
for r in ev.itertuples():
    a0, a1 = (r.P, r.E) if r.strand == "+" else (r.S, r.P)
    if r.cl == 1:
        for _ in range(3):
            s = rng.randint(a0, max(a0, a1 - 30))
            peaks.append((r.chrom, s - 1, s + 20, "AGO2", "HITS-CLIP,CLIPper", "Brain"))
    if rng.random() < 0.3:
        s = rng.randint(r.S, r.E - 20); peaks.append((r.chrom, s - 1, s + 20, "ELAVL1", "PAR-CLIP,PARalyzer", "HEK293"))
    if rng.random() < 0.2:
        s = rng.randint(r.S, r.E - 20); peaks.append((r.chrom, s - 1, s + 20, "PUM2", "eCLIP", "K562"))
    if rng.random() < 0.2:
        s = rng.randint(r.S, r.E - 20); peaks.append((r.chrom, s - 1, s + 20, "FMR1", "PAR-CLIP,PARalyzer", "HEK293"))
    if rng.random() < 0.2:
        s = rng.randint(r.S, r.E - 20); peaks.append((r.chrom, s - 1, s + 20, "PUM1", "eCLIP", "K562"))
with gzip.open(os.path.join(out, "postar_synthetic.txt.gz"), "wt") as g:
    for k, p in enumerate(peaks):
        g.write("%s\t%d\t%d\tpeak%d\t+\t%s\t%s\t%s\tGSE0\t1.0\n" % (p[0], p[1], p[2], k, p[3], p[4], p[5]))

# --- truth table ---
ev.to_csv(os.path.join(out, "truth_events.tsv"), sep="\t", index=False)
print("synthetic inputs written to", out, "events", len(ev), "peaks", len(peaks))
