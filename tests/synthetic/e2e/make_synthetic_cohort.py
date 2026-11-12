#!/usr/bin/env python3
# Write a simulated mini cohort (synthetic only, no real IDs or values) in the input formats of the release pipeline.
# Usage: make_synthetic_cohort.py OUTDIR [--donors N] [--sites N] [--variants N] [--seed S]
# Writes data/raw/ and genotype/ under OUTDIR, plus truth.json with the planted effects.

import argparse
import json
import os

import numpy as np
import pandas as pd

# --- arguments ---
ap = argparse.ArgumentParser()
ap.add_argument("outdir")
ap.add_argument("--donors", type=int, default=120)
ap.add_argument("--sites", type=int, default=240)
ap.add_argument("--variants", type=int, default=30)
ap.add_argument("--seed", type=int, default=20261019)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
os.makedirs(os.path.join(a.outdir, "data", "raw"), exist_ok=True)
os.makedirs(os.path.join(a.outdir, "genotype"), exist_ok=True)

# --- donors and profiles ---
# about 1 in 10 donors has a technical replicate
n = a.donors
ids = [f"D{i:04d}" for i in range(n)]
stage_pre = rng.random(n) < 0.35
age = np.where(stage_pre, rng.uniform(0.02, 0.5, n), rng.uniform(1, 80, n))
sex = np.where(rng.random(n) < 0.5, "Female", "Male")
rin = np.round(rng.normal(7.5, 0.8, n).clip(4, 10), 2)
batch = rng.integers(1, 4, n)
rep = rng.random(n) < 0.10
rows = []
for i in range(n):
    rows.append(dict(sample=ids[i], individual=ids[i], age=round(age[i], 3), sex=sex[i], RIN=rin[i],
                     stage="prenatal" if stage_pre[i] else "postnatal", batch=int(batch[i])))
    if rep[i]:
        rows.append(dict(sample=ids[i] + "_rep", individual=ids[i], age=round(age[i], 3), sex=sex[i],
                         RIN=round(rin[i] + rng.normal(0, 0.1), 2), stage="prenatal" if stage_pre[i] else "postnatal",
                         batch=int(batch[i])))
meta = pd.DataFrame(rows)
meta_out = meta[["sample", "age", "sex", "RIN", "stage", "batch"]]
meta_out.to_csv(os.path.join(a.outdir, "data", "raw", "CHB.all.meta"), sep="\t", index=False)

# --- sites ---
S = a.sites
L = 30
genes = [f"GENE{g:03d}" for g in range(S // 4)]
site_gene = [genes[i % len(genes)] for i in range(S)]
CIS_GENE = genes[10]
GXS_GENE = genes[20]
site_ids = [f"TX{i:04d}|{site_gene[i]}|chr1:{100000 + 500 * i}-{100300 + 500 * i}|+" for i in range(S)]
pd.DataFrame({"site_id": site_ids, "gene_id": site_gene}).to_csv(
    os.path.join(a.outdir, "data", "raw", "phenotype.group"), sep="\t", header=False, index=False)

# --- genotypes: variants with planted cis and G x Stage effects ---
V = a.variants
maf = rng.uniform(0.1, 0.45, V)
dos = np.vstack([rng.binomial(2, p, n) for p in maf]).astype(int)
cis_var = 0
gxs_var = 1
var_gene = [CIS_GENE] * (V // 2) + ["GENEOTHER"] * (V - V // 2)
var_gene[gxs_var] = GXS_GENE
var_ids = [f"chr1:{200000 + 1000 * v}_A_G" for v in range(V)]
# donor-level dosage table, expanded to profile level later
geno_df = pd.DataFrame(dos, index=var_ids, columns=ids)
geno_df.index.name = "variant_id"
geno_df.reset_index().to_csv(os.path.join(a.outdir, "genotype", "geno.dosage.tsv"), sep="\t", index=False)
pd.DataFrame({"variant_id": var_ids, "chr": "chr1", "pos": [200000 + 1000 * v for v in range(V)],
              "ref": "A", "alt": "G", "gene_id": var_gene}).to_csv(
    os.path.join(a.outdir, "genotype", "geno.variants.tsv"), sep="\t", index=False)

# --- PDUI-like matrix over profiles ---
lg = np.log10(meta.age.values + 1.0)
post = (meta.stage.values == "postnatal").astype(float)
donor_pos = {d: k for k, d in enumerate(ids)}
didx = np.array([donor_pos[d] for d in meta.individual.values])
base = rng.uniform(-0.8, 0.8, S)
beta_age = np.zeros(S)
beta_age[:L] = 0.9
beta_age[L:2 * L] = -0.9
beta_cis = np.zeros(S)
beta_gxs = np.zeros(S)
cis_sites = [i for i in range(S) if site_gene[i] == CIS_GENE]
gxs_sites = [i for i in range(S) if site_gene[i] == GXS_GENE]
beta_cis[cis_sites] = 0.8
beta_gxs[gxs_sites] = 1.0
# logit-scale PDUI
eta = (base[:, None] + beta_age[:, None] * lg[None, :]
       + beta_cis[:, None] * dos[cis_var, didx][None, :]
       + beta_gxs[:, None] * (dos[gxs_var, didx] * post)[None, :]
       + rng.normal(0, 0.35, (S, len(meta))))
pdui = 1 / (1 + np.exp(-eta))
pdui = np.clip(pdui, 0.02, 0.98)
prof = pd.DataFrame(pdui, columns=meta["sample"].values)
prof.insert(0, "site_id", site_ids)
prof.to_csv(os.path.join(a.outdir, "data", "raw", "apa_after_combat.tsv"), sep="\t", index=False)

# --- truth ---
# planted: sites 0..29 lengthen and 30..59 shorten with log10(age+1); CIS_GENE sites follow variant 0;
# GXS_GENE sites follow variant 1 in postnatal samples only; other sites are null
truth = dict(
    lengthening_sites=site_ids[:L], shortening_sites=site_ids[L:2 * L],
    cis_gene=CIS_GENE, cis_variant=var_ids[cis_var], cis_sites=[site_ids[i] for i in cis_sites],
    gxs_gene=GXS_GENE, gxs_variant=var_ids[gxs_var], gxs_sites=[site_ids[i] for i in gxs_sites],
    gxs_effect="postnatal only", n_donors=n, n_profiles=int(len(meta)), n_sites=S, seed=a.seed,
)
with open(os.path.join(a.outdir, "truth.json"), "w") as fh:
    json.dump(truth, fh, indent=1)
print(f"wrote {len(meta)} profiles, {n} donors, {S} sites, {V} variants to {a.outdir}")