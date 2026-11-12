#!/usr/bin/env python3
# Methods gate: EUR GWAS with total h2 z (h2 / h2_se, binary baseline run) below 7 are not analysed.
# Bellenguez EUR_AD2 replaces Kunkle EUR_AD. Writes data/06_disease/eur_h2_z_gate.tsv.
import os

import pandas as pd

APA = os.environ["APA_ROOT"]
DATA = os.path.join(APA, "data/06_disease")
Z_CUT = 7.0

# --- EUR h2 z from the binary baseline run ---
base_path = os.path.join(DATA, "sldsc_baseline_h2_tau_sum.tsv")
b = pd.read_csv(base_path, sep="\t")
b["h2"] = pd.to_numeric(b["h2"], errors="coerce")
b["h2_se"] = pd.to_numeric(b["h2_se"], errors="coerce")
b["z"] = b["h2"] / b["h2_se"]
b["source"] = os.path.basename(base_path)
b = b[["population", "gwas", "h2", "h2_se", "z", "source"]]
eur = b[b["population"] == "EUR"].copy()
eur["below_7"] = eur["z"] < Z_CUT
kunkle = eur["gwas"] == "AD"  # replaced by EUR_AD2
eur["excluded_by_gate"] = eur["below_7"] & ~kunkle
eur["analysed_after_substitution"] = ~eur["excluded_by_gate"] & ~kunkle

# --- EUR_AD2 substitute: postnatal row, APOE stratum (apoe_excluded != T) ---
a = pd.read_csv(os.path.join(DATA, "sldsc_ad2_summary.tsv"), sep="\t")
a = a[a["annotation"] == "postnatal_apa_qtl"]
ad2_row = a[a["apoe_excluded"] != "T"].iloc[0]
ad2_z = float(ad2_row["h2"]) / float(ad2_row["h2_se"])
print(f"[AD2 substitute] EUR_AD2 (apoe not excluded) z = {ad2_z:.3f}")

listed = sorted(eur[eur["below_7"] & ~kunkle]["gwas"].tolist())
print(f"[computed list] {listed}")

# --- alternative gate: joint apaQTL+eQTL maxCPP run ---
joint = pd.read_csv(os.path.join(APA, "data/04_qtl/03_sldsc_tau_star.tsv"), sep="\t")
joint = joint.drop_duplicates(["population", "trait"])
joint = joint[joint["population"] == "EUR"].copy()
joint["gwas"] = joint["trait"].str.replace("EUR_", "", regex=False)
joint["z"] = joint["h2"] / joint["h2_se"]
print(f"[alternative: joint maxCPP] EUR z < 7: {sorted(joint.loc[joint['z'] < Z_CUT, 'gwas'].tolist())}")

# --- gate table ---
out = eur[["population", "gwas", "h2", "h2_se", "z", "below_7", "excluded_by_gate",
           "analysed_after_substitution", "source"]]
out.to_csv(os.path.join(DATA, "eur_h2_z_gate.tsv"), sep="\t", index=False, float_format="%.6g")
