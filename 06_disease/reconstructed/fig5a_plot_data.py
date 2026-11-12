#!/usr/bin/env python3
# Fig 5A S-LDSC plot data: apaQTL and eQTL enrichment with 95% intervals and two-stage FDR (EUR).
# Writes data/06_disease/fig5a_sldsc_plot_data.tsv.
import os
import sys

import numpy as np
import pandas as pd

APA = os.environ["APA_ROOT"]
DATA = os.path.join(APA, "data/06_disease")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sldsc_two_stage import bh, load_processed, two_stage_fdr  # noqa: E402

Z975 = 1.96  # 95% interval multiplier used by the consumer table
GWAS_KEY = {"BD": "BIP"}  # output label -> processed-table trait code
ORDER = [("BD", "Postnatal"), ("BD", "Prenatal"), ("PD", "Postnatal"), ("PD", "Prenatal"),
         ("SCZ", "Postnatal"), ("SCZ", "Prenatal"), ("AD", "Prenatal"), ("AD", "Postnatal")]
EQTL_ORDER = [("SCZ", "Prenatal"), ("SCZ", "Postnatal"), ("BD", "Prenatal"), ("BD", "Postnatal"),
              ("PD", "Prenatal"), ("PD", "Postnatal"), ("AD", "Prenatal"), ("AD", "Postnatal")]

# --- apaQTL rows from the binary S-LDSC processed table (EUR, BD/PD/SCZ) ---
df = load_processed(os.path.join(DATA, "sldsc_enrich_processed.tsv"))
two = two_stage_fdr(df)
two = two[two["Population"] == "EUR"].copy()
rows = []
for tr, st in ORDER:
    if tr == "AD":
        continue
    r = two[(two.Trait == GWAS_KEY.get(tr, tr)) & (two.Stage == st)].iloc[0]
    se = float(df.loc[(df.Population == "EUR") & (df.Trait == GWAS_KEY.get(tr, tr)) &
                      (df.Annotation == ("prenatal_apa_qtl" if st == "Prenatal" else "postnatal_apa_qtl")),
                      "Enrichment_std_error"].iloc[0])
    rows.append(dict(Population="EUR", Trait=tr, Stage=st, enrich=r.Enrichment,
                     se=se, p=r.Enrichment_p, type="apaQTL", fdr=r.FDR_two_stage,
                     source="sldsc_two_stage_fdr (processed table)"))

# --- apaQTL rows for AD from the Bellenguez EUR_AD2 stage table (two-stage BH: prenatal, postnatal) ---
ad = pd.read_csv(os.path.join(APA, "data/04_qtl/Fig5a_EUR_AD2_eqtl_sldsc.tsv"), sep="\t")
ad["Stage"] = ad["Stage"].str.capitalize()
ad["fdr2"] = bh(ad["apaQTL_enrichment_p"].to_numpy())
ad = ad.set_index("Stage")
for tr, st in ORDER:
    if tr != "AD":
        continue
    r = ad.loc[st]
    rows.append(dict(Population="EUR", Trait="AD", Stage=st, enrich=r.apaQTL_enrichment,
                     se=r.apaQTL_enrichment_se, p=r.apaQTL_enrichment_p, type="apaQTL",
                     fdr=r.fdr2, source="Fig5a_EUR_AD2 handoff (Bellenguez EUR_AD2)"))
out = pd.DataFrame(rows)

# --- eQTL rows: 1015 eQTL comparison table (AD from the EUR_AD2 table) ---
cmp_ = pd.read_csv(os.path.join(DATA, "sldsc_apa_eqtl_comp.tsv"), sep="\t")
eq_rows = []
for tr, st in EQTL_ORDER:
    if tr == "AD":
        r = ad.loc[st]
        eq_rows.append(dict(Population="EUR", Trait="AD", Stage=st, enrich=r.eQTL_enrichment,
                            se=r.eQTL_enrichment_se, p=r.eQTL_enrichment_p, type="eQTL", fdr=np.nan,
                            source="Fig5a_EUR_AD2 handoff (Bellenguez EUR_AD2)"))
    else:
        r = cmp_[(cmp_.GWAS == f"EUR_{GWAS_KEY.get(tr, tr)}") & (cmp_.Stage == st.lower())].iloc[0]
        eq_rows.append(dict(Population="EUR", Trait=tr, Stage=st, enrich=r.eQTL_enrichment,
                            se=r.eQTL_enrichment_se, p=r.eQTL_enrichment_p, type="eQTL", fdr=np.nan,
                            source="sldsc_apa_eqtl_comparison (1015)"))
out = pd.concat([out, pd.DataFrame(eq_rows)], ignore_index=True)

# --- 95% intervals and output ---
out["lo"] = out["enrich"] - Z975 * out["se"]
out["hi"] = out["enrich"] + Z975 * out["se"]
out.to_csv(os.path.join(DATA, "fig5a_sldsc_plot_data.tsv"), sep="\t", index=False, float_format="%.6g")
