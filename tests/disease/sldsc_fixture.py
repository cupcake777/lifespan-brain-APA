#!/usr/bin/env python3
# Writes a synthetic processed S-LDSC enrichment table (no real GWAS or S-LDSC values).
# Planted: EURSIM1 postnatal apaQTL p = 1e-4, EURSIM1 prenatal p = 0.02; EURSIM2 null-like.
# FDR_p is one BH over every row, as 07_collect_results.py writes it (6 significant digits).
import numpy as np
import pandas as pd

# --- fixture ---

ANN = ["prenatal_apa_qtl", "postnatal_apa_qtl", "lifespan_apa_qtl", "interaction_apa_qtl"]


def make(seed=2026):
    rng = np.random.default_rng(seed)
    rows = []
    fixed = {("EURSIM1", "prenatal_apa_qtl"): 0.02, ("EURSIM1", "postnatal_apa_qtl"): 1e-4,
             ("EURSIM2", "prenatal_apa_qtl"): 0.3, ("EURSIM2", "postnatal_apa_qtl"): 0.6}
    for pop, trait in [("EUR", "EURSIM1"), ("EUR", "EURSIM2"), ("EUR", "EURSIM3"), ("EAS", "EASSIM1")]:
        for ann in ANN:
            p = fixed.get((trait, ann), float(rng.uniform(0.05, 0.95)))
            enr = float(rng.normal(1.0, 0.2)) if (trait, ann) not in fixed else (2.0 if p < 0.05 else 1.0)
            if (trait, ann) == ("EURSIM1", "postnatal_apa_qtl"):
                enr = 2.0
            se = 0.3
            rows.append({"Trait": trait.replace("EURSIM", "SIM").replace("EASSIM", "SIM"), "Population": pop,
                         "Annotation": ann, "Coefficient": 1e-7, "Coefficient_std_error": 1e-7,
                         "Coefficient_z-score": 1.0, "Enrichment": enr, "Enrichment_std_error": se,
                         "Enrichment_p": p})
    df = pd.DataFrame(rows)

    p = np.asarray(df["Enrichment_p"].to_numpy(), float)
    m = len(p)
    o = np.argsort(p, kind="mergesort")
    r = p[o] * m / np.arange(1, m + 1)
    c = np.minimum.accumulate(r[::-1])[::-1]
    q = np.empty(m)
    q[o] = np.minimum(c, 1.0)

    df["FDR_p"] = [f"{x:.6g}" for x in q]
    df["NegLog10P"] = -np.log10(df["Enrichment_p"])
    df["NegLog10FDR"] = -np.log10(q)
    df["Significance"] = np.where(q < 0.05, "sig", "ns")
    df["PopTrait"] = df["Population"] + " " + df["Trait"]
    return df
