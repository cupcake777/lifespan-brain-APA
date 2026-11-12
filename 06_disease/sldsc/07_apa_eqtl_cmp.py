#!/usr/bin/env python3
# Pair stage-specific APA-QTL and eQTL sLDSC estimates (one row per GWAS x stage).
# Reads data/06_disease/regression/*.results under $APA_ROOT.
import csv
import math
import os
from pathlib import Path

APA = Path(os.environ["APA_ROOT"])
REGRESSION_DIR = APA / "data" / "06_disease" / "regression"
OUTPUT_FILE = APA / "data" / "06_disease" / "sldsc_apa_eqtl_comp.tsv"
TRAITS = {
    "EUR": ["MDD", "SCZ", "BIP", "ADHD", "PTSD", "AD", "ASD", "ANX", "INS", "INT", "PD", "EP"],
    "EAS": ["SCZ", "BIP", "AD", "PD", "AIS", "ALS", "EP"],
}
STAGES = ("prenatal", "postnatal")


# --- helpers ---
def read_custom_result(path):
    # the last row of a .results file is the custom annotation
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No result rows in {path}")
    row = rows[-1]
    values = {
        "enrichment": float(row["Enrichment"]),
        "enrichment_se": float(row["Enrichment_std_error"]),
        "enrichment_p": float(row["Enrichment_p"]),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError(f"Non-finite custom result in {path}")
    if not 0 <= values["enrichment_p"] <= 1:
        raise ValueError(f"Invalid enrichment P value in {path}")
    return values


# --- pair APA-QTL and eQTL estimates per GWAS and stage ---
output_rows = []
for population, traits in TRAITS.items():
    for trait in traits:
        for stage in STAGES:
            apa_annotation = f"{stage}_apa_qtl"
            eqtl_annotation = f"{stage}_eqtl"
            apa_path = REGRESSION_DIR / f"{population}_{trait}_{apa_annotation}.results"
            eqtl_path = REGRESSION_DIR / f"{population}_{trait}_{eqtl_annotation}.results"
            apa = read_custom_result(apa_path)
            eqtl = read_custom_result(eqtl_path)
            output_rows.append(
                {
                    "GWAS": f"{population}_{trait}",
                    "Population": population,
                    "Trait": trait,
                    "Stage": stage,
                    "apaQTL_enrichment": apa["enrichment"],
                    "apaQTL_enrichment_se": apa["enrichment_se"],
                    "apaQTL_enrichment_p": apa["enrichment_p"],
                    "eQTL_enrichment": eqtl["enrichment"],
                    "eQTL_enrichment_se": eqtl["enrichment_se"],
                    "eQTL_enrichment_p": eqtl["enrichment_p"],
                    "apaQTL_annotation": apa_annotation,
                    "eQTL_annotation": eqtl_annotation,
                }
            )

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT_FILE.open("w", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=output_rows[0], delimiter="\t")
    writer.writeheader()
    writer.writerows(output_rows)
print(f"Wrote {len(output_rows)} paired GWAS-stage rows to {OUTPUT_FILE}")
