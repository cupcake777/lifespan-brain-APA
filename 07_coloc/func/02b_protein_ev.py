#!/usr/bin/env python3
# MS protein abundance vs candidate genotype dosage, one test per candidate gene.
# Writes data/07_coloc/protein_evidence.tsv.
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# --- helpers ---


def classify(n_tested, slope, pval, detected):
    if not detected:
        return "MS_not_detected"
    if n_tested < 5:
        return "MS_underpowered"
    if pd.notna(pval) and pval < 0.05 and pd.notna(slope) and slope > 0:
        return "MS_higher_in_alt"
    if pd.notna(pval) and pval < 0.05 and pd.notna(slope) and slope < 0:
        return "MS_lower_in_alt"
    return "MS_no_diff"


# --- config ---

APA_ROOT = Path(os.environ["APA_ROOT"])
MS_TABLE = APA_ROOT / "ref/MS_identified_information.txt"
GENO_MATRIX = APA_ROOT / "raw/geno.matrix"
CANDIDATES = APA_ROOT / "data/07_coloc/candidates_prioritized.tsv"
OUTPUT = APA_ROOT / "data/07_coloc/protein_evidence.tsv"

# --- inputs ---

candidates = pd.read_csv(CANDIDATES, sep="\t")

# one row per gene: most unique peptides, then highest coverage
ms = pd.read_csv(MS_TABLE, sep="\t")
ms = ms.loc[
    ms["Gene name"].notna()
    & ms["Gene name"].astype(str).str.strip().ne("")
    & ms["Gene name"].astype(str).str.strip().ne("--")
].copy()
for col in ["Unique peptides", "Coverage [%]"]:
    ms[col] = pd.to_numeric(ms[col], errors="coerce").fillna(0)
ms_best = (
    ms.sort_values(["Gene name", "Unique peptides", "Coverage [%]"], ascending=[True, False, False])
    .drop_duplicates("Gene name", keep="first")
    .copy()
)

# long format: one row per gene and sample
intensity_cols = [col for col in ms_best.columns if col.startswith("Intensity ")]
ms_long = ms_best.melt(
    id_vars=["Gene name", "Protein accession", "Unique peptides"],
    value_vars=intensity_cols,
    var_name="intensity_col",
    value_name="intensity",
)
ms_long["ms_sample_id"] = ms_long["intensity_col"].str.removeprefix("Intensity ")
ms_long["qtl_sample_id"] = ms_long["ms_sample_id"].map(lambda value: str(value).replace("_", "-").strip())
ms_long["intensity"] = pd.to_numeric(ms_long["intensity"], errors="coerce")
ms_long["log2_intensity"] = np.log2(ms_long["intensity"].where(ms_long["intensity"] > 0))

# genotype dosage rows for the candidate variants only
variant_ids = set(candidates["variant_id"].dropna().astype(str))
hits = []
for chunk in pd.read_csv(GENO_MATRIX, sep="\t", chunksize=1000):
    hit = chunk.loc[chunk["variant_id"].isin(variant_ids)].copy()
    if not hit.empty:
        hits.append(hit)
if not hits:
    dosage = pd.DataFrame(columns=["variant_id"])
else:
    dosage = pd.concat(hits, ignore_index=True).drop_duplicates("variant_id", keep="first")
dosage_lookup = dosage.set_index("variant_id") if not dosage.empty else pd.DataFrame()
protein_genes = set(ms_best["Gene name"].astype(str))

# --- tests ---

rows = []
for row in candidates.itertuples(index=False):
    gene = str(row.gene_name)
    candidate_id = str(row.candidate_id)
    variant_id = str(row.variant_id)
    detected = gene in protein_genes
    if not detected or variant_id not in dosage_lookup.index:
        rows.append(
            {
                "candidate_id": candidate_id,
                "gene_name": gene,
                "n_tested": 0,
                "protein_evidence": classify(0, np.nan, np.nan, detected),
                "ms_log2_fold_per_allele": np.nan,
                "ms_pvalue": np.nan,
                "ms_spearman_rho": np.nan,
                "ms_protein_accession": "",
                "ms_unique_peptides": np.nan,
            }
        )
        continue
    gene_ms = ms_long.loc[ms_long["Gene name"].astype(str).eq(gene)].copy()
    dosage_row = pd.to_numeric(dosage_lookup.loc[variant_id], errors="coerce")
    shared = sorted(set(gene_ms["qtl_sample_id"]) & set(dosage_row.index))
    gene_ms = gene_ms.set_index("qtl_sample_id").reindex(shared)
    dosage_vec = dosage_row.reindex(shared)

    # per-gene linear fit of log2 intensity on dosage, Spearman rho
    data = pd.DataFrame({"dosage": dosage_vec, "intensity": gene_ms["log2_intensity"]}).dropna()
    if data.shape[0] < 5 or data["dosage"].nunique() < 2:
        slope, pval, rho = np.nan, np.nan, np.nan
    else:
        rho, pval = spearmanr(data["dosage"], data["intensity"])
        slope = np.polyfit(data["dosage"].astype(float), data["intensity"].astype(float), deg=1)[0]
        slope, pval, rho = float(slope), float(pval), float(rho)

    n_tested = int(pd.DataFrame({"d": dosage_vec, "i": gene_ms["log2_intensity"]}).dropna().shape[0])
    best = ms_best.loc[ms_best["Gene name"].astype(str).eq(gene)].iloc[0]
    rows.append(
        {
            "candidate_id": candidate_id,
            "gene_name": gene,
            "n_tested": n_tested,
            "protein_evidence": classify(n_tested, slope, pval, detected),
            "ms_log2_fold_per_allele": slope,
            "ms_pvalue": pval,
            "ms_spearman_rho": rho,
            "ms_protein_accession": best["Protein accession"],
            "ms_unique_peptides": best["Unique peptides"],
        }
    )

# --- write ---

out = pd.DataFrame(rows)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUTPUT, sep="\t", index=False)
