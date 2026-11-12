#!/usr/bin/env python3
# Fig 5C counts (venn and per-GWAS) from the module 07 gene context table, and the Fig 6A ATP2B1 HC3 stage effects.
# The HC3 refit itself is module 05 (05_gxs) and is not re-run here.
# Writes fig5c_venn_counts.tsv, fig5c_gwas_counts.tsv and fig6a_atp2b1_stage_effects.tsv to data/07_coloc/.
import os

import pandas as pd

root = os.environ["APA_ROOT"]
gl = pd.read_csv(os.path.join(root, "data/07_coloc/gene_gwas_classification_1mb.tsv"), sep="\t")
hc = pd.read_csv(os.path.join(root, "data/05_gxs/intx_hc3_core_groups.tsv"), sep="\t")

# --- fig5c venn ---
e = gl[gl["has_eqtl_data"].astype(str).str.lower() == "true"]
venn = pd.DataFrame({
    "set": ["APA-colocalized genes (with eQTL data)", "eQTL-colocalized genes (PP4>=0.7)",
            "APA-only (no eQTL coloc)"],
    "n": [e["gene"].nunique(),
          e.loc[e["all_eqtl_resources_max_pp4"] >= 0.7, "gene"].nunique(),
          e.loc[e["apa_only_gene_all_gwas"].astype(str).str.lower() == "true", "gene"].nunique()],
})

# --- fig5c per GWAS ---
g = gl.groupby("GWAS").agg(n_genes=("gene", "nunique")).reset_index()
apa = gl[gl["apa_only_gene_gwas"].astype(str).str.lower() == "true"].groupby("GWAS")["gene"].nunique()
eq = gl[gl["all_eqtl_resources_max_pp4"] >= 0.7].groupby("GWAS")["gene"].nunique()
g["n_apa_only"] = g["GWAS"].map(apa).fillna(0).astype(int)
g["n_eqtl_coloc_genes"] = g["GWAS"].map(eq).fillna(0).astype(int)
g = g.sort_values("GWAS")[["GWAS", "n_genes", "n_apa_only", "n_eqtl_coloc_genes"]]

# --- fig6a ATP2B1 stage effects ---
gene = "ATP2B1"
r = hc[hc["group_id"] == gene]
if len(r) != 1:
    raise ValueError(f"expected one HC3 group row for {gene}, found {len(r)}")
r = r.iloc[0]
f6 = pd.DataFrame([{
    "gene": gene, "lead_phenotype_id": r["lead_phenotype_id"], "lead_variant_id": r["lead_variant_id"],
    "beta_pre": r["beta_pre"], "hc3_se_pre": r["hc3_se_pre"], "beta_post": r["beta_post"],
    "hc3_se_post": r["hc3_se_post"], "beta_gxstage": r["beta_gxstage"],
    "hc3_pval_gxstage": r["hc3_pval_gxstage"]}])[["gene", "lead_phenotype_id", "lead_variant_id", "beta_pre", "hc3_se_pre",
    "beta_post", "hc3_se_post", "beta_gxstage", "hc3_pval_gxstage"]]

# --- print and write ---
print(venn.to_string(index=False)); print(g.to_string(index=False)); print(f6.to_string(index=False))
venn.to_csv(os.path.join(root, "data/07_coloc/fig5c_venn_counts.tsv"), sep="\t", index=False)
g.to_csv(os.path.join(root, "data/07_coloc/fig5c_gwas_counts.tsv"), sep="\t", index=False)
f6.to_csv(os.path.join(root, "data/07_coloc/fig6a_atp2b1_stage_effects.tsv"), sep="\t", index=False, float_format="%.17g")
