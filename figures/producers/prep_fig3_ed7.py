# Fig 3 and ED7 inputs: apaQTL gene sets, lead-variant distance, PAS candidates, HID1 betas, coloc PP, ED7 pairs.
# Summary tables only; genotype and per-sample values stay in memory.
import gzip, os
import numpy as np, pandas as pd

ROOT = os.environ["APA_ROOT"]
D4 = f"{ROOT}/data/04_qtl"; OUT3 = f"{ROOT}/data/figures/fig3"; OUT7 = f"{ROOT}/data/figures/ed7"
os.makedirs(OUT3, exist_ok=True); os.makedirs(OUT7, exist_ok=True)

def zread(path, usecols=None):
    return pd.read_csv(gzip.open(path, "rt"), sep="\t", usecols=usecols, low_memory=False)

# Fig 3a / ED7e: apaQTL gene sets and lead-variant distance
cis = {}
for st in ["prenatal", "postnatal"]:
    d = zread(f"{D4}/{st}.cis_qtl.txt.gz", ["phenotype_id", "variant_id", "qval"])
    d = d[d.qval < 0.05].copy()
    d["gene_id"] = d.phenotype_id.str.split("|").str[1]
    p3 = d.phenotype_id.str.split("|").str[2]
    s = p3.str.extract(r"^chr[^:]+:(\d+)-(\d+)$").astype(float)
    pos = d.variant_id.str.split(":").str[1].astype(float)
    d["dist"] = np.minimum(np.abs(pos - s[0]), np.abs(pos - s[1]))
    cis[st] = d
pre_gene = set(cis["prenatal"].gene_id); post_gene = set(cis["postnatal"].gene_id)
venn = pd.DataFrame({"set": ["prenatal", "postnatal", "shared", "prenatal_only", "postnatal_only"],
                     "n": [len(pre_gene), len(post_gene), len(pre_gene & post_gene), len(pre_gene - post_gene), len(post_gene - pre_gene)]})
venn.to_csv(f"{OUT3}/fig3a_venn_counts.tsv", sep="\t", index=False)
dist = pd.concat([cis[s].assign(stage=s.capitalize())[["stage", "phenotype_id", "variant_id", "dist"]] for s in cis])
dist.to_csv(f"{OUT7}/ed7e_lead_distance_by_stage.tsv", sep="\t", index=False)

# Fig 3b: PAS-disrupting apaQTL candidates
pas = zread(f"{D4}/pas.tsv.gz")
pas["stage"] = np.where(pas.candidate_id.str.startswith("prenatal"), "Prenatal", "Postnatal")
pas["phen"] = pas.candidate_id.str.extract(r"^[a-z]+::(.*)::chr[^:]+:\d+:[ACGT]+:[ACGT]+$")[0]
pas["variant_id"] = pas.candidate_id.str.extract(r"::(chr[^:]+:\d+:[ACGT]+:[ACGT]+)$")[0]
h = pas[pas.predicted_pdui_direction.isin(["lengthening", "shortening"]) & (pas.pas_core_change != "not_eligible")]
h = h[["stage", "phen", "variant_id", "pas_label", "pas_core_change", "predicted_pdui_direction"]].drop_duplicates()
sq = []
for st, f in [("Prenatal", "prenatal"), ("Postnatal", "postnatal")]:
    x = zread(f"{D4}/{f}.sig_QTL.txt.gz", ["phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal"])
    x["stage"] = st; sq.append(x)
sq = pd.concat(sq)
m = h.merge(sq, left_on=["stage", "phen", "variant_id"], right_on=["stage", "phenotype_id", "variant_id"], how="left")
m["gene"] = m.phen.str.split("|").str[1]
m["concordant"] = ((m.predicted_pdui_direction == "lengthening") & (m.slope > 0)) | ((m.predicted_pdui_direction == "shortening") & (m.slope < 0))
m["lo"] = m.slope - 1.96 * m.slope_se; m["hi"] = m.slope + 1.96 * m.slope_se
m[["stage", "gene", "variant_id", "pas_label", "pas_core_change", "predicted_pdui_direction", "slope", "slope_se", "lo", "hi", "pval_nominal", "concordant"]].to_csv(f"{OUT3}/fig3b_pas_candidates.tsv", sep="\t", index=False)

# Fig 3c: HID1 variant betas (pipeline sig_QTL)
hid = []
for st, f in [("Prenatal", "prenatal"), ("Postnatal", "postnatal")]:
    x = zread(f"{D4}/{f}.sig_QTL.txt.gz", ["phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal"])
    x = x[x.phenotype_id.str.contains("HID1") & (x.variant_id == "chr17:74961620:T:C")].assign(stage=st)
    hid.append(x)
pd.concat(hid).to_csv(f"{OUT3}/fig3c_hid1_sigqtl.tsv", sep="\t", index=False)

# Fig 3d: coloc PP per gene and stage; mediation calls
pl = []
for st, f in [("Prenatal", "prenatal"), ("Postnatal", "postnatal")]:
    c = pd.read_csv(f"{D4}/coloc_{'pre' if f == 'prenatal' else 'post'}_apa_eqtl.tsv", sep="\t", low_memory=False)
    c = c[["gene_name", "ensg0", "apa_pheno", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4"]].assign(stage=st)
    pl.append(c)
pd.concat(pl).to_csv(f"{OUT3}/fig3d_coloc_pp_by_gene.tsv", sep="\t", index=False)

# ED7a: variant-gene pairs (lead variant per gene) from the two stage sig_QTL files
def lead_vg(path):
    s = zread(path, ["phenotype_id", "variant_id", "slope", "pval_nominal"])
    s["gene"] = s.phenotype_id.str.split("|").str[1]
    s = s.sort_values("pval_nominal").drop_duplicates(["gene", "variant_id"])
    return s[["gene", "variant_id", "slope", "pval_nominal"]]
a = lead_vg(f"{D4}/prenatal.sig_QTL.txt.gz")
b = lead_vg(f"{D4}/postnatal.sig_QTL.txt.gz")
mg = a.merge(b, on=["gene", "variant_id"], suffixes=("_pre", "_post"))
mg.to_csv(f"{OUT7}/ed7a_variant_gene_pairs.tsv", sep="\t", index=False)
