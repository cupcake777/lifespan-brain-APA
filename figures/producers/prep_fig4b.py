# Fig 4b input: corrected PDUI by lead-variant genotype and age for one example gene per G x stage class.
# Individual-level output goes to data/withheld (never upload); the summary sign check goes to data/figures/fig4.
import os, subprocess
import numpy as np, pandas as pd

ROOT = os.environ["APA_ROOT"] + "/"
D4 = ROOT + "data/04_qtl"; OUTW = ROOT + "data/withheld"; OUTF = ROOT + "data/figures/fig4"
os.makedirs(OUTW, exist_ok=True); os.makedirs(OUTF, exist_ok=True)
EX = {"Reversed": "PRDX6", "Amplified": "MIR124-1HG", "Attenuated": "ZEB2"}
gx = pd.read_csv(D4 + "/gxstage_strict_classes_local_checks.tsv", sep="\t")
gx = gx[gx.gene_id.isin(EX.values())].set_index("gene_id")
want = {v.split(":")[0].replace("chr", "") + ":" + v.split(":")[1]: g for g, v in gx.lead_variant_id.items()}
chroms = sorted({int(k.split(":")[0]) for k in want}); last = max(chroms)
p = subprocess.Popen(["zcat", ROOT + "raw/chb.full.vcf.gz"], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
G = {}; hdr = None
for line in p.stdout:
    if line.startswith("##"): continue
    if line.startswith("#CHROM"): hdr = line.rstrip("\n").split("\t"); continue
    c = line[:line.index("\t")]
    if not c.isdigit(): continue
    if int(c) > last: break
    if int(c) not in chroms: continue
    f = line.split("\t", 5); k = c + ":" + f[1]
    if k in want:
        f = line.rstrip("\n").split("\t"); gi = f[8].split(":").index("GT")
        d = []
        for x in f[9:]:
            gt = x.split(":")[gi].replace("|", "/")
            d.append(np.nan if "." in gt else sum(a != "0" for a in gt.split("/")))
        G[want[k]] = dict(zip(hdr[9:], d)); G[want[k] + "_alleles"] = (f[3], f[4])
p.kill()
apa = pd.read_csv(D4 + "/apa_after_combat.tsv", sep="\t", index_col=0)
meta = pd.read_csv(ROOT + "raw/all.912.meta", sep="\t", dtype={"sample": str}).set_index("sample")
st = {s: "Prenatal" for s in open(D4 + "/apa.sample.Prenatal.final").read().split()}
st.update({s: "Postnatal" for s in open(D4 + "/apa.sample.Postnatal.final").read().split()})
rows, chk = [], []
for cls, g in EX.items():
    ph = gx.loc[g, "lead_phenotype_id"]; ref, alt = G[g + "_alleles"]
    for s, v in apa.loc[ph].items():
        if s not in st: continue
        gs = s.replace("_rep", "")
        if gs not in G[g] or pd.isna(G[g][gs]): continue
        rows.append(dict(cls=cls, gene=g, sample=s, stage=st[s], age=meta.loc[s, "age"] if s in meta.index else np.nan,
                         dosage=int(G[g][gs]), ref=ref, alt=alt, pdui=v))
D = pd.DataFrame(rows)
for (g, s), d in D.groupby(["gene", "stage"]):
    b = np.polyfit(d.dosage, d.pdui, 1)[0]
    hc3 = gx.loc[g, "beta_pre" if s == "Prenatal" else "beta_post"]
    chk.append(dict(gene=g, stage=s, n=len(d), ols_slope=b, hc3_beta=hc3, same_sign=np.sign(b) == np.sign(hc3)))
D.to_csv(OUTW + "/fig4b_example_pdui.tsv", sep="\t", index=False)
pd.DataFrame(chk).to_csv(OUTF + "/fig4b_example_signcheck.tsv", sep="\t", index=False)
