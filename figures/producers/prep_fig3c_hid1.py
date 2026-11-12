# Fig 3c (superseded HID1 locus): sequence context and genotype-group PDUI statistics; individual values stay in memory.
import gzip, os
import numpy as np, pandas as pd, statsmodels.api as sm

ROOT = os.environ["APA_ROOT"]
D4 = f"{ROOT}/data/04_qtl"; POS = 74961620; CHROM_VCF = "17"; CHROM_FA = "chr17"
EVENT = "ENST00000528902.1|HID1|chr17:74961597-74961989|-"

# sequence context from hg38 (plus strand, then reverse-complement to the gene strand '-')
fai = {}
for line in open(f"{ROOT}/ref/hg38.fa.fai"):
    p = line.split("\t"); fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
W = 12
L, off, lb, lw = fai[CHROM_FA]
start, end = POS - W, POS + W
n = end - start + 1
s = off + (start - 1)//lb*lw + (start - 1) % lb
with open(f"{ROOT}/ref/hg38.fa", "rb") as f:
    f.seek(s); raw = f.read(n + (n//lb + 2)*(lw - lb) + lw)
plus = raw.replace(b"\n", b"")[:n].decode().upper()
comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
gene = "".join(comp[b] for b in reversed(plus))
ctx = pd.DataFrame({"pos": list(range(POS + W, POS - W - 1, -1)), "gene_strand_base": list(gene)})
ctx["offset"] = ctx.pos - POS
ctx.to_csv(f"{D4}/fig3c_hid1_context.tsv", sep="\t", index=False)

# genotype (ALT copies) at the variant, and PDUI of the event, in memory
vcf_hdr = None; target = None
with gzip.open(f"{ROOT}/raw/chb.full.vcf.gz", "rt") as f:
    for line in f:
        if line.startswith("#CHROM"): vcf_hdr = line.rstrip("\n").split("\t"); continue
        if line.startswith("#"): continue
        c = line[:line.index("\t")]
        if c == CHROM_VCF:
            p = int(line.split("\t", 2)[1])
            if p == POS: target = line.rstrip("\n").split("\t"); break
            if p > POS + 100000: break
assert target is not None, "variant not found"
ref, alt = target[3], target[4]
dos = {}
for i in range(9, len(vcf_hdr)):
    gt = target[i].split(":")[0].replace("|", "/")
    dos[vcf_hdr[i]] = sum(1 for x in gt.split("/") if x == "1")

with open(f"{D4}/apa_after_combat.tsv") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    row = None
    for line in f:
        if line.startswith(EVENT + "\t"):
            row = line.rstrip("\n").split("\t"); break
assert row is not None
pdui = pd.Series(np.array(row[1:], dtype=float), index=hdr[1:])

pre = [l.strip() for l in open(f"{D4}/apa.sample.Prenatal.final") if l.strip()]
post = [l.strip() for l in open(f"{D4}/apa.sample.Postnatal.final") if l.strip()]
stats = []; fits = []
ref_g = comp[ref]; alt_g = comp[alt]
for stage, samples in [("Prenatal", pre), ("Postnatal", post)]:
    y = pdui[samples].values
    x = np.array([dos[s] for s in samples], dtype=float)
    m = sm.OLS(y, sm.add_constant(x)).fit()
    fits.append((stage, len(samples), m.params[1], m.bse[1], m.pvalues[1]))
    for k in [0, 1, 2]:
        v = y[x == k]
        if len(v) == 0: continue
        lab = ref_g + ref_g if k == 0 else (ref_g + alt_g if k == 1 else alt_g + alt_g)
        stats.append(dict(stage=stage, genotype_gene_strand=lab, alt_copies=k, n=len(v),
                          min=v.min(), q1=np.percentile(v, 25), median=np.median(v),
                          q3=np.percentile(v, 75), max=v.max(), mean=v.mean()))
pd.DataFrame(stats).to_csv(f"{D4}/fig3c_hid1_group_stats.tsv", sep="\t", index=False)
pd.DataFrame(fits, columns=["stage", "n", "ols_slope_per_alt_copy", "se", "p"]).to_csv(f"{D4}/fig3c_hid1_ols_check.tsv", sep="\t", index=False)
