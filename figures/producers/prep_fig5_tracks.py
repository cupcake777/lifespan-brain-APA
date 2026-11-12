# Fig 5 d-e regional tracks: GWAS and apaQTL -log10 P from the repaired coloc run, eQTL nominal P, r2 to the apaQTL lead.
# Writes summary tables only (no genotypes) to data/figures/fig5.
import os, json, gzip, subprocess
import numpy as np, pandas as pd

ROOT = os.environ["APA_ROOT"] + "/"
OUT = ROOT + "data/figures/fig5/"; os.makedirs(OUT, exist_ok=True)
LOCI = {"SCAF1": dict(json="EUR_SCZ_pre_chr19_rs7251_SCAF1_locuscompare.json", stage="prenatal", chrom="chr19"),
        "KLC1": dict(json="EUR_SCZ_post_chr14_1_rs_KLC1_locuscompare.json", stage="postnatal", chrom="chr14")}
eq = pd.read_csv(OUT + "fig5_KLC1_SCAF1_eqtl_tracks.tsv.gz", sep="\t")
eq[["chrom", "pos"]] = eq.variant_id.str.split(":", expand=True)[[0, 1]]; eq["pos"] = eq.pos.astype(int)
samples = {s.capitalize(): open(ROOT + f"data/04_qtl/apa.sample.{s.capitalize()}.final").read().split() for s in ("prenatal", "postnatal")}
for g, L in LOCI.items():
    j = json.load(open(OUT + L["json"])); p = pd.DataFrame(j["points"]); p = p[p.gwas_logp.notna() & p.qtl_logp.notna()]
    lead = p.loc[p.qtl_logp.idxmax()]
    L.update(meta=j["meta"], pts=p, lead_snp=lead.snp, lead_pos=int(lead.pos), x0=int(p.pos.min()), x1=int(p.pos.max()))
# one pass over the VCF: dosages for all variants in both windows
want = {L["chrom"].replace("chr", ""): (L["x0"] - 1000, L["x1"] + 1000) for L in LOCI.values()}
proc = subprocess.Popen(["zcat", ROOT + "raw/chb.full.vcf.gz"], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
hdr = None; G = {c: {} for c in want}
for line in proc.stdout:
    if line.startswith("##"): continue
    if line.startswith("#"): hdr = line.rstrip("\n").split("\t"); continue
    c, ps = line.split("\t", 2)[:2]
    if c not in want: continue
    ps = int(ps)
    if not (want[c][0] <= ps <= want[c][1]): continue
    f = line.rstrip("\n").split("\t"); fmt = f[8].split(":"); gi = fmt.index("GT")
    d = []
    for x in f[9:]:
        gt = x.split(":")[gi].replace("|", "/").split("/")
        d.append(np.nan if "." in gt else float(sum(int(a) > 0 for a in gt)))
    G[c][ps] = np.array(d)
proc.stdout.close(); proc.wait()
col = {s: i for i, s in enumerate(hdr[9:])}
for g, L in LOCI.items():
    st = L["stage"].capitalize(); idx = [col[s] for s in samples[st] if s in col]
    D = G[L["chrom"].replace("chr", "")]
    x = D[L["lead_pos"]][idx]
    def r2(pos):
        if pos not in D: return np.nan
        y = D[pos][idx]; m = ~np.isnan(x) & ~np.isnan(y)
        if m.sum() < 20 or np.std(y[m]) == 0: return np.nan
        return float(np.corrcoef(x[m], y[m])[0, 1] ** 2)
    p = L["pts"]
    gw = pd.DataFrame(dict(track="GWAS", pos=p.pos.astype(int), snp=p.snp, logp=p.gwas_logp))
    aq = pd.DataFrame(dict(track="apaQTL", pos=p.pos.astype(int), snp=p.snp, logp=p.qtl_logp))
    e = eq[(eq.gene == g) & (eq.stage == L["stage"])].copy()
    rs = dict(zip(p.pos.astype(int), p.snp))
    e["snp"] = [rs.get(q, v) for q, v in zip(e.pos, e.variant_id)]
    eqt = pd.DataFrame(dict(track="eQTL", pos=e.pos, snp=e.snp, logp=-np.log10(e.pval_nominal)))
    T = pd.concat([gw, aq, eqt]); T = T[(T.pos >= L["x0"]) & (T.pos <= L["x1"])]
    cache = {q: r2(q) for q in T.pos.unique()}
    T["r2"] = T.pos.map(cache); T["gene"] = g; T["lead_pos"] = L["lead_pos"]; T["lead_snp"] = L["lead_snp"]
    T.to_csv(OUT + f"fig5_{g}_tracks.tsv", sep="\t", index=False)
    rows = []
    with gzip.open(ROOT + "ref/gencode.v40.annotation.gtf.gz", "rt") as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.split("\t", 9)
            if f[0] != L["chrom"] or f[2] not in ("gene", "exon"): continue
            s_, e_ = int(f[3]), int(f[4])
            if e_ < L["x0"] or s_ > L["x1"]: continue
            a = f[8]; nm = a.split('gene_name "')[1].split('"')[0]; ty = a.split('gene_type "')[1].split('"')[0]
            if f[2] == "exon" and 'tag "Ensembl_canonical"' not in a: continue
            rows.append((f[0], f[2], s_, e_, f[6], nm, ty))
    pd.DataFrame(rows, columns=["chr", "feat", "start", "end", "strand", "gene", "type"]).drop_duplicates().to_csv(OUT + f"fig5_{g}_genes.tsv", sep="\t", index=False)
