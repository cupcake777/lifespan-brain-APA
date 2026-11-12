# ALT allele frequency in the ARL3 window (chr10:102,200,000-103,300,000, GRCh38) for each QTL sample set.
# Streams raw/chb.full.vcf.gz (no index); writes data/07_coloc/cohort_af_arl3_window.tsv.
import os, subprocess, sys
root = f"{os.environ['APA_ROOT']}/"
sets = {k: [l.strip() for l in open(root + "data/04_qtl/" + f) if l.strip()]
        for k, f in [("apa_prenatal", "apa.sample.Prenatal.final"), ("apa_postnatal", "apa.sample.Postnatal.final"),
                     ("exp_prenatal", "exp.sample.Prenatal.final"), ("exp_postnatal", "exp.sample.Postnatal.final")]}
LO, HI = 102_200_000, 103_300_000
p = subprocess.Popen(["zcat", root + "raw/chb.full.vcf.gz"], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
out = open(root + "data/07_coloc/cohort_af_arl3_window.tsv", "w")
idx = {}; seen10 = False
for line in p.stdout:
    if line.startswith("##"): continue
    if line.startswith("#CHROM"):
        h = line.rstrip("\n").split("\t"); samp = {s: i for i, s in enumerate(h)}
        for k, v in sets.items():
            idx[k] = [samp[s] for s in v if s in samp]
            print(k, "listed", len(v), "in_vcf", len(idx[k]), file=sys.stderr, flush=True)
        out.write("chrom\tpos\tid\tref\talt\t" + "\t".join(f"af_{k}\tn_{k}" for k in sets) + "\n")
        continue
    c = line[:line.index("\t")]
    if c not in ("10", "chr10"):
        if seen10: break
        continue
    seen10 = True
    f = line.rstrip("\n").split("\t")
    pos = int(f[1])
    if pos < LO: continue
    if pos > HI: break
    gi = f[8].split(":").index("GT")
    vals = []
    for k in sets:
        a = n = 0
        for i in idx[k]:
            gt = f[i].split(":")[gi].replace("|", "/")
            if "." in gt: continue
            x, y = gt.split("/")
            a += (x != "0") + (y != "0"); n += 2
        vals += [f"{a/n:.6f}" if n else "NA", str(n // 2)]
    out.write("\t".join([f[0], f[1], f[2], f[3], f[4]] + vals) + "\n")
out.close(); p.kill()
print("done", file=sys.stderr)
