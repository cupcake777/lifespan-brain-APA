# Dumps 0/1/2 ALT-allele dosages for key ARL3-window SNPs from raw/chb.full.vcf.gz (GRCh38), for cohort LD.
# Writes data/07_coloc/cohort_dosage_keysnps.tsv.
import os, subprocess
root = f"{os.environ['APA_ROOT']}/"
keys = {102676884: "rs8354", 102697699: "rs11191369", 102689687: "rs11191368", 102712218: "rs4919669",
        102708826: "rs10883767", 102760520: "rs12571643"}
p = subprocess.Popen(["zcat", root + "raw/chb.full.vcf.gz"], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
out = open(root + "data/07_coloc/cohort_dosage_keysnps.tsv", "w"); seen = False
for line in p.stdout:
    if line.startswith("##"): continue
    if line.startswith("#CHROM"):
        h = line.rstrip("\n").split("\t"); out.write("snp\tid\t" + "\t".join(h[9:]) + "\n"); continue
    c = line[:line.index("\t")]
    if c != "10":
        if seen: break
        continue
    seen = True
    f = line.split("\t", 3); pos = int(f[1])
    if pos > 102_800_000: break
    if pos in keys:
        f = line.rstrip("\n").split("\t"); gi = f[8].split(":").index("GT"); d = []
        for g in f[9:]:
            gt = g.split(":")[gi].replace("|", "/")
            d.append("NA" if "." in gt else str(sum(a != "0" for a in gt.split("/"))))
        out.write(keys[pos] + "\t" + f[2] + "\t" + "\t".join(d) + "\n")
out.close(); p.kill()
