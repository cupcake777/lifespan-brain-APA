# Per-chromosome z-scores of the G x Stage nominal pairs (b_gi / b_gi_se) and phenotype TSS positions.
# Reads data/05_gxs (tensorQTL interaction pairs, GeneRange BED); writes data/04_qtl/interaction_sldsc/.
import gzip
import os
import subprocess
from pathlib import Path

APA_ROOT = Path(os.environ["APA_ROOT"])
PAIRS = APA_ROOT / "data" / "05_gxs" / "interaction.cis_qtl_pairs.txt.gz"
BED = APA_ROOT / "data" / "05_gxs" / "phenotype.GeneRange.bed.gz"
OUT_DIR = APA_ROOT / "data" / "04_qtl" / "interaction_sldsc"
Z_DIR = OUT_DIR / "z"
EXPECTED_HEADER = ["phenotype_id", "variant_id", "start_distance", "end_distance", "af", "ma_samples", "ma_count",
                   "pval_g", "b_g", "b_g_se", "pval_i", "b_i", "b_i_se", "pval_gi", "b_gi", "b_gi_se"]
Z_DIR.mkdir(parents=True, exist_ok=True)

# --- stream pairs, write z = b_gi / b_gi_se per chromosome ---
proc = subprocess.Popen(["zcat", str(PAIRS)], stdout=subprocess.PIPE, bufsize=1024 * 1024)
src = proc.stdout
header = src.readline().decode().rstrip("\n").split("\t")
assert header == EXPECTED_HEADER, f"unexpected pairs header: {header}"

handles, counts = {}, {}
total, bad_z = 0, 0
for line in src:
    total += 1
    f = line.rstrip(b"\n").split(b"\t")
    vid = f[1]
    chrom = vid[3:vid.index(b":")].decode()
    try:
        z = float(f[14]) / float(f[15])
    except (ValueError, ZeroDivisionError):
        bad_z += 1
        continue
    if z != z or z in (float("inf"), float("-inf")):
        bad_z += 1
        continue
    h = handles.get(chrom)
    if h is None:
        h = gzip.open(Z_DIR / f"chr{chrom}.tsv.gz", "wt", compresslevel=5)
        h.write("phenotype_id\tvariant_id\taf\tz\n")
        handles[chrom] = h
        counts[chrom] = 0
    h.write(f"{f[0].decode()}\t{vid.decode()}\t{f[4].decode()}\t{z:.5g}\n")
    counts[chrom] += 1
for h in handles.values():
    h.close()
rc = proc.wait()
assert rc == 0, f"zcat exited {rc}"

with open(OUT_DIR / "extract_counts.tsv", "w") as out:
    out.write("chr\trows\n")
    for chrom in sorted(counts, key=lambda c: (len(c), c)):
        out.write(f"{chrom}\t{counts[chrom]}\n")
    out.write(f"TOTAL\t{sum(counts.values())}\nSOURCE_LINES\t{total}\nBAD_Z\t{bad_z}\n")

# --- phenotype TSS from the BED pid field (transcript|gene|chr:start-end|strand) ---
with gzip.open(BED, "rt") as inp, open(OUT_DIR / "phenotype_pos.tsv", "w") as out:
    out.write("phenotype_id\tchr\ttss\n")
    header = inp.readline()
    assert header.startswith("#Chr"), f"unexpected bed header: {header[:40]}"
    for line in inp:
        pid = line.rstrip("\n").split("\t")[3]
        parts = pid.split("|")
        chrom = parts[2].split(":")[0]
        start, end = parts[2].split(":")[1].split("-")
        tss = start if parts[3] == "+" else end
        out.write(f"{pid}\t{chrom}\t{tss}\n")
