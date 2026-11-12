#!/usr/bin/env python
# Extract nominal eQTL and apaQTL pairs for the eQTL-tested genes of each life stage.
# Writes nominal_eqtl_<stage>.tsv.gz and nominal_apa_<stage>.tsv.gz to data/04_qtl; existing files are reused.
import gzip
import os
import re
import shlex
import subprocess
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

APA_ROOT = Path(os.environ["APA_ROOT"])
D = APA_ROOT / "data" / "04_qtl"
GTF = APA_ROOT / "ref" / "gencode.v40.annotation.gtf.gz"
APA_NOMINAL = Path(os.environ["APA_HPC_HOME"]) / "project" / "work" / "apa" / "01_QTL" / "nominal"
EQTL_NOMINAL = D / "eqtl" / "nominal"
COLS = ["phenotype_id", "variant_id", "af", "ma_count", "pval_nominal", "slope", "slope_se"]


# --- helpers ---
def extract_tsv(src, keys, out, strip_version):
    keys = sorted(keys)
    if not keys:
        raise RuntimeError(f"empty key set for {src}")
    key_file = Path(tempfile.gettempdir()) / f"coloc_keys_{out.name}.txt"
    key_file.write_text("\n".join(keys) + "\n")
    awk_key = 'sub(/\\..*/,"",k);' if strip_version else ""
    awk = f'NR==FNR{{keep[$0]=1;next}} FNR==1{{print;next}} {{k=$1; {awk_key} if(k in keep) print}}'
    cmd = (f"zcat {shlex.quote(str(src))} | awk -F'\\t' -v OFS='\\t' {shlex.quote(awk)} "
           f"{shlex.quote(str(key_file))} - | gzip -c > {shlex.quote(str(out))}")
    try:
        subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    finally:
        key_file.unlink(missing_ok=True)


# --- ENSG symbol map from GENCODE genes ---
sym_to_ensg = {}
with gzip.open(GTF, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        gene = re.search(r'gene_id "([^"]+)"', fields[8])
        name = re.search(r'gene_name "([^"]+)"', fields[8])
        if gene and name and name.group(1) not in sym_to_ensg:
            sym_to_ensg[name.group(1)] = gene.group(1).split(".", 1)[0]

# --- per stage: genes of significant APA loci, then nominal pairs ---
for stage in ("prenatal", "postnatal"):
    with gzip.open(D / f"{stage}.sig_QTL.txt.gz", "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = header.index("phenotype_id")
        apa_keys = {line.rstrip("\n").split("\t")[idx] for line in handle if line.strip()}
    genes = {sym_to_ensg.get(key.split("|", 2)[1], key.split("|", 2)[1])
             for key in apa_keys if "|" in key}
    eqtl_cache = D / f"nominal_eqtl_{stage}.tsv.gz"
    apa_cache = D / f"nominal_apa_{stage}.tsv.gz"
    if not eqtl_cache.exists():
        extract_tsv(EQTL_NOMINAL / f"{stage}.cis_qtl_pairs.txt.gz", genes, eqtl_cache, True)
    if not apa_cache.exists():
        if stage == "prenatal":
            # prenatal nominal pairs were written as parquet chunks by tensorQTL
            keep = pa.array(sorted(apa_keys), type=pa.string())
            files = sorted(APA_NOMINAL.glob(f"{stage}.cis_qtl_pairs.*.parquet"))
            if not files:
                raise RuntimeError(f"no parquet files found in {APA_NOMINAL}")
            with gzip.open(apa_cache, "wt") as handle:
                handle.write("\t".join(COLS) + "\n")
                for path in files:
                    table = pq.read_table(path, columns=COLS)
                    table = table.filter(pc.is_in(table["phenotype_id"], value_set=keep))
                    if table.num_rows == 0:
                        continue
                    rows = table.to_pydict()
                    for i in range(table.num_rows):
                        handle.write("\t".join(str(rows[col][i]) for col in COLS) + "\n")
        else:
            extract_tsv(APA_NOMINAL / f"{stage}.cis_qtl_pairs.txt.gz", apa_keys, apa_cache, False)
    print(f"{stage}: APA phenotypes={len(apa_keys)} genes={len(genes)}")
