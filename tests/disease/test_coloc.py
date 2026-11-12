#!/usr/bin/env python3
# Runs make_coloc_pair.R (coloc.abf and coloc.susie on two synthetic loci), then the PP4 summary step.
# Checks: shared-causal locus PP4 high, distinct-causal locus PP3-dominant, summary reports the shared hit only.
import os
import subprocess
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "07_coloc", "followup")))
import coloc_pp4_summary as prod  # noqa: E402

# --- run coloc on synthetic loci ---

fails = []
rscript = os.environ.get("RSCRIPT", "Rscript")
with tempfile.TemporaryDirectory() as tmp:
    cmd = [rscript, os.path.join(HERE, "make_coloc_pair.R"), tmp]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-2000:]); print(res.stderr[-3000:])
        print("FAIL: make_coloc_pair.R did not run (set RSCRIPT to an Rscript with coloc and susieR)")
        sys.exit(1)
    ac = pd.read_csv(os.path.join(tmp, "synthetic_coloc_results.csv"))
    sus = pd.read_csv(os.path.join(tmp, "synthetic_susie_results.csv"))
    ac["gene"] = ac["Phenotype"].str.split("|").str[1]
    sh = ac[ac.gene == "SIMSHARED"].iloc[0]
    di = ac[ac.gene == "SIMDISTINCT"].iloc[0]
    print(f"[coloc.abf] SIMSHARED PP4={sh.PP4:.4g}; SIMDISTINCT PP3={di.PP3:.4g} PP4={di.PP4:.4g}")
    if not (sh.PP4 >= 0.8):
        fails.append("shared-causal locus PP4 not high (coloc.abf)")
    if not (di.PP3 >= 0.5 and di.PP4 <= 0.1):
        fails.append("distinct-causal locus not PP3-dominant (coloc.abf)")
    s_sh = float(sus.loc[sus.source_file.str.contains("SIMSHARED"), "PP.H4.abf"].iloc[0])
    s_di = float(sus.loc[sus.source_file.str.contains("SIMDISTINCT"), "PP.H4.abf"].iloc[0])
    print(f"[coloc.susie] SIMSHARED PP4={s_sh:.4g}; SIMDISTINCT PP4={s_di:.4g}")
    if not (s_sh >= 0.8):
        fails.append("shared-causal locus PP4 not high (coloc.susie)")
    if not (s_di <= 0.1):
        fails.append("distinct-causal locus PP4 not low (coloc.susie)")

    # --- PP4 summary step ---

    genes = ["SIMSHARED", "SIMDISTINCT"]
    vals = prod.summarize_coloc(ac.drop(columns=["gene"]), pp4_cut=0.7, genes=genes)
    sv = prod.summarize_susie(sus, genes=genes, pp4_cut=0.7)
    print(f"[summary] coloc_hits_raw={vals['coloc_hits_raw']} coloc_hit_genes={vals['coloc_hit_genes']} "
          f"SIMSHARED_hit_qtl={vals['coloc_SIMSHARED_hit_qtl']} SIMDISTINCT_hit_qtl={vals['coloc_SIMDISTINCT_hit_qtl']}; "
          f"susie comparisons={sv['susie_comparisons']} pp4>=0.7={sv['susie_comparisons_pp4ge0.7']}")
    if vals["coloc_hits_raw"] != 1 or vals["coloc_hit_genes"] != 1:
        fails.append("summary did not report exactly the shared-causal hit")
    if vals["coloc_SIMSHARED_hit_qtl"] != "SIM_QTL" or vals["coloc_SIMDISTINCT_hit_qtl"] != "none":
        fails.append("QTL class of hits wrong")
    if sv["susie_comparisons_pp4ge0.7"] != 1:
        fails.append("coloc.susie summary count wrong")

# --- result ---

if fails:
    print("FAIL:", "; ".join(fails))
    sys.exit(1)
print("PASS test_coloc")
