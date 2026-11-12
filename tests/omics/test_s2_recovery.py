#!/usr/bin/env python3
# End-to-end test of s2_element.py on simulated inputs.
# Checks that the planted PUM enrichment in C1 alternative segments and the planted AGO2 CLIP signal are
# recovered, and that C2 shows no PUM enrichment. Exit code 0 on success.
import os
import subprocess
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SCRIPT = os.path.join(REPO, "03_omics", "s2_element.py")

# --- simulate inputs and run s2_element.py ---
tmp = tempfile.mkdtemp(prefix="s2syn_")
syn = os.path.join(tmp, "in")
out = os.path.join(tmp, "out")
subprocess.run([sys.executable, os.path.join(HERE, "make_s2_synthetic.py"), syn], check=True)
env = dict(os.environ)
env.update({
    "APA_ROOT": tmp,
    "APA_S2_DAPARS": os.path.join(syn, "all.apa.count"),
    "APA_S2_SITES": os.path.join(syn, "02_sites_clustered.tsv"),
    "APA_S2_HG38": os.path.join(syn, "hg38_synthetic.fa"),
    "APA_S2_MATURE": os.path.join(syn, "mature.fa"),
    "APA_S2_RPMM": os.path.join(syn, "rpmm.parquet"),
    "APA_S2_POSTAR": os.path.join(syn, "postar_synthetic.txt.gz"),
    "APA_S2_OUT": out,
})
subprocess.run([sys.executable, SCRIPT], env=env, check=True)

# --- read tables ---
t1 = pd.read_csv(os.path.join(out, "S2_aUTR_vs_cUTR_elem.tsv"), sep="\t")
t2 = pd.read_csv(os.path.join(out, "S2_C1C2_vs_nondynamic_elem.tsv"), sep="\t")
t3 = pd.read_csv(os.path.join(out, "S2_C1C2_vs_nondynamic_clip_adj.tsv"), sep="\t")
pum_a = float(t1.loc[t1.motif == "PUM_UGUAHAUA", "rr_AUadj"].iloc[0])
pum_c1 = float(t2.loc[(t2.motif == "PUM_UGUAHAUA") & (t2.group == "C1"), "rate_ratio.Estimate"].iloc[0])
pum_c2 = float(t2.loc[(t2.motif == "PUM_UGUAHAUA") & (t2.group == "C2"), "rate_ratio.Estimate"].iloc[0])
ago_c1 = float(t3.loc[(t3.rbp == "AGO2_brain") & (t3.group == "C1"), "rate_ratio"].iloc[0])
ago_c1_lo = float(t3.loc[(t3.rbp == "AGO2_brain") & (t3.group == "C1"), "lo"].iloc[0])
print("PUM aUTR vs cUTR (AU-adjusted): %.3f" % pum_a)
print("PUM C1 vs non-dynamic: %.3f ; C2: %.3f" % (pum_c1, pum_c2))
print("AGO2 brain C1 vs non-dynamic: %.3f (CI lower %.3f)" % (ago_c1, ago_c1_lo))

# --- checks ---
checks = [
    ("planted PUM enriched in aUTR vs cUTR", pum_a > 1.3),
    ("planted PUM enriched in C1 vs non-dynamic", pum_c1 > 1.3),
    ("planted contrast: C2 PUM rate ratio well below C1 (C2 < C1 / 2)", pum_c2 < pum_c1 / 2),
    ("planted AGO2 signal in C1 (CI excludes 1)", ago_c1_lo > 1.0),
    ("all three tables written", all(os.path.exists(os.path.join(out, f)) for f in
                                    ["S2_aUTR_vs_cUTR_elem.tsv", "S2_C1C2_vs_nondynamic_elem.tsv",
                                     "S2_C1C2_vs_nondynamic_clip_adj.tsv"])),
]
ok = True
for name, passed in checks:
    print(("PASS " if passed else "FAIL ") + name)
    ok = ok and passed
print("S2 synthetic recovery:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
