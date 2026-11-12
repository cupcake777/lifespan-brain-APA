# Synthetic test: reconstruct_S4_gating_eqtl_fixed_RNA.py recovers the planted classes
# and planted Delta-expression on simulated inputs.
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RELEASE = HERE.parents[1]
sys.path.insert(0, str(HERE))
import make_qtl as gen  # noqa: E402

PROD = RELEASE / "05_gxs" / "reconstructed" / "recon_s4_gating.py"

# --- simulate and run producer ---
tmp = Path(tempfile.mkdtemp(prefix="qtl_gate_"))
inp = gen.write_inputs(tmp / "in")
out = tmp / "out"
cmd = [sys.executable, str(PROD),
       "--hc3", str(inp / "hc3_lead_table.tsv"),
       "--expr-nonfixed", str(inp / "expr_nonfixed.tsv"),
       "--expr-fixed", str(inp / "expr_fixed.tsv"),
       "--meta", str(inp / "meta.tsv"),
       "--gtf", str(inp / "synthetic.gtf.gz"),
       "--shared-list", str(inp / "shared_set.txt"),
       "--out-dir", str(out)]
subprocess.run(cmd, check=True, capture_output=True, text=True)
tab = pd.read_csv(out / "S4_gating_eqtl_fixed_RNA.reconstructed.tsv", sep="\t")
truth = json.loads((inp / "truth.json").read_text())

# --- checks ---
fails = []
# classes recovered for every planted interaction gene
cls_rows = tab[tab.cls.notna()]
if len(cls_rows) != 30:
    fails.append(("n_class_rows", len(cls_rows)))
for _, r in cls_rows.iterrows():
    g = int(r.ensg[-7:])
    if r.cls != truth["cls"][str(g)]:
        fails.append(("class", r.ensg, r.cls, truth["cls"][str(g)]))
# planted Delta-expression recovered (both matrices; tolerance ~5 SE)
for _, r in tab.iterrows():
    g = int(r.ensg[-7:])
    planted = truth["delta"][truth["cls"][str(g)]]
    if abs(r.dexpr_v3 - planted) > 0.25:
        fails.append(("dexpr_v3", r.ensg, round(float(r.dexpr_v3), 3), planted))
    if pd.isna(r.dexpr):
        fails.append(("dexpr_missing", r.ensg))
# shared rows = planted shared set; dabs = |post| - |pre|
if int((tab.group == "Shared apaQTL").sum()) != 30:
    fails.append(("n_shared", int((tab.group == "Shared apaQTL").sum())))
dd = (cls_rows.beta_post.abs() - cls_rows.beta_pre.abs() - cls_rows.dabs).abs().max()
if dd > 1e-12:
    fails.append(("dabs", float(dd)))
# statistic direction: interaction genes change expression more than shared
st = json.loads((out / "run_summary.json").read_text())["stats"]
if not (st["median_abs_dexpr_v3_interaction"] > st["median_abs_dexpr_v3_shared_S4"]):
    fails.append(("median_direction", st))
if st["P_wilcoxon_interaction_vs_shared_S4"] > 0.05:
    fails.append(("P_interaction_vs_shared", st["P_wilcoxon_interaction_vs_shared_S4"]))
if st["spearman_rho_dabs_dexpr_v3"] <= 0:
    fails.append(("spearman_sign", st["spearman_rho_dabs_dexpr_v3"]))

if fails:
    print("FAIL", fails[:12])
    sys.exit(1)
print(f"PASS gating_producer: {len(cls_rows)} class rows (planted Reversed/Amplified/Attenuated recovered), "
      f"{int((tab.group == 'Shared apaQTL').sum())} shared rows, planted Delta-expression within 0.25, "
      f"median |Delta| interaction {st['median_abs_dexpr_v3_interaction']:.3f} > shared {st['median_abs_dexpr_v3_shared_S4']:.3f}")
