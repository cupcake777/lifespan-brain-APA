#!/usr/bin/env python3
# Checks magma_gsa_sum.py on synthetic MAGMA .gsa.out files (magma_fixture.py).
# Checks: parsed rows and planted p, summary columns and VARIABLE filter, producer FDR equal to an
# independent BH (within-GWAS and global), and the planted set significant within its GWAS.
import os
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import magma_fixture as fx  # noqa: E402

sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "06_disease", "reconstructed")))
import magma_gsa_sum as prod  # noqa: E402

# --- helpers ---

def bh_loop(p):
    """Independent BH implementation (loop form) for the test."""
    m = len(p)
    order = sorted(range(m), key=lambda i: p[i])
    q = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        val = min(prev, p[i] * m / rank)
        q[i] = val
        prev = val
    return q


# --- run ---

fails = []
with tempfile.TemporaryDirectory() as tmp:
    a, b, sets = fx.make(tmp, seed=2026)

    rows = prod.parse_gsa_out(a)
    names = [r["VARIABLE"] for r in rows]
    if names != ["PLANTED_SET", "NULL_SET_A", "NULL_SET_B"]:
        fails.append(f"parsed variables {names}")
    p_planted = float(rows[0]["P"])
    print(f"[parse] PLANTED_SET NGENES={rows[0]['NGENES']} P={p_planted:g} (planted 0.0003)")
    if abs(p_planted - 3e-4) > 1e-12:
        fails.append("planted p not parsed")

    summ = prod.condition_summary([("synth_cond", b, ["PLANTED_SET"]), ("synth_sets", a, [])])
    cols = list(summ.columns)
    if cols != ["analysis", "variable", "type", "ngenes", "beta", "beta_std", "se", "p", "source"]:
        fails.append(f"summary columns {cols}")
    cond_rows = summ[summ.analysis == "synth_cond"]
    if list(cond_rows.variable) != ["PLANTED_SET"]:
        fails.append("VARIABLE filter not applied")
    print(summ.to_string(index=False))

    # build a 3-set table per GWAS (planted and an independent null GWAS) and check FDR
    tab = pd.read_csv(os.path.join(tmp, "synth_gwas_sets.gsa.out"), sep=r"\s+", comment="#")
    tab = tab.rename(columns={"P": "p"})
    tab["gwas_id"] = "SYNTH_GWAS"
    null_p = np.random.default_rng(7).uniform(0.05, 0.9, 3)
    extra = pd.DataFrame({"VARIABLE": ["A", "B", "C"], "TYPE": ["SET"] * 3, "NGENES": [100, 200, 300],
                          "BETA": [0.1] * 3, "BETA_STD": [0.01] * 3, "SE": [0.02] * 3, "p": null_p,
                          "gwas_id": "SYNTH_GWAS2"})
    tab = pd.concat([tab, extra], ignore_index=True)
    tab["fdr_within_gwas"] = np.nan
    for g, idx in tab.groupby("gwas_id").groups.items():
        tab.loc[idx, "fdr_within_gwas"] = bh_loop(tab.loc[idx, "p"].tolist())
    tab["fdr_global"] = bh_loop(tab["p"].tolist())
    tpath = os.path.join(tmp, "check.tsv")
    tab.to_csv(tpath, sep="\t", index=False)
    res = prod.check_fdr(tab, ["gwas_id"], "p", "fdr_within_gwas", "fdr_global")
    print(f"[fdr] {res}")
    if res["max_abs_diff_within"] > 1e-12 or res["max_abs_diff_global"] > 1e-12:
        fails.append("producer FDR does not match independent BH")

    q_planted = float(tab.loc[tab.VARIABLE == "PLANTED_SET", "fdr_within_gwas"].iloc[0])
    q_null = float(tab.loc[tab.VARIABLE == "NULL_SET_A", "fdr_within_gwas"].iloc[0])
    print(f"[planted] PLANTED_SET within-GWAS q = {q_planted:.3g}; NULL_SET_A q = {q_null:.3g}")
    if not (q_planted < 0.05):
        fails.append("planted set not significant after BH")

if fails:
    print("FAIL:", "; ".join(fails))
    sys.exit(1)
print("PASS test_magma")
