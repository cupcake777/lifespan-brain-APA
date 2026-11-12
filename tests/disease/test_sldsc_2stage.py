#!/usr/bin/env python3
# Checks the two-stage S-LDSC FDR on a synthetic processed table (sldsc_fixture.py).
# Checks: planted postnatal signal recovered, q values equal hand values, null-like traits not significant,
# q independent of lifespan/interaction p values, and the one-BH FDR_p reproduced.
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sldsc_fixture as fx  # noqa: E402

sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "06_disease", "reconstructed")))
import sldsc_two_stage as prod  # noqa: E402

# --- helpers ---

def approx(a, b, tol=1e-12):
    return abs(a - b) <= tol * max(1.0, abs(b))


def q(two, pop, trait, stage):
    r = two[(two.Population == pop) & (two.Trait == trait) & (two.Stage == stage)]
    assert len(r) == 1, (pop, trait, stage, len(r))
    return float(r["FDR_two_stage"].iloc[0])


# --- run ---

fails = []
df_sim = fx.make(seed=2026)
with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, "sim_processed.tsv")
    df_sim.to_csv(path, sep="\t", index=False)
    df = prod.load_processed(path)
    two = prod.two_stage_fdr(df)

    q_post = q(two, "EUR", "SIM1", "Postnatal")
    q_pre = q(two, "EUR", "SIM1", "Prenatal")
    print(f"[planted] EURSIM1 postnatal q = {q_post:.6g} (expected 2e-4); prenatal q = {q_pre:.6g} (expected 0.02)")
    if not (approx(q_post, 2e-4, 1e-9) and q_post < 0.05):
        fails.append("planted postnatal signal not recovered")
    if not approx(q_pre, 0.02, 1e-9):
        fails.append("prenatal two-stage q not equal to the hand value")
    for stage in ["Prenatal", "Postnatal"]:
        if q(two, "EUR", "SIM2", stage) < 0.05:
            fails.append(f"null-like EURSIM2 {stage} called significant")

    # independence from lifespan/interaction rows
    df_mod = df.copy()
    mask = df_mod["Annotation"].isin(["lifespan_apa_qtl", "interaction_apa_qtl"])
    df_mod.loc[mask, "Enrichment_p"] = 1e-9
    two_mod = prod.two_stage_fdr(df_mod)
    q_mod = float(two_mod[(two_mod.Trait == "SIM1") & (two_mod.Population == "EUR") &
                          (two_mod.Stage == "Postnatal")]["FDR_two_stage"].iloc[0])
    print(f"[independence] postnatal q after changing lifespan/interaction p: {q_mod:.6g}")
    if not approx(q_mod, q_post, 1e-12):
        fails.append("two-stage q depends on non-binary annotations")

    # one BH over all rows (what 07_collect_results.py writes as FDR_p)
    g = prod.check_global76(df)
    print(f"[global76] compared {g['n_compared']} rows; max rel diff = {g['max_rel_diff']:.3g}")
    if g["n_compared"] != len(df) or g["max_rel_diff"] > 1e-5:
        fails.append("global BH FDR_p not reproduced")

if fails:
    print("FAIL:", "; ".join(fails))
    sys.exit(1)
print("PASS test_sldsc_2stage")
