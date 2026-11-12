# Synthetic test: HC3 refit functions of 05_gxs/04_hc3_refit.py recover planted
# stage-specific lead effects and assign the planted G x Stage classes.
# The script runs its full pipeline at import time, so only its function definitions
# are extracted (ast) and executed here.
import ast
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
RELEASE = HERE.parents[1]
sys.path.insert(0, str(HERE))
import make_qtl as gen  # noqa: E402

REFIT = RELEASE / "05_gxs" / "04_hc3_refit.py"
WANT = {"residualize_centered", "make_residualizer", "fit_ols", "point_estimate_class",
        "within_stage_standardize"}

# --- load functions ---
ns = {"np": np, "stats": stats}
for node in ast.parse(REFIT.read_text()).body:
    if isinstance(node, ast.FunctionDef) and node.name in WANT:
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(REFIT), "exec"), ns)
missing = WANT - set(ns)
assert not missing, f"functions not found: {missing}"

# --- simulate ---
cls, betas, _ = gen.planted_truth()
inp = gen.write_inputs(Path(tempfile.mkdtemp(prefix="qtl_hc3_")) / "in")
d = np.load(inp / "hc3_inputs.npz")
dos, pheno, cov, stage = d["dos"], d["pheno"], d["cov"], d["stage"]
q = ns["make_residualizer"](cov)

# --- refit and check ---
fails = []
n_checked = 0
for g in range(gen.N_GENES):
    if cls[g] not in ("Reversed", "Amplified", "Attenuated"):
        continue
    n_checked += 1
    raw = [dos[g] * (1 - stage), dos[g] * stage]
    fit = ns["fit_ols"](pheno[g], raw, q, covariate_count=cov.shape[1])
    bh = fit["beta"]
    bp, bq = betas[cls[g]]
    tol_pre = 4 * fit["hc3_se"][0]
    tol_post = 4 * fit["hc3_se"][1]
    if abs(bh[0] - bp) > tol_pre or abs(bh[1] - bq) > tol_post:
        fails.append((g, cls[g], float(bh[0]), float(bh[1])))
    ratio = fit["hc3_se"] / fit["ordinary_se"]
    if not np.all((ratio > 0.8) & (ratio < 1.25)):
        fails.append((g, "hc3_vs_ordinary_se", ratio.tolist()))
    label = ns["point_estimate_class"](bh[0], bh[1])
    label = "Reversed" if label == "Opposing" else label
    if label != cls[g]:
        fails.append((g, "class", cls[g], label))

# within-stage standardisation: mean 0, sd 1 per stage
z = ns["within_stage_standardize"](pheno[0], stage)
for s in (0, 1):
    m = z[stage == s]
    assert abs(m.mean()) < 1e-10 and abs(m.std(ddof=1) - 1) < 1e-10, "within-stage standardisation failed"
assert n_checked == 30, n_checked
if fails:
    print("FAIL", fails[:10])
    sys.exit(1)
print(f"PASS hc3_refit_functions: {n_checked} planted lead pairs refit; betas within 4 HC3 SE; "
      f"planted classes (Reversed/Amplified/Attenuated) assigned; HC3/ordinary SE ratio in (0.8,1.25)")
