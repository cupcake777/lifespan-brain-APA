#!/usr/bin/env python3
# Writes synthetic MAGMA gene-set (.gsa.out) files: a planted set (p = 3e-4) and two null sets.
import os

import numpy as np

# --- helpers ---

HEADER = """# TOTAL_GENES = 18000
# TEST_DIRECTION = one-sided, positive (set), two-sided (covar)
# CONDITIONED_INTERNAL = gene size, gene density, sample size, inverse mac
"""


def write_gsa(path, rows, conditioned=None):
    with open(path, "w") as fh:
        fh.write(HEADER)
        if conditioned:
            fh.write(f"# CONDITIONED_VARIABLES = {conditioned}\n")
            fh.write("VARIABLE       TYPE  MODEL  NGENES         BETA     BETA_STD           SE            P\n")
            for name, ngenes, beta, bstd, se, p in rows:
                fh.write(f"{name:<15}SET   1      {ngenes:<6}  {beta:<10.6g}{bstd:<12.6g}{se:<12.6g}  {p:.6g}\n")
        else:
            fh.write("VARIABLE           TYPE  NGENES         BETA     BETA_STD           SE            P\n")
            for name, ngenes, beta, bstd, se, p in rows:
                fh.write(f"{name:<19}SET   {ngenes:<6}  {beta:<10.6g}{bstd:<12.6g}{se:<12.6g}  {p:.6g}\n")


# --- fixture ---

def make(outdir, seed=2026):
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(seed)
    sets = []
    for name, p in [("PLANTED_SET", 3e-4), ("NULL_SET_A", float(rng.uniform(0.2, 0.9))),
                    ("NULL_SET_B", float(rng.uniform(0.2, 0.9)))]:
        sets.append((name, 300 + int(rng.integers(0, 500)), 0.1 if p < 0.01 else 0.01,
                     0.02 if p < 0.01 else 0.003, 0.03, p))
    a = os.path.join(outdir, "synth_gwas_sets.gsa.out")
    write_gsa(a, sets)
    cond = [("PLANTED_SET", 320, 0.09, 0.02, 0.03, 0.0021), ("NULL_SET_A", 410, 0.001, 0.0002, 0.04, 0.6)]
    b = os.path.join(outdir, "synth_gwas_conditioned.gsa.out")
    write_gsa(b, cond, conditioned="postnatal (set)")
    return a, b, sets
