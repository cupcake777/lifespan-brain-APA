# Synthetic inputs for the QTL tests: fabricated IDs (ENSG9000xxxxx, SYNG_xxx, donors S0001...).
# Planted G x Stage lead effects (b_pre, b_post): Reversed 0-9 (+0.5, -0.7), Amplified 10-19 (+0.3, +1.0),
# Attenuated 20-29 (+0.9, +0.4), Shared 30-59 (+0.4, +0.4), Null 60-199 (0).
# Planted expression change (post minus pre, plus noise): Reversed 0.9, Amplified 1.2, Attenuated 0.4, Shared 0.5, Null 0.
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

# --- constants ---
N_DONORS, N_PRE = 300, 90
N_GENES = 200
SEED = 20261019


# --- helpers ---
def planted_truth():
    cls = {}
    for g in range(N_GENES):
        if g < 10:
            cls[g] = "Reversed"
        elif g < 20:
            cls[g] = "Amplified"
        elif g < 30:
            cls[g] = "Attenuated"
        elif g < 60:
            cls[g] = "Shared"
        else:
            cls[g] = "Null"
    # effects chosen so each post-stage effect is >= 3 HC3 SE from 0 at N = 300 donors
    betas = {"Reversed": (0.5, -0.7), "Amplified": (0.3, 1.0), "Attenuated": (0.9, 0.4),
             "Shared": (0.4, 0.4), "Null": (0.0, 0.0)}
    delta = {"Reversed": 0.9, "Amplified": 1.2, "Attenuated": 0.4, "Shared": 0.5, "Null": 0.0}
    return cls, betas, delta


def gene_ids(g):
    return f"ENSG9000{g:07d}", f"ENST9000{g:07d}.1", f"SYNG_{g:03d}"


def write_inputs(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    cls, betas, delta = planted_truth()
    donors = [f"S{i:04d}" for i in range(N_DONORS)]
    stage = np.array([0] * N_PRE + [1] * (N_DONORS - N_PRE))  # 0 prenatal, 1 postnatal

    pd.DataFrame({"sample": donors, "LifeStage": np.where(stage == 0, "Prenatal", "Postnatal")}).to_csv(
        out / "meta.tsv", sep="\t", index=False)

    # GENCODE-like transcript records
    with gzip.open(out / "synthetic.gtf.gz", "wt") as fh:
        for g in range(N_GENES):
            ensg, enst, sym = gene_ids(g)
            fh.write(f'chrS\tsyn\ttranscript\t1\t100\t.\t+\t.\tgene_id "{ensg}.1"; transcript_id "{enst}"; gene_name "{sym}";\n')

    # expression matrices: fixed and non-fixed share the planted signal, different noise
    for name, seed_off in [("expr_nonfixed.tsv", 1), ("expr_fixed.tsv", 2)]:
        r2 = np.random.default_rng(SEED + seed_off)
        mat = np.zeros((N_GENES, N_DONORS))
        for g in range(N_GENES):
            base = r2.uniform(2, 6)
            mat[g] = base + delta[cls[g]] * stage + r2.normal(0, 0.3, N_DONORS)
        df = pd.DataFrame(mat, columns=donors, index=[gene_ids(g)[0] + ".15" for g in range(N_GENES)])
        df.index.name = "gene_id"
        df.to_csv(out / name, sep="\t")

    # cis dosages and phenotypes for the HC3 refit test (one lead variant per gene)
    cov = rng.normal(size=(N_DONORS, 2))
    dos = rng.binomial(2, 0.3, size=(N_GENES, N_DONORS)).astype(float)
    pheno = np.zeros((N_GENES, N_DONORS))
    for g in range(N_GENES):
        bp, bq = betas[cls[g]]
        pheno[g] = (bp * dos[g] * (1 - stage) + bq * dos[g] * stage
                    + cov @ np.array([0.5, -0.3]) + rng.normal(0, 1.0, N_DONORS))
    np.savez(out / "hc3_inputs.npz", dos=dos, pheno=pheno, cov=cov, stage=stage)

    # HC3-style lead table (betas are the planted values)
    rows = []
    for g in range(N_GENES):
        if cls[g] in ("Reversed", "Amplified", "Attenuated"):
            ensg, enst, sym = gene_ids(g)
            bp, bq = betas[cls[g]]
            rows.append({"group_id": sym, "gene_id": sym, "lead_phenotype_id": f"{enst}|{sym}|chrS:1-100|+",
                         "lead_variant_id": f"chrS:{1000+g}:A:G", "beta_pre": bp, "beta_post": bq})
    pd.DataFrame(rows).to_csv(out / "hc3_lead_table.tsv", sep="\t", index=False)

    shared = [gene_ids(g)[0] for g in range(N_GENES) if cls[g] == "Shared"]
    (out / "shared_set.txt").write_text("\n".join(shared) + "\n")
    (out / "truth.json").write_text(json.dumps({"cls": {str(k): v for k, v in cls.items()},
                                               "delta": delta, "betas": betas}, indent=1))
    return out
