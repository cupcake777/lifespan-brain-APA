#!/usr/bin/env python3
# Reconstructs module 04 tables: ED7 seed recovery, apaQTL totals, mediation call counts.
# Writes *.reconstructed.tsv to data/04_qtl and prints the comparison with the existing tables.

import json
import os
from pathlib import Path

import pandas as pd

# --- paths ---
APA_ROOT = Path(os.environ["APA_ROOT"])
Q4 = APA_ROOT / "data" / "04_qtl"

# --- helpers ---
def gene_key(pid):
    # gene symbol = 2nd '|' field of phenotype_id (ENST|SYMBOL|locus|strand);
    # the full phenotype_id (with transcript version) gives different counts
    parts = str(pid).split("|")
    return parts[1] if len(parts) > 1 else str(pid)

# --- ED7 seed recovery ---
# a seed recovers a prenatal apaQTL gene when its qval < 0.05
pre_sig = pd.read_csv(Q4 / "prenatal.sig_QTL.txt.gz", sep="\t", usecols=["phenotype_id"])
pre_sym = {gene_key(x) for x in pre_sig.phenotype_id.unique()}
counts = pd.Series(0, index=sorted(pre_sym), dtype=int)
rows = []
seeds = sorted(Q4.glob("seed_*_postnatal.cis_qtl.txt.gz"))
if not seeds:
    raise SystemExit("no seed_*_postnatal.cis_qtl.txt.gz in data/04_qtl (files there are named *_post.cis_qtl.txt.gz)")
for f in seeds:
    d = pd.read_csv(f, sep="\t", usecols=["phenotype_id", "qval"])
    rec = {gene_key(x) for x in d.loc[d.qval < 0.05, "phenotype_id"]} & pre_sym
    rows.append({"seed": f.name, "n_recovered": len(rec),
                 "frac": len(rec) / len(pre_sym), "n_prenatal_genes": len(pre_sym)})
    counts[list(rec)] += 1
ed7c = pd.DataFrame(rows)
per_gene = pd.DataFrame({"gene": counts.index, "n_recovered_seeds": counts.values})
ed7c.to_csv(Q4 / "ed7c_seed_recovery.reconstructed.tsv", sep="\t", index=False)
per_gene.to_csv(Q4 / "ed7_per_gene_recovery.n_recovered_seeds.reconstructed.tsv", sep="\t", index=False)

# compare with the existing tables
ed7_res = {}
ref = pd.read_csv(Q4 / "ed7c_seed_recovery.tsv", sep="\t")
m = ed7c.merge(ref, on="seed", suffixes=("", "_ref"))
ed7_res["ed7c_rows_matched"] = int(len(m))
ed7_res["ed7c_max_abs_diff_n_recovered"] = int((m.n_recovered - m.n_recovered_ref).abs().max())
ed7_res["ed7c_max_abs_diff_frac"] = float((m.frac - m.frac_ref).abs().max())
ref = pd.read_csv(Q4 / "ed7_per_gene_recovery.tsv", sep="\t", usecols=["gene", "n_recovered_seeds"])
m = per_gene.merge(ref, on="gene", suffixes=("", "_ref"))
ed7_res["per_gene_rows_matched"] = int(len(m))
ed7_res["per_gene_rows_ref"] = int(len(ref))
ed7_res["per_gene_exact_count"] = int((m.n_recovered_seeds == m.n_recovered_seeds_ref).sum())

# --- apaQTL totals ---
post_sig = pd.read_csv(Q4 / "postnatal.sig_QTL.txt.gz", sep="\t", usecols=["phenotype_id"])
rows = []
for stage, d in [("prenatal", pre_sig), ("postnatal", post_sig)]:
    rows.append({"stage": stage, "apaqtl_genes": int(d.phenotype_id.map(gene_key).nunique()),
                 "significant_pairs": int(len(d))})
t = pd.DataFrame(rows)
t.to_csv(Q4 / "apaqtl_totals.reconstructed.tsv", sep="\t", index=False)
totals_res = {"totals": t.to_dict("records"),
              "expected": {"prenatal_genes": 669, "postnatal_genes": 2249,
                           "prenatal_pairs": 59761, "postnatal_pairs": 261336}}

# --- mediation call counts ---
# registry R028/R029 "mediated in 162": 'Partial mediation' + 'APA -> EXP' (rule inferred from the labels)
med = pd.read_csv(Q4 / "04_mediation_results_module06.tsv", sep="\t")
cnt = med.groupby(["stage", "mediation_call"]).size().rename("N").reset_index()
mediated = med[med.mediation_call.isin(["Partial mediation", "APA -> EXP"])]
med_res = {"n_pairs": int(len(med)), "n_mediated_rule": int(len(mediated)),
           "frac_mediated_rule": float(len(mediated) / len(med))}
cnt.to_csv(Q4 / "mediation_call_counts.reconstructed.tsv", sep="\t", index=False)
ref = pd.read_csv(Q4 / "04_mediation_summary.tsv", sep="\t")
m = cnt.merge(ref[["stage", "mediation_call", "N"]], on=["stage", "mediation_call"], how="outer",
              suffixes=("", "_ref"), indicator=True)
med_res["summary_rows_both"] = int((m["_merge"] == "both").sum())
med_res["summary_N_exact"] = int((m.N == m.N_ref).sum())
med_res["summary_rows_total"] = int(len(m))

# --- summary ---
summary = {"ed7": ed7_res, "apaqtl_totals": totals_res, "mediation": med_res}
(Q4 / "run_summary_04.json").write_text(json.dumps(summary, indent=2, default=str))
print(json.dumps(summary, indent=2, default=str))
