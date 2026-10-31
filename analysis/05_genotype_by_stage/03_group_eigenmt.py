# Group-level G x Stage calls with LD-aware eigenMT (tensorQTL) and BH correction.
# Reads the nominal pair table from 01_tensorqtl.sh and writes data/05_gxs/interaction_*.
from pathlib import Path
import gzip
import os

import numpy as np
import pandas as pd
import polars as pl
import torch

torch.set_num_threads(4)
torch.set_num_interop_threads(1)

from tensorqtl import eigenmt, genotypeio
from tensorqtl.core import filter_maf_interaction, impute_mean

# --- paths and settings ---
APA_ROOT = Path(os.environ["APA_ROOT"])
D = APA_ROOT / "data" / "05_gxs"
Q4 = APA_ROOT / "data" / "04_qtl"
PAIR_FILE = D / "interaction.cis_qtl_pairs.txt.gz"
GROUP_FILE = D / "phenotype.group"
PHENOTYPE_FILE = D / "phenotype.GeneRange.bed.gz"
INTERACTION_FILE = D / "LifeStage.txt"
GENOTYPE_PREFIX = str(APA_ROOT / "raw" / "genotype")
D.mkdir(parents=True, exist_ok=True)
GROUP_Q_CUTOFF = 0.05
PAIR_CHUNK_SIZE = 250_000
PAIR_COLUMNS = [
    "phenotype_id", "variant_id", "start_distance", "end_distance", "af",
    "ma_samples", "ma_count", "pval_g", "b_g", "b_g_se", "pval_i", "b_i",
    "b_i_se", "pval_gi", "b_gi", "b_gi_se",
]
PAIR_SCHEMA = {
    "start_distance": pl.Int64, "end_distance": pl.Int64, "af": pl.Float64,
    "ma_samples": pl.Int64, "ma_count": pl.Int64, "pval_g": pl.Float64,
    "b_g": pl.Float64, "b_g_se": pl.Float64, "pval_i": pl.Float64,
    "b_i": pl.Float64, "b_i_se": pl.Float64, "pval_gi": pl.Float64,
    "b_gi": pl.Float64, "b_gi_se": pl.Float64,
}
PANDAS_DTYPES = {
    "phenotype_id": "string", "variant_id": "string", "start_distance": "Int64",
    "end_distance": "Int64", "af": "float64", "ma_samples": "Int64",
    "ma_count": "Int64", "pval_g": "float64", "b_g": "float64",
    "b_g_se": "float64", "pval_i": "float64", "b_i": "float64",
    "b_i_se": "float64", "pval_gi": "float64", "b_gi": "float64",
    "b_gi_se": "float64",
}


# --- helpers ---
def read_group_qvalues(path):
    stage = pd.read_csv(path, sep="\t", compression="gzip", usecols=["group_id", "qval"],
                        dtype={"group_id": "string", "qval": "float64"})
    stage = stage.dropna(subset=["group_id", "qval"]).sort_values(["group_id", "qval"], kind="mergesort")
    return stage.drop_duplicates("group_id", keep="first")


def read_pair_keys(path):
    pairs = pd.read_csv(path, sep="\t", compression="gzip", usecols=["phenotype_id", "variant_id"], dtype="string")
    return set(zip(pairs["phenotype_id"], pairs["variant_id"]))


# --- group map ---
group_map = pd.read_csv(GROUP_FILE, sep="\t", header=None, names=["phenotype_id", "group_id"], dtype="string")
group_map["gene_id"] = group_map["group_id"].str.split("|", n=1).str[0]

# --- eigenMT effective tests per group ---
phenotype_df = pd.read_csv(PHENOTYPE_FILE, sep="\t", index_col=3, dtype={"#Chr": str, "#chr": str})
phenotype_df.rename(columns={i: i.lower().replace("#chr", "chr") for i in phenotype_df.columns[:3]}, inplace=True)
phenotype_df["start"] += 1
phenotype_pos_df = phenotype_df[["chr", "start", "end"]].copy()
phenotype_df.drop(["chr", "start", "end"], axis=1, inplace=True)
interaction_df = pd.read_csv(INTERACTION_FILE, sep="\t", index_col=0)
interaction_df = interaction_df.loc[phenotype_df.columns].astype(np.float32)
stage = interaction_df.iloc[:, 0].to_numpy(dtype=np.float32)
if not np.all(np.isin(stage, [0.0, 1.0])):
    raise ValueError("LifeStage.txt must contain binary labels 0 and 1")
stage_counts = {value: int(np.sum(stage == value)) for value in (0.0, 1.0)}
if stage_counts != {0.0: 189, 1.0: 508}:
    raise ValueError(f"Unexpected LifeStage counts after sample alignment: {stage_counts}")
group_s = group_map.set_index("phenotype_id").loc[phenotype_df.index, "group_id"]
genotype_df, variant_df = genotypeio.load_genotypes(GENOTYPE_PREFIX, select_samples=phenotype_df.columns)
genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
interaction_mask_t = torch.as_tensor(stage == 1.0, dtype=torch.bool)
igc = genotypeio.InputGeneratorCis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                                   group_s=group_s, window=1_000_000)
records = []
for _, genotypes, _, phenotype_ids, group_id in igc.generate_data(verbose=True):
    genotypes_t = torch.tensor(genotypes, dtype=torch.float32)
    genotypes_t = genotypes_t[:, genotype_ix]
    impute_mean(genotypes_t)
    genotypes_t, _ = filter_maf_interaction(genotypes_t, interaction_mask_t=interaction_mask_t,
                                            maf_threshold_interaction=0.05)
    if genotypes_t.shape[0] == 0:
        continue
    records.append({
        "group_id": str(group_id),
        "num_phenotypes": len(phenotype_ids),
        "tests_emt": int(eigenmt.compute_tests(genotypes_t, var_thresh=0.99, variant_window=200)),
    })
top = pd.DataFrame.from_records(records)
top.to_csv(D / "interaction_group_eigenmt_top.tsv", sep="\t", index=False)

# --- group p-values (min nominal p x phenotypes x eigenMT tests) ---
group_lf = pl.from_pandas(group_map[["phenotype_id", "group_id"]].drop_duplicates()).lazy()
pair_lf = pl.scan_csv(str(PAIR_FILE), separator="\t", has_header=True, schema_overrides=PAIR_SCHEMA,
                      null_values=["NA", "NaN", ""]).select(PAIR_COLUMNS).join(group_lf, on="phenotype_id", how="left")
stats = pair_lf.group_by("group_id").agg(
    pl.len().alias("n_tested_pairs"),
    pl.col("phenotype_id").n_unique().alias("n_tested_phenotypes"),
    pl.col("pval_gi").count().alias("n_valid_pval_gi"),
    pl.col("pval_gi").min().alias("min_pval_gi"),
).collect(engine="streaming").to_pandas()
meta = group_map.groupby("group_id", as_index=False, sort=True).agg(
    gene_id=("gene_id", "first"), n_group_phenotypes=("phenotype_id", "nunique"))
summary = meta.merge(stats, on="group_id", how="left").merge(
    top[["group_id", "tests_emt", "num_phenotypes"]], on="group_id", how="left")
summary["n_tested_pairs"] = summary["n_tested_pairs"].fillna(0).astype("int64")
summary["n_tested_phenotypes"] = summary["n_tested_phenotypes"].fillna(0).astype("int64")
summary["n_valid_pval_gi"] = summary["n_valid_pval_gi"].fillna(0).astype("int64")
summary["min_pval_gi"] = pd.to_numeric(summary["min_pval_gi"], errors="coerce")
summary["group_p_eigenmt"] = np.minimum(1.0, summary["min_pval_gi"] * summary["num_phenotypes"] * summary["tests_emt"])

# --- BH across groups ---
values = pd.to_numeric(summary["group_p_eigenmt"], errors="coerce").to_numpy(dtype=float)
values = np.where(np.isfinite(values), np.clip(values, 0.0, 1.0), 1.0)
order = np.argsort(values, kind="mergesort")
ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
bh = np.empty_like(adjusted)
bh[order] = np.clip(adjusted, 0.0, 1.0)
summary["group_q_eigenmt"] = bh
summary["primary_significant"] = summary["group_q_eigenmt"] < GROUP_Q_CUTOFF

# --- stage annotation from the permutation results ---
prenatal = read_group_qvalues(Q4 / "prenatal.cis_qtl.txt.gz").rename(columns={"qval": "prenatal_group_qval"})
postnatal = read_group_qvalues(Q4 / "postnatal.cis_qtl.txt.gz").rename(columns={"qval": "postnatal_group_qval"})
summary = summary.merge(prenatal, on="group_id", how="left").merge(postnatal, on="group_id", how="left")
summary["prenatal_group_significant"] = summary["prenatal_group_qval"].lt(GROUP_Q_CUTOFF).fillna(False)
summary["postnatal_group_significant"] = summary["postnatal_group_qval"].lt(GROUP_Q_CUTOFF).fillna(False)
summary["stage_group_status"] = np.select([
    summary["prenatal_group_significant"] & summary["postnatal_group_significant"],
    summary["prenatal_group_significant"], summary["postnatal_group_significant"],
], ["Shared group", "Prenatal detected group", "Postnatal detected group"], default="Interaction-only group")

# --- lead pair per group and significant-group pair table (chunked scan) ---
lookup = group_map.set_index("phenotype_id")["group_id"].to_dict()
significant_ids = set(summary.loc[summary["primary_significant"], "group_id"].astype(str))
leads = {}
total_rows = 0
unmapped_rows = 0
pair_output = D / "interaction_pair_summary.tsv.gz"
with gzip.open(pair_output, "wt", encoding="utf-8", newline="") as handle:
    first_output = True
    for chunk in pd.read_csv(PAIR_FILE, sep="\t", compression="gzip", usecols=PAIR_COLUMNS,
                             dtype=PANDAS_DTYPES, chunksize=PAIR_CHUNK_SIZE):
        total_rows += len(chunk)
        chunk["group_id"] = chunk["phenotype_id"].map(lookup)
        unmapped_rows += int(chunk["group_id"].isna().sum())
        valid = chunk.dropna(subset=["group_id", "pval_gi"])
        candidates = valid.sort_values(["group_id", "pval_gi", "pval_g", "phenotype_id", "variant_id"],
                                       kind="mergesort").drop_duplicates("group_id", keep="first")
        for row in candidates.to_dict("records"):
            group_id = str(row["group_id"])
            old = leads.get(group_id)
            new_key = (float(row["pval_gi"]), float(row["pval_g"]), str(row["phenotype_id"]), str(row["variant_id"]))
            old_key = None if old is None else (float(old["pval_gi"]), float(old["pval_g"]),
                                                str(old["phenotype_id"]), str(old["variant_id"]))
            if old is None or new_key < old_key:
                leads[group_id] = row
        sig_chunk = valid[valid["group_id"].isin(significant_ids)][["group_id"] + PAIR_COLUMNS]
        if not sig_chunk.empty:
            sig_chunk.to_csv(handle, sep="\t", index=False, header=first_output, na_rep="NA")
            first_output = False
    if first_output:
        pd.DataFrame(columns=["group_id"] + PAIR_COLUMNS).to_csv(handle, sep="\t", index=False)
if unmapped_rows:
    raise RuntimeError(f"{unmapped_rows} nominal pair rows did not map to phenotype.group")
lead_df = pd.DataFrame.from_records(list(leads.values()))
summary = summary.merge(lead_df[["group_id"] + PAIR_COLUMNS].rename(columns={c: f"lead_{c}" for c in PAIR_COLUMNS}),
                        on="group_id", how="left").sort_values(["group_q_eigenmt", "group_id"], kind="mergesort")
significant = summary[summary["primary_significant"]].copy()
summary.to_csv(D / "interaction_group_summary.tsv", sep="\t", index=False)
significant.to_csv(D / "interaction_significant_groups.tsv", sep="\t", index=False)

# --- lead-pair stage status (exact phenotype_id + variant_id match) ---
prenatal_pairs = read_pair_keys(Q4 / "prenatal.sig_QTL.txt.gz")
postnatal_pairs = read_pair_keys(Q4 / "postnatal.sig_QTL.txt.gz")
status = significant[["group_id", "gene_id", "group_q_eigenmt", "stage_group_status",
                      "lead_phenotype_id", "lead_variant_id"]].copy()
status["prenatal_lead_pair_significant"] = [(str(p), str(v)) in prenatal_pairs
                                           for p, v in zip(status["lead_phenotype_id"], status["lead_variant_id"])]
status["postnatal_lead_pair_significant"] = [(str(p), str(v)) in postnatal_pairs
                                            for p, v in zip(status["lead_phenotype_id"], status["lead_variant_id"])]
status["lead_pair_status"] = np.select([
    status["prenatal_lead_pair_significant"] & status["postnatal_lead_pair_significant"],
    status["prenatal_lead_pair_significant"], status["postnatal_lead_pair_significant"],
], ["Shared lead pair", "Prenatal detected lead pair", "Postnatal detected lead pair"],
    default="Interaction-only lead pair")
status["pair_key_definition"] = "phenotype_id + variant_id; exact string match"
status.to_csv(D / "interaction_stage_mapping_status.tsv", sep="\t", index=False)
print(f"[done] {len(summary)} tested groups; {len(significant)} LD-aware primary groups; {len(lead_df)} groups with leads")
