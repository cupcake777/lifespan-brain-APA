#!/usr/bin/env bash
# APA-QTL x GWAS SMR: BESD build, SMR + HEIDI, aggregation and locus plots.
# Submit from $APA_ROOT: sbatch release/07_coloc/smr/run_smr.sh
#SBATCH --job-name=apa_smr
#SBATCH --output=logs/smr_%j.out
#SBATCH --error=logs/smr_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=CU
: "${APA_ROOT:?source release/config/paths.local.sh}"
: "${APA_HPC_HOME:?source release/config/paths.local.sh}"
set -euo pipefail

# --- paths ---
SMR_HELPERS="${APA_ROOT}/release/07_coloc/smr"
MAKE_BESD_SCRIPT="${SMR_HELPERS}/01_make_besd.py"
HARMONY_TO_MA_SCRIPT="${SMR_HELPERS}/02_harmony_to_ma.py"
DATA_DIR="${APA_ROOT}/data/07_coloc"
PLOT_DIR="${APA_ROOT}/figures/07_coloc"
QTL_INPUT="${QTL_INPUT:-${APA_HPC_HOME}/project/work/apa/04_coloc/input}"
HARMONY_DIR="${HARMONY_DIR:-${APA_GROUP_SHARE}/.ref/GWAS/harmony}"
LD_SRC="${LD_SRC:-${APA_GROUP_SHARE}/.ref/1KG/plink/hg38/1kgp_hg38_bed_biallelic}"
LD_KEEP_DIR="${APA_HPC_HOME}/project/work/apa/04_coloc/refs/plink_keep_runtime"
BESD_DIR="${DATA_DIR}/besd"
GWAS_MA_DIR="${DATA_DIR}/gwas_ma"
SMR_OUT_DIR="${DATA_DIR}/smr_output"
LD_REF_DIR="${LD_REF_DIR:-${DATA_DIR}/refs/plink_hg38_uniqueid_runtime}"
LD_REF="${LD_REF_DIR}/1kgp_hg38_chrpos_ref_alt"
SMR="${SMR:-$(command -v smr || true)}"
CLEAN_QFILE_AFTER_BESD="${CLEAN_QFILE_AFTER_BESD:-1}"
QTL_TYPES=(prenatal postnatal interaction)
AGGREGATE_ONLY="${AGGREGATE_ONLY:-0}"
mkdir -p "${BESD_DIR}" "${GWAS_MA_DIR}" "${SMR_OUT_DIR}" "${PLOT_DIR}" "${APA_ROOT}/logs"

if [[ "${AGGREGATE_ONLY}" != "1" ]]; then

    if [[ -z "${SMR}" || ! -x "${SMR}" ]]; then
        echo "ERROR: SMR executable not found. Set SMR to an executable path or add smr to PATH." >&2
        exit 1
    fi

    # --- unique-ID LD reference ---
    # QTL variant_id uses chr:pos:ref:alt, so the PLINK .bim IDs must use the same ID space.
    # Keep the large genotype payload shared; only rewrite the text .bim and mapping locally.
    mkdir -p "${LD_REF_DIR}"
    ln -sfn "${LD_SRC}.bed" "${LD_REF}.bed"
    ln -sfn "${LD_SRC}.fam" "${LD_REF}.fam"
    if [[ ! -s "${LD_REF}.bim" || "${LD_REF}.bim" -ot "${LD_SRC}.bim" ]]; then
        echo "[$(date)] Preparing unique-ID LD reference: ${LD_REF}"

        DUP_IDS="${LD_REF}.duplicate_base_ids.tmp"
        awk 'BEGIN{OFS="\t"} {chrom=$1; sub(/^chr/,"",chrom); print chrom ":" $4 ":" $6 ":" $5}' "${LD_SRC}.bim" \
            | LC_ALL=C sort -T "${LD_REF_DIR}" -S "${SORT_MEM:-16G}" --parallel="${N_CORES:-8}" \
            | uniq -d > "${DUP_IDS}"

        awk -v dup_ids="${DUP_IDS}" 'BEGIN{
                OFS="\t"
                while ((getline id < dup_ids) > 0) dup[id]=1
                close(dup_ids)
            }
            {
                chrom=$1
                sub(/^chr/,"",chrom)
                base=chrom ":" $4 ":" $6 ":" $5
                id=base
                if (base in dup) {
                    dup_seen[base]++
                    if (dup_seen[base] > 1) id=base "__dup" dup_seen[base]
                }
                print $1,id,$3,$4,$5,$6
            }' "${LD_SRC}.bim" > "${LD_REF}.bim.tmp"

        mv "${LD_REF}.bim.tmp" "${LD_REF}.bim"
        {
            echo "source=${LD_SRC}.bim"
            echo "shared_bed=${LD_SRC}.bed"
            echo "shared_fam=${LD_SRC}.fam"
            echo "duplicate_base_ids=$(wc -l < "${DUP_IDS}")"
            echo "wrote_variants=$(wc -l < "${LD_REF}.bim")"
        } > "${LD_REF}.complete"
    else
        echo "[$(date)] Reusing unique-ID LD reference: ${LD_REF}"
    fi

    N_CORES=${SLURM_CPUS_PER_TASK:-8}
    REBUILD_FREQ_INPUTS="${REBUILD_FREQ_INPUTS:-0}"

    # --- step 1 & 2: QTL tabix -> BESD ---
    echo "[$(date)] Steps 1-2: Converting QTL tabix -> BESD"

    for QTL_TYPE in "${QTL_TYPES[@]}"; do
        if [[ "${REBUILD_FREQ_INPUTS}" == "1" && -f "${BESD_DIR}/${QTL_TYPE}.besd" ]]; then
            echo "  REBUILD_FREQ_INPUTS=1: rebuilding ${QTL_TYPE} qfile input and BESD"
            rm -f "${BESD_DIR}/${QTL_TYPE}.besd" "${BESD_DIR}/${QTL_TYPE}.esi" "${BESD_DIR}/${QTL_TYPE}.epi"
            rm -rf "${BESD_DIR}/${QTL_TYPE}_qfile"
        fi

        TYPE_DIR="${BESD_DIR}/${QTL_TYPE}_qfile"
        if [[ -f "${BESD_DIR}/${QTL_TYPE}.besd" ]]; then
            echo "  skip ${QTL_TYPE} BESD (exists)"
            if [[ "${CLEAN_QFILE_AFTER_BESD}" == "1" && -d "${TYPE_DIR}" ]]; then
                echo "  cleaning qfile intermediate for ${QTL_TYPE}"
                rm -rf "${TYPE_DIR}"
            fi
            continue
        fi

        TABIX="${QTL_INPUT}/${QTL_TYPE}.all_pairs.txt.gz"
        mkdir -p "${TYPE_DIR}"

        echo "  [1] Grouped qfile for ${QTL_TYPE} ..."
        python3 "${MAKE_BESD_SCRIPT}" \
            --tabix "${TABIX}" \
            --out-dir "${TYPE_DIR}" \
            --qtl-type "${QTL_TYPE}"

        echo "  [2] smr --make-besd for ${QTL_TYPE} ..."
        "${SMR}" \
            --qfile "${TYPE_DIR}/${QTL_TYPE}.qfile.txt" \
            --make-besd \
            --out "${BESD_DIR}/${QTL_TYPE}" \
            --thread-num "${N_CORES}"

        echo "  BESD ready: ${BESD_DIR}/${QTL_TYPE}.{besd,epi,esi}"
        if [[ "${CLEAN_QFILE_AFTER_BESD}" == "1" ]]; then
            echo "  cleaning qfile intermediate for ${QTL_TYPE}"
            rm -rf "${TYPE_DIR}"
        fi
    done

    # --- step 3: harmony TSV -> SMR .ma ---
    echo "[$(date)] Step 3: Converting EasyColoc harmony -> SMR .ma"
    if [[ "${REBUILD_FREQ_INPUTS}" == "1" ]]; then
        echo "  REBUILD_FREQ_INPUTS=1: rebuilding GWAS .ma files"
        rm -f "${GWAS_MA_DIR}"/EAS_*.ma "${GWAS_MA_DIR}"/EUR_*.ma
    fi
    python3 "${HARMONY_TO_MA_SCRIPT}" \
        --harmony-dir "${HARMONY_DIR}" \
        --out-dir "${GWAS_MA_DIR}" \
        --include-prefix EAS_ \
        --include-prefix EUR_

    # Collect available .ma files
    mapfile -t MA_FILES < <(find "${GWAS_MA_DIR}" -maxdepth 1 -type f \( -name "EAS_*.ma" -o -name "EUR_*.ma" \) | sort)
    if [[ ${#MA_FILES[@]} -eq 0 ]]; then
        echo "ERROR: No .ma GWAS files found in ${GWAS_MA_DIR}" >&2
        exit 1
    fi
    echo "  Available GWAS: ${#MA_FILES[@]} datasets"

    # --- step 4: SMR + HEIDI ---
    echo "[$(date)] Step 4: Running SMR + HEIDI"

    for QTL_TYPE in "${QTL_TYPES[@]}"; do
        BESD="${BESD_DIR}/${QTL_TYPE}"
        [[ ! -f "${BESD}.besd" ]] && { echo "  skip ${QTL_TYPE}: BESD missing"; continue; }

        for MA in "${MA_FILES[@]}"; do
            GWAS_ID="$(basename "${MA}" .ma)"
            OUT_PREFIX="${SMR_OUT_DIR}/${QTL_TYPE}_${GWAS_ID}"
            GWAS_POP="${GWAS_ID%%_*}"
            LD_KEEP="${LD_KEEP_DIR}/${GWAS_POP}.keep"

            if [[ "${GWAS_POP}" != "EAS" && "${GWAS_POP}" != "EUR" ]]; then
                echo "  WARNING: skip ${QTL_TYPE} x ${GWAS_ID}: unsupported ancestry prefix ${GWAS_POP}"
                continue
            fi

            if [[ -f "${OUT_PREFIX}.smr" || -f "${OUT_PREFIX}.msmr" ]]; then
                echo "  skip ${QTL_TYPE} x ${GWAS_ID} (exists)"
                continue
            fi

            echo "  SMR: ${QTL_TYPE} x ${GWAS_ID} using ${GWAS_POP} LD keep"
            if [[ "${GWAS_POP}" == "EUR" ]]; then
                echo "       EUR cross-ancestry sensitivity: --diff-freq-prop ${EUR_DIFF_FREQ_PROP:-1.0}"
            fi

            SMR_CMD=(
                "${SMR}"
                --bfile   "${LD_REF}" \
                --keep    "${LD_KEEP}" \
                --gwas-summary "${MA}" \
                --beqtl-summary "${BESD}" \
                --smr-multi \
                --heidi-mtd 1 \
                --heidi-min-m 3 \
                --heidi-max-m 20 \
                --cis-wind 2000 \
                --peqtl-smr 1e-3 \
                --peqtl-heidi 1.5e-3 \
                --out "${OUT_PREFIX}" \
                --thread-num "${N_CORES}"
            )
            if [[ "${GWAS_POP}" == "EUR" ]]; then
                SMR_CMD+=(--diff-freq-prop "${EUR_DIFF_FREQ_PROP:-1.0}")
            fi

            "${SMR_CMD[@]}" || echo "  WARNING: SMR failed for ${QTL_TYPE} x ${GWAS_ID}"
        done
    done

else
    echo "[$(date)] AGGREGATE_ONLY=1: reusing existing .msmr outputs"
fi

# --- step 5: aggregate + plot ---
echo "[$(date)] Step 5: Aggregating results"

NEFF_TABLE="${DATA_DIR}/gwas_neff/gwas_effective_tests.tsv"
if [[ ! -s "${NEFF_TABLE}" ]]; then
    echo "ERROR: Missing per-GWAS effective-test table: ${NEFF_TABLE}" >&2
    exit 1
fi

RESULTS_PATH="${DATA_DIR}" NEFF_TABLE="${NEFF_TABLE}" python3 - <<'PYEOF'
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

results_dir = Path(os.environ["RESULTS_PATH"])
neff_path = Path(os.environ["NEFF_TABLE"])
smr_dir  = results_dir / "smr_output"

frames = []
result_files = sorted(smr_dir.glob("*.smr")) + sorted(smr_dir.glob("*.msmr"))
for f in result_files:
    try:
        df = pd.read_csv(f, sep="\t")
        stem  = Path(f).stem               # e.g. prenatal_EAS_SCZ
        parts = stem.split("_", 1)
        df["qtl_type"] = parts[0]
        df["gwas_id"]  = parts[1] if len(parts) > 1 else "unknown"
        frames.append(df)
    except Exception as e:
        print(f"  skip {f}: {e}")

if not frames:
    print("No .smr/.msmr result files found yet.")
    sys.exit(0)

all_smr = pd.concat(frames, ignore_index=True)
neff = pd.read_csv(neff_path, sep="\t")
neff_cols = [
    "gwas_id", "population", "n_gwas_rows", "n_ld_overlap",
    "n_effective_snp", "gwas_alpha", "gwas_threshold", "prune_window_variants",
    "prune_step_variants", "prune_r2",
]
neff = neff[neff_cols]
missing_neff = sorted(set(all_smr["gwas_id"]) - set(neff["gwas_id"]))
if missing_neff:
    raise RuntimeError(f"Missing N_eff values for GWAS IDs: {missing_neff}")
all_smr = all_smr.merge(neff, on="gwas_id", how="left", validate="many_to_one")

p_smr = pd.to_numeric(all_smr.get("p_SMR"), errors="coerce")
p_gwas = pd.to_numeric(all_smr.get("p_GWAS"), errors="coerce")
gwas_threshold = pd.to_numeric(all_smr["gwas_threshold"], errors="coerce")
if "p_HEIDI" in all_smr.columns:
    p_heidi = pd.to_numeric(all_smr["p_HEIDI"], errors="coerce")
    heidi_pass_strict = p_heidi >= 0.05
else:
    p_heidi = pd.Series([pd.NA] * len(all_smr))
    heidi_pass_strict = pd.Series([False] * len(all_smr))

all_smr["p_SMR_num"] = p_smr
all_smr["p_GWAS_num"] = p_gwas
all_smr["p_HEIDI_num"] = p_heidi

group_cols = ["qtl_type", "gwas_id"]
all_smr["n_smr_tests"] = all_smr.groupby(group_cols)["p_SMR_num"].transform("count")
all_smr["p_SMR_BH"] = np.nan
for _, grp in all_smr.groupby(group_cols):
    valid = grp.loc[grp["p_SMR_num"].notna()]
    values = valid["p_SMR_num"].to_numpy(dtype=float)
    if values.size == 0:
        continue
    order = np.argsort(values)
    ranked = values[order]
    q_ranked = ranked * values.size / np.arange(1, values.size + 1)
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_values = np.empty(values.size, dtype=float)
    q_values[order] = np.clip(q_ranked, 0, 1)
    all_smr.loc[valid.index, "p_SMR_BH"] = q_values
all_smr["p_SMR_bonferroni"] = (
    all_smr["p_SMR_num"] * all_smr["n_smr_tests"]
).clip(upper=1)

smr_fdr_pass = all_smr["p_SMR_BH"] < 0.05
smr_bonferroni_pass = all_smr["p_SMR_bonferroni"] < 0.05
gwas_threshold_pass = p_gwas < gwas_threshold
robust_sig = smr_fdr_pass & gwas_threshold_pass & heidi_pass_strict
robust_bonferroni = smr_bonferroni_pass & gwas_threshold_pass & heidi_pass_strict
heidi_failed = smr_fdr_pass & gwas_threshold_pass & p_heidi.notna() & ~heidi_pass_strict
heidi_unavailable = smr_fdr_pass & gwas_threshold_pass & p_heidi.isna()
suggestive = (
    (all_smr["p_SMR_BH"] >= 0.05)
    & (all_smr["p_SMR_BH"] < 0.10)
    & gwas_threshold_pass
    & heidi_pass_strict
)

all_smr["smr_fdr_pass"] = smr_fdr_pass
all_smr["smr_bonferroni_pass"] = smr_bonferroni_pass
all_smr["gwas_threshold_pass"] = gwas_threshold_pass
all_smr["heidi_pass_strict"] = heidi_pass_strict
all_smr["robust_sig"] = robust_sig
all_smr["robust_bonferroni"] = robust_bonferroni
all_smr["suggestive_fdr_0_10"] = suggestive

sig = all_smr[robust_sig]
bonf_sig = all_smr[robust_bonferroni]
suggestive_rows = all_smr[suggestive]
failed = all_smr[heidi_failed]
unavailable = all_smr[heidi_unavailable]
all_smr.to_csv(results_dir / "smr_all_results.tsv", sep="\t", index=False)
sig.to_csv(results_dir / "smr_significant.tsv",   sep="\t", index=False)
bonf_sig.to_csv(results_dir / "smr_bonferroni_significant.tsv", sep="\t", index=False)
suggestive_rows.to_csv(results_dir / "smr_suggestive.tsv", sep="\t", index=False)
failed.to_csv(results_dir / "smr_heidi_failed.tsv", sep="\t", index=False)
unavailable.to_csv(results_dir / "smr_heidi_unavailable.tsv", sep="\t", index=False)
print(
    f"Total SMR tests: {len(all_smr):,} | FDR robust: {len(sig):,} | "
    f"Bonferroni robust: {len(bonf_sig):,} | FDR suggestive: {len(suggestive_rows):,} | "
    f"FDR-significant HEIDI failed: {len(failed):,} | HEIDI unavailable: {len(unavailable):,}"
)

# Summary table
if "p_SMR" in all_smr.columns:
    work = all_smr.copy()

    summary = (
        work.groupby(["qtl_type","gwas_id"], dropna=False)
        .agg(
            n_tested=("p_SMR_num", "size"),
            n_smr_fdr_pass=("smr_fdr_pass", "sum"),
            n_smr_bonferroni_pass=("smr_bonferroni_pass", "sum"),
            n_gwas_threshold_pass=("gwas_threshold_pass", "sum"),
            n_heidi_pass=("heidi_pass_strict", "sum"),
            n_robust_fdr=("robust_sig", "sum"),
            n_robust_bonferroni=("robust_bonferroni", "sum"),
            n_suggestive_fdr_0_10=("suggestive_fdr_0_10", "sum"),
            min_p_smr=("p_SMR_num", "min"),
            min_p_smr_bh=("p_SMR_BH", "min"),
            gwas_threshold=("gwas_threshold", "first"),
            n_effective_snp=("n_effective_snp", "first"),
        )
        .reset_index()
    )
    summary.to_csv(results_dir / "smr_summary.tsv", sep="\t", index=False)
    print("\nSummary:")
    print(summary.to_string(index=False))

    # Gene-level representatives: keep the best already-adjusted probe per gene.
    gene_key = "Gene" if "Gene" in work.columns else "probeID"
    gene_work = work.dropna(subset=["p_SMR_num"]).copy()
    gene_work = gene_work.sort_values(["qtl_type", "gwas_id", gene_key, "p_SMR_BH", "p_SMR_num"])
    gene_min = gene_work.groupby(["qtl_type", "gwas_id", gene_key], as_index=False).first()
    gene_min = gene_min.rename(columns={gene_key: "Gene", "p_SMR_num": "gene_min_p_SMR"})
    n_genes = (
        gene_min.groupby(["qtl_type", "gwas_id"], dropna=False)["Gene"]
        .transform("count")
    )
    gene_min["n_genes_tested"] = n_genes
    gene_cols = [
        "qtl_type", "gwas_id", "Gene", "probeID", "ProbeChr", "Probe_bp",
        "topSNP", "gene_min_p_SMR", "p_SMR_BH", "p_SMR_bonferroni",
        "smr_fdr_pass", "smr_bonferroni_pass", "gwas_threshold",
        "gwas_threshold_pass", "p_HEIDI_num", "heidi_pass_strict", "robust_sig",
        "robust_bonferroni", "n_smr_tests", "n_effective_snp", "n_genes_tested", "b_SMR", "se_SMR",
        "p_GWAS", "p_eQTL", "p_SMR_multi", "nsnp_HEIDI",
    ]
    gene_cols = [c for c in gene_cols if c in gene_min.columns]
    gene_min[gene_cols].sort_values(
        ["qtl_type", "gwas_id", "gene_min_p_SMR"]
    ).to_csv(results_dir / "smr_gene_min.tsv", sep="\t", index=False)
    gene_counts = (
        gene_min.groupby(["qtl_type", "gwas_id"])["Gene"]
        .size()
        .rename("n_genes_tested")
        .reset_index()
    )
    gene_sig = work[work["robust_sig"]].copy()
    gene_sig = gene_sig.sort_values(["qtl_type", "gwas_id", gene_key, "p_SMR_BH", "p_SMR_num"])
    gene_sig = gene_sig.groupby(["qtl_type", "gwas_id", gene_key], as_index=False).first()
    gene_sig = gene_sig.rename(columns={gene_key: "Gene", "p_SMR_num": "gene_min_p_SMR"})
    gene_sig = gene_sig.merge(gene_counts, on=["qtl_type", "gwas_id"], how="left")
    gene_sig[gene_cols].sort_values(
        ["qtl_type", "gwas_id", "p_SMR_BH", "gene_min_p_SMR"]
    ).to_csv(results_dir / "smr_gene_significant.tsv", sep="\t", index=False)

    gene_bonf = work[work["robust_bonferroni"]].copy()
    gene_bonf = gene_bonf.sort_values(["qtl_type", "gwas_id", gene_key, "p_SMR_bonferroni", "p_SMR_num"])
    gene_bonf = gene_bonf.groupby(["qtl_type", "gwas_id", gene_key], as_index=False).first()
    gene_bonf = gene_bonf.rename(columns={gene_key: "Gene", "p_SMR_num": "gene_min_p_SMR"})
    gene_bonf = gene_bonf.merge(gene_counts, on=["qtl_type", "gwas_id"], how="left")
    gene_bonf[gene_cols].sort_values(
        ["qtl_type", "gwas_id", "p_SMR_bonferroni", "gene_min_p_SMR"]
    ).to_csv(results_dir / "smr_gene_bonferroni_significant.tsv", sep="\t", index=False)

    gene_suggestive = work[work["suggestive_fdr_0_10"]].copy()
    gene_suggestive = gene_suggestive.sort_values(["qtl_type", "gwas_id", gene_key, "p_SMR_BH", "p_SMR_num"])
    gene_suggestive = gene_suggestive.groupby(["qtl_type", "gwas_id", gene_key], as_index=False).first()
    gene_suggestive = gene_suggestive.rename(columns={gene_key: "Gene", "p_SMR_num": "gene_min_p_SMR"})
    gene_suggestive = gene_suggestive.merge(gene_counts, on=["qtl_type", "gwas_id"], how="left")
    gene_suggestive[gene_cols].sort_values(
        ["qtl_type", "gwas_id", "p_SMR_BH", "gene_min_p_SMR"]
    ).to_csv(results_dir / "smr_gene_suggestive.tsv", sep="\t", index=False)

    prenatal_candidates = work[
        (work["qtl_type"] == "prenatal")
        & work["smr_fdr_pass"]
        & work["gwas_threshold_pass"]
    ].copy()
    prenatal_candidates = prenatal_candidates.sort_values(
        ["gwas_id", gene_key, "p_SMR_BH", "p_SMR_num"]
    ).groupby(["qtl_type", "gwas_id", gene_key], as_index=False).first()
    prenatal_candidates["candidate_status"] = np.select(
        [
            prenatal_candidates["heidi_pass_strict"],
            prenatal_candidates["p_HEIDI_num"].notna(),
        ],
        ["FDR + HEIDI pass", "FDR + HEIDI fail"],
        default="FDR + HEIDI unavailable",
    )
    status_order = {
        "FDR + HEIDI pass": 0,
        "FDR + HEIDI fail": 1,
        "FDR + HEIDI unavailable": 2,
    }
    prenatal_candidates["status_order"] = prenatal_candidates["candidate_status"].map(status_order)
    prenatal_candidates = prenatal_candidates.sort_values(
        ["status_order", "p_SMR_BH", "p_SMR_num", "gwas_id", gene_key]
    ).drop(columns="status_order")
    prenatal_candidates.insert(
        0, "prenatal_candidate_rank", range(1, len(prenatal_candidates) + 1)
    )
    prenatal_candidates.to_csv(
        results_dir / "smr_prenatal_candidates.tsv", sep="\t", index=False
    )

    gene_sig_counts = (
        gene_sig.groupby(["qtl_type", "gwas_id"])
        .size()
        .rename("n_gene_level_sig")
        .reset_index()
    )
    gene_bonf_counts = (
        gene_bonf.groupby(["qtl_type", "gwas_id"])
        .size()
        .rename("n_gene_bonferroni_sig")
        .reset_index()
    )
    gene_suggestive_counts = (
        gene_suggestive.groupby(["qtl_type", "gwas_id"])
        .size()
        .rename("n_gene_suggestive_fdr_0_10")
        .reset_index()
    )

    gene_summary = (
        gene_min.groupby(["qtl_type", "gwas_id"], dropna=False)
        .agg(
            n_genes_tested=("Gene", "size"),
            min_gene_p_smr=("gene_min_p_SMR", "min"),
            min_p_smr_bh=("p_SMR_BH", "min"),
            min_p_smr_bonferroni=("p_SMR_bonferroni", "min"),
        )
        .reset_index()
    )
    gene_summary = gene_summary.merge(
        gene_sig_counts, on=["qtl_type", "gwas_id"], how="left"
    ).merge(
        gene_bonf_counts, on=["qtl_type", "gwas_id"], how="left"
    ).merge(
        gene_suggestive_counts, on=["qtl_type", "gwas_id"], how="left"
    )
    gene_summary[["n_gene_level_sig", "n_gene_bonferroni_sig", "n_gene_suggestive_fdr_0_10"]] = gene_summary[
        ["n_gene_level_sig", "n_gene_bonferroni_sig", "n_gene_suggestive_fdr_0_10"]
    ].fillna(0).astype(int)
    gene_summary.to_csv(results_dir / "smr_gene_summary.tsv", sep="\t", index=False)
    print("\nGene-level summary:")
    print(gene_summary.to_string(index=False))
PYEOF

# R plot
if command -v Rscript &>/dev/null && [[ -f "${DATA_DIR}/smr_significant.tsv" ]]; then
    RESULTS_PATH="${DATA_DIR}" PLOT_PATH="${PLOT_DIR}" Rscript - <<'REOF'
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
})
results_dir <- Sys.getenv("RESULTS_PATH")
plot_dir <- Sys.getenv("PLOT_PATH")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

gene_sig_path <- file.path(results_dir, "smr_gene_significant.tsv")
if (!file.exists(gene_sig_path)) {
    cat("Gene-level significant table missing\n")
    quit()
}

gene_all <- fread(gene_sig_path)
if (nrow(gene_all) == 0 || !"p_SMR_BH" %in% names(gene_all)) {
    cat("No gene-level hits to plot\n")
    quit()
}

stage_levels <- c("prenatal", "postnatal", "interaction")
stage_labels <- c("Prenatal", "Postnatal", "Interaction")
gene <- copy(gene_all)
gene[, p_SMR_BH_num := as.numeric(p_SMR_BH)]
gene[, gene_min_p_SMR_num := as.numeric(gene_min_p_SMR)]
gene[, b_SMR_num := as.numeric(b_SMR)]
gene <- gene[is.finite(p_SMR_BH_num) & p_SMR_BH_num > 0]
gene[, stage_order := match(qtl_type, stage_levels)]
setorder(gene, stage_order, p_SMR_BH_num, gene_min_p_SMR_num)
gene[, ancestry := sub("_.*$", "", gwas_id)]
gene[, trait := sub("^[A-Z]+_", "", gwas_id)]
gene[, stage := factor(qtl_type, levels=stage_levels, labels=stage_labels)]
gene[, direction := fifelse(b_SMR_num >= 0,
                            "Positive",
                            "Negative")]
gene[, logq := -log10(p_SMR_BH_num)]
gene[, label := sprintf("%s | %s", Gene, trait)]
gene <- gene[order(stage_order, logq, gene_min_p_SMR_num, Gene)]
gene[, label := factor(label, levels=unique(label))]

fwrite(
    gene[, .(
        qtl_type, stage, gwas_id, ancestry, trait, Gene, probeID,
        gene_min_p_SMR, p_SMR_BH, p_SMR_bonferroni, p_HEIDI_num,
        b_SMR, se_SMR, logq, direction
    )],
    file.path(results_dir, "smr_gene_hits_source.tsv"),
    sep="\t"
)

pal <- c(EAS="#0072B2", EUR="#D55E00")
effect_shapes <- c("Positive"=24, "Negative"=25)

base_theme <- theme_classic(base_size=8.5, base_family="sans") +
    theme(
        plot.title=element_text(size=9.5, face="bold", hjust=0),
        plot.subtitle=element_blank(),
        axis.title=element_text(size=8.6),
        axis.text=element_text(size=8.0, color="black"),
        legend.title=element_text(size=8.2),
        legend.text=element_text(size=7.8),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.box="vertical",
        legend.key.height=unit(0.34, "cm"),
        legend.key.width=unit(0.36, "cm"),
        panel.grid.major.x=element_line(color="grey88", linewidth=0.25),
        panel.grid.major.y=element_line(color="grey94", linewidth=0.25),
        panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=8.4, face="bold"),
        strip.text.y.left=element_text(angle=0),
        strip.placement="outside",
        plot.margin=margin(7, 6, 10, 6)
    )

p_gene <- ggplot(gene, aes(x=logq, y=label)) +
    geom_vline(xintercept=-log10(0.05), linetype="22", color="grey55", linewidth=0.35) +
    geom_segment(aes(x=-log10(0.05), xend=logq, y=label, yend=label),
                 color="grey83", linewidth=0.35) +
    geom_point(aes(color=ancestry, shape=direction),
               size=2.8, alpha=0.95) +
    scale_color_manual(values=pal, name="Ancestry") +
    scale_shape_manual(values=effect_shapes, name=expression(beta[SMR])) +
    scale_x_continuous(expand=expansion(mult=c(0.01, 0.08))) +
    labs(
        title="Within-analysis FDR SMR genes",
        x=expression(-log[10](q[SMR])),
        y=NULL
    ) +
    guides(
        color=guide_legend(order=1, nrow=1, override.aes=list(size=3.0)),
        shape=guide_legend(order=2, nrow=1, override.aes=list(size=3.0, color="grey20"))
    ) +
    base_theme +
    facet_grid(stage ~ ., scales="free_y", space="free_y", switch="y") +
    theme(
        axis.text.y=element_text(size=8.0)
    )

figure_height <- max(4.8, 2.2 + 0.28 * nrow(gene))
summary_pdf <- file.path(plot_dir, "smr_gene_hits.pdf")
summary_png <- file.path(plot_dir, "smr_gene_hits.png")
ggsave(summary_pdf, p_gene, width=6.2, height=figure_height, device=cairo_pdf)
ggsave(summary_png, p_gene, width=6.2, height=figure_height, dpi=320, bg="white")

main_summary <- gene[, .(
    n_risk_genes=uniqueN(Gene),
    genes=paste(sort(unique(Gene)), collapse=", "),
    min_p_SMR_BH=min(p_SMR_BH_num, na.rm=TRUE),
    max_logq=max(logq, na.rm=TRUE),
    min_gene_p_SMR=min(gene_min_p_SMR_num, na.rm=TRUE),
    n_positive=sum(direction == "Positive", na.rm=TRUE),
    n_negative=sum(direction == "Negative", na.rm=TRUE)
), by=.(qtl_type, stage, ancestry, trait)]
main_summary[, stage_order := match(qtl_type, stage_levels)]
trait_levels <- main_summary[, .(
    strongest=max(max_logq, na.rm=TRUE),
    total_genes=sum(n_risk_genes)
), by=trait][order(strongest, total_genes, trait), trait]
main_summary[, trait := factor(trait, levels=trait_levels)]
main_summary[, stage := factor(as.character(stage), levels=stage_labels)]
main_summary[, ancestry := factor(ancestry, levels=names(pal))]
setorder(main_summary, stage_order, trait, ancestry)
fwrite(
    main_summary[, .(
        qtl_type, stage, ancestry, trait, n_risk_genes, genes,
        min_p_SMR_BH, max_logq, min_gene_p_SMR, n_positive, n_negative
    )],
    file.path(results_dir, "smr_main_figure_source.tsv"),
    sep="\t"
)

size_breaks <- sort(unique(main_summary$n_risk_genes))
p_main <- ggplot(main_summary, aes(x=stage, y=trait, group=ancestry)) +
    geom_point(aes(color=ancestry, size=n_risk_genes),
               alpha=0.94, position=position_dodge(width=0.55)) +
    scale_color_manual(values=pal, name="Ancestry", drop=FALSE) +
    scale_size_continuous(name="Risk genes", range=c(3.2, 7.2), breaks=size_breaks) +
    scale_x_discrete(drop=FALSE) +
    labs(
        title="SMR risk genes collapsed by disease",
        x=NULL,
        y=NULL
    ) +
    guides(
        color=guide_legend(order=1, nrow=1, override.aes=list(size=3.4)),
        size=guide_legend(order=2, nrow=1)
    ) +
    base_theme +
    theme(
        panel.grid.major.x=element_line(color="grey92", linewidth=0.25),
        panel.grid.major.y=element_line(color="grey90", linewidth=0.25),
        axis.text.x=element_text(size=8.6),
        axis.text.y=element_text(size=8.2),
        legend.box="horizontal"
    )

main_height <- max(3.6, 1.8 + 0.40 * uniqueN(main_summary$trait))
main_pdf <- file.path(plot_dir, "smr_main_figure.pdf")
main_png <- file.path(plot_dir, "smr_main_figure.png")
ggsave(main_pdf, p_main, width=5.6, height=main_height, device=cairo_pdf)
ggsave(main_png, p_main, width=5.6, height=main_height, dpi=320, bg="white")

n_labels <- uniqueN(gene$label)
max_label_chars <- max(nchar(as.character(gene$label)))
n_legend_items <- uniqueN(gene$ancestry) + uniqueN(gene$direction)
n_facets <- uniqueN(gene$stage)
layout_status <- if (n_labels <= 30 && max_label_chars <= 28 && n_legend_items <= 6 && n_facets <= 3) "PASS" else "REVIEW"
layout_validation <- data.table(
    metric=c("visible_labels", "max_label_characters", "legend_items", "stage_facets", "overall"),
    value=c(n_labels, max_label_chars, n_legend_items, n_facets, layout_status),
    status=c(
        ifelse(n_labels <= 30, "PASS", "REVIEW"),
        ifelse(max_label_chars <= 28, "PASS", "REVIEW"),
        ifelse(n_legend_items <= 6, "PASS", "REVIEW"),
        ifelse(n_facets <= 3, "PASS", "REVIEW"),
        layout_status
    )
)
fwrite(layout_validation, file.path(results_dir, "smr_gene_hits_layout_validation.tsv"), sep="\t")

main_labels <- uniqueN(main_summary$trait)
main_max_trait_chars <- max(nchar(as.character(main_summary$trait)))
main_legend_items <- uniqueN(main_summary$ancestry) + length(size_breaks)
main_layout_status <- if (main_labels <= 18 && main_max_trait_chars <= 18 && main_legend_items <= 8) "PASS" else "REVIEW"
main_layout_validation <- data.table(
    metric=c("collapsed_points", "visible_traits", "max_trait_characters", "legend_items", "max_genes_per_point", "overall"),
    value=c(nrow(main_summary), main_labels, main_max_trait_chars, main_legend_items, max(main_summary$n_risk_genes), main_layout_status),
    status=c(
        ifelse(nrow(main_summary) <= 30, "PASS", "REVIEW"),
        ifelse(main_labels <= 18, "PASS", "REVIEW"),
        ifelse(main_max_trait_chars <= 18, "PASS", "REVIEW"),
        ifelse(main_legend_items <= 8, "PASS", "REVIEW"),
        "PASS",
        main_layout_status
    )
)
fwrite(main_layout_validation, file.path(results_dir, "smr_main_figure_layout_validation.tsv"), sep="\t")

main_index <- data.table(
    rank=1L,
    figure_type="collapsed_trait_stage_gene_count",
    pdf=main_pdf,
    png=main_png,
    source=file.path(results_dir, "smr_main_figure_source.tsv"),
    detailed_gene_source=file.path(results_dir, "smr_gene_hits_source.tsv"),
    n_gene_level_sig=nrow(gene),
    n_collapsed_points=nrow(main_summary),
    n_traits=uniqueN(main_summary$trait),
    n_prenatal_sig=sum(gene$qtl_type == "prenatal"),
    n_postnatal_sig=sum(gene$qtl_type == "postnatal"),
    n_interaction_sig=sum(gene$qtl_type == "interaction")
)
fwrite(main_index, file.path(results_dir, "smr_main_figure_index.tsv"), sep="\t")

cat("Plot:", file.path(plot_dir, "smr_gene_hits.png"), "\n")
REOF
fi

echo "[$(date)] SMR pipeline complete."
echo "  FDR hits          : ${DATA_DIR}/smr_significant.tsv"
echo "  Bonferroni hits   : ${DATA_DIR}/smr_bonferroni_significant.tsv"
echo "  Suggestive q<0.10 : ${DATA_DIR}/smr_suggestive.tsv"
echo "  Prenatal ranking  : ${DATA_DIR}/smr_prenatal_candidates.tsv"
echo "  Summary table     : ${DATA_DIR}/smr_summary.tsv"
echo "  Plot              : ${PLOT_DIR}/smr_gene_hits.png"

