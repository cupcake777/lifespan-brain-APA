#!/usr/bin/env bash
# Extract mash inputs: shared phenotypes, strong keys and nominal pairs per QTL class and stage,
# plus the APA isoform-switch table. Writes data/04_qtl/mash/ and data/04_qtl/05_apa_isoform_switch.tsv.
#SBATCH -J apa_eqtl_mash
#SBATCH -p cu
#SBATCH -c 16
#SBATCH --mem=96G
set -euo pipefail
: "${APA_ROOT:?source config/paths.sh}"
: "${APA_HPC_HOME:?source config/paths.sh}"

D="$APA_ROOT/data/04_qtl"
OUT="$D/mash"
APA_DIR="$APA_HPC_HOME/project/work/apa/01_QTL"
EQTL_DIR="$D/eqtl"
STAGES=(prenatal postnatal)
QTL_CLASSES=(apa eqtl)
declare -A SAMPLE_MOD=([apa]=211 [eqtl]=307)

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  mkdir -p "$APA_ROOT/logs"
  sbatch -o "$APA_ROOT/logs/05_mashr_%j.out" -e "$APA_ROOT/logs/05_mashr_%j.err" "$0"
  exit 0
fi
mkdir -p "$OUT"

# --- phenotype universe ---
for qtl in "${QTL_CLASSES[@]}"; do
  if [[ "$qtl" == apa ]]; then perm="$D"; else perm="$EQTL_DIR/permutation"; fi
  out="$OUT/${qtl}.shared_phenotypes.txt"
  if [[ ! -s "$out" ]]; then
    tmp_pre="$(mktemp)"; tmp_post="$(mktemp)"
    gzip -cd "$perm/prenatal.cis_qtl.txt.gz" | awk 'NR > 1 {print $1}' | sort -u > "$tmp_pre"
    gzip -cd "$perm/postnatal.cis_qtl.txt.gz" | awk 'NR > 1 {print $1}' | sort -u > "$tmp_post"
    comm -12 "$tmp_pre" "$tmp_post" > "$out"
    rm -f "$tmp_pre" "$tmp_post"
  fi
done

# --- APA isoform switch ---
isoform_out="$D/05_apa_isoform_switch.tsv"
if [[ ! -s "$isoform_out" ]]; then
  tmp_pre="$(mktemp)"; tmp_post="$(mktemp)"
  gzip -cd "$APA_DIR/permutation/prenatal.cis_qtl.txt.gz" |
    awk -v OFS='\t' 'NR > 1 {split($1, f, "|"); print f[2], $1}' |
    sort -u -k1,1 > "$tmp_pre"
  gzip -cd "$APA_DIR/permutation/postnatal.cis_qtl.txt.gz" |
    awk -v OFS='\t' 'NR > 1 {split($1, f, "|"); print f[2], $1}' |
    sort -u -k1,1 > "$tmp_post"
  {
    printf 'gene\tprenatal_phenotype_id\tpostnatal_phenotype_id\n'
    join -t $'\t' -j 1 "$tmp_pre" "$tmp_post" |
      awk -F'\t' -v OFS='\t' '$2 != $3 {print $1, $2, $3}'
  } > "$isoform_out"
  rm -f "$tmp_pre" "$tmp_post"
fi

# --- strong keys ---
for qtl in "${QTL_CLASSES[@]}"; do
  out="$OUT/${qtl}.strong_keys.txt"
  if [[ ! -s "$out" ]]; then
    if [[ "$qtl" == apa ]]; then perm="$D"; else perm="$EQTL_DIR/permutation"; fi
    {
      for stage in "${STAGES[@]}"; do
        gzip -cd "$perm/${stage}.cis_qtl.txt.gz" |
          awk -v OFS='\t' -v shared="$OUT/${qtl}.shared_phenotypes.txt" '
            BEGIN {while ((getline line < shared) > 0) keep[line] = 1}
            NR > 1 && ($1 in keep) {print $1, $7}'
      done
    } | sort -u > "$out"
  fi
done

# --- nominal pairs: apa prenatal is parquet chunks; everything else is a gzipped table ---
pids=()
labels=()
for qtl in "${QTL_CLASSES[@]}"; do
  for stage in "${STAGES[@]}"; do
    out="$OUT/${qtl}.${stage}.pairs.tsv"
    if [[ ! -s "$out" ]]; then
      (
        if [[ "$qtl" == apa ]]; then nominal="$APA_DIR/nominal"; else nominal="$EQTL_DIR/nominal"; fi
        if [[ "$qtl" == apa && "$stage" == prenatal ]]; then
          PARQUET_DIR="$nominal" PARQUET_STAGE="$stage" python3 - <<'PY'
import glob
import os
import sys

import pyarrow.parquet as pq

pattern = os.path.join(os.environ["PARQUET_DIR"], "%s.cis_qtl_pairs.*.parquet" % os.environ["PARQUET_STAGE"])
cols = ["phenotype_id", "variant_id", "start_distance", "end_distance",
        "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se"]
out = sys.stdout
for path in sorted(glob.glob(pattern)):
    pf = pq.ParquetFile(path)
    for batch in pf.iter_batches(batch_size=500000, columns=cols):
        d = batch.to_pydict()
        ph, vid = d["phenotype_id"], d["variant_id"]
        sd, ed, af = d["start_distance"], d["end_distance"], d["af"]
        ms, mc, pv = d["ma_samples"], d["ma_count"], d["pval_nominal"]
        sl, se = d["slope"], d["slope_se"]
        for i in range(batch.num_rows):
            if sl[i] is None or se[i] is None or pv[i] is None or af[i] is None:
                continue
            out.write("%s\t%s\t%d\t%d\t%g\t%d\t%d\t%g\t%g\t%g\n" % (
                ph[i], vid[i], sd[i], ed[i], af[i], ms[i], mc[i], pv[i], sl[i], se[i]))
PY
        else
          gzip -cd "$nominal/${stage}.cis_qtl_pairs.txt.gz"
        fi |
        awk -v OFS='\t' \
            -v shared="$OUT/${qtl}.shared_phenotypes.txt" \
            -v strong="$OUT/${qtl}.strong_keys.txt" \
            -v mod="${SAMPLE_MOD[$qtl]}" '
          BEGIN {
            while ((getline line < shared) > 0) keep[line] = 1
            while ((getline line < strong) > 0) {
              split(line, a, "\t"); strongkey[a[1] SUBSEP a[2]] = 1
            }
            print "phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal", "set"
          }
          $1 == "phenotype_id" {next}
          !($1 in keep) {next}
          {
            is_strong = (($1 SUBSEP $2) in strongkey)
            if (!is_strong) {
              split($2, v, ":")
              if (v[2] % mod != 0) next
            }
            print $1, $2, $9, $10, $8, (is_strong ? "strong" : "random")
          }' > "$out.tmp"
        mv -f "$out.tmp" "$out"
      ) &
      pids+=("$!")
      labels+=("${qtl} ${stage}")
    fi
  done
done
n_failed=0
for i in "${!pids[@]}"; do
  if ! wait "${pids[$i]}"; then
    echo "ERROR: extraction failed: ${labels[$i]}" >&2
    n_failed=$((n_failed + 1))
  fi
done
(( n_failed == 0 )) || exit 1

for qtl in "${QTL_CLASSES[@]}"; do
  for stage in "${STAGES[@]}"; do
    f="$OUT/${qtl}.${stage}.pairs.tsv"
    [[ -s "$f" ]] || { echo "ERROR: empty extraction $f" >&2; exit 1; }
  done
done
