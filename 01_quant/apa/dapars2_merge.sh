#!/usr/bin/env bash
# Merge per-chromosome DaPars2 results (chr1-22) into one matrix; strip '_PDUI' from sample headers.
set -euo pipefail

WORKDIR="${WORKDIR:?}"
FINAL_OUTPUT="${FINAL_OUTPUT:-${WORKDIR}/all.apa.count}"
rm -f "${FINAL_OUTPUT}"

first=1
for i in {1..22}; do
  chr_file="${WORKDIR}/Dapars2_out_chr${i}/Dapars2_result_temp.chr${i}.txt"
  if [[ "${first}" -eq 1 ]]; then
    read -r gene fit_value predicted loci rest < <(head -n 1 "${chr_file}")
    header="Gene\tfit_value\tPredicted_Proximal_APA\tLoci"
    for sample_path in ${rest}; do
      sample_name=$(basename "${sample_path}")
      header="${header}\t${sample_name%_PDUI}"
    done
    echo -e "${header}" > "${FINAL_OUTPUT}"
    first=0
  fi
  tail -n +2 "${chr_file}" >> "${FINAL_OUTPUT}"
done
echo "Merged DaPars2 results: ${FINAL_OUTPUT}"