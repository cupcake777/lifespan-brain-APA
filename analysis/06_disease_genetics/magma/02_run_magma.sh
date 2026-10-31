#!/usr/bin/env bash
# MAGMA gene analysis per GWAS and APA-QTL gene-set enrichment, with within-GWAS and global BH FDR.
: "${APA_ROOT:?source config/paths.sh}"

set -euo pipefail

# --- config ---

project_dir="${APA_ROOT}/data/06_disease/magma_run"
manifest="${APA_ROOT}/data/06_disease/gwas_manifest.tsv"
gene_loc="${APA_ROOT}/ref/magma/NCBI38.gene.loc"

eas_plink_dir="${APA_ROOT}/ref/sLDSC/EAS_hg38/plink_files"
eur_plink_dir="${APA_ROOT}/ref/sLDSC/EUR_hg38/plink_files"

gene_sets="${project_dir}/gene_sets/apa_qtl_gene_sets.sets"
enrichment_dir="${project_dir}/enrichment"
logs_dir="${APA_ROOT}/logs"

# --- directories ---

mkdir -p "${project_dir}/work/EAS" "${project_dir}/work/EUR"
mkdir -p "${project_dir}/annotations/EAS" "${project_dir}/annotations/EUR"
mkdir -p "${project_dir}/gene_sets" "${enrichment_dir}" "${logs_dir}"

# --- merge reference PLINK files and baseline annotations ---

for population in EAS EUR; do
  if [[ "${population}" == "EAS" ]]; then
    plink_dir="${eas_plink_dir}"
  else
    plink_dir="${eur_plink_dir}"
  fi

  work_dir="${project_dir}/work/${population}"
  bfile_prefix="${work_dir}/${population}_hg38_genome"
  merge_list="${work_dir}/plink_merge_list.txt"
  baseline_annot="${project_dir}/annotations/${population}/baseline.genes.annot"

  if [[ ! -s "${bfile_prefix}.bed" || ! -s "${bfile_prefix}.bim" || ! -s "${bfile_prefix}.fam" ]]; then
    : > "${merge_list}"
    for chr in $(seq 2 22); do
      printf "%s\n" "${plink_dir}/1000G.${population}.hg38.${chr}" >> "${merge_list}"
    done

    plink \
      --bfile "${plink_dir}/1000G.${population}.hg38.1" \
      --merge-list "${merge_list}" \
      --make-bed \
      --out "${bfile_prefix}" \
      > "${logs_dir}/plink_merge_${population}.stdout.log" 2>&1
  fi

  if [[ ! -s "${baseline_annot}" ]]; then
    magma \
      --annotate \
      --snp-loc "${bfile_prefix}.bim" \
      --gene-loc "${gene_loc}" \
      --out "${project_dir}/annotations/${population}/baseline" \
      > "${logs_dir}/magma_annotate_baseline_${population}.stdout.log" 2>&1
  fi
done

# --- build APA-QTL gene sets ---

Rscript "${project_dir}/build_qtl_gene_sets.R" \
  > "${logs_dir}/build_qtl_gene_sets.stdout.log" 2>&1

# --- run genome-wide baseline gene analysis ---

tail -n +2 "${manifest}" | while IFS=$'\t' read -r gwas_id population gwas_path; do
  bfile_prefix="${project_dir}/work/${population}/${population}_hg38_genome"
  baseline_annot="${project_dir}/annotations/${population}/baseline.genes.annot"
  out_dir="${enrichment_dir}/${gwas_id}"
  baseline_prefix="${out_dir}/baseline"

  mkdir -p "${out_dir}"

  if [[ ! -s "${baseline_prefix}.genes.raw" ]]; then
    magma \
      --bfile "${bfile_prefix}" \
      --pval "${gwas_path}" use=SNPID,P ncol=N \
      --gene-annot "${baseline_annot}" \
      --out "${baseline_prefix}" \
      > "${logs_dir}/magma_${gwas_id}_baseline.stdout.log" 2>&1
  fi
done

# --- run APA-QTL gene-set enrichment ---

tail -n +2 "${manifest}" | while IFS=$'\t' read -r gwas_id population gwas_path; do
  out_dir="${enrichment_dir}/${gwas_id}"
  baseline_prefix="${out_dir}/baseline"
  gsa_prefix="${out_dir}/apa_qtl_sets"

  if [[ ! -s "${gsa_prefix}.gsa.out" ]]; then
    magma \
      --gene-results "${baseline_prefix}.genes.raw" \
      --set-annot "${gene_sets}" \
      --out "${gsa_prefix}" \
      > "${logs_dir}/magma_${gwas_id}_apa_qtl_sets.stdout.log" 2>&1
  fi
done

# --- summarize enrichment results ---

summary="${enrichment_dir}/magma_gsa_summary.tsv"
printf "gwas_id\tpopulation\tset_name\ttype\tngenes\tbeta\tbeta_std\tse\tp\tfdr_within_gwas\tfdr_global\tsignificant_within_gwas\tsource_gsa\n" > "${summary}"

tail -n +2 "${manifest}" | while IFS=$'\t' read -r gwas_id population gwas_path; do
  gsa_out="${enrichment_dir}/${gwas_id}/apa_qtl_sets.gsa.out"
  awk -v gwas_id="${gwas_id}" -v population="${population}" -v source_gsa="${gsa_out}" '
    /^VARIABLE[[:space:]]+/ {
      for (i = 1; i <= NF; i++) col[$i] = i
      next
    }
    NF > 0 && $1 !~ /^#/ && ("VARIABLE" in col) && ($col["TYPE"] == "SET") {
      print gwas_id "\t" population "\t" $col["VARIABLE"] "\t" $col["TYPE"] "\t" \
            $col["NGENES"] "\t" $col["BETA"] "\t" $col["BETA_STD"] "\t" \
            $col["SE"] "\t" $col["P"] "\tNA\tNA\tNA\t" source_gsa
    }
  ' "${gsa_out}" >> "${summary}"
done

Rscript -e '
library(data.table)
summary_path <- "${APA_ROOT}/data/06_disease/magma_run/enrichment/magma_gsa_summary.tsv"
x <- fread(summary_path)
x[, p := as.numeric(p)]
x[, fdr_within_gwas := as.numeric(fdr_within_gwas)]
x[, fdr_global := as.numeric(fdr_global)]
x[, fdr_within_gwas := p.adjust(p, method = "BH"), by = gwas_id]
x[, fdr_global := p.adjust(p, method = "BH")]
x[, significant_within_gwas := fdr_within_gwas < 0.05]
fwrite(x, summary_path, sep = "\t", na = "NA")
' > "${logs_dir}/summarize_magma_gsa.stdout.log" 2>&1

# --- summarize warnings ---

warning_register="${enrichment_dir}/magma_warning_register.tsv"
printf "source\twarning\n" > "${warning_register}"

shopt -s nullglob
warning_sources=(
  "${enrichment_dir}"/*/*.log
  "${logs_dir}"/magma_*_baseline.stdout.log
  "${logs_dir}"/magma_*_apa_qtl_sets.stdout.log
  "${logs_dir}"/build_qtl_gene_sets.stdout.log
  "${logs_dir}"/summarize_magma_gsa.stdout.log
)

if [[ "${#warning_sources[@]}" -gt 0 ]]; then
  printf "%s\0" "${warning_sources[@]}" |
    xargs -0 awk '
      /WARNING/ {
        warning = $0
        sub(/^[[:space:]]+/, "", warning)
        print FILENAME "\t" warning
      }
    ' >> "${warning_register}"
fi
