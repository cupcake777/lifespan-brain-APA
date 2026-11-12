#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"

tmp_dir="$(mktemp -d)"
trap 'if [[ -n "${server_pid:-}" ]]; then kill "${server_pid}" 2>/dev/null || true; fi; rm -rf "${tmp_dir}"' EXIT

mkdir -p "${tmp_dir}/results/abf" "${tmp_dir}/results/plots" "${tmp_dir}/results/plot_data"
cat > "${tmp_dir}/results/abf/GWAS_A_rs100_locus_results.csv" <<'CSV'
GWAS_ID,QTL_ID,Locus,Phenotype,PP4,n_snps,plot_file,locuscompare_plot_file,locuscompare_json_file
GWAS_A,postnatal,rs100,GENE1,0.92,41,plots/GWAS_A_postnatal_rs100_GENE1_coloc.png,plots/GWAS_A_postnatal_rs100_GENE1_locuscompare.png,plot_data/GWAS_A_postnatal_rs100_GENE1_locuscompare.json
CSV
printf 'plot\n' >"${tmp_dir}/results/plots/GWAS_A_postnatal_rs100_GENE1_coloc.png"
printf 'locuscompare\n' >"${tmp_dir}/results/plots/GWAS_A_postnatal_rs100_GENE1_locuscompare.png"
printf 'pdf\n' >"${tmp_dir}/results/plots/GWAS_A_postnatal_rs100_GENE1_locuscompare.pdf"
printf '{"meta":{"lead_snp":"rs100"},"points":[{"snp":"rs100"}]}\n' >"${tmp_dir}/results/plot_data/GWAS_A_postnatal_rs100_GENE1_locuscompare.json"

./easycoloc report "${tmp_dir}/results" --host 127.0.0.1 --port 0 --no-open >"${tmp_dir}/report.log" 2>&1 &
server_pid=$!

server_url=""
for _ in $(seq 1 15); do
  if ! kill -0 "${server_pid}" 2>/dev/null; then
    cat "${tmp_dir}/report.log" >&2 || true
    echo "[SMOKE] report-web server exited early" >&2
    exit 1
  fi

  if grep '\[REPORT-WEB\] url:' "${tmp_dir}/report.log" >/dev/null 2>&1; then
    server_url="$(sed -n 's/.*\[REPORT-WEB\] url: \(http:\/\/[^[:space:]]*\).*/\1/p' "${tmp_dir}/report.log" | tail -n 1)"
    if [[ -n "${server_url}" ]]; then
      break
    fi
  fi
  sleep 1
done

if [[ -z "${server_url}" ]]; then
  cat "${tmp_dir}/report.log" >&2 || true
  echo "[SMOKE] failed to capture server url" >&2
  exit 1
fi

for _ in $(seq 1 15); do
  if ! kill -0 "${server_pid}" 2>/dev/null; then
    cat "${tmp_dir}/report.log" >&2 || true
    echo "[SMOKE] report-web server exited before API became reachable" >&2
    exit 1
  fi
  if curl -fsS "${server_url}/api/report-data" >/dev/null 2>&1; then
    break
  fi
  sleep 1
done

curl -fsS "${server_url}/api/report-data" | grep '"project_name"' >/dev/null
curl -fsS "${server_url}/api/report-data" | grep '"locuscompare_plot_file"' >/dev/null
curl -fsS "${server_url}/plot_data/GWAS_A_postnatal_rs100_GENE1_locuscompare.json" | grep '"lead_snp"' >/dev/null
curl -fsSI "${server_url}/plots/GWAS_A_postnatal_rs100_GENE1_locuscompare.pdf" | grep -i '^Content-Type: application/pdf' >/dev/null
curl -fsS "${server_url}/does-not-exist" | grep -i "<!doctype html>" >/dev/null
grep '\[REPORT-WEB\] url:' "${tmp_dir}/report.log" >/dev/null

echo "[SMOKE] report-web cli smoke test passed"
