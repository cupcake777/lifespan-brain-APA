#!/usr/bin/env Rscript

files_to_check <- c(
  "run_coloc.R",
  "tools/rerun_plots.R",
  "src/utils_config.R",
  "src/utils_download.R",
  "src/utils_reference.R",
  "src/utils_bootstrap.R",
  "src/utils_output.R",
  "src/utils_runtime.R",
  "src/utils_helpers.R",
  "src/utils_analysis.R",
  "src/utils_report.R",
  "src/utils_sensitivity.R",
  "src/utils_format.R",
  "src/utils_plot.R",
  "src/utils_plot_locuscompare.R",
  "tools/doctor_easycoloc.R",
  "tools/bootstrap_references.R",
  "tools/list_reference_requirements.R",
  "tools/monitor_easycoloc.R",
  "tools/check_run_completion.R",
  "tools/build_output_manifest.R",
  "tools/build_harmony_qc_report.R",
  "tools/clean_intermediates.R",
  "tools/finalize_run.R",
  "tools/harmonize_gwas.R",
  "tools/summarize_run_status.R",
  "tools/checks/smoke_project_fixture.R",
  "tools/checks/smoke_test_plotting.R",
  "tools/checks/smoke_test_bootstrap_refs.R",
  "tools/checks/smoke_test_bootstrap_demo.R",
  "tools/checks/smoke_test_gtex_bootstrap.R",
  "tools/checks/smoke_test_1kg_setup.R",
  "tools/checks/smoke_test_ld_cache.R",
  "tools/checks/smoke_test_ld_timeout.R",
  "tools/checks/smoke_test_chunked_harmonization.R",
  "tools/checks/smoke_test_harmonization_fallback.R",
  "tools/checks/smoke_test_doctor.R",
  "tools/checks/smoke_test_refs.R",
  "tools/checks/smoke_test_runtime_monitor.R",
  "tools/checks/smoke_test_shard_runtime.R",
  "tools/checks/smoke_test_output_dir_tools.R",
  "examples/synthetic_plot_data.R",
  "examples/synthetic_plot_data_v2.R",
  "examples/synthetic_plot_demo.R"
)

for (file in files_to_check) {
  parse(file = file)
}

cat("[PARSE] checked", length(files_to_check), "files\n")
