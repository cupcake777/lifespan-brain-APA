export type ReportRow = {
  gwas_id: string;
  qtl_id: string;
  locus: string;
  phenotype: string;
  pp4: number;
  n_snps: number;
  plot_file: string;
  locuscompare_plot_file: string;
  locuscompare_json_file: string;
  source_file: string;
};

export type ReportMeta = {
  project_name: string;
  results_dir: string;
  generated_at: string;
  report_version: string;
};

export type ReportSummary = {
  total_tests: number;
  significant_pp4_08: number;
  significant_pp4_07: number;
  mean_pp4: number;
  unique_genes: number;
  unique_loci: number;
  [key: string]: number | string | null | undefined;
};

export type ReportAsset = {
  name: string;
  rel_path: string;
};

export type RuntimeHeartbeat = Record<string, unknown> | null;

export type ReportPayload = {
  meta: ReportMeta;
  summary: ReportSummary;
  results: ReportRow[];
  assets: {
    plots: ReportAsset[];
    plot_data: ReportAsset[];
    downloads: ReportAsset[];
  };
  runtime: {
    heartbeat?: RuntimeHeartbeat;
    monitor_snapshot?: Record<string, unknown> | null;
    [key: string]: unknown;
  };
  warnings: string[];
};
