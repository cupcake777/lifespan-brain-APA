import type { ReportAsset, ReportPayload, ReportRow, RuntimeHeartbeat } from "../types";

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null && !Array.isArray(value);
}

function asString(value: unknown, fallback = ""): string {
  return typeof value === "string" ? value : fallback;
}

function asNumber(value: unknown, fallback = 0): number {
  if (typeof value === "number" && Number.isFinite(value)) {
    return value;
  }

  if (typeof value === "string" && value.trim().length > 0) {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : fallback;
  }

  return fallback;
}

function normalizeAssetList(value: unknown): ReportAsset[] {
  if (!Array.isArray(value)) {
    return [];
  }

  return value
    .filter(isRecord)
    .map((item) => ({
      name: asString(item.name),
      rel_path: asString(item.rel_path)
    }));
}

function normalizeReportRow(value: unknown): ReportRow | null {
  if (!isRecord(value)) {
    return null;
  }

  return {
    gwas_id: asString(value.gwas_id),
    qtl_id: asString(value.qtl_id),
    locus: asString(value.locus),
    phenotype: asString(value.phenotype),
    pp4: asNumber(value.pp4),
    n_snps: asNumber(value.n_snps),
    plot_file: asString(value.plot_file),
    locuscompare_plot_file: asString(value.locuscompare_plot_file),
    locuscompare_json_file: asString(value.locuscompare_json_file),
    source_file: asString(value.source_file)
  };
}

function appendIssue(issues: string[], condition: boolean, message: string): void {
  if (condition) {
    issues.push(message);
  }
}

export function normalizeReportPayload(raw: unknown): ReportPayload {
  const issues: string[] = [];
  const root = isRecord(raw) ? raw : {};

  appendIssue(issues, !isRecord(raw), "Report payload was not a JSON object.");

  const meta = isRecord(root.meta) ? root.meta : {};
  appendIssue(issues, !isRecord(root.meta), "meta was missing or invalid.");

  const summary = isRecord(root.summary) ? root.summary : {};
  appendIssue(issues, !isRecord(root.summary), "summary was missing or invalid.");

  const results = Array.isArray(root.results)
    ? root.results.map(normalizeReportRow).filter((row): row is ReportRow => row !== null)
    : [];
  appendIssue(issues, !Array.isArray(root.results), "results was missing or invalid.");

  const assetsRoot = isRecord(root.assets) ? root.assets : {};
  appendIssue(issues, !isRecord(root.assets), "assets was missing or invalid.");

  const runtimeRoot = isRecord(root.runtime) ? root.runtime : {};
  appendIssue(issues, !isRecord(root.runtime), "runtime was missing or invalid.");

  const runtime: ReportPayload["runtime"] = {
    heartbeat: isRecord(runtimeRoot.heartbeat) ? runtimeRoot.heartbeat : null,
    monitor_snapshot: isRecord(runtimeRoot.monitor_snapshot) ? runtimeRoot.monitor_snapshot : null
  };

  const warnings = Array.isArray(root.warnings)
    ? root.warnings.filter((warning): warning is string => typeof warning === "string")
    : [];
  appendIssue(issues, !Array.isArray(root.warnings), "warnings was missing or invalid.");

  return {
    meta: {
      project_name: asString(meta.project_name, "EasyColoc Report"),
      results_dir: asString(meta.results_dir),
      generated_at: asString(meta.generated_at),
      report_version: asString(meta.report_version)
    },
    summary: {
      total_tests: asNumber(summary.total_tests),
      significant_pp4_08: asNumber(summary.significant_pp4_08),
      significant_pp4_07: asNumber(summary.significant_pp4_07),
      mean_pp4: asNumber(summary.mean_pp4),
      unique_genes: asNumber(summary.unique_genes),
      unique_loci: asNumber(summary.unique_loci)
    },
    results,
    assets: {
      plots: normalizeAssetList(assetsRoot.plots),
      plot_data: normalizeAssetList(assetsRoot.plot_data),
      downloads: normalizeAssetList(assetsRoot.downloads)
    },
    runtime,
    warnings: [...warnings, ...issues]
  };
}

export function filterResults(
  rows: ReportRow[],
  filters: { query: string; minPp4: number }
): ReportRow[] {
  const needle = filters.query.trim().toLowerCase();

  return rows.filter((row) => {
    const matchesQuery =
      needle.length === 0 ||
      [row.gwas_id, row.qtl_id, row.locus, row.phenotype, row.source_file].some((value) =>
        String(value ?? "").toLowerCase().includes(needle)
      );

    return matchesQuery && row.pp4 >= filters.minPp4;
  });
}

export function topHits(rows: ReportRow[], limit: number): ReportRow[] {
  if (limit <= 0) {
    return [];
  }

  return [...rows].sort((left, right) => right.pp4 - left.pp4).slice(0, limit);
}

export function rowKey(row: ReportRow): string {
  return [
    row.gwas_id,
    row.qtl_id,
    row.locus,
    row.phenotype,
    row.source_file,
    row.pp4,
    row.n_snps
  ].join("::");
}

export function formatPp4(value: number): string {
  return Number.isFinite(value) ? value.toFixed(2) : "NA";
}

export function formatHeartbeatStage(heartbeat: RuntimeHeartbeat | undefined): string | null {
  if (!heartbeat || typeof heartbeat !== "object") {
    return null;
  }

  const stage = heartbeat.stage;
  if (typeof stage === "string" && stage.trim().length > 0) {
    return stage;
  }
  if (Array.isArray(stage)) {
    const parts = stage.filter((value): value is string => typeof value === "string" && value.trim().length > 0);
    return parts.length > 0 ? parts.join(" / ") : null;
  }

  if (heartbeat.missing === true) {
    return "heartbeat missing";
  }
  if (heartbeat.invalid_json === true) {
    return "heartbeat invalid";
  }

  return null;
}

export function formatHeartbeatMessage(heartbeat: RuntimeHeartbeat | undefined): string | null {
  if (!heartbeat || typeof heartbeat !== "object") {
    return null;
  }

  const message = heartbeat.message_text;
  return typeof message === "string" && message.trim().length > 0 ? message : null;
}
