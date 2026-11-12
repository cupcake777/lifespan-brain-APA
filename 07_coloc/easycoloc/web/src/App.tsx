import { useEffect, useState } from "react";
import { AssetsPanel } from "./components/AssetsPanel";
import { LocusDetail } from "./components/LocusDetail";
import { MetricCard } from "./components/MetricCard";
import { ResultsTable } from "./components/ResultsTable";
import {
  filterResults,
  formatHeartbeatMessage,
  formatHeartbeatStage,
  formatPp4,
  normalizeReportPayload,
  rowKey,
  topHits
} from "./lib/report";
import type { ReportPayload, ReportRow } from "./types";

export default function App() {
  const [payload, setPayload] = useState<ReportPayload | null>(null);
  const [query, setQuery] = useState("");
  const [minPp4, setMinPp4] = useState(0.5);
  const [selected, setSelected] = useState<ReportRow | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;

    async function loadReport() {
      try {
        const response = await fetch("/api/report-data");
        if (!response.ok) {
          throw new Error(`HTTP ${response.status}`);
        }

        const data = normalizeReportPayload(await response.json());
        if (!cancelled) {
          setPayload(data);
          setSelected(data.results[0] ?? null);
        }
      } catch (err) {
        if (!cancelled) {
          setError(err instanceof Error ? err.message : "Failed to load report data");
        }
      }
    }

    void loadReport();

    return () => {
      cancelled = true;
    };
  }, []);

  const filteredRows = payload ? filterResults(payload.results, { query, minPp4 }) : [];
  const selectedKey = selected ? rowKey(selected) : null;
  const topRows = topHits(filteredRows, 5);
  const heartbeatStage = formatHeartbeatStage(payload?.runtime.heartbeat);
  const heartbeatMessage = formatHeartbeatMessage(payload?.runtime.heartbeat);
  const plotCount = payload?.assets.plots.length ?? 0;
  const downloadCount = payload?.assets.downloads.length ?? 0;

  useEffect(() => {
    if (!payload) {
      return;
    }

    if (filteredRows.length === 0) {
      if (selected !== null) {
        setSelected(null);
      }
      return;
    }

    if (!selected) {
      setSelected(filteredRows[0]);
      return;
    }

    const stillVisible = filteredRows.some((row) => rowKey(row) === selectedKey);
    if (!stillVisible) {
      setSelected(filteredRows[0]);
    }
  }, [filteredRows, payload, selected, selectedKey]);

  if (error) {
    return (
      <main className="app-shell">
        <section className="hero hero--compact">
          <p className="eyebrow">EasyColoc Report</p>
          <h1>Unable to load report data</h1>
          <p className="subtitle">{error}</p>
        </section>
      </main>
    );
  }

  if (!payload) {
    return (
      <main className="app-shell">
        <section className="hero hero--compact">
          <p className="eyebrow">EasyColoc Report</p>
          <h1>Loading interactive report</h1>
          <p className="subtitle">Fetching `/api/report-data` from the local report server.</p>
        </section>
      </main>
    );
  }

  return (
    <main className="app-shell">
      <section className="hero">
        <div className="hero__content">
          <div>
            <p className="eyebrow">EasyColoc Report</p>
            <h1>{payload.meta.project_name}</h1>
            <p className="subtitle">
              Publication-ready colocalization review with ranked PP4 signals, locus drill-downs,
              run provenance, and indexed outputs in one browser view.
            </p>
            <div className="hero__actions" aria-label="Report highlights">
              <span>ABF + SuSiE</span>
              <span>Build-aware references</span>
              <span>Interactive QC</span>
            </div>
          </div>
          <dl className="hero__meta">
            <div>
              <dt>Results</dt>
              <dd>{payload.meta.results_dir}</dd>
            </div>
            <div>
              <dt>Generated</dt>
              <dd>{payload.meta.generated_at}</dd>
            </div>
            <div>
              <dt>Version</dt>
              <dd>{payload.meta.report_version}</dd>
            </div>
          </dl>
        </div>
        <div className="hero__visual" aria-label="Report overview">
          <div className="signal-card signal-card--primary">
            <span>Mean PP4</span>
            <strong>{formatPp4(payload.summary.mean_pp4)}</strong>
            <small>{payload.summary.significant_pp4_08} signals at PP4 &gt;= 0.80</small>
          </div>
          <div className="signal-card">
            <span>Indexed Plots</span>
            <strong>{plotCount}</strong>
            <small>{downloadCount} downloadable outputs</small>
          </div>
          <div className="signal-bars" aria-hidden="true">
            <i style={{ height: "42%" }} />
            <i style={{ height: "68%" }} />
            <i style={{ height: "36%" }} />
            <i style={{ height: "84%" }} />
            <i style={{ height: "58%" }} />
            <i style={{ height: "92%" }} />
          </div>
        </div>
        <div className="hero__status">
          {heartbeatStage ? <p className="runtime-pill">Run status: {heartbeatStage}</p> : null}
          {heartbeatMessage ? <p className="hero__note">{heartbeatMessage}</p> : null}
          {payload.warnings.length > 0 ? (
            <ul className="warning-list">
              {payload.warnings.map((warning) => (
                <li key={warning}>{warning}</li>
              ))}
            </ul>
          ) : null}
        </div>
      </section>

      <section className="metric-grid">
        <MetricCard label="Total tests" value={payload.summary.total_tests} tone="accent" />
        <MetricCard label="PP4 >= 0.80" value={payload.summary.significant_pp4_08} />
        <MetricCard label="PP4 >= 0.70" value={payload.summary.significant_pp4_07} />
        <MetricCard label="Mean PP4" value={formatPp4(payload.summary.mean_pp4)} />
        <MetricCard label="Unique genes" value={payload.summary.unique_genes} />
        <MetricCard label="Unique loci" value={payload.summary.unique_loci} />
      </section>

      <section className="workspace">
        <div className="explorer">
          <div className="explorer__header">
            <div>
              <p className="eyebrow">Explore Results</p>
              <h2>Screen coloc signals interactively</h2>
            </div>
            <p className="explorer__count">
              {filteredRows.length} of {payload.results.length} rows shown
            </p>
          </div>

          <div className="explorer-toolbar">
            <label className="field">
              <span>Search</span>
              <input
                type="search"
                value={query}
                onChange={(event) => setQuery(event.target.value)}
                placeholder="GWAS, QTL, locus, phenotype, source file"
              />
            </label>

            <label className="field field--range">
              <span>Minimum PP4</span>
              <div className="range-control">
                <input
                  type="range"
                  min="0"
                  max="1"
                  step="0.05"
                  value={minPp4}
                  onChange={(event) => setMinPp4(Number(event.target.value))}
                />
                <strong>{formatPp4(minPp4)}</strong>
              </div>
            </label>
          </div>

          <div className="top-hit-strip">
            {topRows.length === 0 ? (
              <p className="top-hit-strip__empty">No top hits at the current threshold.</p>
            ) : (
              topRows.map((row) => {
                const key = rowKey(row);

                return (
                  <button
                    key={key}
                    type="button"
                    className="top-hit"
                    data-active={selectedKey === key}
                    onClick={() => setSelected(row)}
                  >
                    <span>{row.locus}</span>
                    <strong>{formatPp4(row.pp4)}</strong>
                    <small>
                      {row.phenotype} · {row.qtl_id}
                    </small>
                  </button>
                );
              })
            )}
          </div>

          <ResultsTable rows={filteredRows} selectedKey={selectedKey} onSelect={setSelected} />
        </div>

        <LocusDetail row={selected} />
      </section>

      <AssetsPanel plots={payload.assets.plots} downloads={payload.assets.downloads} />
    </main>
  );
}
