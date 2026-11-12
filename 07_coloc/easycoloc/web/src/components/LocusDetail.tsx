import { formatPp4 } from "../lib/report";
import type { ReportRow } from "../types";

type LocusDetailProps = {
  row: ReportRow | null;
};

export function LocusDetail({ row }: LocusDetailProps) {
  if (!row) {
    return (
      <aside className="detail-panel detail-panel--empty">
        <p className="eyebrow">Locus Detail</p>
        <h2>Select a locus</h2>
        <p className="detail-panel__copy">
          Use the table or top-hit strip to inspect one coloc signal in detail.
        </p>
      </aside>
    );
  }

  const previewPath = row.locuscompare_plot_file || row.plot_file;
  const previewLabel = row.locuscompare_plot_file ? "Locuscompare preview" : "Coloc preview";
  const previewIsImage = /\.(png|svg|jpe?g)$/i.test(previewPath);
  const links = [
    { label: "PNG", path: previewPath && /\.png$/i.test(previewPath) ? previewPath : "" },
    { label: "PDF", path: previewPath.replace(/\.(png|pdf)$/i, ".pdf") },
    { label: "JSON", path: row.locuscompare_json_file }
  ].filter((item) => item.path);

  return (
    <aside className="detail-panel">
      <p className="eyebrow">Locus Detail</p>
      <h2>{row.locus}</h2>
      <p className="detail-panel__copy">
        {row.phenotype} in {row.qtl_id} against {row.gwas_id}
      </p>
      <dl className="detail-grid">
        <div>
          <dt>GWAS</dt>
          <dd>{row.gwas_id}</dd>
        </div>
        <div>
          <dt>QTL</dt>
          <dd>{row.qtl_id}</dd>
        </div>
        <div>
          <dt>Phenotype</dt>
          <dd>{row.phenotype}</dd>
        </div>
        <div>
          <dt>PP4</dt>
          <dd>{formatPp4(row.pp4)}</dd>
        </div>
        <div>
          <dt>n_snps</dt>
          <dd>{row.n_snps}</dd>
        </div>
        <div>
          <dt>Plot</dt>
          <dd>{previewPath || "NA"}</dd>
        </div>
        <div>
          <dt>Source</dt>
          <dd>{row.source_file || "NA"}</dd>
        </div>
      </dl>
      <section className="detail-preview" aria-label={previewLabel}>
        {previewIsImage ? (
          <a href={previewPath} target="_blank" rel="noreferrer">
            <img src={previewPath} alt={`${previewLabel} for ${row.locus}`} />
          </a>
        ) : previewPath ? (
          <a className="detail-preview__placeholder" href={previewPath} target="_blank" rel="noreferrer">
            Open available plot
          </a>
        ) : (
          <div className="detail-preview__placeholder">No plot asset matched this result row.</div>
        )}
      </section>
      <div className="detail-links" aria-label="Plot downloads">
        {links.length === 0 ? (
          <span>No plot downloads available</span>
        ) : (
          links.map((item) => (
            <a key={`${item.label}-${item.path}`} href={item.path} target="_blank" rel="noreferrer">
              {item.label}
            </a>
          ))
        )}
      </div>
    </aside>
  );
}
