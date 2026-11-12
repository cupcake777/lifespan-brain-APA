import type { ReportAsset } from "../types";

type AssetsPanelProps = {
  plots: ReportAsset[];
  downloads: ReportAsset[];
};

function AssetList({
  title,
  items,
  emptyLabel
}: {
  title: string;
  items: ReportAsset[];
  emptyLabel: string;
}) {
  return (
    <section className="asset-panel__section">
      <div className="asset-panel__header">
        <h3>{title}</h3>
        <span>{items.length}</span>
      </div>
      {items.length === 0 ? (
        <p className="asset-panel__empty">{emptyLabel}</p>
      ) : (
        <ul className="asset-list">
          {items.map((item) => (
            <li key={`${title}-${item.rel_path}`} className="asset-list__item">
              <div>
                <strong>{item.name}</strong>
                <p>{item.rel_path}</p>
              </div>
              <a href={item.rel_path} target="_blank" rel="noreferrer">Open</a>
            </li>
          ))}
        </ul>
      )}
    </section>
  );
}

export function AssetsPanel({ plots, downloads }: AssetsPanelProps) {
  return (
    <section className="asset-panel">
      <div className="asset-panel__intro">
        <p className="eyebrow">Run Assets</p>
        <h2>Plots and downloadable outputs</h2>
        <p className="asset-panel__copy">
          Indexed from the EasyColoc results directory so you can cross-check figures and exported tables.
          Paths are shown relative to that directory.
        </p>
      </div>
      <div className="asset-panel__grid">
        <AssetList
          title="Plots"
          items={plots}
          emptyLabel="No plots were indexed for this run."
        />
        <AssetList
          title="Downloads"
          items={downloads}
          emptyLabel="No downloadable tables were indexed for this run."
        />
      </div>
    </section>
  );
}
