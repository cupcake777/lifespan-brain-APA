import { formatPp4, rowKey } from "../lib/report";
import type { ReportRow } from "../types";

type ResultsTableProps = {
  rows: ReportRow[];
  selectedKey: string | null;
  onSelect: (row: ReportRow) => void;
};

export function ResultsTable({ rows, selectedKey, onSelect }: ResultsTableProps) {
  return (
    <div className="table-panel">
      <table className="results-table">
        <thead>
          <tr>
            <th>GWAS</th>
            <th>QTL</th>
            <th>Locus</th>
            <th>Phenotype</th>
            <th>PP4</th>
            <th>n_snps</th>
          </tr>
        </thead>
        <tbody>
          {rows.length === 0 ? (
            <tr>
              <td className="results-table__empty" colSpan={6}>
                No loci match the current query and PP4 threshold.
              </td>
            </tr>
          ) : (
            rows.map((row) => {
              const key = rowKey(row);

              return (
                <tr
                  key={key}
                  data-selected={selectedKey === key}
                  onClick={() => onSelect(row)}
                >
                  <td className="results-table__text-cell">{row.gwas_id}</td>
                  <td className="results-table__text-cell">{row.qtl_id}</td>
                  <td className="results-table__text-cell results-table__text-cell--compact">
                    {row.locus}
                  </td>
                  <td className="results-table__text-cell">{row.phenotype}</td>
                  <td>{formatPp4(row.pp4)}</td>
                  <td>{row.n_snps}</td>
                </tr>
              );
            })
          )}
        </tbody>
      </table>
    </div>
  );
}
