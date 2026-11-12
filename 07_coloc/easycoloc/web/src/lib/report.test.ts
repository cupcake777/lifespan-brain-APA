import { describe, expect, it } from "vitest";
import { filterResults, normalizeReportPayload, rowKey, topHits } from "./report";

const rows = [
  {
    gwas_id: "GWAS_A",
    qtl_id: "postnatal",
    locus: "rs100",
    phenotype: "GENE1",
    pp4: 0.92,
    n_snps: 41,
    plot_file: "plots/GWAS_A_rs100_coloc.png",
    locuscompare_plot_file: "plots/GWAS_A_rs100_locuscompare.png",
    locuscompare_json_file: "plot_data/GWAS_A_rs100_locuscompare.json",
    source_file: "GWAS_A_rs100_locus_results.csv"
  },
  {
    gwas_id: "GWAS_B",
    qtl_id: "prenatal",
    locus: "rs101",
    phenotype: "GENE2",
    pp4: 0.41,
    n_snps: 25,
    plot_file: "",
    locuscompare_plot_file: "",
    locuscompare_json_file: "",
    source_file: "GWAS_B_rs101_locus_results.csv"
  }
];

describe("report helpers", () => {
  it("filters by free-text query and minimum pp4", () => {
    const filtered = filterResults(rows, { query: "gene1", minPp4: 0.8 });

    expect(filtered).toHaveLength(1);
    expect(filtered[0].locus).toBe("rs100");
  });

  it("returns top hits in descending PP4 order", () => {
    const hits = topHits(rows, 1);

    expect(hits).toHaveLength(1);
    expect(hits[0].pp4).toBe(0.92);
  });

  it("builds stable keys that distinguish rows from different source files", () => {
    const variantRow = {
      ...rows[0],
      source_file: "GWAS_A_rs100_locus_results_recomputed.csv"
    };

    expect(rowKey(rows[0])).not.toBe(rowKey(variantRow));
  });

  it("normalizes malformed payloads into safe defaults", () => {
    const normalized = normalizeReportPayload({
      meta: { project_name: "demo" },
      results: [{ gwas_id: "GWAS_A", pp4: "0.91" }],
      assets: { plots: "nope" },
      warnings: null
    });

    expect(normalized.meta.project_name).toBe("demo");
    expect(normalized.meta.results_dir).toBe("");
    expect(normalized.summary.total_tests).toBe(0);
    expect(normalized.results).toHaveLength(1);
    expect(normalized.results[0].qtl_id).toBe("");
    expect(normalized.results[0].pp4).toBe(0.91);
    expect(normalized.results[0].plot_file).toBe("");
    expect(normalized.results[0].locuscompare_plot_file).toBe("");
    expect(normalized.results[0].locuscompare_json_file).toBe("");
    expect(normalized.assets.plots).toEqual([]);
    expect(normalized.assets.plot_data).toEqual([]);
    expect(normalized.assets.downloads).toEqual([]);
    expect(normalized.warnings.length).toBeGreaterThan(0);
  });
});
