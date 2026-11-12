# Disease track (modules 06-07): bugs and broken links in existing files

Reconstructed 2026-10-19. Existing release files were not edited. Each item gives the file, the line, what breaks, and a proposed one-line (or short) fix for the owner to apply.

1. release/06_disease_heritability_geneset/followup/registry_checks/check_m06_07.py, lines 13-14 and 23.
   The binary S-LDSC FDR (fdr4) is BH across the four annotations, labelled as within-GWAS. The manuscript text uses BH across the two stages (prenatal, postnatal). The file has no two-stage column.
   Fix: before line 14, restrict to the two stages, e.g. b2 = b[b.Annotation.isin(['prenatal_apa_qtl','postnatal_apa_qtl'])], compute fdr2 on b2 grouped by (Population, Trait), and emit a _fdr2 value.

2. release/06_disease_heritability_geneset/followup/registry_checks/check_m06_07.py, line 7.
   OUT = A + '1016/registry_1017/m06_07/' hard-wires an output folder inside the working tree hand-off (1016), so running the release script overwrites a file outside the release.
   Fix: OUT = os.environ.get('RERUN_OUT', A + '1016/registry_1017/m06_07/') or an argument; the verifier verify_rerun_values.py already redirects it to a temp folder.

3. release/06_disease_heritability_geneset/followup/registry_checks/build_registry.py, line 2 (docstring) and line 294.
   The docstring says the script writes gwas_manifest_registered_m06_07.tsv. The code writes only registry.tsv (line 294). The manifest file has no producer in release.
   Fix: correct the docstring, or add the manifest writer if the file is wanted.

4. release/06_disease_heritability_geneset/sLDSC/07_collect_results.py, line 114 (FDR loop) and its docstring.
   FDR_p is one BH over every row it reads (76 rows: 19 GWAS x 4 annotations). The docstring does not say so. The manuscript text uses a two-stage family, which this file does not produce.
   Fix: document the 76-row family in the docstring; optionally add an FDR_2stage column computed per (Population, Trait) over the prenatal and postnatal rows (see sldsc_two_stage_fdr.py, function two_stage_fdr).

5. Registry notes (docs/number_registry.tsv rows for R092 and R058; m06_07 registry README mismatch table).
   The notes say the processed FDR reproduces only with BH across 4 annotations, and that the family is 4-annotation. Checked here: the processed FDR_p is the 76-row BH; the text values follow the two-stage family.
   Fix: change the family description in the notes to "two stages within GWAS (sldsc_two_stage_fdr.py)".

6. Registry row R026 and its note: "master summary 438 is the outlier". Not supported. twas_ranked_all.tsv.gz covers 17 GWAS; EAS_EP and EUR_EP are absent. The master summary's extra hit is EUR_EP postnatal (p = 3.8e-10, below threshold).
   Also: no release script writes twas_master_summary.tsv or twas_ranked_all.tsv.gz (04_summarize_spredixcan_association.py is per GWAS). The aggregation step is missing from release.
   Fix: decide whether EP is analysed; either add EP to the ranked aggregate (438) or state 17 GWAS in Results (437). Add the aggregation script to release.

7. Two S-LDSC runs disagree on EUR PD h2: 1013 sldsc_baseline_h2_tau_summary.tsv gives 0.0153 (z 8.05); 1015 03_sldsc_tau_star.tsv gives 0.0147 (z 6.39). The z < 7 gate statement (R057) is correct only for the baseline run.
   Fix: name the run in the Methods gate sentence.

8. results_final/RESULTS_MANIFEST.tsv rows for 06 magma/eur_scz_*.gsa.out say "unlocated (... registry R117 says server run)". The files are present locally at 1013/collect/eur_scz_*.gsa.out.
   Fix: update those rows to point at 1013/collect.
   Also 1013/collect/run_eur_scz_followup.sh (driver) is not readable in this session (permission denied); the owner should confirm the summary writer.

9. release/06_disease_heritability_geneset/MAGMA/run_enrichment.R, line 227, writes enrichment_summary.tsv. The registry and RESULTS_MANIFEST cite it as the producer of magma_gsa_summary.tsv.
   Fix: confirm the producer of magma_gsa_summary.tsv; if none, record it as unlocated.

10. release/07_locus_colocalization_smr/coloc/audit_susie_convergence.R, line 591, writes all_susie_results.csv. The registry says "lists source only" and the RESULTS_MANIFEST says unlocated.
    Fix: update the manifest; the producer is identified (inputs not local).

11. release/07_locus_colocalization_smr/followup/handoff_module07/ (module 07 context tables) and the RESULTS_MANIFEST rows for 1016/handoff_module07_1016/context/*.tsv cite "upstream 14_context/*.R on HPC not in release".
    Fix: add the 14_context producers to release or mark those rows server-only.

12. 1014/figure/round1_panels/fig6c_arl3_eqtl_track.R, line 2 (PENDING_SERVER) and line 5 (hard-coded absolute path into the working tree). Outside release; not edited.
    Fix: produce Fig6c_ARL3_eqtl_regional.tsv on the server; read the path from an environment variable.

13. Coloc copies: 1013/collect/coloc_prenatal_apa_eqtl.tsv and coloc_postnatal_apa_eqtl.tsv are an older run than the frozen 1015 copies (662 vs 661 and 2,227 vs 2,225 rows; PP4 >= 0.7 counts 165 vs 186 and 576 vs 655). The RESULTS_MANIFEST lists the 1015 copies as producer-identified; the 1013 copies are listed as unlocated.
    Fix: mark the 1013 copies superseded in the manifest.

14. Fig 5B input (results_final/figures/fig5/fig5b_magma_plot_data.tsv): no producer in release or working tree; EAS_BIP prenatal p = 0.696310 appears in no local MAGMA table.
    Fix: record the producer or drop the file from the figure input list.
