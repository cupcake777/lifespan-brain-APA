# v19 数字溯源：未解决项（2026-10-17）

来源：[number_registry.tsv](number_registry.tsv)，共 420 行，由三个分模块 registry（m01_03、m04_05、m06_07）合并。
去重规则：段落号和数值都相同的行，§2.4 疾病分析以 06–07 的登记为准；同一数字只要有一行复现成功，就不再保留 table-only 那一行。

状态计数：reproduced 266，table-only 39，no-producer 36，server-only 29，mismatch 21，not-rerun 19，arithmetic 10

## A. 文本与结果表不一致（mismatch，需要改正文或改源表）

| reg_id | 段落 | 正文值 | 复核值 | 来源表 | 说明 |
|---|---|---|---|---|---|
| N0008 | Abstract P3 | 40% | 38.6 (27/70) | 1016/handoff_module07_1016/context/gene_gwas_classification_1mb.tsv | Results and Discussion say 27 of 70 = 39%. Abstract says 40%. Use 39%. |
| N0009 | Introduction P9 | 72; 40% | 38.6 (27/70) | 1016/handoff_module07_1016/context/gene_gwas_classification_1mb.tsv | Same 40% vs 39% issue as the Abstract. |
| N0036 | §2.1 P15 | 3.7 × 10⁻⁵ | 3.1 × 10⁻⁵ (refit from adjusted data, β −0.020022) | 1014/outputs/11_loeuf_adjusted_models.tsv (3.10e-5); 1014/11_loeuf_adjusted_models.tsv (3.7e-5, β −0.019825) | The 3.7e-5 value comes from a superseded root-level model whose β (−0.019825) is not reproduced from the current adjusted data. |
| N0042 | §2.2 P17 | ≤4% | C1 3.5%; C2 4.8% | 1014/figure/data/fig2a_event_coupling.tsv | C2 exceeds 4% on the Pearson r² definition; Spearman ρ² gives 2.8% for C2. Author must state the definition or change the bound. |
| N0118 | 2.3 P23 | 0.001 | 0.003362 | 1016/r4_nan_pvalues/r4_nan_adjusted_p.tsv | R4: naive logit P=0.0034 (OR 0.38); gene-clustered SE 0.046; GEE exchangeable 0.025; pair-level Fisher 1.3e-6; gene-level Fisher 0.0034 (18/69 vs 22/39). No specification gives 0.001. Unresolved; see unresolved.md. |
| N0165 | Results 2.4 P28 | FDR = 0.013 | 0.017897970578517733 | 1013/collect/13_sldsc_binary/sldsc_enrichment_processed.tsv | Correct value 0.018. 0.013 is the nominal P (0.0134). The 76-test FDR is 0.107. Fix to 0.018 or state the family. |
| N0167 | Results 2.4 P28 | FDR = 0.004 | 0.0050209054485088 | 1013/collect/13_sldsc_binary/sldsc_enrichment_processed.tsv | Enrichment 1.41 reproduces (1.4124). FDR 0.00502 rounds to 0.005, not 0.004. |
| N0198 | Results 2.4 P30 | 437 | 438 | 1013/collect/twas_master_summary.tsv | Master summary sums to 438 (postnatal 262 vs 261 in ranked table). Fix at source. |
| N0201 | Results 2.4 P30 | 50 of 76 | 64 of 104 | 1013/collect/all_susie_results.csv | all_susie_results.csv: 64 of 104. susie_rerun_summary.tsv: 66 of 106. Neither gives 50/76. Text likely from an earlier run. |
| N0206 | 2.4 P31 | rs7047 | rs7251 | 1013/collect/all_colocalization_results.csv | The coloc table lead for SCAF1 (EUR SCZ prenatal) is rs7251, not rs7047. Handoff R3 flags the same issue. rs7047 may be the external reference 3′aQTL lead; check before the figure text is fixed. |
| N0224 | Results 2.4 P32 | PP4 = 0.99 | 0.986 (interaction) | 1013/collect/all_colocalization_results.csv | PP4 reproduces, but the colocalizing QTL is the interaction term, not the prenatal apaQTL. Prenatal and postnatal ATP2B1 hits are below 0.7. Relabel as G x Stage. |
| N0225 | Results 2.4 P32 | PP4 = 0.86 | 0.857 (interaction) | 1013/collect/all_colocalization_results.csv | Value reproduces. Same label issue (interaction term, not prenatal apaQTL). ATP2B1 coloc.susie not evaluable (0 comparisons). |
| N0231 | Results 2.4 P32 | qualitative | no OGFOD2 x EUR_INT comparison; EAS SCZ only (max 0.15) | 1013/collect/all_susie_results.csv | Unsupported. Also omitted: a prenatal CHB eQTL colocalizes with this signal (PP4 0.899), so OGFOD2 is APA + eQTL, not APA-only. |
| N0245 | Figure 6 legend P39 | qualitative (label) | interaction | 1013/collect/all_colocalization_results.csv | Same relabel as Results 2.4. |
| N0260 | Discussion P53 | ATP2B1 + 24 = 25 | 26 in total (ATP2B1 + 25 others) | Results para 25 (2.3) | Results para 25 reads "25 stage-dependent genes" beyond ATP2B1. Total 26, not 25. Outside module 06-07; cross-check only. |
| N0321 | Methods: external validation P91 | 25,178 | 25,190 (table) | 1014/outputs/08_brainspan_liftover_qc.tsv (n_strictly_matched_dynamic_sites 25,190; n_matched_raw_events 25,190) | Text 25,178 vs table 25,190 (both columns). |
| N0331 | Methods: protein P94 | 20,334; 5,252 | 20,450 pairs; 5,264 genes (table); 20,447 correlations | 1014/results/03_site_protein_pairs.tsv (20,450 rows); 1014/results/01_sample_qc.tsv (primary_universe 20450 / 5264); 1014/results/09_global_triad_correlations.tsv (20,447 rows) | Three different counts in local tables (20,450 / 20,447) vs text 20,334. Gene count 5,264 vs 5,252. |
| N0349 | Methods P101 | 0.001 (Fig 4e P, Methods model) | P 0.0034 (naive); 0.046 (gene-clustered); 0.025 (GEE) | 1016/r4_nan_pvalues/r4_nan_adjusted_p.tsv | The Methods unit (gene × cell type) and clustering are not specified; the three variants bracket 0.001 from above. |
| N0385 | Methods: partitioned heritability P111 | EUR z<7 (3 listed) | EUR PD z=6.39; EUR AD (Kunkle) z=6.55; ANX 1.96; EP 5.64; INS 5.97 | 1015/handoff_eqtl_1015/03_sldsc_tau_star.tsv | EUR PD (6.39) and AD Kunkle (6.55) are below 7 and were analysed, but the EUR list omits them. The EUR PD result is in Results. Resolve. |
| N0386 | Methods: partitioned heritability P111 | two stages | across 4 annotations | 1013/collect/13_sldsc_binary/sldsc_enrichment_processed.tsv | Reported FDRs reproduce only with BH across the 4 annotations. Correct the Methods wording. |
| N0396 | Methods: TWAS P115 | 6,007; 7,326 | unique model_id: prenatal 5981 ; postnatal 7292 | 1013/collect/twas_ranked_all.tsv.gz | Ranked table has 5,981 and 7,292. Training outputs are server-only, so the gap cannot be resolved locally. |

## B. 子代理报告的未解决项（包括无 producer、版本冲突、口径未定义）

| 模块 | 编号 | 数值 | 问题 | 建议 |
|---|---|---|---|---|
| 01-03 | M03-027 | 3.7 × 10⁻⁵ | mismatch: BH-adjusted p for LOEUF membership coefficient; rerun/table gives 3.1 × 10⁻⁵ (refit from adjusted data, β −0.020022) | The 3.7e-5 value comes from a superseded root-level model whose β (−0.019825) is not reproduced from the current adjusted data. |
| 01-03 | M03-033 | ≤4% | mismatch: variance explained, PDUI change vs expression change; rerun/table gives C1 3.5%; C2 4.8% | C2 exceeds 4% on the Pearson r² definition; Spearman ρ² gives 2.8% for C2. Author must state the definition or change the bound. |
| 01-03 | M03-119 | 25,178 | mismatch: matched events; rerun/table gives 25,190 (table) | Text 25,178 vs table 25,190 (both columns). |
| 01-03 | M03-129 | 20,334; 5,252 | mismatch: APA–protein pairs; genes; rerun/table gives 20,450 pairs; 5,264 genes (table); 20,447 correlations | Three different counts in local tables (20,450 / 20,447) vs text 20,334. Gene count 5,264 vs 5,252. |
| 01-03 | M03-022 | 6.56 | version conflict between accepted outputs/ table and root-level table (outputs 6.52; root 6.56) | Decide canonical module 02 run; rerun 05_func.R or update text. |
| 01-03 | M03-023 | 2.85 | version conflict between accepted outputs/ table and root-level table (outputs 2.84; root 2.85) | Decide canonical module 02 run; rerun 05_func.R or update text. |
| 01-03 | M03-025 | 5.15 | version conflict between accepted outputs/ table and root-level table (outputs 5.14; root 5.15) | Decide canonical module 02 run; rerun 05_func.R or update text. |
| 01-03 | M03-029 | 0.002 | version conflict between accepted outputs/ table and root-level table (outputs 4.2e-4; root 2.1e-3) | Decide canonical module 02 run; rerun 05_func.R or update text. |
| 01-03 | M03-PIK3CA | At PIK3CA, six of seven events | producer not located; not checked | Locate Fig 2e producer and verify. |
| 01-03 | M03-read-count | rho 0.19/-0.17; pi1 0.67 vs 0.64; protein 51.7% vs 49.4% | read-count length artifact (1016/session_handoff_1017.md): values reproduced but provisional | Complete Codex R1 (UTR-free expression) and rerun §2.2 coupling before submission. |
| 01-03 | M03-timing | 92%, 96%, 52%, 37% | reproduce only with median-curve estimator; Methods does not define it | Add the estimator definition to Methods. |
| 01-03 | M03-VERIFY | [[VERIFY: justify 40 surrogate variables ...]] | placeholder in Methods protein paragraph | Resolve placeholder with num.sv estimate and robustness check. |
| 01-03 | M03-module02-status | 63,097 / 39,580 etc. | results.md (10-12) holds module 02 pending donor-level rerun; acceptance file (10-14) accepts; GAMs not refit in this pass | Confirm governing status; rerun donor-level GAM if required. |
| 01-03 | M03-report-39555 | 39,580 dynamic sites | 1014/outputs/report.md says 39,555 | Regenerate report.md; per-event table gives 39,580. |
| 01-03 | M03-repeats | 66 repeat libraries vs 65 bridge pairs | one repeat without retained original not explained | Add sentence to Methods. |
| 01-03 | M03-839 | 839 donors | no sample-to-donor map in local files; donor counts unverified | Provide donor ID map or verify with cohort metadata. |
| 04-05 | R068 / R170 | 0.001 (Fig 4e, Nan reversed vs other pairs) | Not reproduced under any specified model: naive logit P=0.0034, gene-clustered SE 0.046, GEE exchangeable 0.025, pair-level Fisher 1.3e-6, gene-level Fisher 0.0034. Methods leave the unit and clustering unspecified. | Name the exact model and unit in Methods and quote the matching P (e.g. 0.003 naive, or report the gene-clustered 0.046); or supply the script that yields 0.001. |
| 04-05 | R005 | 40% (Abstract) | Results §2.4 gives 27 of 70 (39%); Abstract rounds to 40%. | Use 39% in the Abstract. |
| 04-05 | R127 | 50 of the 76 evaluable signals (coloc.susie) | all_susie_results.csv has 104 loci; 64 with max PP4>=0.7. The 76 evaluable and 50 corroborated counts are not reproducible locally. | Identify the 76 evaluable comparisons (credible-set availability) and export the 50 from the Furong susie run. |
| 04-05 | R132 | lead variant rs7047 (SCAF1) | Coloc lead for SCAF1 with EUR SCZ prenatal is rs7251 (PP4 0.915). rs7047 is a reported external 3′aQTL lead, not the coloc lead. | Say "rs7251 is the colocalized lead; rs7047 is the reported 3′aQTL" or align the locus lead. |
| 04-05 | R117 | P = 8.5e-3 (MAGMA prenatal conditional on postnatal) | Only condition_hide (all-tested background) is local; the --model condition=postnatal run is on Furong. | Export magma condition=postnatal output for EUR SCZ into 1013/collect/15_magma. |
| 04-05 | TWAS | 437 gene–trait associations; 108 prenatal, 140 postnatal genes | Ranked TWAS table is local, but the Bonferroni-per-GWAS-and-stage counting rule was not reproduced in this run. | Export the significant-association table from the Furong TWAS run with its threshold column. |
| 04-05 | R151 | 12 reversed loci | The HC3 table has 193 rows; no panel-d selection column. | Add the locus-selection rule (e.g. largest \|β\| reversed, or CI excludes 0 both stages) to the figure source table. |
| 04-05 | R137 | 0.89 (ATP2B1 ADHD × cerebellar eQTL) | GTEx coloc table is keyed by signal_id without gene or GWAS labels; ATP2B1 signal not mapped. | Add gene and GWAS columns to gtex_brain_eqtl_coloc_all_tissues.tsv. |
| 04-05 | R128 | 27 (39%) APA-only at 26 loci (of 70 genes) | Module 07 eqtl_coloc.tsv present but counts not recomputed here. | Recompute APA-only classification in module 07 and record the rule. |
| 04-05 | R198/R202 | 1e-5 instrument P; 2 of 72 genes unclassified | Not checked against SMR instrument columns. | Add apaQTL instrument P column to smr_significant.tsv. |
| 04-05 | Methods 98 | 45,392 pairs with estimates in 20 subsets | Per-subset pair files are not local. | Export the per-seed pair counts to 1014/figure/data. |
| 04-05 | Methods 122 | cell-composition accuracy (0.77, 0.98; 129/774 mixtures; 43/258 donors) | No producer or output table located locally. | Locate the BisqueRNA validation output; otherwise cite as server-only. |
| 04-05 | R069-R076 | S4 gating table (Fig 4f,g) | S4_gating_eqtl_fixed_RNA.tsv has no producer script; it carries both dexpr and dexpr_v3; text values match dexpr_v3 only. | Add producer for S4 gating table; name dexpr_v3 in Methods. |
| 04-05 | R091-R100 | FDR column (sldsc_enrichment_processed.tsv) | FDR_p is BH across 4 annotations and has an older AD row (postnatal 1.40); the text uses 2-stage BH within GWAS and AD2. | Regenerate the processed table with 2-stage BH and AD2 (Fig5a_EUR_AD2) so the table matches the text. |
| 04-05 | R015/R016 | 16 of 23 apaQTL variants | 23 are stage × event × variant candidate rows (17 unique variants). | Write "16 of 23 PAS-disrupting candidates (17 variants)". |
| 04-05 | R032/R033 | 97.8% of 49,957 pairs | Gene × variant key gives 49,957 (97.8%). Phenotype × variant key gives 42,692 (99.9%). | State the key in Methods (gene × variant). |
| 04-05 | framework.md Module 05 | 559/593 ATP2B1 saturation (framework) | Superseded: table gives 502/508 (manuscript value). | Update framework.md Module 05 README. |
| 04-05 | §2.4 placeholders | Supplementary Table [[maxCPP]], [[GWAS]], [[VERIFY]] | Unfilled placeholders in v19. | Fill before submission. |
| 06-07 | R006 | EUR BD prenatal FDR = 0.013 | Table within-GWAS BH FDR is 0.0179 (76-test 0.107); 0.013 is the nominal P 0.0134. | Replace with FDR = 0.018, or state the FDR family explicitly. |
| 06-07 | R008 | EUR SCZ postnatal FDR = 0.004 | Table FDR is 0.00502, which rounds to 0.005. | Replace with FDR = 0.005. |
| 06-07 | R029 | coloc.susie corroborated 50 of the 76 evaluable signals | all_susie_results.csv gives 64 of 104 comparisons with PP4>=0.7; susie_rerun_summary.tsv gives 66 of 106. Neither matches. | Replace with current counts and state the denominator rule. |
| 06-07 | R002 | 40% missed by eQTL maps (Abstract; also Introduction R003) | 27/70 = 38.6%; Results and Discussion say 39%. | Use 39% in Abstract and Introduction. |
| 06-07 | R036 | prenatally restricted apaQTL colocalized with EUR ADHD (PP4 = 0.99) and EAS SCZ  | PP4 values reproduce, but the colocalizing QTL is the G x Stage interaction term (QTL_ID = interaction). Prenatal and postnatal ATP2B1 PP4 < 0.7. coloc.susie not evaluable (0 comparisons). | Relabel as G x Stage/interaction; do not describe as a prenatal-specific apaQTL colocalization. |
| 06-07 | R043 | OGFOD2 prenatal apaQTL colocalized with intelligence ... with coloc.susie suppor | No coloc.susie comparison exists for OGFOD2 x EUR_INT (only EAS SCZ, max PP4 0.15). Also omitted: a prenatal CHB eQTL colocalizes with this signal (PP4 0.899), so OGFOD2 is APA + eQTL, not APA-only. | Remove the coloc.susie claim or run it; state the eQTL colocalization. |
| 06-07 | R047 | At ATP2B1 and 24 other stage-dependent genes (Discussion) | Results para 25 says 25 stage-dependent genes beyond ATP2B1 (26 total). Discussion gives 25 total. Outside module 06-07 scope; cross-check. | Reconcile the saturated-gene count between Results 2.3 and the Discussion. |
| 06-07 | R057 | EUR GWAS with h2 z < 7 not analysed: anxiety, epilepsy, insomnia | EUR PD (h2 z 6.39) and EUR AD Kunkle (6.55) are below 7 and were analysed. The EUR PD binary result (1.89, FDR 0.046) is reported in Results. The EAS list matches exactly. | Apply the z<7 gate to EUR PD or caveat its S-LDSC result; correct the EUR list. |
| 06-07 | R058 | binary FDR BH-corrected across the two stages within each GWAS | Reported FDRs reproduce only with BH across the 4 binary annotations (prenatal, postnatal, lifespan, interaction). | Correct the Methods wording. |
| 06-07 | R026 | 437 TWAS gene-trait associations | Ranked table gives 437 (correct); twas_master_summary.tsv sums to 438 (postnatal 262 vs 261). | Fix the master summary at source. |
| 06-07 | R063 | retained models 6,007 prenatal and 7,326 postnatal | Unique model_id in twas_ranked_all: 5,981 prenatal and 7,292 postnatal. Training outputs are server-only, so this cannot be resolved locally. | Confirm the count from the server training log and correct the text. |
| 06-07 | R010 | AD (1.79-fold; FDR = 0.022) using Bellenguez EUR_AD2 | Values reproduce from 1014/sldsc_ad2_summary.tsv, but the GWAS manifest has only Kunkle EUR_AD (21,984 / 41,942) and no AD2 row. The '19 GWAS' count holds only if AD2 replaces Kunkle. The AD2 N (85,934 / 401,577; sum 487,511) has no local table. | Add the EUR_AD2 row (registered copy written) and state the substitution in Methods. |
| 06-07 | R077 | EAS BD N used: 4,479 cases / 75,725 controls (PGC bip2024 no23andMe) | User-stated. The manifest N (5,969 / 145,911) is the paper total including 23andMe. HPC EAS_BIP N columns not yet checked locally. | Check the HPC file N columns, then update the 1013 manifest with author approval. |
| 06-07 | R078 | EUR BD N used: 59,287 cases / 781,022 controls (PGC bip2024 no23andMe) | User-stated. Manifest N (132,046 / 2,322,339) is the paper total including 23andMe. | As for EAS BD. |
| 06-07 | R023 | 543 of 669 prenatal, 1,869 of 2,249 postnatal, 154 of 193 interaction genes (MAG | Only 154 is checkable (interaction set). Per-GWAS ngenes are 526 (EUR SCZ prenatal) and 1,824 (postnatal), not the 543 / 1,869 denominators. Gene-location file is server-only. | Report the denominator rule with the per-GWAS ngenes. |
| 06-07 | R052 | Supplementary Table [[GWAS]] and other [[...]] placeholders (Methods 109, 111, R | Placeholders remain in the manuscript text: [[GWAS]], [[maxCPP]], [[MAGMA_trajectory]], [[saturation]]. | Fill from the registered GWAS manifest and the supplementary table plan before submission. |
| 06-07 | GWAS-CITES | GWAS citations: EUR_ANX (Otowa 2016 vs manifest Levey 2020), EUR_INS (Hammerschl | The manifest literature column disagrees with manuscript reference numbers for four GWAS. Several candidate mappings are unverified. | Confirm each reference and update the manifest and bibliography. |

## C. 各模块声明的覆盖限制（复核深度）

- **01–03：** GAM 重拟合、置换零分布、bootstrap、敏感性分析没有重跑（标为 not-rerun 或 table-only）；原始 APA 计数矩阵只核对了行数；服务器步骤（ChinaMAP、GATK、PLINK）标为 server-only。
- **04–05：** tensorQTL、eigenMT/HC3、mashr、mediation、coloc 都是读本地结果表核对的，没有重跑；R 脚本（mashr、coloc.R、mediate.R、fig3b）没有执行，hexamer 16/23 改用 python 复现；Nan 2026 的完整 sumstats 已重新流式读取，结果与工作簿一致。
- **06–07：** v19 没有 §2.5 和 §2.6，疾病结果在 §2.4；GWAS manifest 没有改动，带登记的副本放在 m06_07/；ARL3 本地 coloc 使用队列 AF 作为代理；BIP 的 N 依据 figshare 元数据（no23andMe），HPC 上的文件还没核实。

上述限制意味着：标为 table-only 或 server-only 的数字，只核对了表和正文是否一致，没有核对分析本身。