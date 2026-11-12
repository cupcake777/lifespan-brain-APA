这份分析涵盖了全基因组测序、RNA-seq、蛋白质组、统计遗传学和多层级 GWAS 整合，流程很长。

为了清晰、模块化地整理脚本库（Scripts Repository），下面按 00 至 07 模块梳理分析流水线。现有分析脚本和已确认的实际执行流程共同作为方法与参数来源；由公司执行或必须通过网页完成的步骤明确记录执行边界，不将“本机脚本未保留”误写为“该步骤未执行”。

---

### 建议脚本库目录架构（Repository Structure）

```text
1011/
├── 00_genotype_qc_imputation/       # WGS 质控、ChinaMAP 填充与遗传主成分
├── 01_rnaseq_apa_quant/             # 表达量定量、DaPars2 定量与批次校正
├── 02_lifespan_trajectories/        # GAM 拟合、C1/C2 聚类与功能富集 (Fig. 1)
├── 03_multiomics_uncoupling/        # APA vs RNA 解耦与蛋白质组效应 (Fig. 2)
├── 04_stage_specific_apaqtl/        # 产前/产后 apaQTL、保守性与 eQTL 解耦 (Fig. 3)
├── 05_gxstage_interaction/          # G×Stage 交互、表达门控与表型饱和 (Fig. 4)
├── 06_disease_heritability_geneset/ # S-LDSC、MAGMA 与 TWAS (Fig. 5)
└── 07_locus_colocalization_smr/     # coloc.susie、SMR 与重点位点绘图 (Fig. 5 & 6)

```

仓库中的 `1011/00_*` 至 `1011/07_*` 是 publication-facing 入口；入口内的
链接指向保留的历史 producer 目录，因此不会改变原脚本的工作目录、绝对路径
或下游结果位置。

---

### 详细流水线清单与脚本对照表

#### **Module 00: WGS 质控、填充与遗传主成分分析**

*针对全基因组测序（WGS）的标准化清洗与遗传协变量提取。*

* **00_0_wgs_alignment_gvcf（local scripts + company/manual CLI checkpoint）**：所有样本使用同一套 reference、resource files 和参数完成 WGS 比对、去重复、BQSR 与 GATK HaplotypeCaller per-sample gVCF calling；部分样本由本地运行，部分由测序公司按实验室提供的流程分块批量运行。部分中间命令由 CLI 手动执行，未全部保留为独立 script。
* **00_1_joint_genotyping.sh**（`1010/preprocess/01_joint_genotyping.sh`）：GATK 4.5 `GenomicsDBImport` 与 `GenotypeGVCFs` 按染色体联合分型，再合并为 cohort VCF。
* **00_2_vqsr_and_sample_qc.sh**：GATK VQSR；随后进行样本缺失率、杂合度、KING 亲缘关系及位点缺失率/HWE/MAF 质控。
* **00_3_chinamap_upload_prep.sh**（`1010/preprocess/05_prepare_imputation.sh`）：整理、重命名染色体并拆分 upload-ready VCF；ChinaMAP phase 1 填充在网页端运行，下载结果作为下游 input checkpoint。
* **00_4_stage_variant_filtering.sh**（`1010/preprocess/00_4_stage_variant_filtering.sh`）：在 genotype cohort 的产前（224 人）和产后（603 人）样本中分别过滤 MAF $\ge 0.05$、MAC $\ge 20$、HWE $P \ge 10^{-6}$，取两个阶段共有的 5,014,001 个位点并生成统一 PGEN；该脚本对应保留的 process log，`RAW_VCF` 指向 ChinaMAP 下载后的 genotype input。
* **00_5_stage_genotype_preparation.R**：从上述统一 PGEN 中提取通过 phenotype/sample QC 的产前（189 人）和产后（508 人）stage-specific PGEN，供 apaQTL/eQTL mapping 使用。
* **00_6_genotype_pca.R**：PLINK 2.0 LD 剪枝（`--indep-pairwise 50 5 0.2`） $\to$ 计算前 20 个阶段特异 Genotype PCs。

---

#### **Module 01: RNA-seq 表达量与 APA 定量及校正**

*转录组清洗、3′UTR 事件定量以及核心的 Bridge-batch 批次校正。*

* **External checkpoint: RNA-seq read filtering**：测序公司使用 SOAPnuke v1.5.6 完成 reads filtering 并生成 HTML 质控报告；本机未保留原始 FASTQ、报告和 exact command，因此只在 Methods 中记录该步骤，不补写独立 SOAPnuke script。
* **01_1_rnaseq_star_rsem.sh**（`1010/preprocess/rna_quant_v2.sh`）：FastQC $\to$ STAR 2.7/RSEM 1.3 基因和转录本定量；论文主分析使用 GENCODE v40 reference，并从配对的 `.genes.results`/`.isoforms.results` 一次导出 sample-aligned gene/isoform count 与 TPM 四张矩阵。merge 同时写 source-to-matrix sample manifest、统一 5 个经 RSEM count identity 验证的 aliases，并拒绝 duplicate sample IDs 或不一致的 feature order；retained-matrix readers 使用同一组 verified mappings，避免旧 header 静默丢样本。
* **01_2_expression_qc_normalization.R**（`1010/preprocess/rna_pre_step01_v2.r`）：读入 retained gene matrices 时将已验证的第二个 `XYA20171125` column 区分为 `20171125C`，并验证 count/TPM 的 sample 和 feature order；合并两个 LifeStage 后做 RIN $\ge 5$ 过滤，使用 `set.seed(2026)` 和 *XIST*/*UTY* 的 log1p 两中心 k-means 自适应阈值进行性别一致性检查，随后做 TMM 归一化（edgeR）和分 LifeStage 的 WGCNA 网络连通性离群样本剔除；ComBat 模型保护 LifeStage，并在存在时以 batch9 为 reference。
* **01_3_dapars2_quantification.py/sh**（`1010/preprocess/apa_quant/run_dapars2.sh`）：串联 coordinate-sorted RNA-seq BAM 的四列 coverage 和 `flagstat` read depth、3′UTR/config preparation、DaPars2 常染色体并行运行与 chromosome merge。源码来自 `${APA_HPC_HOME}/project/work/ref_repo/DaPars2` commit `fb81c6c`（README 标记 v2.1，`git describe` 为 `v2.1-14-gfb81c6c`）；事件要求前 99 个位置平均 coverage 严格 `>15` 的 profiles 超过 30%，且 3′UTR span $\ge 150$ bp。
* **01_4_pdui_imputation_qc.R**：低质量事件过滤 $\to$ 软阈值 SVD（softImpute） $\to$ 经验贝叶斯矩阵分解（flashier + ebnm）填补缺失值 $\to$ 马氏距离排除多变量离群样本（保留 815 profiles / 750 人）。
* **01_5_batch_correction_evaluation.R**（`1010/batch/00_batch_correction_v2.R`）：测试 8 种批次校正策略（ComBat、limma、RUV-III 等）；不使用 age-derived metrics，仅在 corrected candidates 中按 4 个等权 domains 选择：bridge concordance（Pearson、CCC）、bridge error（MAE、RMSE）、batch removal（residual coefficient、batch variance、batch ASW）和 biological preservation（LifeStage variance、LifeStage ASW）。相关 metrics 先在 domain/subcomponent 内合并，Raw 仅作 baseline；综合结果 prioritize **ComBat_REF**，并将选中矩阵显式导出到 downstream APA TSV。*(对应 ED Fig. 2)*
* **01_6_bridge_method_selection.R**（`1010/batch/01_bridge_method_selection.R`）：在 65 个 bridge donor pairs 上评估所有 correction methods，并以 10,000 次 donor-paired bootstrap 直接比较 ComBat_REF 与 limma。Pearson、CCC 和 MAE 支持 ComBat_REF，RMSE 小幅支持 limma；该 trade-off 在 diagnostics 中披露，不再要求完整复制 limma downstream pipeline。

---

#### **Module 02: 生命轨迹建模、聚类与生物学功能（服务于 Fig. 1）**

* **02_0_donor_level_input.R**（`1011/02_lifespan_trajectories/pipeline/01_cov.R`）：在 ComBat_REF 批次校正后，将 65 对 technical-replicate profiles 在 logit-PDUI scale 上按 donor 取均值，同时对 RIN 取 donor mean；生成 750 个 independent donors 的 APA matrix 和 covariate table。Replicate profiles 仅保留给上游 batch/bridge QC，不作为 lifespan inference 的独立观测。
* **02_1_supported_age_mask.R**：基于滑动窗口密度评估年龄采样充分区间（Supported age range：GW14 至产后 2 个月、11 至 106 岁）。*(对应 Fig. 1b)*
* **02_2_gam_trajectory_fitting.R**：在 donor-level matrix 上使用 mgcv 广义相加模型（薄板样条拟合 transformed age，校正性别、donor-mean RIN），并对年龄效应做 BH-FDR；修订后的事件数待主流程重跑。
* **02_3_permutation_directional_events.R**：在 donor level 做年龄标签置换（1,000 次），标定端点位移阈值并筛选定向事件；修订后的阈值和事件数待主流程重跑。*(对应 Fig. 1c)*
* **02_4_trajectory_clustering.R**：Z-score 标准化轨迹 PCA 降维 $\to$ k-means 聚类识别 **C1（远端伸长）** 与 **C2（近端缩短）**。*(对应 Fig. 1c/d)*
* **02_5_external_validation_brainspan.R**：BrainSpan (n=26) 队列的数据比对、liftOver 坐标转换与轨迹拟合验证。*(对应 ED Fig. 3)*
* **02_6_gene_classification_enrichment.R**：基因三分类（C1-specific, C2-specific, Shared） $\to$ 控制长度/表达量后的通路富集（GO/KEGG/Autism）与 gnomAD LOEUF 约束分析。*(对应 Fig. 1d/e/f)*

---

#### **Module 03: 多组学解耦与蛋白质水平功能验证（服务于 Fig. 2）**

* **03_1_apa_mrna_residual_uncoupling.R**：校正年龄、性别、RIN 后，计算 APA 残差与基因/转录本表达量残差的加权相关性（证实 $\rho \approx 0$）。*(对应 Fig. 2a)*
* **03_2_proteomics_preprocessing.R**：DIA 蛋白组 MaxQuant 结果处理、离群过滤、缺失值填补、sva 校正（40 个 SVs）。*(对应 ED Fig. 5)*
* **03_3_apa_protein_association.R**：在 366 例匹配样本中建模：`Protein ~ APA + Age + mRNA`，筛选 10 个基因中 13 对显著相关 pairs，计算 $\Delta R^2$ 增量。*(对应 Fig. 2b/c/d)*
* **03_4_mirna_seed_matching_inpp4a.py/R**：*INPP4A* 延伸 3′UTR 区域提取 $\to$ miRBase 种子匹配 $\to$ miRNATissueAtlas2 脑高表达 miRNA 排序与位点标注。*(对应 Fig. 2e)*

---

#### **Module 04: 阶段特异 apaQTL 绘图与保守性分析（服务于 Fig. 3）**

* **04_1_celltype_deconvolution.R**：BisqueRNA 结合单核参考集计算 7 类脑细胞比例 $\to$ ILR 转换 $\to$ 偏 F 检验排除细胞混杂。*(对应 ED Fig. 6)*
* **04_2_phenotype_pc_selection.R**：PCAForQTL (runBE) 自动选择表型主成分（产前 21 PCs，产后 37 PCs）。
* **04_3_stage_specific_tensorqtl.py**：tensorQTL 运行产前（n=189）与产后（n=508）顺式 apaQTL（10,000 次置换，Storey $q < 0.05$）。*(对应 Fig. 3a)*
* **04_4_eqtl_mapping_coloc_mediation.R**：同批样本 eQTL 计算 $\to$ coloc.abf 检验 apaQTL 与 eQTL 共定位；对 PP4 $\ge 0.8$ pairs 运行自定义 OLS 正向/反向中介模型和 2,000 次 donor bootstrap。*(对应 Fig. 3b)*
* **04_5_motif_disruption_hid1.R**：PAS 核心基序（AAUAAA 及其变体）扫描与突变对齐，绘制 *HID1* 位点基因型-PDUI 关联。*(对应 Fig. 3c)*
* **04_6_cross_stage_sharing_mashr.R**：提取共享显著 pairs 效应相关性（$\rho = 0.88$） $\to$ 运行 mashr 评估效应量跨期波动幅度。*(对应 Fig. 3d & ED Fig. 7)*
* **04_7_downsampling_power_test.R**：产后 189 样本（20 次随机重抽样）下采样测试，评估产前特异位点恢复率。*(对应 ED Fig. 7)*
* **04_8_braingvex_replication.R**：BrainGVEX (n=375 欧洲人) apaQTL 计算与 $\pi_1$ 复现率评估。*(对应 ED Fig. 6)*

---

#### **Module 05: 基因型-发育时序交互与机制解析（服务于 Fig. 4）**

* **05_1_gxstage_interaction_tensorqtl.py**：全样本（n=697）拟合交互模型 `APA ~ G + Stage + G:Stage + Covs`（44 个表型 PCs，eigenMT 多重检验校正）。
* **05_2_hc3_robust_refitting.R**：对 lead pairs 采用 OLS + HC3 异方差稳健标准误重拟合，锁定 193 个稳健交互基因。*(对应 Fig. 4a/b)*
* **05_3_expression_gating_test.R**：计算跨出生表达量变化（$\Delta\text{expr}$）与效应变化（$\Delta\text{effect}$）相关性，比较放大组与衰减组表达差异。*(对应 Fig. 4c)*
* **05_4_representative_gating_loci.R**：提取并绘制 *TESPA1*、*SLCO1A2* 产前不表达、产后开启的双层表达与 QTL 图。*(对应 Fig. 4d)*
* **05_5_saturation_analysis_atp2b1.R**：*ATP2B1* 产后截断饱和统计（559/593 人打满），绘制产前充分分布 vs 产后饱和平线的对比箱线图。*(对应 Fig. 4e)*

---

#### **Module 06: 全基因组与基因集疾病关联分析（服务于 Fig. 5）**

* **06_1_gwas_harmonization.sh/R**：19 个 GWAS 摘要统计量统一 liftOver 到 GRCh38，对齐效应等位基因与链反转。
* **06_2_sldsc_heritability_enrichment.sh**：
* 并行保留 APA-only 分支（prenatal/postnatal/interaction/lifespan APA 注释）与 joint apaQTL+eQTL maxCPP 分支（prenatal/postnatal、各自独立注释列）；
* 两个分支均基于 1000G 种族匹配面板计算 LD scores，并结合 Baseline-LD v2.2 模型计算疾病遗传度富集倍数及显著性（$z \ge 7$ 疾病）。*(对应 Fig. 5a)*


* **06_3_magma_geneset_analysis.sh**：MAGMA 1.10 进行产前、产后、交互基因集的竞争性富集检验（条件分析排除背景基因影响）。*(对应 Fig. 5b)*
* **06_4_twas_spredixcan.py/R**：弹性网（glmnet）训练产前与产后 APA 预测模型 $\to$ S-PrediXcan 计算各性状关联。*(对应 Fig. 5c)*

---

#### **Module 07: 位点级共定位、SMR 与关键致病位点绘图（服务于 Fig. 5 & 6）**

* **07_1_coloc_screening_susie.R**：
* PLINK 剪枝定义 GWAS risk loci（$P < 10^{-7}, r^2 < 0.2$）；
* 运行 coloc.abf 与 coloc.susie（考虑多因果位点），筛选 104 个信号与 72 个基因；
* 筛查 GTEx 13 个脑区及 CHB eQTL，标记出 40% 的 APA-only 基因。*(对应 Fig. 5c 概览)*


* **07_2_smr_heidi_analysis.sh**：SMR 1.3.2 测试 GWAS 与 apaQTL 的多效性关联并执行 HEIDI 异质性检验。
* **07_3_plot_prenatal_scz_loci.R**：绘制产前独家、无 eQTL 支持的硬核致病位点（*SCAF1* 与 *NRBP2*）的区域关联与对比图。*(对应 Fig. 5d/e)*
* **07_4_plot_stage_dependent_atp2b1.R**：绘制 *ATP2B1* 产前特异 apaQTL 与 ADHD/SCZ 共定位图，并附皮层 eQTL 阴性对照。*(对应 Fig. 6a)*
* **07_5_plot_pleiotropy_12q24.R**：绘制 12q24.31 位点 *MTRFR* (SCZ) 与 *OGFOD2* (智力) 的双轨解绑 LocusZoom。*(对应 Fig. 6b)*
* **07_6_plot_shared_locus_arl3.R**：绘制 10q24.32 位点 *ARL3* 的 APA 与 eQTL 协同共定位图。*(对应 Fig. 6c)*

---

### 整理建议：

您可以直接在您的分析机里建好这些 `00_` 到 `07_` 的文件夹。在拼图时，需要哪张图的源数据（Source Data）或重绘某个 panel，直接去对应的模块里找脚本，整套代码库就会井井有条，完全符合顶刊接收后开源代码（Code Availability）的标准！
