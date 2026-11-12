# Fig 1e enrichment hierarchy

GO/KEGG terms use Fisher `fdr_by_gene_set < 0.05` for the main clusters. Clustering uses the 11,322 genes with complete covariates, foreground-intersected gene-set Jaccard distance, average linkage, and a 0.7 tree cut. The TSVs also report cuts at 0.6 and 0.8; separate FDR < 0.10 tables are provided. `n_fg` is the original Fisher count; `n_fg_covariate_bg` is the clustering overlap count.

GO slim mapping uses GO.db ancestors and the downloaded generic slim; terms without a generic slim ancestor remain explicitly unmapped. KEGG mapping uses the downloaded BRITE br08901 hierarchy. The 2024 MSigDB possessive disease names are normalized to current KEGG titles.

## C1spec: 3 terms, 2 clusters

The C1-specific terms separate translation/ribosome overlap from a GPCR activity singleton.

- `C1spec_01` (2 terms): These 2 overlapping terms center on cytoplasmic translation, with frequent title words: cytoplasmic, ribosome, translation.
- `C1spec_02` (1 term): This singleton is the g protein coupled receptor activity term.

## C2spec: 27 terms, 12 clusters

The C2-specific clusters include respiration/energy production, nucleotide-related metabolism, and separate proteostasis terms.

- `C2spec_01` (9 terms): These 9 overlapping terms center on generation of precursor metabolites and energy, with frequent title words: disease, energy, oxidative.
- `C2spec_02` (7 terms): These 7 overlapping terms center on nucleobase containing small molecule metabolic process, with frequent title words: atp, containing, nucleoside.
- `C2spec_03` (1 term): This singleton is the metallo sulfur cluster assembly term.
- `C2spec_04` (1 term): This singleton is the mitochondrial translation term.
- `C2spec_05` (1 term): This singleton is the mitochondrion organization term.
- `C2spec_06` (1 term): This singleton is the peptidyl proline modification term.
- `C2spec_07` (2 terms): These 2 overlapping terms center on isomerase activity, with frequent title words: intramolecular, isomerase, oxidoreductase.
- `C2spec_08` (1 term): This singleton is the ligase activity term.
- `C2spec_09` (1 term): This singleton is the unfolded protein binding term.
- `C2spec_10` (1 term): This singleton is the amino sugar and nucleotide sugar metabolism term.
- `C2spec_11` (1 term): This singleton is the glycolysis gluconeogenesis term.
- `C2spec_12` (1 term): This singleton is the proteasome term.

## Shared: 71 terms, 36 clusters

The Shared clusters span cytoskeleton/junction organization, transport and hormone response, with several neural-development terms.

- `Shared_01` (8 terms): These 8 overlapping terms center on regulation of cytoskeleton organization, with frequent title words: organization, actin, filament.
- `Shared_02` (3 terms): These 3 overlapping terms center on l amino acid transport, with frequent title words: acid, amino, transport.
- `Shared_03` (6 terms): These 6 overlapping terms center on response to hormone, with frequent title words: response, hormone, peptide.
- `Shared_04` (2 terms): These 2 overlapping terms center on cell junction organization, with frequent title words: cell, junction, assembly.
- `Shared_05` (4 terms): These 4 overlapping terms center on cell projection morphogenesis, with frequent title words: cell, development, morphogenesis.
- `Shared_06` (1 term): This singleton is the cerebral cortex development term.
- `Shared_07` (1 term): This singleton is the commissural neuron axon guidance term.
- `Shared_08` (1 term): This singleton is the cytokinesis term.
- `Shared_09` (1 term): This singleton is the cytosolic transport term.
- `Shared_10` (4 terms): These 4 overlapping terms center on post golgi vesicle mediated transport, with frequent title words: golgi, membrane, plasma.
- `Shared_11` (1 term): This singleton is the heart development term.
- `Shared_12` (1 term): This singleton is the maintenance of location term.
- `Shared_13` (3 terms): These 3 overlapping terms center on mitotic cell cycle process, with frequent title words: cell, cycle, mitotic.
- `Shared_14` (5 terms): These 5 overlapping terms center on kinase activity, with frequent title words: kinase, serine, modification.
- `Shared_15` (1 term): This singleton is the positive regulation of mitotic cell cycle term.
- `Shared_16` (2 terms): These 2 overlapping terms center on regulation of transmembrane transport, with frequent title words: transmembrane, transport.
- `Shared_17` (1 term): This singleton is the positive regulation of transport term.
- `Shared_18` (1 term): This singleton is the process utilizing autophagic mechanism term.
- `Shared_19` (2 terms): These 2 overlapping terms center on protein localization to cell periphery, with frequent title words: localization, cell, membrane.
- `Shared_20` (1 term): This singleton is the regulation of dendritic spine morphogenesis term.
- `Shared_21` (1 term): This singleton is the regulation of microtubule based process term.
- `Shared_22` (1 term): This singleton is the regulation of pinocytosis term.
- `Shared_23` (2 terms): These 2 overlapping terms center on regulation of protein containing complex disassembly, with frequent title words: complex, containing, depolymerization.
- `Shared_24` (1 term): This singleton is the response to growth factor term.
- `Shared_25` (1 term): This singleton is the stem cell differentiation term.
- `Shared_26` (1 term): This singleton is the actin binding term.
- `Shared_27` (2 terms): These 2 overlapping terms center on cell adhesion molecule binding, with frequent title words: binding, adhesion, cadherin.
- `Shared_28` (1 term): This singleton is the calmodulin binding term.
- `Shared_29` (2 terms): These 2 overlapping terms center on transcription factor binding, with frequent title words: binding, factor, transcription.
- `Shared_30` (1 term): This singleton is the histone h3 methyltransferase activity term.
- `Shared_31` (3 terms): These 3 overlapping terms center on phospholipid binding, with frequent title words: binding, phosphatidylinositol, phosphate.
- `Shared_32` (1 term): This singleton is the protein domain specific binding term.
- `Shared_33` (1 term): This singleton is the adherens junction term.
- `Shared_34` (2 terms): These 2 overlapping terms center on focal adhesion, with frequent title words: actin, adhesion, cytoskeleton.
- `Shared_35` (1 term): This singleton is the pancreatic cancer term.
- `Shared_36` (1 term): This singleton is the pathways in cancer term.

The clusters are descriptive overlap groups; singleton terms and GO/KEGG pathway titles are not independent mechanistic confirmations.

## Broad gene sets

T3 uses the 11,322-gene complete-covariate universe for one-sided Fisher tests and logistic models with `scale(log1p(n_sites))`, `scale(log1p(gene_length))`, and `scale(mean_expression)`. Fisher P values are one-sided; logistic P values and Wald 95% CIs are two-sided. BH FDR is calculated separately for all 18 Fisher tests and all 18 adjusted tests.

SynGO input is the supplied version 1.3 archive (the original request named 1.2). The all-synapse set is the union of symbols in `genes.xlsx` and `annotations.xlsx`; pre/post sets use annotations for the specified GO roots and all descendants in `ontologies.xlsx`. The supplied HGNC S file duplicated the L group, so RPS uses a fresh HGNC group 728 snapshot in this handoff.

- C1spec: Fisher enrichment in HGNC_RPL_RPS; adjusted enrichment in HGNC_RPL_RPS; adjusted depletion in MitoCarta3, HGNC_Proteasome (FDR < 0.05).
- C2spec: Fisher enrichment in MitoCarta3, HGNC_Proteasome; adjusted enrichment in MitoCarta3, HGNC_Proteasome; adjusted depletion in HGNC_RPL_RPS (FDR < 0.05).
- Shared: Fisher enrichment in SynGO_all, SynGO_presynapse, SynGO_postsynapse; adjusted enrichment in none; adjusted depletion in SynGO_all, SynGO_postsynapse, HGNC_RPL_RPS (FDR < 0.05).

Shared SynGO_all changes from Fisher enrichment to adjusted depletion; these are different estimands and should not be presented as concordant evidence.

