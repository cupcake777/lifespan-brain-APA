# Interaction-gene and APA-only risk-gene context

## A. Interaction-gene features

The tested background contained 13,566 unique genes from the interaction model. The HC3-core foreground contained 193 genes: 77 interaction-only, 15 prenatal-only, 59 postnatal-only, and 42 detected in both stage-stratified analyses. The stage-stratified apaQTL sets contained 669 prenatal and 2,249 postnatal genes; 563 were significant in both stages, leaving 521 shared apaQTL comparator genes after excluding the 42 interaction genes.

LOEUF used the gnomAD v4.1.1 MANE-preferred, otherwise canonical-preferred, within-gene median. Models adjusted for `log1p(number of directional events)`, gene length, and weighted mean expression, with continuous covariates scaled as in the Gene-level architecture implementation.

| Comparison | Complete target/reference | Adjusted beta | 95% CI | P |
|---|---:|---:|---:|---:|
| 193 interaction genes vs shared apaQTL comparator | 151 / 417 | -0.0186 | -0.0548 to 0.0176 | 0.314 |
| 77 interaction-only vs 116 other interaction genes | 66 / 85 | -0.00115 | -0.0609 to 0.0586 | 0.970 |

No GO BP, GO MF, KEGG, GWAS Catalog, or HPO term passed the within-library Fisher FDR < 0.10 candidate screen, so no term entered the covariate-adjusted logistic confirmation step and the requested FDR < 0.05 table is empty. No prenatal or postnatal seven-class cell-marker set passed within-stage BH FDR < 0.05 for either foreground.

Conclusion: interaction genes, including the interaction-only subset, did not show detectable differential constraint, functional enrichment, or broad cell-class marker enrichment under the prespecified models and multiple-testing families.

## B. APA-only risk genes

The 151 PP4 >= 0.7 rows collapsed to the prespecified 104 unique `GWAS x stage x APA event` signals spanning 72 genes. For each signal, the harmonized GWAS was restricted to the index position +/-500 kb and tested against the same gene's stage-matched CHB nominal eQTLs and GTEx v11 nominal eQTLs in 13 brain tissues using `coloc.abf` with `p1 = p2 = 1e-4` and `p12 = 1e-5`. CHB used prenatal N = 182 and postnatal N = 497; the five interaction signals were evaluated in both CHB stages and assigned the larger PP4.

- 40/104 signal rows were APA-only and 64/104 were APA + eQTL after gene-within-GWAS classification.
- 37/87 gene-GWAS pairs were APA-only and 50/87 were APA + eQTL.
- 29/72 genes were APA-only across all of their matched GWAS; these involved 27 distinct GWAS index loci.
- The 40 APA-only signal rows represented 29 distinct GWAS index loci before requiring the gene to be APA-only across every matched GWAS.
- CHB was evaluable for 102/104 signals. `Metazoa_SRP` and `RN7SL331P` were absent from the CHB nominal tables and retain missing CHB PP4 rather than negative evidence.
- 36/104 signals had CHB PP4 >= 0.7. CHB added qualifying evidence for 13 genes that were APA-only under GTEx alone: ARL3, C2orf69, CKMT1A, DENND4A, ENSG00000262202, ENSG00000277701, H4C12, LNPEP, PBRM1, RBL2, SBNO1, SNORD3C, and WDR46.
- GTEx was evaluable for 96/104 signals. Eight signals had no corresponding gene rows in any of the 13 tissues and retain missing GTEx PP4.

Conclusion: after adding stage-matched CHB and 13-tissue GTEx brain eQTL evidence, 29 of 72 APA-colocalized genes remained nominated only by APA across all matched GWAS, involving 27 distinct GWAS index loci.

The optional genome-wide CHB eQTL-GWAS comparison was not run in this targeted update.
