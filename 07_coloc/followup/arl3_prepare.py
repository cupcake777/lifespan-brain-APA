# Builds the ARL3 x EAS SCZ coloc inputs from the Fig 6c locus JSONs and ARL3 coloc tables.
# Writes in_apaqtl_eassc.tsv and in_eqtl_regional_maf.tsv to data/07_coloc/.
import json, os
import numpy as np, pandas as pd
D = f'{os.environ["APA_ROOT"]}/data/07_coloc/'
# --- apaQTL points ---
rows = []
for stage, fn in [('postnatal', 'Figure6c_ARL3_locuscompare.json'),
                  ('prenatal', 'Fig6c_ARL3_pre_locuscompare.json')]:
    for p in json.load(open(D + fn))['points']:
        rows.append(dict(stage=stage, snp=p['snp'], pos=p['pos'], gwas_logp=p['gwas_logp'], apa_logp=p['qtl_logp']))
pts = pd.DataFrame(rows)
sv = pd.read_csv(D + 'arl3_scz_validation_input.tsv', sep='\t')
bd = pd.read_csv(D + 'arl3_bd_coloc_input.tsv', sep='\t')
apa = pts.merge(sv[['snp', 'stage', 'af_apa']], on=['snp', 'stage'], how='inner')
# --- eQTL MAF ---
eq = pd.read_csv(D + 'Fig6c_ARL3_eqtl_reg.tsv', sep='\t')
eq = eq.merge(bd[['stage', 'rsid', 'af_exp']], on=['stage', 'rsid'], how='left')
af_map = apa.drop_duplicates(['snp', 'stage']).set_index(['snp', 'stage']).af_apa
eq['maf_eq'] = eq.af_exp
miss = eq.maf_eq.isna()
eq.loc[miss, 'maf_eq'] = [af_map.get((s, t), np.nan) for s, t in zip(eq.loc[miss, 'rsid'], eq.loc[miss, 'stage'])]
apa.to_csv(D + 'in_apaqtl_eassc.tsv', sep='\t', index=False)
eq.to_csv(D + 'in_eqtl_regional_maf.tsv', sep='\t', index=False)
print(apa.groupby('stage').size().to_dict(), eq.groupby('stage').size().to_dict(), eq.maf_eq.isna().sum())
