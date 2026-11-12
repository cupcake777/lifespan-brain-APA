#!/usr/bin/env python
# Recompute coloc-dependent manuscript numbers from one coloc run.
# usage: ms_numbers.py <all_colocalization_results.csv> <eqtl_classification.tsv> <out_prefix>
# eqtl_classification.tsv: one row per GWAS x stage x gene signal; columns gene, GWAS, stage, chb_eqtl_max_pp4, gtex_max_pp4 (NaN = not testable).
import sys, json, numpy as np, pandas as pd

# --- inputs ---
coloc, cls_path, out = sys.argv[1:4]
THR = 0.7
co = pd.read_csv(coloc, low_memory=False)
co['gene'] = co.Phenotype.str.split('|').str[1]
x = co.Phenotype.str.extract(r'\|chr(\w+):(\d+)-'); co['chr'] = x[0]; co['pos'] = pd.to_numeric(x[1])

# --- signals: unique GWAS x QTL stage x gene with apaQTL PP4 >= 0.7 ---
sig = co[co.PP4 >= THR].copy()
xi = sig.plot_file.astype(str).str.extract(r'_chr(\w+?)_(\d+)_rs'); sig['ichr'] = xi[0]; sig['ipos'] = pd.to_numeric(xi[1])
S = (sig.sort_values('PP4', ascending=False).drop_duplicates(['GWAS_ID', 'QTL_ID', 'gene'])
        .rename(columns={'GWAS_ID': 'GWAS', 'QTL_ID': 'stage', 'PP4': 'apa_pp4'}))
res = {'signals': len(S), 'genes': S.gene.nunique(), 'gwas': S.GWAS.nunique(), 'rows_pp4_ge_0.7': len(sig)}

# --- join eQTL classification ---
cls = pd.read_csv(cls_path, sep='\t')
m = S.merge(cls[['gene', 'GWAS', 'stage', 'chb_eqtl_max_pp4', 'gtex_max_pp4']], on=['gene', 'GWAS', 'stage'], how='left')
res['signals_without_classification_row'] = int(m.chb_eqtl_max_pp4.isna().__and__(m.gtex_max_pp4.isna()).sum())
m['eqtl_max'] = m[['chb_eqtl_max_pp4', 'gtex_max_pp4']].max(axis=1)
m['signal_class'] = np.where(m.eqtl_max.isna(), 'NA', np.where(m.eqtl_max >= THR, 'Shared', 'apaQTL-only'))

# --- APA-only genes: have eQTL data and max eQTL PP4 < 0.7 ---
G = m.groupby('gene').agg(has=('eqtl_max', lambda v: v.notna().any()), emax=('eqtl_max', 'max'))
apa_only = G[G.has & (G.emax < THR)].index
res.update({'genes_with_eqtl_data': int(G.has.sum()), 'apa_only_genes': len(apa_only),
            'apa_only_pct': round(100 * len(apa_only) / max(G.has.sum(), 1), 1)})

# --- 1-Mb loci of APA-only gene signals, merged within GWAS ---
res['signals_missing_index_pos'] = int(m.ipos.isna().sum())
L = m[m.gene.isin(apa_only)][['GWAS', 'ichr', 'ipos']].drop_duplicates().rename(columns={'ichr': 'chr', 'ipos': 'pos'}).sort_values(['GWAS', 'chr', 'pos'])
n = 0
for (_, _), d in L.groupby(['GWAS', 'chr']):
    last = None
    for p in d.pos:
        if last is None or p - last > 1e6: n += 1
        last = p
res['apa_only_loci_1mb'] = n

# --- named loci ---
NAMED = [('SCAF1','EUR_SCZ'),('NRBP2','EUR_SCZ'),('ATP2B1','EUR_ADHD'),('ATP2B1','EAS_SCZ'),('MTRFR','EAS_SCZ'),('MTRFR','EUR_SCZ'),
         ('OGFOD2','EUR_INT'),('ARL3','EAS_SCZ'),('ARL3','EUR_SCZ'),('NDUFS3','EUR_AD'),('NPAS1','EUR_INT')]
named = {f'{g}x{t}': (co[(co.gene == g) & (co.GWAS_ID == t)].groupby('QTL_ID').PP4.max().round(3).to_dict() or 'not tested') for g, t in NAMED}
res['named'] = named
res['per_gwas_genes'] = S.groupby('GWAS').gene.nunique().to_dict()

# --- outputs ---
m.to_csv(out + '_signals.tsv', sep='\t', index=False)
json.dump(res, open(out + '_numbers.json', 'w'), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
