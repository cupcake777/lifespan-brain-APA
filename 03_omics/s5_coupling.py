import os
import numpy as np
import pandas as pd
import patsy
from scipy.stats import spearmanr

# APA vs expression coupling per site cluster, along the trajectory stages.
# Reads data/; writes the cluster summary to data/03_omics/.

# --- metadata and clustered sites ---
root = os.environ['APA_ROOT'] + '/data/'
meta = pd.read_csv(root + '04_qtl/CHB.all.meta', sep='\t').set_index('sample')
sc = pd.read_csv(root + '02_traj/02_sites_clustered.tsv', sep='\t', usecols=['site_id', 'cluster_id', 'delta_logit'], low_memory=False).dropna(subset=['cluster_id']).set_index('site_id')

# --- stage segments ---
sup = pd.read_csv(root + '02_traj/01_support.tsv', sep='\t')
am = sup.age_model.values
ok = sup.is_supported.astype(bool).values
segs = np.split(np.arange(len(ok)), np.where(np.diff(ok.astype(int)) != 0)[0] + 1)
segs = [s for s in segs if ok[s[0]]]
e, l = segs[0], segs[1]

# --- APA model ---
Apa = pd.read_csv(root + '04_qtl/apa_after_combat.tsv', sep='\t', index_col=0)
Apa = Apa.loc[sc.index]
ma = meta.reindex(Apa.columns)
ka = (ma.age.notna() & ma.RIN.notna()).values
Apa = Apa.loc[:, ka]
ma = ma[ka]

xa = np.log10(np.maximum(ma.age.values + 1, 1e-6))
sa = (ma.sex.values == 'Male').astype(float)
dmA = patsy.dmatrix("cr(x, df=5) + sexb + RIN", {'x': xa, 'sexb': sa, 'RIN': ma.RIN.values.astype(float)}, return_type='dataframe')
Y = Apa.values.T.astype(float)
bA = np.linalg.lstsq(dmA.values, Y, rcond=None)[0]
gA = np.asarray(patsy.build_design_matrices([dmA.design_info], {'x': am, 'sexb': np.full(len(am), sa.mean()), 'RIN': np.full(len(am), ma.RIN.astype(float).mean())})[0])
PA = gA @ bA
rA = (Y - dmA.values @ bA).std(axis=0, ddof=dmA.shape[1])

# --- expression model ---
E = pd.read_csv(root + '01_quant/expr_after_combat.tsv', sep='\t', index_col=0)
E.index = E.index.str.replace(r'\..*$', '', regex=True)
m = meta.reindex(E.columns)
keep = m.age.notna() & m.RIN.notna()
E = E.loc[:, keep]
m = m[keep]

am_s = np.log10(np.maximum(m.age.values + 1, 1e-6))
sexb = (m.sex.values == 'Male').astype(float)
dm = patsy.dmatrix("cr(x, df=5) + sexb + RIN", {'x': am_s, 'sexb': sexb, 'RIN': m.RIN.values.astype(float)}, return_type='dataframe')
beta = np.linalg.lstsq(dm.values, E.values.T, rcond=None)[0]
grid = np.asarray(patsy.build_design_matrices([dm.design_info], {'x': am, 'sexb': np.full(len(am), sexb.mean()), 'RIN': np.full(len(am), m.RIN.astype(float).mean())})[0])
P = grid @ beta
rsd = (E.values.T - dm.values @ beta).std(axis=0, ddof=dm.shape[1])

# --- gene symbols ---
gagg = pd.read_csv(root + '02_traj/05_gene_aggregation.tsv', sep='\t')
sym2ens = dict(zip(gagg.gene_symbol, gagg.gene_id_base))

# --- event table and gene join ---
ev = pd.DataFrame({'site': sc.index, 'cl': sc.cluster_id.values,
    'apa_early': (PA[e[-1]] - PA[e[0]]) / rA,
    'apa_gap': (PA[l[0]] - PA[e[-1]]) / rA,
    'apa_late': (PA[l[-1]] - PA[l[0]]) / rA})
ev['gene'] = ev.site.str.split('|').str[1].map(sym2ens)

exg = pd.DataFrame({'gene': E.index,
    'ex_early': (P[e[-1]] - P[e[0]]) / rsd,
    'ex_gap': (P[l[0]] - P[e[-1]]) / rsd,
    'ex_late': (P[l[-1]] - P[l[0]]) / rsd})

j = ev.merge(exg, on='gene', how='inner')
j['apa_dev'] = j.apa_early + j.apa_gap
j['ex_dev'] = j.ex_early + j.ex_gap

# --- per-cluster summary ---
cp = pd.DataFrame([dict(
    cluster=f'C{int(c)}',
    n_events=int((j.cl == c).sum()),
    rho_dev_dPDUI_vs_dExpr=round(spearmanr(j.apa_dev[j.cl == c], j.ex_dev[j.cl == c]).correlation, 3),
    median_abs_late_change_APA_resSD=round(j.apa_late[j.cl == c].abs().median(), 3),
    median_abs_late_change_expr_resSD=round(j.ex_late[j.cl == c].abs().median(), 3),
    frac_events_APA_more_stable_late=round((j.apa_late[j.cl == c].abs() < j.ex_late[j.cl == c].abs()).mean(), 3))
    for c in [1.0, 2.0]])

cp.to_csv(root + '03_omics/S5_apa_vs_expression_coupling.tsv', sep='\t', index=False)