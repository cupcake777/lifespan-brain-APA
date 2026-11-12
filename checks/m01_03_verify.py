# Recomputes the module 02 numbers cited in the manuscript text from local tables.
# Writes rerun_values.json to data/02_traj/.

import os, re, gzip, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# --- helpers ---
def to_bool(x): return x.astype(str).str.lower().isin(['true','1'])
def rd_art(p): return pd.read_csv(p,sep='\t')
def z(x): return (x-x.mean())/x.std()
def kmeans1d(x,k=2,nstart=25,seed=0):
    rng=np.random.default_rng(seed); best=None
    for _ in range(nstart):
        cen=np.sort(rng.choice(x,k,replace=False)).astype(float)
        for _ in range(200):
            lab=np.argmin(np.abs(x[:,None]-cen[None,:]),axis=1)
            new=np.array([x[lab==j].mean() if np.any(lab==j) else cen[j] for j in range(k)])
            if np.allclose(new,cen): break
            cen=new
        w=np.sum((x-cen[lab])**2)
        if best is None or w<best[0]: best=(w,cen.copy())
    return best[1]

# --- paths ---
BASE = os.environ['APA_ROOT']
D02 = BASE + '/data/02_traj/'
D03 = BASE + '/data/03_omics/'
RAW = BASE + '/raw/'
res = {}

# --- trajectory counts ---
s=pd.read_csv(D02+'02_sites.tsv',sep='\t',low_memory=False)
res['A_n_sites']=len(s); res['A_n_genes']=int(s.gene_id.nunique())
res['A_n_fdr05']=int((s.age_qvalue<0.05).sum()); res['A_n_q_na']=int(s.age_qvalue.isna().sum())
res['A_n_dynamic']=int(to_bool(s.dynamic_pass).sum()); res['A_n_union']=int(to_bool(s.dynamic_any_pass).sum())
res['A_n_transient']=int(to_bool(s.transient_pass).sum())
res['A_delta_thr']=float(s.delta_logit_threshold.iloc[0]); res['A_amp_thr']=float(s.amplitude_logit_threshold.iloc[0])
sc=pd.read_csv(D02+'02_sites_clustered.tsv',sep='\t',low_memory=False)
dy=sc[to_bool(sc.dynamic_pass)].copy()
res['A_cluster_sites']={int(k):int(v) for k,v in dy.cluster_id.value_counts().items()}
gc=dy.groupby('gene_id').cluster_id.agg(lambda x: tuple(sorted(set(int(v) for v in x))))
res['A_gene_classes']={'C1_specific':int((gc==(1,)).sum()),'C2_specific':int((gc==(2,)).sum()),'Shared':int((gc==(1,2)).sum()),'total_genes_directional':int(len(gc))}
res['A_loess_check']=None

# --- timing fractions from fitted curves ---
traj=pd.read_csv(D02+'traj_logit.csv').set_index('site_id')
sup=pd.read_csv(D02+'01_support.tsv',sep='\t').sort_values('grid_index')
gcols=[c for c in traj.columns if c.startswith('grid_')]
S=sup.is_supported.values.astype(bool); GA=sup.grid_age.values
dy_idx=dy.set_index('site_id')
M=traj.loc[dy_idx.index, gcols].values[:,S]; ages=GA[S]
start=M[:,0]; end=M[:,-1]; disp=end-start
res['B_endpoint_vs_delta_maxabs']=float(np.nanmax(np.abs(disp-dy_idx.delta_logit.values)))
cl=dy_idx.cluster_id.values
i11=int(np.argmin(np.abs(ages-11.0))); i2m=int(np.argmin(np.abs(ages-2/12)))
for c in [1,2]:
    m=cl==c; mc=np.median(M[m],axis=0); D=mc[-1]-mc[0]
    res[f'B_C{c}_frac_by_11y_medcurve']=float((mc[i11]-mc[0])/D)
    res[f'B_C{c}_frac_by_postnatal_m2_medcurve']=float((mc[i2m]-mc[0])/D)
    res[f'B_C{c}_grid_age_11y']=float(ages[i11])

# --- Fig 2a Spearman / Pearson R2 ---
f2=pd.read_csv(D02+'fig2a_event_coupling.tsv',sep='\t')
res['C_fig2a_rows']=len(f2)
for c in [1,2]:
    sub=f2[f2.cl==c]
    rho=spearmanr(sub.apa_dev,sub.ex_dev)[0]; r=pearsonr(sub.apa_dev,sub.ex_dev)[0]
    res[f'C_C{c}_n']=len(sub); res[f'C_C{c}_spearman']=float(rho); res[f'C_C{c}_pearson_R2']=float(r*r)

# --- LOEUF membership refit, BH over 5 coefficients ---
ld=pd.read_csv(D02+'ed4c_loeuf.tsv',sep='\t').dropna(subset=['loeuf'])
ld['y']=np.log1p(ld.loeuf); ld['shared']=(ld.membership=='Shared').astype(int)
ld['zS']=z(np.log1p(ld.n_sites)); ld['zL']=z(np.log1p(ld.gene_length)); ld['zE']=z(ld.mean_expression)
m=smf.ols('y ~ shared + zS + zL + zE',data=ld).fit()
pv=m.pvalues.values; q=multipletests(pv,method='fdr_bh')[1]
res['D_loeuf_n']=len(ld); res['D_shared_beta']=float(m.params['shared']); res['D_shared_p']=float(m.pvalues['shared'])
res['D_shared_fdr_bh5']=float(q[list(m.params.index).index('shared')])

# --- enrichment values, root vs outputs ---
eo=pd.read_csv(D02+'08_enrich_adjusted_outputs.tsv',sep='\t'); er=pd.read_csv(D02+'08_enrich_adjusted.tsv',sep='\t')
keys=[('GOBP_CYTOPLASMIC_TRANSLATION','C1-spec','Five-state'),('KEGG_RIBOSOME','C1-spec','Five-state'),
      ('GOBP_OXIDATIVE_PHOSPHORYLATION','C2-spec','Five-state'),('KEGG_PROTEASOME','C2-spec','Five-state'),
      ('Alzheimers Disease','Shared','Membership')]
res['E_enrich']={}
for t,gs,lv in keys:
    ro=eo[(eo.term==t)&(eo.gene_set==gs)&(eo.level==lv)]; rr=er[(er.term==t)&(er.gene_set==gs)&(er.level==lv)]
    res['E_enrich'][f'{t}|{gs}|{lv}']={'outputs_OR':float(ro.adjusted_odds_ratio.iloc[0]),'outputs_fdr':float(ro.adjusted_fdr.iloc[0]),
                                       'root_OR':float(rr.adjusted_odds_ratio.iloc[0]),'root_fdr':float(rr.adjusted_fdr.iloc[0])}

# --- S1/S2 artifact tables (pi1, nulls, element rate ratios) ---
s1r=rd_art(D03+'S1_apa_rna_cis_vs_trans_371.tsv').set_index('metric')
s1p=rd_art(D03+'S1_apa_prot_cis_vs_trans_371.tsv').set_index('metric')
res['F_pi1_rna']={'obs':float(s1r.loc['pi1','observed']),'null':float(s1r.loc['pi1','null_mean'])}
res['F_pi1_protein']={'obs':float(s1p.loc['pi1','observed']),'null':float(s1p.loc['pi1','null_mean'])}
res['F_protein_n_q05']={'obs':float(s1p.loc['n_q05','observed']),'null':float(s1p.loc['n_q05','null_mean'])}
res['F_protein_frac_neg']={'obs':float(s1p.loc['frac_neg','observed']),'null':float(s1p.loc['frac_neg','null_mean'])}
e1=rd_art(D03+'S2_aUTR_vs_cUTR_elem.tsv').set_index('motif')
res['F_aUTR_rr_AUadj']={k:float(e1.loc[k,'rr_AUadj']) for k in ['CPE_UUUUAU','PUM_UGUAHAUA','ARE_AUUUA']}
e2=rd_art(D03+'S2_C1C2_vs_nondynamic_elem.tsv')
res['F_C1C2_PUM']={g:float(e2[(e2.motif=='PUM_UGUAHAUA')&(e2.group==g)]['rate_ratio.Estimate'].iloc[0]) for g in ['C1','C2']}
e3=rd_art(D03+'S2_C1C2_vs_nondynamic_clip_adj.tsv')
res['F_C1C2_AGO2']={g:float(e3[(e3.rbp=='AGO2_brain')&(e3.group==g)]['rate_ratio'].iloc[0]) for g in ['C1','C2']}

# --- batch-9 and sample metadata ---
meta=pd.read_csv(RAW+'all.912.meta',sep='\t',dtype={'sample':str})
res['G_meta_samples']=len(meta); res['G_RIN_ge5']=int((meta.RIN>=5).sum())
b9=meta[(meta.batch=='batch9')&(meta.RIN>=5)]
res['G_batch9_RIN5']={'n':len(b9),'prenatal':int((b9.age<0).sum()),'postnatal':int((b9.age>=0).sum()),'repeat_libraries':int(b9['sample'].str.contains('_rep').sum())}
res['G_prenatal_RIN5']=int(((meta.age<0)&(meta.RIN>=5)).sum())

# --- expression filter (producer rna_pre_step01_v2.r logic) ---
gid_chr={}; gid_name={}
with gzip.open(BASE+'/ref/gencode.v40.annotation.gtf.gz','rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9 or f[2]!='gene': continue
        gu=re.search(r'gene_id "([^"]+)"',f[8]).group(1).split('.')[0]
        gid_chr[gu]=f[0]; nm=re.search(r'gene_name "([^"]+)"',f[8]); gid_name[gu]=nm.group(1) if nm else ''
ec=pd.read_csv(RAW+'all.912.expr.count',sep='\t',index_col=0).astype(np.float32)
et=pd.read_csv(RAW+'all.912.expr.tpm',sep='\t',index_col=0).astype(np.float32)
alias={'XYA20180516B':'XYA20180516','XYA20210813B':'XYA20210813','XYA20230626B':'XYA20230626','XYA20231026B':'XYA20231026','XYA20240412A':'XYA20240412'}
cols=[('20171125C' if c=='XYA20171125.1' else alias.get(c,c)) for c in ec.columns]
ec.columns=cols; et.columns=cols
mset=meta.set_index('sample')
common=[c for c in cols if c in mset.index]
res['H_expr_cols']=len(cols); res['H_cols_in_meta']=len(common); res['H_cols_not_in_meta']=[c for c in cols if c not in mset.index]
mm=mset.loc[common]; mm=mm[mm.RIN>=5]
res['H_RIN5']=len(mm)
bare=[g.split('.')[0] for g in ec.index]
xi=[i for i,g in enumerate(bare) if gid_name.get(g)=='XIST']; ui=[i for i,g in enumerate(bare) if gid_name.get(g)=='UTY']
cnt=ec[mm.index]
xe=cnt.iloc[xi].sum().values; ue=cnt.iloc[ui].sum().values
thx=np.expm1(kmeans1d(np.log1p(xe)).mean()); thu=np.expm1(kmeans1d(np.log1p(ue)).mean())
xpos=xe>thx; upos=ue>thu
imp=np.where(xpos&~upos,'female',np.where(upos,'male',''))
amb=(~xpos)&(~upos)
rep=mm.sex.map({1:'female',0:'male'}).values
mism=np.array([(r!='' and i!='' and r!=i) for r,i in zip(rep,imp)])
drop=set(mm.index[mism|amb])
keep=[s for s in mm.index if s not in drop]
res['H_sex_mismatch']=int(mism.sum()); res['H_ambiguous']=int(amb.sum()); res['H_after_sex_check']=len(keep)
c2=ec[keep].values; t2=et[keep].values
autos=np.array([gid_chr.get(g,'NA') in [f'chr{i}' for i in range(1,23)] for g in bare])
thr=0.2*len(keep)
keepg=((t2>0.1).sum(1)>=thr)&((c2>=6).sum(1)>=thr)&autos
res['H_genes_retained_producer_filter']=int(keepg.sum())
res['H_genes_retained_827_variant']=int((((et[mm.index].values>0.1).sum(1)>=0.2*len(mm.index))&((ec[mm.index].values>=6).sum(1)>=0.2*len(mm.index))&autos).sum())
with open(D02+'rerun_values.json','w') as fh:
    json.dump(res,fh,indent=1,default=str)
