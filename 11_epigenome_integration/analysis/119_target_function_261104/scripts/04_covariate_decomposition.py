#!/usr/bin/env python
"""119: which covariate removes the naive enrichment - GC, gene length, or both?
Gene-level CMH under four stratifications: none (= Fisher), GC quintile only,
length quintile only, GC x length (primary)."""
import os, resource
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
RS='/Users/okaban/bioinfo/rna-seq'
OUT=f'{RS}/11_epigenome_integration/analysis/119_target_function_261104/tables'
SETS=['prom_4mC','prom_6mA','body_4mC','body_6mA','bodycall_4mC','bodycall_6mA']
exec(open(os.path.join(os.path.dirname(__file__),'02_enrichment.py')).read().split('if __name__')[0])

ann=pd.read_csv(f'{OUT}/genes_universe_261104.tsv',sep='\t')
cats=pd.read_csv(f'{OUT}/categories_261104.tsv',sep='\t')
U=make_gene_units(); n=len(U['ids'])
idx={g:i for i,g in enumerate(U['ids'])}
keep,Cl=[],[]
for r in cats.itertuples():
    v=np.zeros(n)
    for g in r.genes.split(';'):
        if g in idx: v[idx[g]]=1.0
    if 5<=v.sum()<=n-5: keep.append(r.category); Cl.append(v)
C=np.array(Cl).T
qgc=quintile(U['gc']); qln=quintile(U['size'])
SCHEMES={'none':(np.zeros(n,int),1),'GC_only':(qgc,5),'length_only':(qln,5),
         'GC_x_length':(qgc*5+qln,25)}
rows=[]
for sname,(st,k) in SCHEMES.items():
    for s in SETS:
        t=U['target'][s].astype(float)
        orm,chi2=cmh_vec(C,st,t,k)
        p=stats.chi2.sf(chi2,1); q=multipletests(p,method='fdr_bh')[1]
        for j,cat in enumerate(keep):
            rows.append(dict(stratification=sname,set=s,category=cat,n_cat=int(C[:,j].sum()),
                             or_cmh=orm[j],p=p[j],q=q[j]))
R=pd.DataFrame(rows); R.to_csv(f'{OUT}/covariate_decomposition_261104.tsv',sep='\t',index=False)
piv=R[R.q<.05].groupby(['stratification','set']).size().unstack(fill_value=0).reindex(SCHEMES.keys())
piv=piv.reindex(columns=SETS,fill_value=0)
piv.to_csv(f'{OUT}/covariate_decomposition_summary_261104.tsv',sep='\t')
print(piv.to_string())
print('peak_rss_gb',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9)
