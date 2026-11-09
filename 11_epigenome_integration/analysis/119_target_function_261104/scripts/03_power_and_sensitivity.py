#!/usr/bin/env python
"""119: (i) minimum detectable OR by category size K; (ii) sensitivity analysis with the
universe restricted to genes that carry an experimental Jeong TSS (promoter sets only)."""
import os, resource
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

RS  = '/Users/okaban/bioinfo/rna-seq'
OUT = f'{RS}/11_epigenome_integration/analysis/119_target_function_261104/tables'
SETS = ['prom_4mC','prom_6mA','body_4mC','body_6mA','bodycall_4mC','bodycall_6mA']
MIN_CAT = 5
rng = np.random.default_rng(20261104)
exec(open(os.path.join(os.path.dirname(__file__), '02_enrichment.py')).read().split('if __name__')[0])

# ---------- (i) minimum detectable OR ----------
def mdor(N, K, p, alpha, power=0.80, side='up'):
    za = stats.norm.isf(alpha/2); zb = stats.norm.ppf(power)
    lo, hi = (1.0001, 200.0) if side == 'up' else (1/200.0, 0.9999)
    def z(theta):
        p1 = theta*p/(1-p+theta*p)
        a = K*p1; b = K*(1-p1); c = (N-K)*p; d = (N-K)*(1-p)
        if min(a,b,c,d) <= 0: return 0.0
        v = 1/a+1/b+1/c+1/d
        return abs(np.log(theta))/np.sqrt(v)
    f = lambda th: z(th) - (za+zb)
    if side == 'up':
        if f(hi) < 0: return np.nan
        for _ in range(200):
            mid = np.sqrt(lo*hi)
            if f(mid) < 0: lo = mid
            else: hi = mid
        return hi
    else:
        if f(lo) < 0: return np.nan
        for _ in range(200):
            mid = np.sqrt(lo*hi)
            if f(mid) < 0: hi = mid
            else: lo = mid
        return lo

ann = pd.read_csv(f'{OUT}/genes_universe_261104.tsv', sep='\t')
res = pd.read_csv(f'{OUT}/enrichment_full_261104.tsv', sep='\t')
rows = []
for unit in ['gene','tu','tile']:
    sub = res[res.unit == unit]
    if sub.empty: continue
    m = sub.category.nunique()
    N = int(sub.n_unit.iloc[0])
    for s in SETS:
        p = float(sub[sub.set == s].n_target.iloc[0]) / N
        for K in [5,10,20,50,100,200,500,1000]:
            if K >= N: continue
            for aname, alpha in (('BH-floor(0.05/m)', 0.05/m), ('nominal 0.05', 0.05)):
                rows.append(dict(unit=unit, set=s, n_unit=N, n_target=int(p*N),
                                 prevalence=p, n_cat_tested=m, K=K, alpha_name=aname,
                                 alpha=alpha,
                                 min_OR_up=mdor(N,K,p,alpha,0.80,'up'),
                                 max_OR_down=mdor(N,K,p,alpha,0.80,'down')))
pw = pd.DataFrame(rows)
pw.to_csv(f'{OUT}/power_floor_261104.tsv', sep='\t', index=False)
print('power floor rows', len(pw), flush=True)

# ---------- (ii) TSS-restricted sensitivity (gene level, promoter sets) ----------
cats = pd.read_csv(f'{OUT}/categories_261104.tsv', sep='\t')
sub = ann[ann.has_tss].reset_index(drop=True)
globals()['ann'] = sub
U = make_gene_units()
globals()['SETS'] = ['prom_4mC','prom_6mA']
globals()['NPERM'] = 1000
df, perm, nc, strat = run_unit(U, cats, 'gene_TSSonly')
df.to_csv(f'{OUT}/sensitivity_tssonly_261104.tsv', sep='\t', index=False)
print('TSS-restricted universe', len(sub), 'cats', nc,
      'CMH FDR5%', {s: int((df[df.set==s].q_cmh<.05).sum()) for s in ['prom_4mC','prom_6mA']},
      'Fisher FDR5%', {s: int((df[df.set==s].q_fisher<.05).sum()) for s in ['prom_4mC','prom_6mA']},
      'perm FWER', perm, flush=True)
print('peak_rss_gb', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9, flush=True)
