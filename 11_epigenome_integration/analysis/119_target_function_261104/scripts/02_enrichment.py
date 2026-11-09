#!/usr/bin/env python
"""119: functional bias of methylation-TARGET genes.

Primary estimator: Cochran-Mantel-Haenszel OR stratified by GC quintile x length
quintile (25 strata).  Naive Fisher reported as contrast only.
Statistical unit ladder: gene -> transcription unit -> 20 kb tile.
"""
import os, sys, json, resource
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

RS  = '/Users/okaban/bioinfo/rna-seq'
D   = f'{RS}/11_epigenome_integration/analysis/119_target_function_261104'
OUT = f'{D}/tables'
GENOME_LEN = 8667507
TILE = 20000
NPERM = 1000
MIN_CAT = 5
rng = np.random.default_rng(20261104)

ann  = pd.read_csv(f'{OUT}/genes_universe_261104.tsv', sep='\t')
cats = pd.read_csv(f'{OUT}/categories_261104.tsv', sep='\t')
SETS = ['prom_4mC','prom_6mA','body_4mC','body_6mA','bodycall_4mC','bodycall_6mA']

# ---------------- unit construction ----------------
def quintile(x):
    return pd.qcut(pd.Series(x).rank(method='first'), 5, labels=False).to_numpy()

def make_gene_units():
    gid = ann.gene_id.to_numpy()
    tgt = {s: ann[s].to_numpy(bool) for s in SETS}
    return dict(name='gene', ids=gid, target=tgt,
                gc=ann.gc.to_numpy(), size=ann.length.to_numpy(),
                member=[{g} for g in gid])

def make_tu_units():
    a = ann.sort_values(['strand','start']).reset_index(drop=True)
    tu_id = np.empty(len(a), int); cur = -1; prev = None
    for i, r in enumerate(a.itertuples()):
        if prev is None or r.strand != prev.strand or (r.start - prev.end - 1) > 50:
            cur += 1
        tu_id[i] = cur; prev = r
    a['tu'] = tu_id
    grp = a.groupby('tu')
    span_gc = grp.apply(lambda d: np.average(d.gc, weights=d.length), include_groups=False).to_numpy()
    size    = grp.length.sum().to_numpy()
    tgt = {s: grp[s].any().to_numpy() for s in SETS}
    member = [set(d) for d in grp.gene_id.apply(list)]
    return dict(name='tu', ids=np.array([f'TU{i}' for i in grp.size().index]),
                target=tgt, gc=span_gc, size=size, member=member)

def make_tile_units():
    a = ann.copy()
    a['tile'] = ((a.start + a.end)//2) // TILE
    grp = a.groupby('tile')
    gc   = grp.apply(lambda d: np.average(d.gc, weights=d.length), include_groups=False).to_numpy()
    size = grp.size().to_numpy().astype(float)        # gene count (tiles are equal width)
    tgt = {s: grp[s].any().to_numpy() for s in SETS}
    member = [set(d) for d in grp.gene_id.apply(list)]
    return dict(name='tile', ids=np.array([f'tile{i}' for i in grp.size().index]),
                target=tgt, gc=gc, size=size, member=member)

# ---------------- vectorised CMH ----------------
def cmh_vec(C, strat, t, nstrat):
    """C: (nunit, ncat) float membership; strat: (nunit,) int; t: (nunit,) float target.
    returns (or_mh, chi2) arrays of length ncat."""
    S = np.zeros((len(strat), nstrat)); S[np.arange(len(strat)), strat] = 1.0
    Ns = S.sum(0)                               # (k,)
    m1 = t @ S                                  # targets per stratum
    m0 = Ns - m1
    n1 = C.T @ S                                # (ncat,k) category size per stratum
    n0 = Ns[None, :] - n1
    a  = C.T @ (S * t[:, None])                 # (ncat,k)
    b  = n1 - a
    c  = m1[None, :] - a
    d  = n0 - c
    with np.errstate(divide='ignore', invalid='ignore'):
        ok = Ns > 1
        E  = np.where(ok, n1 * m1[None, :] / np.where(Ns == 0, 1, Ns)[None, :], 0.0)
        V  = np.where(ok, n1 * n0 * m1[None, :] * m0[None, :] /
                      (np.where(Ns == 0, 1, Ns)[None, :]**2 *
                       np.where(Ns > 1, Ns - 1, 1)[None, :]), 0.0)
        R  = np.where(ok, a * d / np.where(Ns == 0, 1, Ns)[None, :], 0.0).sum(1)
        Sx = np.where(ok, b * c / np.where(Ns == 0, 1, Ns)[None, :], 0.0).sum(1)
        num = np.abs(a.sum(1) - E.sum(1)) - 0.5
        num = np.maximum(num, 0.0)
        chi2 = np.where(V.sum(1) > 0, num**2 / np.where(V.sum(1) > 0, V.sum(1), 1), 0.0)
        orm = np.where(Sx > 0, R / np.where(Sx > 0, Sx, 1), np.nan)
    return orm, chi2

def run_unit(U, cat_df, label):
    n = len(U['ids'])
    idx = {g: i for i, g in enumerate(U['ids'])}
    gene2unit = {}
    for i, mem in enumerate(U['member']):
        for g in mem: gene2unit[g] = i
    strat = quintile(U['gc']) * 5 + quintile(U['size'])
    nstrat = 25
    # category membership matrix
    keep, Cl = [], []
    for r in cat_df.itertuples():
        gs = set(r.genes.split(';'))
        v = np.zeros(n)
        for g in gs:
            u = gene2unit.get(g)
            if u is not None: v[u] = 1.0
        if v.sum() >= MIN_CAT and v.sum() <= n - MIN_CAT:
            keep.append((r.source, r.category)); Cl.append(v)
    C = np.array(Cl).T if Cl else np.zeros((n, 0))
    res = []
    for s in SETS:
        t = U['target'][s].astype(float)
        orm, chi2 = cmh_vec(C, strat, t, nstrat)
        p = stats.chi2.sf(chi2, 1)
        # naive Fisher
        pf, orf = np.empty(len(keep)), np.empty(len(keep))
        for j in range(len(keep)):
            cj = C[:, j].astype(bool); tj = t.astype(bool)
            tab = [[int((cj & tj).sum()), int((cj & ~tj).sum())],
                   [int((~cj & tj).sum()), int((~cj & ~tj).sum())]]
            orf[j], pf[j] = stats.fisher_exact(tab)
        q  = multipletests(p,  method='fdr_bh')[1] if len(p)  else np.array([])
        qf = multipletests(pf, method='fdr_bh')[1] if len(pf) else np.array([])
        for j, (src, cat) in enumerate(keep):
            res.append(dict(unit=label, set=s, source=src, category=cat,
                            n_cat=int(C[:, j].sum()), n_target=int(t.sum()),
                            n_unit=n, a=int((C[:, j].astype(bool) & t.astype(bool)).sum()),
                            or_cmh=orm[j], p_cmh=p[j], q_cmh=q[j],
                            or_fisher=orf[j], p_fisher=pf[j], q_fisher=qf[j]))
    # ---- permutation calibration (CMH, FDR5%) ----
    perm = {}
    for s in SETS:
        t = U['target'][s].astype(float)
        hits = 0
        for _ in range(NPERM):
            tp = rng.permutation(t)
            _, chi2 = cmh_vec(C, strat, tp, nstrat)
            pp = stats.chi2.sf(chi2, 1)
            if len(pp) and multipletests(pp, method='fdr_bh')[0].any(): hits += 1
        perm[s] = hits / NPERM
    return pd.DataFrame(res), perm, C.shape[1], strat

# ---------------- validate vectorised CMH against statsmodels ----------------
def validate(U, cat_df):
    strat = quintile(U['gc']) * 5 + quintile(U['size'])
    gene2unit = {}
    for i, mem in enumerate(U['member']):
        for g in mem: gene2unit[g] = i
    n = len(U['ids'])
    r = cat_df.iloc[0]
    v = np.zeros(n)
    for g in r.genes.split(';'):
        u = gene2unit.get(g)
        if u is not None: v[u] = 1.0
    t = U['target']['body_4mC'].astype(float)
    orm, chi2 = cmh_vec(v[:, None], strat, t, 25)
    tabs = []
    for k in range(25):
        m = strat == k
        cj = v[m].astype(bool); tj = t[m].astype(bool)
        tabs.append([[int((cj & tj).sum()), int((cj & ~tj).sum())],
                     [int((~cj & tj).sum()), int((~cj & ~tj).sum())]])
    st = StratifiedTable(np.array(tabs).transpose(1, 2, 0))
    print('VALIDATE cmh  mine OR=%.6f chi2=%.6f | statsmodels OR=%.6f stat=%.6f' %
          (orm[0], chi2[0], st.oddsratio_pooled, st.test_null_odds(correction=True).statistic), flush=True)

if __name__ == '__main__':
    units = dict(gene=make_gene_units(), tu=make_tu_units(), tile=make_tile_units())
    validate(units['gene'], cats)
    allres, allperm, ncats = [], {}, {}
    for lab, U in units.items():
        df, perm, nc, strat = run_unit(U, cats, lab)
        allres.append(df); allperm[lab] = perm; ncats[lab] = nc
        print(lab, 'units', len(U['ids']), 'cats tested', nc,
              'FDR5% hits CMH', {s: int((df[df.set == s].q_cmh < .05).sum()) for s in SETS},
              'Fisher', {s: int((df[df.set == s].q_fisher < .05).sum()) for s in SETS}, flush=True)
    R = pd.concat(allres, ignore_index=True)
    R.to_csv(f'{OUT}/enrichment_full_261104.tsv', sep='\t', index=False)
    pd.DataFrame([dict(unit=u, set=s, perm_fwer=v, n_cat_tested=ncats[u])
                  for u, d in allperm.items() for s, v in d.items()]).to_csv(
        f'{OUT}/permutation_calibration_261104.tsv', sep='\t', index=False)
    print('peak_rss_gb', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9, flush=True)
