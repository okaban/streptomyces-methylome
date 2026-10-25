import pandas as pd, numpy as np, json, sys
from collections import defaultdict
libs = {'T1': ['1-1','1-2','1-3'], 'T2': ['2-1','2-3','2-4'], 'T3': ['3-2','3-3','3-4']}

def pair_stats(df, maxd=60_000, min_reads=20):
    d = df[['read','pos','state']].sort_values(['read','pos'])
    acc = defaultdict(lambda: [0,0,0,0])
    for _rd, g in d.groupby('read', sort=False):
        pos = g.pos.values; st = g.state.values
        if len(pos) < 2: continue
        for a in range(len(pos)-1):
            for b in range(a+1, len(pos)):
                if pos[b]-pos[a] > maxd: break
                k = (pos[a], pos[b]); e = acc[k]
                e[0]+=1; e[1]+=int(st[a]==st[b]); e[2]+=st[a]; e[3]+=st[b]
    rows=[]
    for (i,jx),(n,same,si,sj) in acc.items():
        if n < min_reads: continue
        pi,pj = si/n, sj/n
        rows.append((jx-i, n, same/n, pi*pj+(1-pi)*(1-pj)))
    return pd.DataFrame(rows, columns=['dist','n_reads','obs_conc','exp_conc'])

PR={}
for s in ['a01','tgg']:
    for r in sum(libs.values(), []):
        PR[(s,r)] = pd.read_csv(f'pr2b/{s}_{r}.tsv', sep='\t', header=None,
                                names=['read','pos','strand','prob','state'])
res=[]
for s, mark in [('tgg','GCCGGC 4mC'), ('a01','AAGCCCG 6mA')]:
    for tp, rs in libs.items():
        for r in rs:
            d = PR[(s,r)]
            P_ = pair_stats(d); P_ = P_[P_.dist > 10]
            if not len(P_): continue
            P_['excess'] = P_.obs_conc - P_.exp_conc
            w = P_.n_reads; occ = d.state.mean()
            bist = 1.0 - (occ**2 + (1-occ)**2)
            ex = np.average(P_.excess, weights=w)
            res.append((mark, tp, r, len(P_), int(w.sum()), occ, ex,
                        np.average(P_.obs_conc, weights=w), bist, 100*ex/bist))
H = pd.DataFrame(res, columns=['mark','tp','lib','n_pairs','n_obs','occupancy',
                               'excess_conc','obs_conc','bistable_pred','pct_of_bistable'])
pc=[]
for tp, rs in libs.items():
    for r in rs:
        P_ = pair_stats(PR[('a01',r)], maxd=200, min_reads=20)
        P_['excess'] = P_.obs_conc - P_.exp_conc
        same = P_[P_.dist <= 10]
        pc.append((tp, r, len(same), int(same.n_reads.sum()),
                   np.average(same.obs_conc, weights=same.n_reads),
                   np.average(same.excess, weights=same.n_reads)))
PCd = pd.DataFrame(pc, columns=['tp','lib','n_pairs','n_obs','obs_conc','excess'])
DIST=[]
for s, mark in [('tgg','GCCGGC 4mC'), ('a01','AAGCCCG 6mA')]:
    for tp, rs in libs.items():
        for r in rs:
            P_ = pair_stats(PR[(s,r)])
            P_['excess'] = P_.obs_conc - P_.exp_conc
            P_ = P_[P_.dist > 10]
            for iv, g in P_.groupby(pd.cut(P_.dist,[10,200,500,2000,5000,15000,60000]), observed=True):
                if len(g) < 5: continue
                DIST.append((mark, r, iv.right, len(g), np.average(g.excess, weights=g.n_reads)))
D26 = pd.DataFrame(DIST, columns=['mark','lib','dist_max','n_pairs','excess'])
H.to_csv('H26b.csv', index=False); PCd.to_csv('PC26b.csv', index=False); D26.to_csv('D26b.csv', index=False)
six = H[H['mark'].str.contains('6mA')]
print("bistable", round(six.bistable_pred.mean(),3),
      "| within-instance", round(PCd.excess.mean(),3),
      "| between", round(six.excess_conc.mean(),3))
for mk in ['GCCGGC 4mC','AAGCCCG 6mA']:
    print(mk, "pct_of_bistable", round(H[H['mark']==mk].pct_of_bistable.mean(),1))
