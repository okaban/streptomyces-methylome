#!/usr/bin/env python
"""119: build gene universe, modifiable-site universe, 6 target gene sets, covariates.

Outputs tables/genes_universe_261104.tsv (one row per chromosomal gene) and
tables/sites_universe_261104.tsv (one row per modifiable position).
"""
import os, re, sys, resource
import numpy as np, pandas as pd

RS   = '/Users/okaban/bioinfo/rna-seq'
OUT  = f'{RS}/11_epigenome_integration/analysis/119_target_function_261104/tables'
FNA  = f'{RS}/11_epigenome_integration/data/NC_003888.3.fna'
ANN  = f'{RS}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
COG  = f'{RS}/12_supplementary_figures/analysis/12_supplementary_260202_v1/tables/gene_COG_classification.tsv'
KEGG = f'{RS}/15_paper_figures/scripts/newfig_261101/tables/kegg_br_sco00001_261102.keg'
TSS  = f'{RS}/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv'
PUP  = '/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015'
T1   = ['1-1','1-2','1-3']
os.makedirs(OUT, exist_ok=True)

# ---------- genome ----------
seq = []
with open(FNA) as fh:
    for ln in fh:
        if not ln.startswith('>'): seq.append(ln.strip())
G = ''.join(seq).upper()
L = len(G)
print('genome', L, flush=True)

COMP = str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(COMP)[::-1]

def find_all(hay, needle):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i); i = hay.find(needle, i+1)
    return out

# modifiable positions.  0-based genomic coord of the modified base; strand = strand
# the motif reads on.  TGGCCGGC C5 -> offset 4; AAGCCCG C5 -> offset 4; AAGCCCG A1 -> offset 0.
SPEC = [('TGGCCGGC', 4, '4mC', 'TGGCCGGC_C5'),
        ('AAGCCCG',  4, '4mC', 'AAGCCCG_C5'),
        ('AAGCCCG',  0, '6mA', 'AAGCCCG_A1')]
rows = []
for motif, off, mark, label in SPEC:
    k = len(motif)
    for st, pat in (('+', motif), ('-', rc(motif))):
        for s in find_all(G, pat):
            pos = s + off if st == '+' else s + (k-1-off)
            rows.append((pos, st, motif, off, mark, label))
sites = pd.DataFrame(rows, columns=['pos','strand','motif','offset','mark','label'])
sites = sites.sort_values(['label','pos']).reset_index(drop=True)
print(sites.groupby(['label','strand']).size(), flush=True)
# sanity: base identity
bad = 0
for r in sites.itertuples():
    b = G[r.pos]
    want = 'C' if r.mark == '4mC' else 'A'
    if r.strand == '-': want = {'C':'G','A':'T'}[want]
    if b != want: bad += 1
assert bad == 0, bad

# ---------- T1 calls (pooled over 3 replicates, project rule: dep>=30, lvl>=0.5) ----------
cols = ['chrom','start','end','code','score','strand','ts','te','col','Nvalid','pct',
        'Nmod','Ncan','Noth','Ndel','Nfail','Ndiff','Nnocall']
acc = {}
for lib in T1:
    d = pd.read_csv(f'{PUP}/{lib}_nofilt.bed', sep='\t', header=None, names=cols,
                    usecols=['start','code','strand','Nvalid','Nmod','Nfail'])
    d = d[d.code.str.startswith('21839') | d.code.str.startswith('a,')]
    d['motif'] = d.code.str.split(',').str[1]
    d['off']   = d.code.str.split(',').str[2].astype(int)
    for r in d.itertuples():
        k = (r.start, r.strand)
        v = acc.get(k)
        if v is None: acc[k] = [r.Nvalid, r.Nmod, r.Nfail]
        else: v[0]+=r.Nvalid; v[1]+=r.Nmod; v[2]+=r.Nfail
    del d
cal = pd.DataFrame([(k[0],k[1],v[0],v[1],v[2]) for k,v in acc.items()],
                   columns=['pos','strand','Nvalid','Nmod','Nfail'])
del acc
cal['dep'] = cal.Nvalid + cal.Nfail
cal['lvl'] = np.where(cal.Nvalid>0, cal.Nmod/cal.Nvalid.replace(0,np.nan), np.nan)
cal['eligible'] = (cal.dep>=30) & (cal.Nvalid>0)
cal['called']   = cal.eligible & (cal.lvl>=0.5)
sites = sites.merge(cal[['pos','strand','dep','lvl','eligible','called']],
                    on=['pos','strand'], how='left')
sites['eligible'] = sites.eligible.fillna(False)
sites['called']   = sites.called.fillna(False)
print(sites.groupby(['label']).agg(n=('pos','size'), elig=('eligible','sum'),
                                   called=('called','sum'), occ=('called','mean')), flush=True)
sites.to_csv(f'{OUT}/sites_universe_261104.tsv', sep='\t', index=False)

# ---------- genes ----------
ann = pd.read_csv(ANN, sep='\t', dtype=str)
ann = ann[ann.contig == 'NC_003888.3'].copy()
ann['start'] = ann.start.astype(int); ann['end'] = ann.end.astype(int)
ann = ann.sort_values('start').reset_index(drop=True)
print('chromosomal genes', len(ann), flush=True)
assert len(ann) == 7863

gc, glen = [], []
for r in ann.itertuples():
    s = G[r.start-1:r.end]            # annotation is 1-based inclusive
    glen.append(len(s))
    gc.append((s.count('G')+s.count('C'))/max(len(s),1))
ann['length'] = glen; ann['gc'] = gc

# ---------- TSS -> promoter windows ----------
tss = pd.read_csv(TSS)
tss = tss[tss.gene_id.notna()].copy()
gidx = {g:i for i,g in enumerate(ann.gene_id)}
tss = tss[tss.gene_id.isin(gidx)].copy()
print('TSS with chromosomal gene', len(tss), 'unique genes', tss.gene_id.nunique(), flush=True)
prom = {}   # gene_id -> list of (lo, hi) 0-based half-open
for r in tss.itertuples():
    p = int(r.tss_position)           # 1-based
    if r.strand == '+': lo, hi = p-1-300, p-1
    else:               lo, hi = p, p+300
    prom.setdefault(r.gene_id, []).append((max(lo,0), min(hi,L)))

# ---------- assign sites to genes ----------
gs = ann.start.to_numpy()-1; ge = ann.end.to_numpy()   # 0-based half-open
order = np.argsort(gs); gs_s = gs[order]; ge_s = ge[order]
gene_ids = ann.gene_id.to_numpy()

def genes_overlapping(lo, hi):
    """indices (into ann) of genes overlapping [lo,hi)."""
    j = np.searchsorted(gs_s, hi, side='left')
    out = []
    for k in range(j-1, -1, -1):
        if ge_s[k] > lo: out.append(order[k])
        if gs_s[k] < lo - 60000: break
    return out

body_has  = {m: np.zeros(len(ann), bool) for m in ['4mC','6mA']}
body_call = {m: np.zeros(len(ann), bool) for m in ['4mC','6mA']}
for r in sites.itertuples():
    for gi in genes_overlapping(r.pos, r.pos+1):
        body_has[r.mark][gi] = True
        if r.called: body_call[r.mark][gi] = True

prom_has = {m: np.zeros(len(ann), bool) for m in ['4mC','6mA']}
spos = {m: np.sort(sites.loc[sites.mark==m,'pos'].to_numpy()) for m in ['4mC','6mA']}
for g, wins in prom.items():
    gi = gidx[g]
    for lo, hi in wins:
        for m in ['4mC','6mA']:
            a = np.searchsorted(spos[m], lo); b = np.searchsorted(spos[m], hi)
            if b > a: prom_has[m][gi] = True

ann['has_tss'] = ann.gene_id.isin(prom).to_numpy()
for m in ['4mC','6mA']:
    ann[f'prom_{m}']      = prom_has[m]
    ann[f'body_{m}']      = body_has[m]
    ann[f'bodycall_{m}']  = body_call[m]

# ---------- KEGG ----------
sco2gene = {}
for r in ann.itertuples():
    if isinstance(r.old_locus_tag, str) and r.old_locus_tag.strip():
        sco2gene[r.old_locus_tag.strip()] = r.gene_id
kegg = {}
cur = None
with open(KEGG) as fh:
    for ln in fh:
        if ln.startswith('C'):
            mm = re.match(r'C\s+(\d{5})\s+(.*?)\s*$', ln)
            if mm: cur = f'ko{mm.group(1)}|{mm.group(2).split(" [")[0]}'
            else: cur = None
        elif ln.startswith('D') and cur:
            mm = re.search(r'\bSCO\d+\b', ln)
            if mm:
                g = sco2gene.get(mm.group(0))
                if g: kegg.setdefault(cur, set()).add(g)
        elif ln.startswith('B'): cur = None
print('KEGG pathways', len(kegg), 'genes mapped', len(set().union(*kegg.values())), flush=True)

# ---------- COG ----------
cg = pd.read_csv(COG, sep='\t', dtype=str)
cg = cg[cg.gene_id.isin(gidx)]
cogs = {}
for r in cg.itertuples():
    v = r.COG_category
    if not isinstance(v, str) or not v.strip(): continue
    for part in re.split(r'[;,]', v):
        part = part.strip()
        if not part: continue
        code = part.split(' - ')[0].strip()
        if len(code) != 1: continue
        cogs.setdefault(f'COG_{code}|{part}', set()).add(r.gene_id)
print('COG cats', len(cogs), flush=True)

# serialise category membership
rows = []
for src, dd in (('KEGG', kegg), ('COG', cogs)):
    for cat, gset in dd.items():
        rows.append((src, cat, len(gset), ';'.join(sorted(gset))))
pd.DataFrame(rows, columns=['source','category','n','genes']).to_csv(
    f'{OUT}/categories_261104.tsv', sep='\t', index=False)

ann.to_csv(f'{OUT}/genes_universe_261104.tsv', sep='\t', index=False)
summ = {c: int(ann[c].sum()) for c in ann.columns if c.startswith(('prom_','body_','bodycall_'))}
summ['has_tss'] = int(ann.has_tss.sum())
print('SET SIZES', summ, flush=True)
print('peak_rss_gb', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9, flush=True)
