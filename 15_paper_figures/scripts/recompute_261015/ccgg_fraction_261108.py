"""L65 の主張「4mC was called at X% of CCGG / GGCCGG positions」を正準規則で再計算する。
正準規則 = 3 ライブラリ合算（和集合インデックス）・真深度 Nvalid+Nfail>=30・Nvalid>0・level>=0.5・refok。
分母 = その規則で「測定可能」なシトシン位置のうち、当該文脈に属するもの。
文脈の取り方を 2 通り出す:
  all_C  : CCGG / GGCCGG 出現内の全シトシン（両鎖）
  modC   : 認識配列で実際に修飾される C のみ（CCGG の 2 番目の C = TGGCCGGC の C5 に対応）
"""
import numpy as np, pandas as pd, pathlib, re, sys
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014')
from _gccggc import genome
HERE = pathlib.Path('/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/recompute_261015')
C = HERE / "pileup_cache_261015"
g = genome(); L = len(g)
TP = {"T1": ["1-1","1-2","1-3"], "T2": ["2-1","2-3","2-4"], "T3": ["3-2","3-3","3-4"]}

def ctx_sets(core, rc):
    """core = plus-strand word, rc = its reverse complement.
    returns (all_C_set, modC_set) of (pos, plus) cytosine positions."""
    allC, modC = set(), set()
    # modified C index within `core` (0-based): second C of CCGG
    for m in re.finditer('(?=%s)' % core, g):
        i = m.start()
        for k, ch in enumerate(core):
            if ch == 'C': allC.add((i+k, 1))
        modC.add((i + MODIDX[core], 1))
    for m in re.finditer('(?=%s)' % rc, g):
        j = m.start()
        for k, ch in enumerate(rc):
            if ch == 'G': allC.add((j+k, 0))      # G on plus == C on minus
        modC.add((j + len(rc) - 1 - MODIDX[core], 0))
    return allC, modC

MODIDX = {'CCGG': 1, 'GGCCGG': 3}   # 0-based index of the modified C on the plus strand
def rcomp(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]

def pooled(libs, code):
    frames = []
    for lib in libs:
        z = np.load(C / f"{lib}.npz")
        m = (z["code"] == code)
        frames.append(pd.DataFrame({"start": z["start"][m], "plus": z["plus"][m],
            "Nvalid": z["Nvalid"][m], "Nmod": z["Nmod"][m], "Nfail": z["Nfail"][m],
            "refok": z["refok"][m]}).set_index(["start","plus"]))
    uni = frames[0].index
    for f in frames[1:]: uni = uni.union(f.index)
    s = sum(f[["Nvalid","Nmod","Nfail"]].reindex(uni).fillna(0) for f in frames)
    s["refok"] = frames[0]["refok"].reindex(uni).groupby(level=[0,1]).first()
    bad = s.refok.isna()
    if bad.any():
        for f in frames[1:]:
            s.loc[bad.values,"refok"] = s.loc[bad.values,"refok"].fillna(
                f["refok"].reindex(uni[bad.values]))
    return s.reset_index()

CTX = {}
for core in ('CCGG','GGCCGG'):
    CTX[core] = ctx_sets(core, rcomp(core))
    print(core, 'all_C', len(CTX[core][0]), 'modC', len(CTX[core][1]), flush=True)

rows = []
for tp, libs in TP.items():
    s = pooled(libs, 1)                      # code 1 = 4mC
    s["lvl"] = s.Nmod / s.Nvalid.replace(0, np.nan)
    meas = s[((s.Nvalid + s.Nfail) >= 30) & (s.Nvalid > 0) & (s.refok.astype(bool))]
    called = meas[meas.lvl >= 0.5]
    mk = set(zip(meas.start, meas.plus)); ck = set(zip(called.start, called.plus))
    for core in ('CCGG','GGCCGG'):
        for lab, S in zip(('all_C','modC'), CTX[core]):
            d = len(mk & S); n = len(ck & S)
            rows.append(dict(timepoint=tp, context=core, c_definition=lab,
                genomic_positions=len(S), measurable=d, called_4mC=n,
                pct_of_measurable=round(100*n/d, 4) if d else float('nan'),
                pct_of_genomic=round(100*n/len(S), 4)))
    print(tp, 'done', flush=True)
r = pd.DataFrame(rows)
r.to_csv(HERE/"ccgg_fraction_261108.tsv", sep="\t", index=False)
print(r.to_string(index=False))
