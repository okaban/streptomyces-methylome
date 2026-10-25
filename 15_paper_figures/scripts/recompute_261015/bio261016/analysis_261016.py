"""未帰属 6mA（T1 1,531 部位）の生物学的意味の検証。2026-10-16。
仮説・予測・判定規則は PREDICTIONS_261016.md に測定前に固定済み。
実行: cd .../recompute_261015 && PYTHONPATH=../supp_261014 python bio261016/analysis_261016.py
"""
import numpy as np, pandas as pd, pathlib, re, sys
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import StratifiedTable
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014')
from _gccggc import genome

HERE = pathlib.Path('/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/recompute_261015')
C = HERE / "pileup_cache_261015"; OUT = HERE / "bio261016"; OUT.mkdir(exist_ok=True)
g = genome(); L = len(g)
TP = {"T1": ["1-1","1-2","1-3"], "T2": ["2-1","2-3","2-4"], "T3": ["3-2","3-3","3-4"]}
MOT = "AAGCCCG"; W = 12

aag_A = set()
for m in re.finditer("AAGCCCG", g): aag_A |= {(m.start(), True), (m.start()+1, True)}
for m in re.finditer("CGGGCTT", g): aag_A |= {(m.start()+6, False), (m.start()+5, False)}

def pooled6(tp):
    """時点のプールされた 6mA パイルアップ（参照塩基照合済み、pre-filter 深度 >= 30）。"""
    fr = []
    for lib in TP[tp]:
        z = np.load(C / f"{lib}.npz"); m = (z["code"] == 0)
        fr.append(pd.DataFrame({"start": z["start"][m], "plus": z["plus"][m], "Nvalid": z["Nvalid"][m],
            "Nmod": z["Nmod"][m], "Nfail": z["Nfail"][m], "refok": z["refok"][m]}).set_index(["start","plus"]))
    uni = fr[0].index
    for f in fr[1:]: uni = uni.union(f.index)
    s = sum(f[["Nvalid","Nmod","Nfail"]].reindex(uni).fillna(0) for f in fr)
    ok = fr[0]["refok"].reindex(uni)
    for f in fr[1:]: ok = ok.fillna(f["refok"].reindex(uni))
    s["refok"] = ok.fillna(False).astype(bool); s = s.reset_index()
    s["depth"] = s.Nvalid + s.Nfail; s["lvl"] = s.Nmod / s.Nvalid.replace(0, np.nan)
    return s[s.refok & (s.depth >= 30) & (s.Nvalid > 0)]

e = pooled6("T1")
e["called"] = e.lvl >= 0.5
e["in_motif"] = [(p,b) in aag_A for p,b in zip(e.start, e.plus)]
e = e[~e.in_motif].copy()           # 2,412,520 適格、うち called 1,531

# ---------- ゲノム注釈 ----------
GFF = ("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/"
       "data/GCF_000203835.1_ASM20383v1_genomic.gff")
cds = np.zeros(L, bool); rna = np.zeros(L, bool); mge = np.zeros(L, bool); genes = []
for line in open(GFF):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != "NC_003888.3": continue
    s0, e0, ty, att, strand = int(f[3])-1, int(f[4]), f[2], f[8], f[6]
    if ty == "CDS": cds[s0:e0] = True
    if ty in ("rRNA","tRNA","ncRNA","tmRNA"): rna[s0:e0] = True
    if ty == "CDS" and re.search("transposase|integrase|recombinase|IS[0-9]|mobile", att, re.I): mge[s0:e0] = True
    if ty == "gene":
        gid = re.search(r"locus_tag=([^;]+)", att)
        genes.append((s0, e0, strand, gid.group(1) if gid else ""))
genes.sort()
gap_id = np.full(L, -1, np.int32); gaps = []
for i in range(len(genes)-1):
    a, b = genes[i], genes[i+1]
    if b[0] <= a[1]: continue
    cls = {("+","+"):"tandem", ("-","-"):"tandem", ("+","-"):"convergent", ("-","+"):"divergent"}[(a[2], b[2])]
    gap_id[a[1]:b[0]] = len(gaps)
    gaps.append(dict(gap=len(gaps), start=a[1], end=b[0], length=b[0]-a[1], cls=cls, up=a[3], dn=b[3]))
gapdf = pd.DataFrame(gaps)

# ---------- 配列特徴（すべてアルゴリズム定義。経験的な構造パラメータ表は使わない） ----------
def runs_of(mask):
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)
ga = np.frombuffer(g.encode(), np.uint8)
isgc = ((ga == ord("G")) | (ga == ord("C"))).astype(np.int32); cumgc = np.concatenate([[0], np.cumsum(isgc)])
hp_any = np.zeros(L, bool); hp_base = np.zeros(L, np.uint8)
for b in "ACGT":                      # b*5 が [p-5, p+6) に出現（deepdive_261015 と同一定義）
    eqb = (ga == ord(b)).astype(np.int8); c = np.concatenate([[0], np.cumsum(eqb)])
    s = np.flatnonzero((c[5:] - c[:-5]) == 5)
    d = np.zeros(L+2, np.int32); np.add.at(d, np.maximum(s-1, 0), 1); np.add.at(d, np.minimum(s+6, L+1), -1)
    m = np.cumsum(d)[:L] > 0; new = m & ~hp_any; hp_base[new] = ord(b); hp_any |= m
g4base = np.zeros(L, np.int16)        # G4Hunter: G 連長 n に +min(n,4)、C 連長に -min(n,4)
for base, sign in (("G",1), ("C",-1)):
    st, en = runs_of(ga == ord(base)); ln = (np.minimum(en-st, 4)*sign).astype(np.int16); w = en-st
    idx = np.repeat(np.arange(len(st)), w)
    offs = np.arange(w.sum()) - np.repeat(np.cumsum(np.concatenate([[0], w[:-1]])), w)
    g4base[st[idx] + offs] = ln[idx]
cg4 = np.concatenate([[0], np.cumsum(g4base.astype(np.int64))])
g4mot = np.zeros(L, bool)             # G{3,}N{1,7} を 4 回（両鎖）
for base in ("G","C"):
    st, en = runs_of(ga == ord(base)); keep = (en-st) >= 3; st, en = st[keep], en[keep]
    gp = st[1:] - en[:-1]; ok = (gp >= 1) & (gp <= 7)
    for i in np.flatnonzero(ok[:-2] & ok[1:-1] & ok[2:]): g4mot[st[i]:en[i+3]] = True
pur = ((ga == ord("A")) | (ga == ord("G"))); change = pur[1:] != pur[:-1]
brk = np.concatenate([[-1], np.flatnonzero(~change), [L-1]]); seglen = np.diff(brk)
cov = np.repeat(seglen, seglen).astype(np.int32)   # 交代 RY 連長（Z-DNA 代理）

TSSF = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv"
t = pd.read_csv(TSSF)
sel = t[t.category.astype(str).str.upper().isin(["P","S"]) & t.gene_id.notna()
        & (t.gene_id.astype(str).str.strip().str.len() > 0)]
tss = np.sort(sel.tss_position.astype(int).values)
pos = e.start.to_numpy()
i = np.searchsorted(tss, pos)
lo = np.where(i > 0, pos - tss[np.clip(i-1, 0, None)], 10**9)
hi = np.where(i < len(tss), tss[np.clip(i, None, len(tss)-1)] - pos, 10**9)
e["d_tss"] = np.minimum(lo, hi); e["promoter"] = e.d_tss <= 300
e["in_cds"] = cds[pos]; e["rna_gene"] = rna[pos]; e["mge"] = mge[pos]
e["intergenic"] = (~e.in_cds) & (~e.promoter)
e["gc50"] = ((cumgc[np.minimum(pos+26, L)] - cumgc[np.maximum(pos-25, 0)])
             / (np.minimum(pos+26, L) - np.maximum(pos-25, 0)))
e["homopolymer"] = hp_any[pos]; e["hp_base"] = [chr(x) if x else "" for x in hp_base[pos]]
e["g4_motif"] = g4mot[pos]
e["g4h"] = ((cg4[np.minimum(pos+26, L)] - cg4[np.maximum(pos-25, 0)])
            / (np.minimum(pos+26, L) - np.maximum(pos-25, 0)))
e["altRY"] = cov[pos]; e["gap_id"] = gap_id[pos]
K = 21; code = np.full(256, -1, np.int64)
for ii, b in enumerate("ACGT"): code[ord(b)] = ii
a = code[ga]; km = np.zeros(L-K+1, np.int64); bad = np.zeros(L-K+1, bool)
for ii in range(K):
    seg = a[ii:L-K+1+ii]; km = km*4 + np.maximum(seg, 0); bad |= (seg < 0)
km[bad] = -1
u, cnt = np.unique(km, return_counts=True); occ = cnt[np.searchsorted(u, km)]
e["kmer21_occ"] = occ[np.clip(pos - K//2, 0, L-K)]; e["repeat"] = e.kmer21_occ > 1

U = {}
for tp in ("T2","T3"):
    s = pooled6(tp); s = s[s.lvl >= 0.5]
    U[tp] = set(zip(s.start.tolist(), s.plus.tolist())) - aag_A
key = list(zip(e.start.tolist(), e.plus.tolist()))
e["in_T2"] = [k in U["T2"] for k in key]; e["in_T3"] = [k in U["T3"] for k in key]
e["nrep"] = e.in_T2.astype(int) + e.in_T3.astype(int)

# 修飾鎖の読み枠での配列窓
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
raw = [g[max(0, p-W):p+W+1] for p in pos]
seqs = np.array([s if pl else rc(s) for s, pl in zip(raw, e.plus.to_numpy())], dtype=object)
arr = np.frombuffer("".join(seqs).encode(), np.uint8).reshape(len(seqs), 2*W+1)
assert (arr[:, W] == ord("A")).all()
MOTB = np.frombuffer(MOT.encode(), np.uint8)
k0 = arr[:, W:W+7]; k1 = arr[:, W-1:W+6]          # 部位が A0 / A1 を演じる 2 登録
d0 = (k0 != MOTB).sum(1); d1 = (k1 != MOTB).sum(1)
e["reg"] = np.where(d0 <= d1, 0, 1); e["hd_aagcccg"] = np.minimum(d0, d1)
best = np.where(e.reg.values[:, None] == 0, k0, k1)
def maxrun(b):
    bestr = np.zeros(len(b), np.int8); cur = np.zeros(len(b), np.int8)
    for j in range(b.shape[1]):
        cur = np.where(b[:, j], cur+1, 0); bestr = np.maximum(bestr, cur)
    return bestr
e["crun3"] = maxrun(arr[:, W+1:W+9] == ord("C"))   # 3' 側 +1..+8 の最長 C 連続
e["grun3"] = maxrun(arr[:, W+1:W+9] == ord("G"))   # 同一窓・同一 GC の内部対照
e["b5p"] = [chr(x) for x in arr[:, W-1]]
cl = e[e.called]; nc = e[~e.called]

# ---------- 統計ヘルパ ----------
def mh_or(x, y, col, strat_col="gc50", q=10):
    edges = np.unique(np.quantile(y[strat_col], np.linspace(0, 1, q+1)))
    xb = np.clip(np.searchsorted(edges, x[strat_col], side="right")-1, 0, len(edges)-2)
    yb = np.clip(np.searchsorted(edges, y[strat_col], side="right")-1, 0, len(edges)-2)
    num = den = 0.0
    for k in range(len(edges)-1):
        xs = x[xb == k]; ys = y[yb == k]
        if len(xs) == 0 or len(ys) == 0: continue
        A = int(xs[col].sum()); B = len(xs)-A; Cc = int(ys[col].sum()); D = len(ys)-Cc
        n = A+B+Cc+D
        if n: num += A*D/n; den += B*Cc/n
    return num/den if den > 0 else np.nan

def cmh(x, y, col, q=10):
    edges = np.unique(np.quantile(y["gc50"], np.linspace(0, 1, q+1)))
    xb = np.clip(np.searchsorted(edges, x["gc50"], side="right")-1, 0, len(edges)-2)
    yb = np.clip(np.searchsorted(edges, y["gc50"], side="right")-1, 0, len(edges)-2)
    tabs = []
    for k in range(len(edges)-1):
        xs = x[xb == k]; ys = y[yb == k]
        if len(xs) == 0 or len(ys) == 0: continue
        tabs.append([[int(xs[col].sum()), len(xs)-int(xs[col].sum())],
                     [int(ys[col].sum()), len(ys)-int(ys[col].sum())]])
    T = StratifiedTable(np.dstack([np.array(t) for t in tabs]))
    return T.oddsratio_pooled, T.test_null_odds().pvalue, T.oddsratio_pooled_confint()

rng = np.random.default_rng(0)
def strat_diff(x, y, col, q=10, nboot=1000):
    """GC 十分位層別の重み付き平均差とブートストラップ 95 % CI。"""
    edges = np.unique(np.quantile(y["gc50"], np.linspace(0, 1, q+1)))
    xb = np.clip(np.searchsorted(edges, x["gc50"], side="right")-1, 0, len(edges)-2)
    yb = np.clip(np.searchsorted(edges, y["gc50"], side="right")-1, 0, len(edges)-2)
    yv = y[col].values
    ymean = np.array([yv[yb == k].mean() if (yb == k).sum() else np.nan for k in range(len(edges)-1)])
    xv = x[col].values
    def wd(idx):
        bb = xb[idx]; vals = xv[idx]; num = den = 0.0
        for k in range(len(edges)-1):
            m = bb == k
            if m.sum() and np.isfinite(ymean[k]): num += m.sum()*(vals[m].mean()-ymean[k]); den += m.sum()
        return num/den
    base = wd(np.arange(len(xv)))
    bs = np.array([wd(rng.integers(0, len(xv), len(xv))) for _ in range(nboot)])
    return base, np.percentile(bs, 2.5), np.percentile(bs, 97.5)

# ---------- B3 再現 vs 一過性（12 検定、BH） ----------
cal = cl.copy(); cal["repro"] = cal.nrep == 2
R = cal[cal.repro]; Tr = cal[~cal.repro]
rows = []
for col in ["promoter","in_cds","intergenic","rna_gene","mge","homopolymer","repeat","g4_motif"]:
    A = int(R[col].sum()); B = len(R)-A; Cc = int(Tr[col].sum()); D = len(Tr)-Cc
    orr, p = stats.fisher_exact([[A,B],[Cc,D]])
    rows.append(dict(feature=col, repro_pct=100*R[col].mean(), trans_pct=100*Tr[col].mean(),
                     OR=orr, p=p, MH_gc_repro=mh_or(R, nc, col), MH_gc_trans=mh_or(Tr, nc, col)))
for col in ["gc50","g4h","altRY","lvl","depth","Nvalid","d_tss","kmer21_occ"]:
    rows.append(dict(feature=col+"_median", repro_pct=float(R[col].median()), trans_pct=float(Tr[col].median()),
                     OR=np.nan, p=stats.mannwhitneyu(R[col], Tr[col])[1], MH_gc_repro=np.nan, MH_gc_trans=np.nan))
B3 = pd.DataFrame(rows); B3["q_BH"] = multipletests(B3.p, method="fdr_bh")[1]
B3.to_csv(OUT/"B3_repro_vs_transient_261016.tsv", sep="\t", index=False)

# ---------- B1 遺伝子間間隙の向き（4 検定） ----------
gcls = gapdf.set_index("gap")["cls"]; glen = gapdf.set_index("gap")["length"]
ing = e[e.gap_id >= 0].copy(); ing["cls"] = gcls.reindex(ing.gap_id).values
ing["gaplen"] = glen.reindex(ing.gap_id).values
cg = ing[ing.called]; ng = ing[~ing.called]
chi2, pchi, dof, _ = stats.chi2_contingency(pd.crosstab(ing.called, ing.cls).values)
b1 = []
for k in ["convergent","divergent","tandem"]:
    cgk = cg.assign(_f=(cg.cls == k)); ngk = ng.assign(_f=(ng.cls == k))
    A = int(cgk._f.sum()); B = len(cgk)-A; Cc = int(ngk._f.sum()); D = len(ngk)-Cc
    orr, p = stats.fisher_exact([[A,B],[Cc,D]])
    b1.append(dict(klass=k, called_n=A, called_pct=100*A/len(cgk), nc_pct=100*Cc/len(ngk),
                   OR=orr, p=p, MH_gc=mh_or(cgk, ngk, "_f")))
B1 = pd.DataFrame(b1); B1["q_BH"] = multipletests(B1.p, method="fdr_bh")[1]
B1.attrs["chi2"] = (chi2, dof, pchi)
X = sm.add_constant(pd.DataFrame({"gc50": ing.gc50.values, "log_depth": np.log10(ing.depth.values),
    "convergent": (ing.cls == "convergent").astype(int).values,
    "divergent": (ing.cls == "divergent").astype(int).values,
    "log_gaplen": np.log10(ing.gaplen.values.astype(float)+1)}))
m1 = sm.Logit(ing.called.astype(int).values, X.values).fit(disp=0)
B1_logit = pd.DataFrame({"term": X.columns, "OR": np.exp(m1.params), "p": m1.pvalues})
B1.to_csv(OUT/"B1_gap_class_261016.tsv", sep="\t", index=False)
B1_logit.to_csv(OUT/"B1_gap_logit_261016.tsv", sep="\t", index=False)

# ---------- B2 構造（4 検定） ----------
b2 = []
for col in ["g4_motif","homopolymer"]:
    A = int(cl[col].sum()); B = len(cl)-A; Cc = int(nc[col].sum()); D = len(nc)-Cc
    orr, p = stats.fisher_exact([[A,B],[Cc,D]]); orp, pc, ci = cmh(cl, nc, col)
    b2.append(dict(metric=col, called_pct=100*cl[col].mean(), nc_pct=100*nc[col].mean(),
                   crude_OR=orr, crude_p=p, CMH_OR=orp, CMH_lo=ci[0], CMH_hi=ci[1], CMH_p=pc))
for col in ["g4h","altRY"]:
    base, lo, hi = strat_diff(cl, nc, col)
    b2.append(dict(metric=col+"_GCstrat_diff", called_pct=float(cl[col].mean()), nc_pct=float(nc[col].mean()),
                   crude_OR=np.nan, crude_p=np.nan, CMH_OR=base, CMH_lo=lo, CMH_hi=hi, CMH_p=np.nan))
B2 = pd.DataFrame(b2); B2.to_csv(OUT/"B2_structure_261016.tsv", sep="\t", index=False)

# ---------- B5 cognate への近さ（4 検定）+ 単一ミスマッチ 21 変異体 ----------
b5 = []
for thr in (1,2,3):
    e["_f"] = e.hd_aagcccg <= thr; c2 = e[e.called]; n2 = e[~e.called]
    A = int(c2._f.sum()); B = len(c2)-A; Cc = int(n2._f.sum()); D = len(n2)-Cc
    orr, p = stats.fisher_exact([[A,B],[Cc,D]]); orp, pc, ci = cmh(c2, n2, "_f")
    b5.append(dict(test=f"hd<={thr}", called_n=A, called_pct=100*A/len(c2), nc_pct=100*Cc/len(n2),
                   crude_OR=orr, crude_p=p, CMH_OR=orp, CMH_lo=ci[0], CMH_hi=ci[1], CMH_p=pc))
b5.append(dict(test="hd median", called_n=len(cl), called_pct=float(cl.hd_aagcccg.median()),
               nc_pct=float(nc.hd_aagcccg.median()), crude_OR=np.nan,
               crude_p=stats.mannwhitneyu(cl.hd_aagcccg, nc.hd_aagcccg)[1],
               CMH_OR=np.nan, CMH_lo=np.nan, CMH_hi=np.nan, CMH_p=np.nan))
pd.DataFrame(b5).to_csv(OUT/"B5_near_cognate_261016.tsv", sep="\t", index=False)
m = (e.hd_aagcccg == 1).values
bb = best[m]; pm = np.argmax(bb != MOTB, axis=1); pbase = bb[np.arange(m.sum()), pm]
sub = pd.DataFrame({"pos1": pm+1, "base": [chr(x) for x in pbase], "called": e.called.values[m]})
tab = sub.groupby(["pos1","base"]).agg(n=("called","size"), k=("called","sum")).reset_index()
tab["occupancy_pct"] = 100*tab.k/tab.n; tab["cognate_base"] = [MOT[p-1] for p in tab.pos1]
KK = int(sub.called.sum()); NN = len(sub)
tab["OR_vs_other_hd1"] = [stats.fisher_exact([[int(r.k), int(r.n-r.k)],
    [KK-int(r.k), (NN-int(r.n))-(KK-int(r.k))]])[0] for r in tab.itertuples()]
tab["p"] = [stats.fisher_exact([[int(r.k), int(r.n-r.k)],
    [KK-int(r.k), (NN-int(r.n))-(KK-int(r.k))]])[1] for r in tab.itertuples()]
tab["q_BH"] = multipletests(tab.p, method="fdr_bh")[1]
tab.sort_values("occupancy_pct", ascending=False).to_csv(
    OUT/"single_mismatch_variant_occupancy_261016.tsv", sep="\t", index=False)

# ---------- 3' C 連続の用量反応と G 連続の内部対照（2 検定） ----------
dose = []
for col in ("crun3","grun3"):
    for v in range(0, 7):
        mm = e[col].values >= v
        dose.append(dict(feature=f"3' {col[0].upper()}-run >= {v}", n=int(mm.sum()),
                         k=int(e.called.values[mm].sum()), rate_pct=100*e.called.values[mm].mean()))
    e["_f"] = e[col] >= 3; orp, pc, ci = cmh(e[e.called], e[~e.called], "_f")
    dose.append(dict(feature=f"CMH(GC decile) {col}>=3", n=np.nan, k=np.nan, rate_pct=np.nan,
                     CMH_OR=orp, CMH_lo=ci[0], CMH_hi=ci[1], CMH_p=pc))
pd.DataFrame(dose).to_csv(OUT/"crun_dose_response_261016.tsv", sep="\t", index=False)

# ---------- 多変量（全適格アデニン） ----------
Xg = sm.add_constant(pd.DataFrame({"gc50": e.gc50.values, "log_depth": np.log10(e.depth.values),
    "homopolymer": e.homopolymer.astype(int).values, "g4_motif": e.g4_motif.astype(int).values,
    "intergenic": e.intergenic.astype(int).values, "promoter": e.promoter.astype(int).values,
    "repeat": e.repeat.astype(int).values, "altRY": e.altRY.values,
    "hd_le2": (e.hd_aagcccg <= 2).astype(int).values, "crun3_ge3": (e.crun3 >= 3).astype(int).values}))
mg = sm.Logit(e.called.astype(int).values, Xg.values).fit(disp=0)
pd.DataFrame({"term": Xg.columns, "OR": np.exp(mg.params), "p": mg.pvalues}).to_csv(
    OUT/"joint_logit_261016.tsv", sep="\t", index=False)

# ---------- B6 5' 塩基規則（1 検定） ----------
selc = e[e.crun3 >= 3]
t6 = selc.groupby("b5p").agg(n=("called","size"), k=("called","sum"))
t6["rate_pct"] = 100*t6.k/t6.n
e["_f"] = e.b5p.isin(list("ACT")); s2 = e[e.crun3 >= 3]
orp, pc, ci = cmh(s2[s2.called], s2[~s2.called], "_f")
t6.loc["H_vs_G_CMH"] = [np.nan, np.nan, np.nan]
t6.to_csv(OUT/"B6_5prime_base_261016.tsv", sep="\t")
print("B6 H vs G CMH OR=%.3f CI=[%.3f,%.3f] p=%.3g" % (orp, ci[0], ci[1], pc))

# ---------- B4 発現（2 検定） ----------
fc = pd.read_csv("/Users/okaban/bioinfo/rna-seq/03_quantification/data/featureCounts/featureCounts_M145_s2.txt",
                 sep="\t", comment="#")
cnt = fc.iloc[:, 6:].values; lenkb = fc.Length.values/1000
rpk = cnt/lenkb[:, None]; tpm = 1e6*rpk/rpk.sum(0)
gx = pd.DataFrame({"gene": fc.Geneid, "chrom": fc.Chr.astype(str), "start": fc.Start.astype(str),
                   "end": fc.End.astype(str), "length": fc.Length, "tpm_T1": tpm[:, :3].mean(1)})
gx = gx[gx.chrom.str.contains("NC_003888.3")].copy()
gx["s0"] = gx.start.str.split(";").str[0].astype(int)-1
gx["e0"] = gx.end.str.split(";").str[-1].astype(int)
gidx = np.full(L, -1, np.int32); gx = gx.reset_index(drop=True)
for j, (aa, bb2) in enumerate(zip(gx.s0.values, gx.e0.values)): gidx[aa:bb2] = j
e["gene_j"] = gidx[e.start.to_numpy()]
gg = e[e.gene_j >= 0].groupby("gene_j").agg(n_elig=("called","size"), n_site=("called","sum"),
                                            gc_gene=("gc50","mean"))
G = gx.join(gg, how="inner"); G["has_site"] = G.n_site > 0; G["logtpm"] = np.log10(G.tpm_T1+1)
G["d_oriC"] = np.abs(G.s0-4270000); cent = ((G.s0+G.e0)/2).values
G["gene_density"] = [int(((cent > c-25000) & (cent < c+25000)).sum()) for c in cent]
G["dens_site"] = G.n_site/G.n_elig
def partial_spearman(df, x, y, covs):
    rk = df[[x,y]+covs].rank()
    rx = sm.OLS(rk[x], sm.add_constant(rk[covs])).fit().resid
    ry = sm.OLS(rk[y], sm.add_constant(rk[covs])).fit().resid
    return stats.spearmanr(rx, ry)
r_raw, p_raw = stats.spearmanr(G.dens_site, G.logtpm)
r_p, p_p = partial_spearman(G, "dens_site", "logtpm", ["gc_gene","d_oriC","gene_density","length"])
pd.DataFrame([dict(test="MWU logTPM has_site vs not", stat=np.nan,
                   p=stats.mannwhitneyu(G[G.has_site].logtpm, G[~G.has_site].logtpm)[1]),
              dict(test="Spearman dens_site~logTPM (raw)", stat=r_raw, p=p_raw),
              dict(test="partial (gene GC, oriC dist, gene density, length)", stat=r_p, p=p_p)]
            ).to_csv(OUT/"B4_expression_261016.tsv", sep="\t", index=False)

cl = e[e.called]
cl.drop(columns=[c for c in ["_f","in_motif"] if c in cl]).to_csv(
    OUT/"unassigned_T1_sites_annotated_261016.tsv", sep="\t", index=False)
print("done")
