#!/usr/bin/env python3
"""125 - 6mA deep dive: B1 replication dilution, B2 bipartite/long/asymmetric motif
search, B5 strand x expression interaction, B6 4mC-conditioned 6mA occupancy.

2026-11-04.  Companion: b4_orphan_mtase_261104.py (B4), and the report
Writing/sixmA_open_questions_261104.md.

Nothing upstream is recomputed.  Inputs:
  adenine cache  117_minor_motifs_261104/cache/{pos_all,str_all,Nv,Nm,Nf,ctx,
                 fg_unassigned6mA_*,bg_unassigned6mA_*,is_major}.npy
  cytosine rows  cache/<lib>_c.tsv, pulled once from the canonical pileup by:
      for lib in 1-1 1-2 1-3 2-1 2-3 2-4 3-2 3-3 3-4; do
        awk -v OFS='\t' 'NR==FNR{k[$1"_"$2]=1; next}
          $4=="21839" && (($2"_"$6) in k){print $2,$6,$10,$12,$16}' \
          cache/c_targets.tsv $PILEUP/${lib}_pileup.bed > cache/${lib}_c.tsv &
      done; wait
    (the 1.3 GB/library pileup BEDs are never read into a DataFrame)
  genome         /Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa
  genes / TPM    05_annotation/.../gene_annotation_basic.tsv, featureCounts_M145_s2.txt

Conventions
  occupancy  pooled Nmod/Nvalid >= 0.5 at pooled true depth (Nvalid+Nfail) >= 30
  local GC   win_excl: +-25 bp minus the strand-oriented motif offsets and the 5'
             determinant -> 43 bp retained for all three mark-motif pairs
             (6mA A1 of AAGCCCG -1..+6; 4mC C5 of AAGCCCG -5..+2; 4mC C5 of TGGCCGGC -4..+3)
  oriC       4,270,310 (as in 112_integration_blanks_261019/scripts/core.py)
  B2 fg/bg   identical to 117_: fg = unassigned 6mA reproducible in >= 2 timepoints
             (1,492), bg = 2,412,468 callable adenines with AAGCCCG A1/A2 removed
"""
import os, re, gc, time, resource, numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy import stats
from scipy.stats import hypergeom
from statsmodels.stats.proportion import proportion_confint
from statsmodels.stats.power import TTestIndPower

RNA = "/Users/okaban/bioinfo/rna-seq"
AN  = f"{RNA}/11_epigenome_integration/analysis"
OUT = f"{AN}/125_sixmA_deepdive_261104"
CACHE, TAB = f"{OUT}/cache", f"{OUT}/tables"
C117 = f"{AN}/117_minor_motifs_261104/cache"
REF_FA = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
LIBS = ["1-1","1-2","1-3","2-1","2-3","2-4","3-2","3-3","3-4"]
TPS  = {"T1":[0,1,2], "T2":[3,4,5], "T3":[6,7,8]}
ORIC, WINB, WIN = 4_270_310, 200_000, 12
os.makedirs(TAB, exist_ok=True)

genome = "".join(l.strip() for l in open(REF_FA) if not l.startswith(">")).upper()
L = len(genome); assert L == 8_667_507
grc = genome.translate(str.maketrans("ACGT","TGCA"))[::-1]
gcode = np.frombuffer(genome.encode(), dtype=np.uint8)
isGC  = ((gcode == ord("G")) | (gcode == ord("C"))).astype(np.int8)

def motif_positions(motif, offsets):
    """0-based strand-aware positions of `offsets` within every occurrence."""
    pp, ss = [], []
    for m in re.finditer("(?=" + motif + ")", genome):
        for o in offsets: pp.append(m.start()+o); ss.append("+")
    for m in re.finditer("(?=" + motif + ")", grc):
        for o in offsets: pp.append(L-1-(m.start()+o)); ss.append("-")
    return np.array(pp, np.int64), np.array(ss, "<U1")

def local_gc_excl(pos, strand, excl_rel, half=25):
    off  = np.arange(-half, half+1)
    keep = np.array([o for o in off if o not in set(excl_rel)])
    plus = (strand == "+")
    idx  = np.clip(pos[:,None] + np.where(plus[:,None], keep[None,:], -keep[None,:]), 0, L-1)
    return isGC[idx].mean(axis=1), len(keep)

# =============================== load caches ================================
pos_all = np.load(f"{C117}/pos_all.npy"); str_all = np.load(f"{C117}/str_all.npy")
Nv = np.load(f"{C117}/Nv.npy"); Nm = np.load(f"{C117}/Nm.npy"); Nf = np.load(f"{C117}/Nf.npy")
key_a = {(p, "+" if s else "-"): i for i, (p, s) in enumerate(zip(pos_all, str_all))}
mp, ms = motif_positions("AAGCCCG", [0])
ia = np.array([key_a[(p, s)] for p, s in zip(mp, ms)])          # all 1,334 present

# replication proxy: mean true depth per 200 kb window / timepoint median
wid = (pos_all // WINB).astype(np.int32); nwin = wid.max()+1
cnt = np.bincount(wid, minlength=nwin)
depth_win = np.vstack([np.bincount(wid, weights=(Nv[i].astype(float)+Nf[i]), minlength=nwin)
                       / np.maximum(cnt,1) for i in range(9)])
log2cn = np.log2(depth_win / np.median(depth_win, axis=1, keepdims=True))
oriter = np.array([depth_win[i, ORIC//WINB] /
                   np.mean([depth_win[i,0], depth_win[i,nwin-1]]) for i in range(9)])
print("ori/ter per library:", np.round(oriter,2))        # 4.20 / 1.00 / 1.17 per timepoint

cfr = {lib: pd.read_csv(f"{CACHE}/{lib}_c.tsv", sep="\t", header=None,
                        names=["pos0","strand","Nvalid","Nmod","Nfail"]
                        ).set_index(["pos0","strand"]) for lib in LIBS}
def cyt_arrays(pos, strand):
    idx = pd.MultiIndex.from_arrays([pos, strand])
    v = np.zeros((9,len(pos)), np.int32); m = v.copy(); f = v.copy()
    for i, lib in enumerate(LIBS):
        r = cfr[lib].reindex(idx)
        v[i] = r.Nvalid.fillna(0); m[i] = r.Nmod.fillna(0); f[i] = r.Nfail.fillna(0)
    return v, m, f

SITES = {}
SITES["6mA_AAGCCCG_A1"]  = dict(pos=mp, strand=ms, Nv=Nv[:,ia], Nm=Nm[:,ia], Nf=Nf[:,ia],
                                excl=list(range(-1,7)))
for name, motif, off, excl in [("4mC_TGGCCGGC_C5","TGGCCGGC",4,range(-4,4)),
                               ("4mC_AAGCCCG_C5","AAGCCCG",4,range(-5,3))]:
    p, s = motif_positions(motif, [off]); v, m_, f = cyt_arrays(p, s)
    SITES[name] = dict(pos=p, strand=s, Nv=v, Nm=m_, Nf=f, excl=list(excl))

# ============================ per-site table (B1) ===========================
rows = []
for name, d in SITES.items():
    pos, strand = d["pos"], d["strand"]
    gc43, _ = local_gc_excl(pos, strand, d["excl"])
    plus = (strand == "+")
    i5 = np.clip(np.where(plus, pos-1, pos+1), 0, L-1)
    b5 = np.array([chr(x) for x in np.where(plus, gcode[i5],
                   [{65:84,67:71,71:67,84:65}.get(x,78) for x in gcode[i5]])])
    for tp, ix in TPS.items():
        nv = d["Nv"][ix].sum(0).astype(float); nm = d["Nm"][ix].sum(0).astype(float)
        td = (d["Nv"][ix] + d["Nf"][ix]).sum(0).astype(float)
        rows.append(pd.DataFrame(dict(
            mark=name, tp=tp, pos=pos, strand=strand,
            level=np.where(nv>0, nm/np.maximum(nv,1), np.nan),
            true_depth=td, nvalid=nv, nmod=nm, gc43=gc43, b5=b5,
            oriC_dist_Mb=np.abs(pos-ORIC)/1e6,
            replichore=np.where(pos>=ORIC,"right","left"),
            log2cn=log2cn[ix][:, pos//WINB].mean(0))))
S = pd.concat(rows, ignore_index=True)
S["callable"] = S.true_depth >= 30
S["occ"] = ((S.level >= 0.5) & S["callable"]).astype(int)
S["leading"] = np.where(((S.strand=="+")&(S.replichore=="right")) |
                        ((S.strand=="-")&(S.replichore=="left")), "leading", "lagging")
S.to_csv(f"{TAB}/b1_site_table_261104.tsv.gz", sep="\t", index=False, compression="gzip")

# B1-1 per-timepoint gradient models (depth-adjusted fit kept as a sensitivity:
# at T1 corr(oriC distance, log2 depth) = -0.971, so depth is the mediator)
res = []
for mark in S["mark"].unique():
    for tp in ["T1","T2","T3"]:
        g = S[(S["mark"]==mark)&(S.tp==tp)&S["callable"]].copy()
        g["logd"] = np.log2(g.true_depth); g["gcz"] = (g.gc43-g.gc43.mean())/g.gc43.std()
        base = "oriC_dist_Mb + gcz + logd" + (" + C(b5)" if mark.startswith("6mA") else "")
        for resp, form, w in [("occ_binary", f"occ ~ {base}", None),
                              ("level_binom", f"level ~ {base}", g.nvalid)]:
            m_ = smf.glm(form, data=g, family=sm.families.Binomial(),
                         **({"freq_weights": w} if w is not None else {})).fit()
            b, se = m_.params["oriC_dist_Mb"], m_.bse["oriC_dist_Mb"]
            res.append(dict(mark=mark, tp=tp, response=resp, predictor="oriC_dist_Mb", n=len(g),
                            beta=b, se=se, ci_lo=b-1.96*se, ci_hi=b+1.96*se,
                            p=m_.pvalues["oriC_dist_Mb"], gc_beta=m_.params["gcz"]))
        m_ = smf.glm("level ~ log2cn + gcz" + (" + C(b5)" if mark.startswith("6mA") else ""),
                     data=g, family=sm.families.Binomial(), freq_weights=g.nvalid).fit()
        b, se = m_.params["log2cn"], m_.bse["log2cn"]
        res.append(dict(mark=mark, tp=tp, response="level_binom", predictor="log2cn", n=len(g),
                        beta=b, se=se, ci_lo=b-1.96*se, ci_hi=b+1.96*se, p=m_.pvalues["log2cn"]))
pd.DataFrame(res).to_csv(f"{TAB}/b1_replication_gradient_models_261104.tsv", sep="\t", index=False)

# B1-2 the decisive test: is the gradient T1-specific?  (no depth adjustment)
rows = []
for mark in S["mark"].unique():
    g = S[(S["mark"]==mark)&S["callable"]].copy()
    g["gcz"] = (g.gc43-g.gc43.mean())/g.gc43.std(); g["isT1"] = (g.tp=="T1").astype(int)
    m_ = smf.glm("level ~ oriC_dist_Mb * isT1 + gcz" + (" + C(b5)" if mark.startswith("6mA") else ""),
                 data=g, family=sm.families.Binomial(), freq_weights=g.nvalid).fit()
    k = "oriC_dist_Mb:isT1"
    rows.append(dict(mark=mark, n=len(g), main_beta=m_.params["oriC_dist_Mb"],
                     inter_beta=m_.params[k], inter_se=m_.bse[k], inter_p=m_.pvalues[k],
                     inter_ci_lo=m_.params[k]-1.96*m_.bse[k],
                     inter_ci_hi=m_.params[k]+1.96*m_.bse[k]))
IT = pd.DataFrame(rows)
IT.to_csv(f"{TAB}/b1_gradient_timepoint_interaction_261104.tsv", sep="\t", index=False)

# B1-3 culture unit (n=9), and B1-4 leading vs lagging
cr = []
for name, d in SITES.items():
    td = d["Nv"] + d["Nf"]; lv = np.where(d["Nv"]>0, d["Nm"]/np.maximum(d["Nv"],1), np.nan)
    for i, lib in enumerate(LIBS):
        ok = td[i] >= 10                      # canonical per-replicate floor
        cr.append(dict(mark=name, lib=lib, tp="T"+str(i//3+1), n=int(ok.sum()),
                       occ_pct=100*(lv[i][ok]>=0.5).mean(),
                       mean_level=100*np.nanmean(lv[i][ok]), oriter=oriter[i]))
CL = pd.DataFrame(cr); CL.to_csv(f"{TAB}/b1_per_culture_occupancy_261104.tsv", sep="\t", index=False)
ll = []
for mark in S["mark"].unique():
    for tp in ["T1","T2","T3"]:
        g = S[(S["mark"]==mark)&(S.tp==tp)&S["callable"]].copy()
        g["gcz"] = (g.gc43-g.gc43.mean())/g.gc43.std(); g["lead"] = (g.leading=="leading").astype(int)
        m_ = smf.glm("level ~ lead + gcz" + (" + C(b5)" if mark.startswith("6mA") else ""),
                     data=g, family=sm.families.Binomial(), freq_weights=g.nvalid).fit()
        ll.append(dict(mark=mark, tp=tp, n_lead=int((g.lead==1).sum()), n_lag=int((g.lead==0).sum()),
                       occ_lead=100*g[g.lead==1].occ.mean(), occ_lag=100*g[g.lead==0].occ.mean(),
                       beta=m_.params["lead"], se=m_.bse["lead"],
                       OR=np.exp(m_.params["lead"]), p=m_.pvalues["lead"]))
pd.DataFrame(ll).to_csv(f"{TAB}/b1_leading_lagging_261104.tsv", sep="\t", index=False)

# B1-5 minimum detectable effects
pw = TTestIndPower(); rows = []
for name in SITES:
    g = CL[CL["mark"]==name]
    sd = np.sqrt((g[g.tp=="T1"].occ_pct.std()**2*2 + g[g.tp!="T1"].occ_pct.std()**2*5)/7)
    dd = pw.solve_power(nobs1=3, ratio=2.0, alpha=0.05, power=0.8, alternative="two-sided")
    rows.append(dict(test="culture-level T1 vs T2/T3", mark=name, pooled_sd_pp=sd, mde_pp=dd*sd,
                     observed_diff_pp=g[g.tp=="T1"].occ_pct.mean()-g[g.tp!="T1"].occ_pct.mean()))
for _, r in IT.iterrows():
    p0 = S[(S["mark"]==r["mark"])&(S.tp=="T1")&S["callable"]].occ.mean()
    lo = np.log(p0/(1-p0))
    f = lambda b: 100*(1/(1+np.exp(-(lo+b*4.27))) - p0)
    rows.append(dict(test="T1-specific oriC gradient (interaction)", mark=r["mark"],
                     pooled_sd_pp=np.nan, mde_pp=f(1.96*r.inter_se),
                     observed_diff_pp=f(r.inter_beta)))
pd.DataFrame(rows).to_csv(f"{TAB}/b1_power_261104.tsv", sep="\t", index=False)

del Nv, Nm, Nf; gc.collect()

# ===================== B2 bipartite / long / asymmetric =====================
ctx = np.load(f"{C117}/ctx.npy")                       # (N,25) int8, offsets -12..+12
FG = np.load(f"{C117}/fg_unassigned6mA_repro.npy")     # 1,499 -> 1,492 inside bg
BG = np.load(f"{C117}/bg_unassigned6mA_repro.npy")     # 2,412,468
cs = np.concatenate([[0], np.cumsum(isGC.astype(np.int64))])
def win_gc(pos, plusmask, half=25, excl=range(-1,7)):
    lo = np.clip(pos-half,0,L); hi = np.clip(pos+half+1,0,L)
    tot = cs[hi]-cs[lo]; n = hi-lo
    for o in excl:
        i = np.where(plusmask, pos+o, pos-o); ok = (i>=0)&(i<L); i = np.clip(i,0,L-1)
        tot = tot - np.where(ok, isGC[i], 0); n = n - ok.astype(np.int64)
    return tot/np.maximum(n,1)
gc_all = win_gc(pos_all, str_all==1)
depth_all = (np.load(f"{C117}/Nv.npy") + np.load(f"{C117}/Nf.npy")).sum(0)

def word_codes(cols):
    sub = ctx[:, [c+WIN for c in cols]]
    valid = (sub < 4).all(axis=1)
    code = np.zeros(len(sub), np.int64)
    for i in range(sub.shape[1]): code = code*4 + np.where(sub[:,i]<4, sub[:,i], 0)
    return code, valid

ARCH = []                                   # pre-registered architecture space
for b1 in (3,4):
    for b2 in (3,4):
        for gapn in range(1,11):
            for j in range(b1+b2):
                offs = [-j+i for i in range(b1)] + [-j+b1+gapn+i for i in range(b2)]
                if min(offs) >= -WIN and max(offs) <= WIN and 0 in offs:
                    ARCH.append(("bipartite", b1, b2, gapn, j, tuple(offs)))
for k in (9,10,11):
    for j in range(k):
        offs = [-j+i for i in range(k)]
        if min(offs) >= -WIN and max(offs) <= WIN:
            ARCH.append(("contiguous", k, 0, 0, j, tuple(offs)))
n_words = sum(4**(len(a[5])-1) for a in ARCH)       # 15,543,296

fg_idx = np.flatnonzero(FG & BG); bg_idx = np.flatnonzero(BG)
gq = np.searchsorted(np.quantile(gc_all[bg_idx], [.2,.4,.6,.8]), gc_all, side="right")
dq = np.searchsorted(np.quantile(depth_all[bg_idx], [.2,.4,.6,.8]), depth_all, side="right")
strat = gq*5 + dq
bg_by_s = {s: bg_idx[strat[bg_idx]==s] for s in range(25)}
fg_cnt  = {s: int((strat[fg_idx]==s).sum()) for s in range(25)}
rng = np.random.default_rng(20261104); NPERM = 1000
PERM = np.vstack([np.concatenate([rng.choice(bg_by_s[s], fg_cnt[s], replace=False)
                                  for s in range(25) if fg_cnt[s]>0]) for _ in range(NPERM)])

I2B = np.array(list("ACGT"))
decode = lambda c,k: "".join(I2B[(c//(4**i))%4] for i in range(k-1,-1,-1))
MIN_BG, MIN_FG = 20, 5
perm_max = np.zeros(NPERM); obs = []
for kind, p1, p2, gapn, j, offs in ARCH:
    k = len(offs); code, valid = word_codes(offs)
    bgi = bg_idx[valid[bg_idx]]; fgi = fg_idx[valid[fg_idx]]
    bgc = np.bincount(code[bgi], minlength=4**k); NB, NF = len(bgi), len(fgi)
    u, n = np.unique(code[fgi], return_counts=True); bb = bgc[u]
    keep = (bb>=MIN_BG) & (n>=MIN_FG)
    if keep.any():
        pv = hypergeom.sf(n[keep]-1, NB, bb[keep], NF)
        for w, nf, nbw, pp in zip(u[keep], n[keep], bb[keep], pv):
            obs.append(dict(kind=kind, b1=p1, b2=p2, gap=gapn, a_offset=j, k=k,
                            offsets=",".join(map(str,offs)), word=decode(int(w),k),
                            n_fg=int(nf), n_bg=int(nbw), occ_pct=100*nf/nbw, p=pp,
                            neglog10p=-np.log10(max(pp,1e-300))))
    for b in range(NPERM):                      # FWER by max statistic over the whole space
        uu, nn = np.unique(code[PERM[b]], return_counts=True); b2c = bgc[uu]
        kk = (b2c>=MIN_BG) & (nn>=MIN_FG)
        if kk.any():
            v = -np.log10(np.maximum(hypergeom.sf(nn[kk]-1, NB, b2c[kk], NF), 1e-300)).max()
            if v > perm_max[b]: perm_max[b] = v
OB = pd.DataFrame(obs).sort_values("neglog10p", ascending=False)
OB.to_csv(f"{TAB}/b2_bipartite_scan_261104.tsv", sep="\t", index=False)
np.save(f"{CACHE}/b2_perm_max.npy", perm_max)
thr = float(np.quantile(perm_max, 0.95))                                    # 5.98
pas = OB[OB.neglog10p > thr]
tt = pas.sort_values("occ_pct", ascending=False).head(9).copy()
lo, hi = proportion_confint(tt.n_fg, tt.n_bg, method="wilson")
tt["ci_lo_pct"], tt["ci_hi_pct"] = 100*lo, 100*hi
tt.to_csv(f"{TAB}/b2_top_candidates_261104.tsv", sep="\t", index=False)

fl = []                                   # detection floor at the FWER threshold
for nb in (20,50,100,200,500,1000,2000,5000):
    nf = next((c for c in range(1, min(nb,len(fg_idx))+1)
               if -np.log10(max(hypergeom.sf(c-1,len(bg_idx),nb,len(fg_idx)),1e-300)) > thr), 0)
    fl.append(dict(motif_instances_n_bg=nb, min_detectable_n_fg=nf,
                   min_detectable_occupancy_pct=round(100*nf/nb,2)))
pd.DataFrame(fl).to_csv(f"{TAB}/b2_detection_floor_261104.tsv", sep="\t", index=False)

def occ_of(pattern, a_rel):               # near-cognate panel, same fg/bg rule
    offs = [i-a_rel for i in range(len(pattern))]
    if min(offs) < -WIN or max(offs) > WIN: return None
    sub = ctx[:, [o+WIN for o in offs]]; ok = np.ones(len(sub), bool)
    for i, ch in enumerate(pattern):
        if ch == "N": continue
        allowed = [{"A":0,"C":1,"G":2,"T":3}[b] for b in
                   (ch if ch in "ACGT" else {"S":"CG","W":"AT","R":"AG","Y":"CT","H":"ACT"}[ch])]
        ok &= np.isin(sub[:,i], allowed)
    nb = int((ok&BG).sum()); nf = int((ok&FG&BG).sum())
    return nf, nb, (100*nf/nb if nb else np.nan)
pd.DataFrame([dict(pattern=p, a_index=j, note=t, n_fg=r[0], n_bg=r[1], occ_pct=round(r[2],2))
              for p,j,t in [("AACCCG",0,"near-cognate, A1"),("AACCCG",1,"near-cognate, A2"),
                            ("CGGCAACCCG",5,"top contiguous 10-mer hit"),
                            ("AANCCCG",0,"A1, N at +2"),("AASCCCG",0,"A1, S at +2"),
                            ("AAGCCC",0,"motif minus 3' G")]
              for r in [occ_of(p,j)] if r]).to_csv(
    f"{TAB}/b2_near_cognate_panel_261104.tsv", sep="\t", index=False)

sa = []                                   # genomic-strand asymmetry
g6 = S[(S["mark"]=="6mA_AAGCCCG_A1") & S["callable"]]
for tp in ["T1","T2","T3"]:
    h = g6[g6.tp==tp]; a = h[h.strand=="+"]; b = h[h.strand=="-"]
    orr, p = stats.fisher_exact([[a.occ.sum(), len(a)-a.occ.sum()],
                                 [b.occ.sum(), len(b)-b.occ.sum()]])
    sa.append(dict(set="AAGCCCG A1 (major motif)", tp=tp, n_plus=len(a), n_minus=len(b),
                   occ_plus=100*a.occ.mean(), occ_minus=100*b.occ.mean(), OR=orr, p=p))
pl = (str_all==1); mm = FG & BG
nfp, nbp = int((mm&pl).sum()), int((BG&pl).sum())
nfm, nbm = int((mm&~pl).sum()), int((BG&~pl).sum())
orr, p = stats.fisher_exact([[nfp, nbp-nfp], [nfm, nbm-nfm]])
sa.append(dict(set="unassigned 6mA (repro)", tp="pooled", n_plus=nbp, n_minus=nbm,
               occ_plus=100*nfp/nbp, occ_minus=100*nfm/nbm, OR=orr, p=p))
pd.DataFrame(sa).to_csv(f"{TAB}/b2_strand_asymmetry_261104.tsv", sep="\t", index=False)

# ============ B6 4mC-conditioned 6mA occupancy (site and culture units) ============
inst = ([(m.start(),"+",m.start(),m.start()+4) for m in re.finditer("(?=AAGCCCG)", genome)] +
        [(m.start(),"-",L-1-m.start(),L-1-(m.start()+4)) for m in re.finditer("(?=AAGCCCG)", grc)])
IN = pd.DataFrame(inst, columns=["mstart","strand","posA1","posC5"])
a6 = S[S["mark"]=="6mA_AAGCCCG_A1"][["pos","strand","tp","occ","level","callable","gc43","b5"]]
c4 = S[S["mark"]=="4mC_AAGCCCG_C5"][["pos","strand","tp","occ","level","callable"]]
J = (IN.merge(a6, left_on=["posA1","strand"], right_on=["pos","strand"])
       .merge(c4, left_on=["posC5","strand","tp"], right_on=["pos","strand","tp"],
              suffixes=("_A","_C")))
J = J[J.callable_A & J.callable_C]
rows, adj = [], []
for tp in ["T1","T2","T3"]:
    g = J[J.tp==tp]; p1 = g[g.occ_C==1]; p0 = g[g.occ_C==0]
    orr, p = stats.fisher_exact([[p1.occ_A.sum(), len(p1)-p1.occ_A.sum()],
                                 [p0.occ_A.sum(), len(p0)-p0.occ_A.sum()]])
    rows.append(dict(tp=tp, n_4mC_pos=len(p1), n_4mC_neg=len(p0),
                     occ6mA_given_4mC_pos=100*p1.occ_A.mean(),
                     occ6mA_given_4mC_neg=100*p0.occ_A.mean(),
                     diff_pp=100*(p1.occ_A.mean()-p0.occ_A.mean()), OR=orr, p=p,
                     lvl6mA_given_4mC_pos=100*p1.level_A.mean(),
                     lvl6mA_given_4mC_neg=100*p0.level_A.mean()))
    h = g.copy(); h["gcz"] = (h.gc43-h.gc43.mean())/h.gc43.std()
    m_ = smf.glm("occ_A ~ occ_C + gcz + C(b5)", data=h, family=sm.families.Binomial()).fit()
    b, se = m_.params["occ_C"], m_.bse["occ_C"]
    adj.append(dict(tp=tp, n=len(h), beta_4mC=b, se=se, OR_adj=np.exp(b),
                    ci_lo=np.exp(b-1.96*se), ci_hi=np.exp(b+1.96*se), p=m_.pvalues["occ_C"]))
pd.DataFrame(rows).to_csv(f"{TAB}/b6_conditional_occupancy_site_261104.tsv", sep="\t", index=False)
pd.DataFrame(adj).to_csv(f"{TAB}/b6_conditional_adjusted_261104.tsv", sep="\t", index=False)
d6, d4 = SITES["6mA_AAGCCCG_A1"], SITES["4mC_AAGCCCG_C5"]
k6 = {(p,s):i for i,(p,s) in enumerate(zip(d6["pos"], d6["strand"]))}
k4 = {(p,s):i for i,(p,s) in enumerate(zip(d4["pos"], d4["strand"]))}
i6 = np.array([k6[(p,s)] for p,s in zip(IN.posA1, IN.strand)])
i4 = np.array([k4[(p,s)] for p,s in zip(IN.posC5, IN.strand)])
rows = []
for li, lib in enumerate(LIBS):
    tdA = (d6["Nv"][li]+d6["Nf"][li])[i6]; tdC = (d4["Nv"][li]+d4["Nf"][li])[i4]
    lvA = np.where(d6["Nv"][li][i6]>0, d6["Nm"][li][i6]/np.maximum(d6["Nv"][li][i6],1), np.nan)
    lvC = np.where(d4["Nv"][li][i4]>0, d4["Nm"][li][i4]/np.maximum(d4["Nv"][li][i4],1), np.nan)
    ok = (tdA>=10) & (tdC>=10)
    rows.append(dict(lib=lib, tp="T"+str(li//3+1), n_pos=int((ok&(lvC>=0.5)).sum()),
                     n_neg=int((ok&(lvC<0.5)).sum()),
                     occ_given_pos=100*np.mean(lvA[ok&(lvC>=0.5)]>=0.5),
                     occ_given_neg=100*np.mean(lvA[ok&(lvC<0.5)]>=0.5)))
B6c = pd.DataFrame(rows); B6c["diff_pp"] = B6c.occ_given_pos - B6c.occ_given_neg
B6c.to_csv(f"{TAB}/b6_conditional_occupancy_culture_261104.tsv", sep="\t", index=False)

# ============== B5 strand x expression (single-strand access) ===============
GA = pd.read_csv(f"{RNA}/05_annotation/analysis/05_annotation_260128_v1/tables/"
                 "gene_annotation_basic.tsv", sep="\t")
GA = GA[GA.contig=="NC_003888.3"].copy()
fc = pd.read_csv(f"{RNA}/03_quantification/data/featureCounts/featureCounts_M145_s2.txt",
                 sep="\t", comment="#")
fc.columns = list(fc.columns[:6]) + LIBS
rpk = fc[LIBS].to_numpy(float) / (fc["Length"].to_numpy(float)[:,None]/1000)
TPM = pd.DataFrame(rpk/rpk.sum(0)*1e6, columns=LIBS); TPM["gene_id"] = fc["Geneid"].values
for t, c in zip(["T1","T2","T3"], [LIBS[:3], LIBS[3:6], LIBS[6:]]): TPM[t] = TPM[c].mean(1)
GA = GA.merge(TPM[["gene_id","T1","T2","T3"]], on="gene_id", how="left")
GA["glen"] = GA.end - GA.start + 1
order = np.argsort(GA.start.to_numpy()); gs_s = GA.start.to_numpy()[order]; ge_s = GA.end.to_numpy()[order]
def host_gene(pos):
    j = np.searchsorted(gs_s, pos, side="right") - 1
    out = np.full(len(pos), -1)
    ok = (j>=0) & (pos <= np.where(j>=0, ge_s[np.clip(j,0,len(ge_s)-1)], -1))
    out[ok] = order[j[ok]]; return out
g6 = S[S["mark"]=="6mA_AAGCCCG_A1"].copy(); g6["gene_idx"] = host_gene(g6.pos.to_numpy())
g6 = g6[(g6.gene_idx>=0) & g6["callable"]].copy()
g6["gene_strand"] = GA.strand.to_numpy()[g6.gene_idx]
g6["glen"] = GA.glen.to_numpy()[g6.gene_idx]
g6["tpm"] = [GA[t].to_numpy()[i] for t,i in zip(g6.tp, g6.gene_idx)]
g6["site_orient"] = np.where(g6.strand==g6.gene_strand, "coding", "template")
rows, inter, quint = [], [], []
for tp in ["T1","T2","T3"]:
    for orient in ["coding","template"]:
        h = g6[(g6.tp==tp)&(g6.site_orient==orient)].copy()
        h["ltpm"]=np.log2(h.tpm+1); h["gcz"]=(h.gc43-h.gc43.mean())/h.gc43.std()
        h["logd"]=np.log2(h.true_depth); h["lglen"]=np.log2(h.glen)
        m_ = smf.glm("occ ~ ltpm + gcz + logd + lglen + C(b5)", data=h,
                     family=sm.families.Binomial()).fit()
        b, se = m_.params["ltpm"], m_.bse["ltpm"]
        rows.append(dict(tp=tp, stratum=orient, n=len(h), beta_logTPM=b, se=se, OR=np.exp(b),
                         ci_lo=np.exp(b-1.96*se), ci_hi=np.exp(b+1.96*se), p=m_.pvalues["ltpm"]))
    h = g6[g6.tp==tp].copy(); h["q"] = pd.qcut(h.tpm, 5, labels=False, duplicates="drop")
    for q in sorted(h.q.dropna().unique()):
        k = h[h.q==q]; a = k[k.site_orient=="coding"]; b_ = k[k.site_orient=="template"]
        if len(a)<10 or len(b_)<10: continue
        tab = [[a.occ.sum(), len(a)-a.occ.sum()], [b_.occ.sum(), len(b_)-b_.occ.sum()]]
        ci = stats.contingency.odds_ratio(tab).confidence_interval()
        quint.append(dict(tp=tp, tpm_quintile=int(q)+1, tpm_median=k.tpm.median(),
                          n_coding=len(a), n_template=len(b_), occ_coding=100*a.occ.mean(),
                          occ_template=100*b_.occ.mean(), OR=stats.fisher_exact(tab)[0],
                          ci_lo=ci.low, ci_hi=ci.high, p=stats.fisher_exact(tab)[1]))
    h["coding"]=(h.site_orient=="coding").astype(int); h["ltpm"]=np.log2(h.tpm+1)
    h["gcz"]=(h.gc43-h.gc43.mean())/h.gc43.std(); h["logd"]=np.log2(h.true_depth)
    h["lglen"]=np.log2(h.glen)
    m_ = smf.glm("occ ~ coding*ltpm + gcz + logd + lglen + C(b5)", data=h,
                 family=sm.families.Binomial()).fit()
    b, se = m_.params["coding:ltpm"], m_.bse["coding:ltpm"]
    inter.append(dict(tp=tp, n=len(h), inter_beta=b, se=se, OR=np.exp(b),
                      ci_lo=np.exp(b-1.96*se), ci_hi=np.exp(b+1.96*se),
                      p=m_.pvalues["coding:ltpm"], coding_main_OR=np.exp(m_.params["coding"])))
pd.DataFrame(rows).to_csv(f"{TAB}/b5_expression_within_strand_261104.tsv", sep="\t", index=False)
pd.DataFrame(quint).to_csv(f"{TAB}/b5_strand_OR_by_expression_quintile_261104.tsv", sep="\t", index=False)
pd.DataFrame(inter).to_csv(f"{TAB}/b5_strand_expression_interaction_261104.tsv", sep="\t", index=False)

# leading/lagging adjusted for host-gene orientation (the obvious confounder)
g6["lead"] = np.where(g6.leading=="leading", 1, 0)
rows = []
for tp in ["T1","T2","T3"]:
    h = g6[g6.tp==tp].copy(); h["coding"]=(h.site_orient=="coding").astype(int)
    h["gcz"]=(h.gc43-h.gc43.mean())/h.gc43.std()
    m1 = smf.glm("occ ~ lead + gcz + C(b5)", data=h, family=sm.families.Binomial()).fit()
    m2 = smf.glm("occ ~ lead + coding + gcz + C(b5)", data=h, family=sm.families.Binomial()).fit()
    rows.append(dict(tp=tp, n=len(h), OR_unadj=np.exp(m1.params["lead"]), p_unadj=m1.pvalues["lead"],
                     OR_adj_geneorient=np.exp(m2.params["lead"]), p_adj=m2.pvalues["lead"],
                     coding_OR=np.exp(m2.params["coding"])))
pd.DataFrame(rows).to_csv(f"{TAB}/b1_leading_lagging_geneadj_261104.tsv", sep="\t", index=False)

print("peak RSS GB", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9, 3))
