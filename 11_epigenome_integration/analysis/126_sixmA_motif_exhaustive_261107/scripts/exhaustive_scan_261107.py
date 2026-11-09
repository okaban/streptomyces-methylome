#!/usr/bin/env python3
"""126 - exhaustive 6mA recognition-sequence search over the spaces 117_/125_ did
not cover.  2026-11-07.

Stages
  arch   contiguous k = 12..16; bipartite 3-4 bp blocks with spacer N11..N20;
         three-block 3-3-3 with two spacers N1..N10.  Same hypergeometric test,
         same MIN_BG/MIN_FG filter and the same GC x depth matched permutation
         null (1,000 draws, seed 20261104, byte-identical construction to 125_)
         so the FWER threshold is directly comparable.
  iupac  full enumeration of IUPAC-degenerate contiguous patterns, k = 4..8,
         every offset covering the modified adenine, up to 3 degenerate letters
         from {R,Y,S,W,K,M,B,D,H,V,N}.  Counts are obtained by folding the full
         concrete k-mer contingency tensor, so every pattern in the space is
         evaluated exactly, not sampled.  A 1,000-draw permutation max over
         1.2 M template instances is not affordable; the significance threshold
         for this stage is Bonferroni over the tested cells (conservative
         relative to the permutation FWER) and the primary criterion, as in
         125_, is OCCUPANCY, not significance.

Inputs (nothing upstream is recomputed)
  117_minor_motifs_261104/cache/{pos_all,str_all,Nv,Nf,fg_unassigned6mA_repro,
      bg_unassigned6mA_repro}.npy
  126_.../cache/ctxw30.npy   +-30 bp strand-oriented context, verified to agree
      with 117_/cache/ctx.npy on the shared +-12 bp window
  fg = 1,492 unassigned 6mA sites reproducible in >= 2 timepoints
  bg = 2,412,468 callable adenines with AAGCCCG A1/A2 removed
"""
import os, sys, time, itertools, numpy as np, pandas as pd
from scipy.stats import hypergeom
from multiprocessing import Pool

AN   = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
OUT  = f"{AN}/126_sixmA_motif_exhaustive_261107"
C117 = f"{AN}/117_minor_motifs_261104/cache"
CACHE, TAB = f"{OUT}/cache", f"{OUT}/tables"
WIDE, MIN_BG, MIN_FG, NPERM, SEED = 30, 20, 5, 1000, 20261104
I2B = np.array(list("ACGT"))

def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# ------------------------------------------------------------------ setup ---
ctxw = np.load(f"{CACHE}/ctxw30.npy", mmap_mode="r")
pos_all = np.load(f"{C117}/pos_all.npy"); str_all = np.load(f"{C117}/str_all.npy")
FG = np.load(f"{C117}/fg_unassigned6mA_repro.npy")
BG = np.load(f"{C117}/bg_unassigned6mA_repro.npy")
fg_idx = np.flatnonzero(FG & BG); bg_idx = np.flatnonzero(BG)
NF_TOT, NB_TOT = len(fg_idx), len(bg_idx)

def _perm_matrix():
    """GC x depth stratified draw, construction identical to 125_ B2."""
    g = "".join(l.strip() for l in open("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa")
                if not l.startswith(">")).upper()
    L = len(g); gc8 = np.frombuffer(g.encode(), np.uint8)
    isGC = ((gc8 == ord("G")) | (gc8 == ord("C"))).astype(np.int8)
    csum = np.concatenate([[0], np.cumsum(isGC.astype(np.int64))])
    pos, pl = pos_all, (str_all == 1)
    lo = np.clip(pos - 25, 0, L); hi = np.clip(pos + 26, 0, L)
    tot = csum[hi] - csum[lo]; n = (hi - lo).astype(np.int64)
    for o in range(-1, 7):
        i = np.where(pl, pos + o, pos - o); ok = (i >= 0) & (i < L); i = np.clip(i, 0, L - 1)
        tot = tot - np.where(ok, isGC[i], 0); n = n - ok.astype(np.int64)
    gc_all = tot / np.maximum(n, 1)
    depth_all = (np.load(f"{C117}/Nv.npy") + np.load(f"{C117}/Nf.npy")).sum(0)
    gq = np.searchsorted(np.quantile(gc_all[bg_idx], [.2, .4, .6, .8]), gc_all, side="right")
    dq = np.searchsorted(np.quantile(depth_all[bg_idx], [.2, .4, .6, .8]), depth_all, side="right")
    strat = gq * 5 + dq
    by_s = {s: bg_idx[strat[bg_idx] == s] for s in range(25)}
    cnt  = {s: int((strat[fg_idx] == s).sum()) for s in range(25)}
    rng = np.random.default_rng(SEED)
    return np.vstack([np.concatenate([rng.choice(by_s[s], cnt[s], replace=False)
                                      for s in range(25) if cnt[s] > 0]) for _ in range(NPERM)])

PERM = None
def _init(perm):
    global PERM; PERM = perm

# -------------------------------------------------------------- stage arch ---
def build_arch():
    A = []
    for k in range(12, 17):                       # (2) contiguous 12-16
        for j in range(k):
            o = tuple(range(-j, -j + k))
            if min(o) >= -WIDE and max(o) <= WIDE: A.append(("contiguous", k, 0, 0, 0, j, o))
    for b1 in (3, 4):                             # (3) bipartite, spacer N11-N20
        for b2 in (3, 4):
            for g in range(11, 21):
                for j in range(b1 + b2):
                    o = tuple([-j + i for i in range(b1)] +
                              [-j + b1 + g + i for i in range(b2)])
                    if min(o) >= -WIDE and max(o) <= WIDE: A.append(("bipartite", b1, b2, g, 0, j, o))
    for g1 in range(1, 11):                       # (4) three-block 3-N-3-N-3
        for g2 in range(1, 11):
            for j in range(9):
                o = tuple([-j + i for i in range(3)] +
                          [-j + 3 + g1 + i for i in range(3)] +
                          [-j + 6 + g1 + g2 + i for i in range(3)])
                if min(o) >= -WIDE and max(o) <= WIDE: A.append(("three_block", 3, 3, g1, g2, j, o))
    return A

def codes_of(offs, rows):
    sub = ctxw[np.ix_(rows, [o + WIDE for o in offs])]
    valid = (sub < 4).all(1)
    c = np.zeros(len(sub), np.int64)
    for i in range(sub.shape[1]):
        c = c * 4 + np.where(sub[:, i] < 4, sub[:, i], 0)
    return c, valid

def decode(c, k):
    return "".join(I2B[(c // (4 ** i)) % 4] for i in range(k - 1, -1, -1))

def scan_arch(a):
    kind, b1, b2, g1, g2, j, offs = a
    k = len(offs)
    cb, vb = codes_of(offs, bg_idx); cf, vf = codes_of(offs, fg_idx)
    cb = cb[vb]; cf = cf[vf]; NB, NF = len(cb), len(cf)
    ub, nb = np.unique(cb, return_counts=True)
    uf, nf = np.unique(cf, return_counts=True)
    p = np.searchsorted(ub, uf); p = np.clip(p, 0, len(ub) - 1)
    bb = np.where(ub[p] == uf, nb[p], 0)
    keep = (bb >= MIN_BG) & (nf >= MIN_FG)
    hits = []
    if keep.any():
        pv = hypergeom.sf(nf[keep] - 1, NB, bb[keep], NF)
        for w, a_, b_, pp in zip(uf[keep], nf[keep], bb[keep], pv):
            hits.append(dict(kind=kind, b1=b1, b2=b2, gap1=g1, gap2=g2, a_offset=j, k=k,
                             offsets=",".join(map(str, offs)), word=decode(int(w), k),
                             n_fg=int(a_), n_bg=int(b_), occ_pct=100 * a_ / b_, p=pp,
                             neglog10p=-np.log10(max(pp, 1e-300))))
    # permutation max over this architecture
    pm = np.zeros(NPERM)
    code_by_row = np.full(len(BG), -1, np.int64); code_by_row[bg_idx[vb]] = cb
    pc = code_by_row[PERM]                              # (NPERM, nfg)
    for r in range(NPERM):
        v = pc[r]; v = v[v >= 0]
        if len(v) == 0: continue
        uu, nn = np.unique(v, return_counts=True)
        q = np.clip(np.searchsorted(ub, uu), 0, len(ub) - 1)
        b2c = np.where(ub[q] == uu, nb[q], 0)
        kk = (b2c >= MIN_BG) & (nn >= MIN_FG)
        if kk.any():
            pm[r] = -np.log10(np.maximum(hypergeom.sf(nn[kk] - 1, NB, b2c[kk], len(v)), 1e-300)).max()
    return hits, pm

# ------------------------------------------------------------- stage iupac ---
DEG = {"R": (0, 2), "Y": (1, 3), "S": (1, 2), "W": (0, 3), "K": (2, 3), "M": (0, 1),
       "B": (1, 2, 3), "D": (0, 2, 3), "H": (0, 1, 3), "V": (0, 1, 2), "N": (0, 1, 2, 3)}
DEGL = list(DEG)

KEEP_OCC, KEEP_NLP = 3.0, 5.98     # keep a row if it is occupancy-relevant or FWER-relevant

def scan_iupac(task):
    """Exact enumeration of every IUPAC pattern of length k with the modified A at
    index j and up to `maxdeg` degenerate letters.  Counts come from folding the
    full concrete k-mer contingency tensor, so no pattern is sampled or skipped."""
    k, j, maxdeg = task
    offs = tuple(range(-j, -j + k))
    cb, vb = codes_of(offs, bg_idx); cf, vf = codes_of(offs, fg_idx)
    cb = cb[vb]; cf = cf[vf]; NB, NF = len(cb), len(cf)
    BGT = np.bincount(cb, minlength=4 ** k).reshape((4,) * k).astype(np.int64)
    FGT = np.bincount(cf, minlength=4 ** k).reshape((4,) * k).astype(np.int64)
    keep = []; n_cells = 0; n_patterns = 0; max_occ = 0.0
    for d in range(0, maxdeg + 1):
        for posset in itertools.combinations(range(k), d):
            for letters in itertools.product(DEGL, repeat=d):
                bt, ft = BGT, FGT
                for pi, lt in zip(posset, letters):
                    sl = [slice(None)] * k; sl[pi] = list(DEG[lt])
                    bt = bt[tuple(sl)].sum(axis=pi, keepdims=True)
                    ft = ft[tuple(sl)].sum(axis=pi, keepdims=True)
                n_patterns += bt.size
                m = (ft >= MIN_FG) & (bt >= MIN_BG)
                n_cells += int(m.sum())
                if not m.any(): continue
                fv = ft[m].astype(np.int64); bv = bt[m].astype(np.int64)
                occ = 100.0 * fv / bv
                max_occ = max(max_occ, float(occ.max()))
                pre = (occ >= KEEP_OCC)
                pv = hypergeom.sf(fv - 1, NB, bv, NF)
                nlp = -np.log10(np.maximum(pv, 1e-300))
                sel = pre | (nlp >= KEEP_NLP)
                if not sel.any(): continue
                ix = np.argwhere(m)[sel]
                for row, a_, b_, o_, n_ in zip(ix, fv[sel], bv[sel], occ[sel], nlp[sel]):
                    pat = list(I2B[list(row)])
                    for pi, lt in zip(posset, letters): pat[pi] = lt
                    keep.append((k, j, d, "".join(pat), int(a_), int(b_), float(o_), float(n_)))
    return keep, n_cells, n_patterns, max_occ, (k, j)

# ------------------------------------------------------------------- main ---
if __name__ == "__main__":
    stage = sys.argv[1]
    os.makedirs(TAB, exist_ok=True)
    log("fg", NF_TOT, "bg", NB_TOT)
    if stage == "arch":
        perm = _perm_matrix(); log("PERM", perm.shape)
        A = build_arch(); nw = sum(4 ** len(a[6]) for a in A)
        log("architectures", len(A), "words", nw)
        with Pool(10, initializer=_init, initargs=(perm,)) as pool:
            allhits = []; pmax = np.zeros(NPERM); done = 0
            for h, pm in pool.imap_unordered(scan_arch, A, chunksize=4):
                allhits += h; pmax = np.maximum(pmax, pm); done += 1
                if done % 100 == 0: log("arch", done, "/", len(A), "hits", len(allhits))
        H = pd.DataFrame(allhits).sort_values("neglog10p", ascending=False)
        H.to_csv(f"{TAB}/a_arch_scan_261107.tsv", sep="\t", index=False)
        np.save(f"{CACHE}/a_arch_perm_max.npy", pmax)
        thr = float(np.quantile(pmax, 0.95))
        log("FWER thr", thr, "passing", int((H.neglog10p > thr).sum()),
            "max occ among passing", float(H[H.neglog10p > thr].occ_pct.max()) if (H.neglog10p > thr).any() else None)
        fl = [dict(motif_instances_n_bg=nb,
                   min_detectable_n_fg=next((c for c in range(1, min(nb, NF_TOT) + 1)
                       if -np.log10(max(hypergeom.sf(c - 1, NB_TOT, nb, NF_TOT), 1e-300)) > thr), 0))
              for nb in (20, 50, 100, 200, 500, 1000, 2000, 5000, 20000)]
        for r in fl: r["min_detectable_occupancy_pct"] = round(100 * r["min_detectable_n_fg"] / r["motif_instances_n_bg"], 3)
        pd.DataFrame(fl).assign(space="arch", fwer_thr=thr).to_csv(
            f"{TAB}/a_arch_detection_floor_261107.tsv", sep="\t", index=False)
        pd.DataFrame([dict(space="contiguous_12_16", n_arch=sum(1 for a in A if a[0] == "contiguous"),
                           n_words=sum(4 ** len(a[6]) for a in A if a[0] == "contiguous")),
                      dict(space="bipartite_N11_N20", n_arch=sum(1 for a in A if a[0] == "bipartite"),
                           n_words=sum(4 ** len(a[6]) for a in A if a[0] == "bipartite")),
                      dict(space="three_block_3N3N3", n_arch=sum(1 for a in A if a[0] == "three_block"),
                           n_words=sum(4 ** len(a[6]) for a in A if a[0] == "three_block"))]).to_csv(
            f"{TAB}/a_arch_space_261107.tsv", sep="\t", index=False)
        log("arch done")
    elif stage == "iupac":
        tasks = [(k, j, 3) for k in range(4, 9) for j in range(k)]
        with Pool(10) as pool:
            allhits = []; cells = 0; pats = 0; mx = 0.0
            for h, nc, npt, mo, who in pool.imap_unordered(scan_iupac, tasks):
                allhits += h; cells += nc; pats += npt; mx = max(mx, mo)
                log("iupac", who, "kept", len(h), "cells", nc, "max_occ", round(mo, 3))
        H = pd.DataFrame(allhits, columns=["k", "a_offset", "n_deg", "pattern",
                                           "n_fg", "n_bg", "occ_pct", "neglog10p"])
        H["bonf_p"] = np.minimum(10.0 ** (-H.neglog10p) * cells, 1.0)
        H = H.sort_values("occ_pct", ascending=False)
        H.to_csv(f"{TAB}/a_iupac_scan_261107.tsv", sep="\t", index=False)
        thr = -np.log10(0.05 / max(cells, 1))
        fl = [dict(motif_instances_n_bg=nb,
                   min_detectable_n_fg=next((c for c in range(1, min(nb, NF_TOT) + 1)
                       if -np.log10(max(hypergeom.sf(c - 1, NB_TOT, nb, NF_TOT), 1e-300)) > thr), 0))
              for nb in (20, 50, 100, 200, 500, 1000, 2000, 5000, 20000)]
        for r in fl: r["min_detectable_occupancy_pct"] = round(100 * r["min_detectable_n_fg"] / r["motif_instances_n_bg"], 3)
        pd.DataFrame(fl).assign(space="iupac", bonf_thr=thr, n_tested_cells=cells,
                                n_patterns_enumerated=pats).to_csv(
            f"{TAB}/a_iupac_detection_floor_261107.tsv", sep="\t", index=False)
        log("iupac done: patterns", pats, "cells", cells, "thr", thr,
            "max occ over ALL tested cells", mx,
            "| passing Bonferroni", int((H.neglog10p > thr).sum()),
            "| max occ among passing",
            float(H[H.neglog10p > thr].occ_pct.max()) if (H.neglog10p > thr).any() else None)
