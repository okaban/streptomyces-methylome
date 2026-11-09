#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Minor / unassigned-motif analysis for the M145 methylome (2026-11-04).

    python3 minor_motifs_261104.py [--stage all|kmer|spatial|expr]

Author question (2026-11-04): apart from the two major recognition sequences,
is there any MINOR motif with (a) sequence enrichment, (b) spatial clustering
on the chromosome, or (c) an association with transcription?

Design decisions
----------------
site set     canonical high-confidence sites, 90_/canonical_sites.load_canonical()
             (per-replicate depth >= 10, 3 replicates, weighted freq >= 50 %).
             "unassigned" = motif column == 'other', i.e. not TGGCCGGC-C5 and
             not AAGCCCG (C5 / A1 / A2) under the canonical strand-aware matcher.
             4mC: 0 unassigned sites at every timepoint -> no 4mC minor-motif
             analysis is possible; recorded and skipped.
             6mA: 2,681 unique (position, strand) pooled over T1-T3;
             1,499 of them called at >= 2 timepoints ("reproducible").
NULL         every test uses CALLABLE positions as the denominator, never the
             whole genome. callable(TP) = adenine position whose per-replicate
             Nvalid >= 10 in all three replicates of that timepoint, from the
             canonical genome-wide pileup (pileup_canonical_261015, 'a' rows).
             Positions of the major 6mA motif (AAGCCCG A1/A2) are removed from
             both foreground and background so no test can be driven by them.
k-mer        k = 4..8, every offset whose k-mer covers the modified adenine.
             One-sided Fisher (hypergeometric) per k-mer, BH over all
             167,424 (k, offset, k-mer) cells, q < 0.05.
spatial      (i) two-sample KS of site coordinate vs callable coordinate;
             (ii) 50 kb window overdispersion against a callable-proportional
                  expectation, with a Monte-Carlo null drawn from callable;
             (iii) core (1.50-7.167 Mb) vs arms;
             (iv) overlap with antiSMASH 29 regions, the six rRNA operons and
                  the 16 IslandPath-DIMOB genomic islands.
expression   procedure of 114_pergene_corr_261102 (same inputs, same depth
             floor of 30 on PRE-FILTER depth, same exact Spearman null at
             n = 9, same 200 circular shifts, seed 20261102).

Outputs: ../tables/*.tsv
"""
import os, re, sys, json
import numpy as np
import pandas as pd
from itertools import permutations
from scipy.stats import hypergeom, ks_2samp, rankdata, chi2

RNA = "/Users/okaban/bioinfo/rna-seq"
AN = f"{RNA}/11_epigenome_integration/analysis"
sys.path.insert(0, f"{AN}/90_per_timepoint_census_audit")
import canonical_sites as cs

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.normpath(os.path.join(HERE, "..", "tables"))
CACHE = os.path.normpath(os.path.join(HERE, "..", "cache"))
REF_FA = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
GFF = f"{AN}/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff"
TSSCSV = f"{AN}/18_tss_analyses/jeong2016_all_tss.csv"
BGC = f"{RNA}/15_paper_figures/scripts/supp_261014/data/bgc_regions_29.csv"
DIMOB = f"{RNA}/16_genomic_islands_261102/islandpath/dimob_NC_003888.3.txt"
FC = f"{RNA}/03_quantification/data/featureCounts/featureCounts_M145_s2.txt"

LIBS = ["1-1", "1-2", "1-3", "2-1", "2-3", "2-4", "3-2", "3-3", "3-4"]
TPS = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"],
       "T3": ["3-2", "3-3", "3-4"]}
CHROM = "NC_003888.3"
CORE_LO, CORE_HI = cs.CORE_LO, cs.CORE_HI
MIN_REP_DEPTH = 10      # canonical callable floor, per replicate
MIN_DEPTH_EXPR = 30     # 114_ floor, per culture, pre-filter depth
MIN_TPM = 1.0
N_SHIFTS = 200
SEED = 20261102
WIN = 12                # +-12 bp strand-oriented context around the modified A
KMIN, KMAX = 4, 8

B2I = {"A": 0, "C": 1, "G": 2, "T": 3}
I2B = np.array(list("ACGT"))


def log(*a):
    print(*a, flush=True)


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


# ----------------------------------------------------------------- genome ---
def load_genome():
    g = "".join(l.strip() for l in open(REF_FA) if not l.startswith(">")).upper()
    assert len(g) == 8_667_507, len(g)
    return g


def motif_positions(genome, grc, L, motif, offsets):
    """0-based positions of `offsets` within every strand-aware occurrence.

    Returns (pos, strand) arrays; strand '+' for the forward match.
    """
    IUPAC = {b: b for b in "ACGT"}
    IUPAC.update({"K": "[GT]", "R": "[AG]", "Y": "[CT]", "S": "[CG]",
                  "W": "[AT]", "M": "[AC]", "N": "."})
    pat = "".join(IUPAC[c] for c in motif)
    pp, ss = [], []
    for m in re.finditer("(?=" + pat + ")", genome):
        for o in offsets:
            pp.append(m.start() + o)
            ss.append("+")
    for m in re.finditer("(?=" + pat + ")", grc):
        for o in offsets:
            pp.append(L - 1 - (m.start() + o))
            ss.append("-")
    return np.array(pp, dtype=np.int64), np.array(ss, dtype=object)


def context_codes(genome, grc, L, pos, strand):
    """Strand-oriented [-WIN, +WIN] context, as int8 codes (0..3, 4 = N)."""
    n = len(pos)
    out = np.full((n, 2 * WIN + 1), 4, dtype=np.int8)
    gfwd = np.frombuffer(genome.encode(), dtype=np.uint8)
    code = np.full(256, 4, dtype=np.int8)
    for b, i in B2I.items():
        code[ord(b)] = i
    comp = np.array([3, 2, 1, 0, 4], dtype=np.int8)      # A<->T, C<->G
    gcode = code[gfwd]
    plus = strand == "+"
    off = np.arange(-WIN, WIN + 1)
    idx = pos[:, None] + np.where(plus[:, None], off[None, :], -off[None, :])
    ok = (idx >= 0) & (idx < L)
    idxc = np.clip(idx, 0, L - 1)
    vals = gcode[idxc]
    vals = np.where(plus[:, None], vals, comp[vals])
    out = np.where(ok, vals, 4).astype(np.int8)
    return out


# ---------------------------------------------------------------- pileup ----
def load_pileup():
    """per-library Nvalid / Nmod / Nfail at every adenine position reported."""
    cols = ["pos0", "strand", "Nvalid", "Nmod", "Nfail"]
    frames = {}
    for lib in LIBS:
        d = pd.read_csv(f"{CACHE}/{lib}_a.tsv", sep="\t", header=None, names=cols,
                        dtype={"pos0": np.int64, "strand": "category",
                               "Nvalid": np.int32, "Nmod": np.int32,
                               "Nfail": np.int32})
        d = d.drop_duplicates(["pos0", "strand"])
        frames[lib] = d
        log(f"  {lib}: {len(d):,} adenine positions")
    base = frames[LIBS[0]][["pos0", "strand"]]
    for lib in LIBS[1:]:
        base = base.merge(frames[lib][["pos0", "strand"]], how="outer",
                          on=["pos0", "strand"])
    base = base.sort_values(["pos0", "strand"]).reset_index(drop=True)
    key = pd.MultiIndex.from_frame(base[["pos0", "strand"]])
    n = len(base)
    Nv = np.zeros((9, n), dtype=np.int32)
    Nm = np.zeros((9, n), dtype=np.int32)
    Nf = np.zeros((9, n), dtype=np.int32)
    for i, lib in enumerate(LIBS):
        r = frames[lib].set_index(["pos0", "strand"]).reindex(key)
        Nv[i] = r.Nvalid.fillna(0).to_numpy()
        Nm[i] = r.Nmod.fillna(0).to_numpy()
        Nf[i] = r.Nfail.fillna(0).to_numpy()
    return base, Nv, Nm, Nf


def main():
    os.makedirs(TAB, exist_ok=True)
    rng = np.random.default_rng(SEED)
    genome = load_genome()
    L = len(genome)
    grc = rc(genome)

    # ---------------------------------------------------------- site sets --
    hc = cs.load_canonical()
    cnt = (hc.groupby(["mod", "timepoint", "motif"]).size()
             .rename("n").reset_index())
    cnt.to_csv(f"{TAB}/site_census_261104.tsv", sep="\t", index=False)
    n4 = int(((hc["mod"] == "4mC") & (hc.motif == "other")).sum())
    log(f"unassigned 4mC sites (all timepoints): {n4}")

    oth = hc[(hc["mod"] == "6mA") & (hc.motif == "other")].copy()
    ntp = oth.groupby(["position", "strand"]).timepoint.nunique()
    uni = ntp.reset_index().rename(columns={"timepoint": "n_tp"})
    uni["strand"] = uni.strand.astype(str)
    sets = {"unassigned6mA_any": uni[uni.n_tp >= 1],
            "unassigned6mA_repro": uni[uni.n_tp >= 2]}
    log({k: len(v) for k, v in sets.items()})

    # ---------------------------------------------------------- callable ---
    log("loading canonical genome-wide adenine pileup ...")
    base, Nv, Nm, Nf = load_pileup()
    pos_all = base.pos0.to_numpy()
    str_all = base.strand.astype(str).to_numpy()
    # keep only positions that really are adenine on their own strand
    gcode = np.frombuffer(genome.encode(), dtype=np.uint8)
    is_A = np.where(str_all == "+", gcode[pos_all] == ord("A"),
                    gcode[pos_all] == ord("T"))
    log(f"  pileup rows {len(pos_all):,}; adenine on own strand {is_A.sum():,}")
    base, Nv, Nm, Nf = base[is_A], Nv[:, is_A], Nm[:, is_A], Nf[:, is_A]
    pos_all, str_all = pos_all[is_A], str_all[is_A]

    callable_tp = {}
    for tp, libs in TPS.items():
        ii = [LIBS.index(l) for l in libs]
        callable_tp[tp] = (Nv[ii] >= MIN_REP_DEPTH).all(axis=0)
    n_tp_callable = np.sum([callable_tp[t] for t in TPS], axis=0)
    callable_any = n_tp_callable >= 1
    callable_repro = n_tp_callable >= 2

    # drop the major 6mA motif positions from foreground and background
    mp, ms = motif_positions(genome, grc, L, "AAGCCCG", [0, 1])
    major = set(zip(mp.tolist(), ms.tolist()))
    is_major = np.array([(p, s) in major for p, s in zip(pos_all, str_all)])
    log(f"  AAGCCCG A1/A2 positions in pileup: {is_major.sum():,}")

    keymap = {(p, s): i for i, (p, s) in enumerate(zip(pos_all, str_all))}

    rows = []
    for name, df in sets.items():
        idx = np.array([keymap.get((int(p), s), -1)
                        for p, s in zip(df.position, df.strand)])
        miss = int((idx < 0).sum())
        idx = idx[idx >= 0]
        bgmask = callable_repro if name.endswith("repro") else callable_any
        bgmask = bgmask & ~is_major
        fg = np.zeros(len(pos_all), bool)
        fg[idx] = True
        fg &= ~is_major
        rows.append(dict(site_set=name, n_sites=len(df), n_not_in_pileup=miss,
                         n_sites_major_motif=int(np.sum(np.isin(
                             idx, np.flatnonzero(is_major)))),
                         n_foreground=int(fg.sum()),
                         n_callable_background=int(bgmask.sum()),
                         pct_fg_in_background=round(
                             100 * float((fg & bgmask).sum() / max(fg.sum(), 1)), 2)))
        sets[name] = dict(df=df, fg=fg, bg=bgmask)
    pd.DataFrame(rows).to_csv(f"{TAB}/site_and_callable_counts_261104.tsv",
                              sep="\t", index=False)
    log(pd.DataFrame(rows).to_string(index=False))

    np.save(f"{CACHE}/pos_all.npy", pos_all)
    np.save(f"{CACHE}/str_all.npy", np.where(str_all == "+", 1, 0).astype(np.int8))
    np.save(f"{CACHE}/Nv.npy", Nv)
    np.save(f"{CACHE}/Nm.npy", Nm)
    np.save(f"{CACHE}/Nf.npy", Nf)
    np.save(f"{CACHE}/callable_any.npy", callable_any)
    np.save(f"{CACHE}/callable_repro.npy", callable_repro)
    np.save(f"{CACHE}/is_major.npy", is_major)
    for name in sets:
        np.save(f"{CACHE}/fg_{name}.npy", sets[name]["fg"])
        np.save(f"{CACHE}/bg_{name}.npy", sets[name]["bg"])
    log("stage 1 (sites + callable) done")
    return dict(genome=genome, grc=grc, L=L, pos_all=pos_all, str_all=str_all,
                Nv=Nv, Nm=Nm, Nf=Nf, sets=sets, is_major=is_major,
                callable_any=callable_any, callable_repro=callable_repro)


# ------------------------------------------------------------ stage 2 ------
def stage_kmer(S):
    """k-mer enrichment at the modified adenine, callable background."""
    genome, grc, L = S["genome"], S["grc"], S["L"]
    pos_all, str_all = S["pos_all"], S["str_all"]
    log("building strand-oriented contexts ...")
    ctx = context_codes(genome, grc, L, pos_all, str_all)
    assert (ctx[:, WIN] == 0).all(), "centre base must be adenine"

    all_rows = []
    n_tests_total = sum(k * 4 ** (k - 1) for k in range(KMIN, KMAX + 1))
    for name in ("unassigned6mA_repro", "unassigned6mA_any"):
        fg = S["sets"][name]["fg"] & S["sets"][name]["bg"]
        bg = S["sets"][name]["bg"]
        Nbg, nfg = int(bg.sum()), int(fg.sum())
        log(f"{name}: foreground {nfg:,} / callable {Nbg:,}")
        fi = np.flatnonzero(fg)
        bi = np.flatnonzero(bg)
        for k in range(KMIN, KMAX + 1):
            for o in range(-(k - 1), 1):          # k-mer start, A at index -o
                cols = np.arange(WIN + o, WIN + o + k)
                Cf = ctx[np.ix_(fi, cols)].astype(np.int64)
                Cb = ctx[np.ix_(bi, cols)].astype(np.int64)
                okf = (Cf < 4).all(axis=1)
                okb = (Cb < 4).all(axis=1)
                w = 4 ** np.arange(k - 1, -1, -1)
                cf = np.bincount((Cf[okf] * w).sum(axis=1), minlength=4 ** k)
                cb = np.bincount((Cb[okb] * w).sum(axis=1), minlength=4 ** k)
                # only k-mers with A at the modified position are testable
                digit = (np.arange(4 ** k) // (4 ** (k - 1 - (-o)))) % 4
                sel = np.flatnonzero(digit == 0)
                K = cb[sel]
                x = cf[sel]
                p = np.ones(len(sel))
                nz = K > 0
                p[nz] = hypergeom.sf(x[nz] - 1, Nbg, K[nz], nfg)
                keep = np.flatnonzero((x >= 5) & (p < 1e-3))
                for j in keep:
                    code = sel[j]
                    seq = "".join(I2B[(code // 4 ** (k - 1 - t)) % 4]
                                  for t in range(k))
                    all_rows.append(dict(
                        site_set=name, k=k, offset=o, kmer=seq,
                        mod_index_in_kmer=-o + 1,
                        n_sites_with_kmer=int(x[j]),
                        n_callable_with_kmer=int(K[j]),
                        frac_sites=round(x[j] / nfg, 6),
                        frac_callable=round(K[j] / Nbg, 8),
                        oe=round((x[j] / nfg) / max(K[j] / Nbg, 1e-12), 3),
                        occupancy_pct=round(100 * x[j] / K[j], 3),
                        p=float(p[j])))
            log(f"  k={k} done")
    df = pd.DataFrame(all_rows)
    # BH over the full test family (all k, offsets, A-centred k-mers), both sets
    out = []
    for name, sub in df.groupby("site_set"):
        pv = np.concatenate([sub.p.to_numpy(),
                             np.ones(n_tests_total - len(sub))])
        q = bh(pv)[:len(sub)]
        sub = sub.assign(q_BH=q, n_tests_in_family=n_tests_total)
        out.append(sub)
    df = pd.concat(out).sort_values(["site_set", "p"])
    df.to_csv(f"{TAB}/kmer_enrichment_261104.tsv", sep="\t", index=False)
    log(f"k-mer cells reported: {len(df)}, q<0.05: {int((df.q_BH < 0.05).sum())}")

    # ---- collapse significant cells into distinct candidate motifs --------
    cand = []
    for name in ("unassigned6mA_repro", "unassigned6mA_any"):
        sub = df[(df.site_set == name) & (df.q_BH < 0.05)].sort_values("p")
        fg = S["sets"][name]["fg"] & S["sets"][name]["bg"]
        fi = np.flatnonzero(fg)
        covered = []
        for _, r in sub.iterrows():
            k, o = int(r.k), int(r.offset)
            cols = np.arange(WIN + o, WIN + o + k)
            C = ctx[np.ix_(fi, cols)].astype(np.int64)
            w = 4 ** np.arange(k - 1, -1, -1)
            code = sum(B2I[b] * 4 ** (k - 1 - t) for t, b in enumerate(r.kmer))
            mem = set(np.flatnonzero(((C * w).sum(axis=1) == code)
                                     & (C < 4).all(axis=1)).tolist())
            red = max([len(mem & c) / len(mem) for c in covered], default=0.0)
            if red < 0.8:
                covered.append(mem)
                cand.append(dict(site_set=name, kmer=r.kmer, offset=o,
                                 mod_index_in_kmer=int(r.mod_index_in_kmer),
                                 n_sites=int(r.n_sites_with_kmer),
                                 n_callable=int(r.n_callable_with_kmer),
                                 oe=r.oe, occupancy_pct=r.occupancy_pct,
                                 p=r.p, q_BH=r.q_BH, redundancy_with_kept=red))
    cd = pd.DataFrame(cand)
    cd.to_csv(f"{TAB}/candidate_minor_motifs_261104.tsv", sep="\t", index=False)
    log(f"distinct candidate motifs: "
        f"{cd.groupby('site_set').size().to_dict() if len(cd) else 0}")
    np.save(f"{CACHE}/ctx.npy", ctx)
    return df, cd, ctx


# ------------------------------------------------------------ features -----
def load_features(L):
    feats = {}
    bgc = pd.read_csv(BGC)
    feats["antiSMASH_29_regions"] = list(zip(bgc.start.astype(int),
                                             bgc.end.astype(int)))
    gi = pd.read_csv(DIMOB, sep="\t", header=None, names=["id", "start", "end"])
    feats["IslandPath_DIMOB_islands"] = list(zip(gi.start.astype(int),
                                                 gi.end.astype(int)))
    rr = []
    for line in open(GFF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 4 and f[2] == "rRNA":
            rr.append((int(f[3]), int(f[4])))
    rr.sort()
    ops, cur = [], None
    for s, e in rr:
        if cur and s - cur[1] < 10_000:
            cur = (cur[0], max(cur[1], e))
        else:
            if cur:
                ops.append(cur)
            cur = (s, e)
    ops.append(cur)
    feats["rRNA_operons"] = [(max(1, s - 1000), min(L, e + 1000)) for s, e in ops]
    return feats


def in_intervals(pos, ivs, L):
    m = np.zeros(L + 1, bool)
    for s, e in ivs:
        m[max(1, s):min(L, e) + 1] = True
    return m[np.clip(pos + 1, 0, L)]


# ------------------------------------------------------------ stage 3 ------
def stage_spatial(S, cd, ctx):
    from scipy.stats import fisher_exact
    L = S["L"]
    pos_all = S["pos_all"]
    feats = load_features(L)
    nfeat = {k: (len(v), sum(e - s + 1 for s, e in v)) for k, v in feats.items()}
    log("features: " + json.dumps(nfeat))
    WINSZ = 50_000
    nwin = int(np.ceil(L / WINSZ))
    rng = np.random.default_rng(SEED)

    testsets = []
    for name in ("unassigned6mA_repro", "unassigned6mA_any"):
        fg = S["sets"][name]["fg"] & S["sets"][name]["bg"]
        testsets.append((name, fg, S["sets"][name]["bg"]))
    for _, r in cd[cd.site_set == "unassigned6mA_repro"].iterrows():
        name = f"motif:{r.kmer}@{int(r.mod_index_in_kmer)}"
        base_fg = S["sets"]["unassigned6mA_repro"]["fg"] & S["sets"]["unassigned6mA_repro"]["bg"]
        bg = S["sets"]["unassigned6mA_repro"]["bg"]
        k, o = len(r.kmer), int(r.offset)
        cols = np.arange(WIN + o, WIN + o + k)
        code = np.array([B2I[b] for b in r.kmer], dtype=np.int8)
        has = (ctx[:, cols] == code[None, :]).all(axis=1)
        testsets.append((name, base_fg & has, bg & has))

    rows = []
    for name, fg, bg in testsets:
        sp = pos_all[fg]
        bp = pos_all[bg]
        n, N = len(sp), len(bp)
        if n < 5:
            rows.append(dict(site_set=name, test="SKIPPED_n<5", n_sites=n,
                             n_callable=N))
            continue
        ks = ks_2samp(sp, bp, method="asymp")
        rows.append(dict(site_set=name, test="KS_coordinate_vs_callable",
                         n_sites=n, n_callable=N, statistic=round(float(ks.statistic), 4),
                         p=float(ks.pvalue)))
        wi_s = np.bincount(sp // WINSZ, minlength=nwin).astype(float)
        wi_b = np.bincount(bp // WINSZ, minlength=nwin).astype(float)
        exp = n * wi_b / wi_b.sum()
        ok = exp > 0
        stat = float((((wi_s - exp) ** 2) / np.maximum(exp, 1e-9))[ok].sum())
        draws = rng.choice(N, size=(2000, n), replace=True)
        mc = np.empty(2000)
        wb = bp // WINSZ
        for t in range(2000):
            c = np.bincount(wb[draws[t]], minlength=nwin).astype(float)
            mc[t] = (((c - exp) ** 2) / np.maximum(exp, 1e-9))[ok].sum()
        rows.append(dict(site_set=name, test="window50kb_overdispersion_MC",
                         n_sites=n, n_callable=N, statistic=round(stat, 1),
                         p=float((1 + (mc >= stat).sum()) / 2001),
                         note=f"df={int(ok.sum()) - 1}; chi2 asymptotic "
                              f"p={chi2.sf(stat, max(int(ok.sum()) - 1, 1)):.3g}"))
        for lab, mask_fn in [
                ("core_vs_arm", lambda p: (p >= CORE_LO) & (p <= CORE_HI))] + [
                (f"overlap:{fk}", (lambda ivs: (lambda p: in_intervals(p, ivs, L)))(iv))
                for fk, iv in feats.items()]:
            a = int(mask_fn(sp).sum())
            b = n - a
            c = int(mask_fn(bp).sum())
            d = N - c
            orr, pv = fisher_exact([[a, b], [c, d]])
            rows.append(dict(site_set=name, test=lab, n_sites=n, n_callable=N,
                             statistic=round(float(orr), 4), p=float(pv),
                             note=f"sites_in={a}/{n} ({100*a/n:.2f}%), "
                                  f"callable_in={c}/{N} ({100*c/N:.2f}%)"))
    sp = pd.DataFrame(rows)
    if "p" in sp:
        sp["q_BH"] = np.nan
        m = sp.p.notna()
        sp.loc[m, "q_BH"] = bh(sp.loc[m, "p"].to_numpy())
    sp.to_csv(f"{TAB}/spatial_tests_261104.tsv", sep="\t", index=False)
    log(sp[["site_set", "test", "n_sites", "statistic", "p", "q_BH"]]
        .to_string(index=False))
    return sp


# ----------------------------------------------------------- stage 3b -----
def local_gc(genome, pos, half=50):
    g = np.frombuffer(genome.encode(), dtype=np.uint8)
    isgc = ((g == ord("G")) | (g == ord("C"))).astype(np.int32)
    cum = np.concatenate([[0], np.cumsum(isgc)])
    L = len(g)
    a = np.clip(pos - half, 0, L)
    b = np.clip(pos + half + 1, 0, L)
    return (cum[b] - cum[a]) / np.maximum(b - a, 1)


def stage_spatial_matched(S, cd, ctx, n_mc=2000):
    """Every spatial statistic re-tested against a GC- and depth-matched null.

    The null draws, for each (local-GC quintile x mean-depth quintile) cell,
    exactly as many callable positions as the observed set has in that cell.
    """
    genome, L, pos_all = S["genome"], S["L"], S["pos_all"]
    Nv = S["Nv"]
    feats = load_features(L)
    WINSZ = 50_000
    nwin = int(np.ceil(L / WINSZ))
    rng = np.random.default_rng(SEED)

    gcv = local_gc(genome, pos_all)
    dep = Nv.mean(axis=0)

    testsets = []
    for name in ("unassigned6mA_repro", "unassigned6mA_any"):
        fg = S["sets"][name]["fg"] & S["sets"][name]["bg"]
        testsets.append((name, fg, S["sets"][name]["bg"]))
    cdr = cd[cd.site_set == "unassigned6mA_repro"].sort_values("p").head(10)
    for _, r in cdr.iterrows():
        k, o = len(r.kmer), int(r.offset)
        cols = np.arange(WIN + o, WIN + o + k)
        code = np.array([B2I[b] for b in r.kmer], dtype=np.int8)
        has = (ctx[:, cols] == code[None, :]).all(axis=1)
        base_fg = (S["sets"]["unassigned6mA_repro"]["fg"]
                   & S["sets"]["unassigned6mA_repro"]["bg"])
        testsets.append((f"motif:{r.kmer}@{int(r.mod_index_in_kmer)}",
                         base_fg & has, S["sets"]["unassigned6mA_repro"]["bg"] & has))

    rows = []
    for name, fg, bg in testsets:
        bi = np.flatnonzero(bg)
        fi = np.flatnonzero(fg)
        n = len(fi)
        if n < 5:
            continue
        gq = np.searchsorted(np.quantile(gcv[bi], [.2, .4, .6, .8]), gcv)
        dq = np.searchsorted(np.quantile(dep[bi], [.2, .4, .6, .8]), dep)
        strat = gq * 5 + dq
        pools = {s: bi[strat[bi] == s] for s in range(25)}
        need = np.bincount(strat[fi], minlength=25)
        draws = np.empty((n_mc, n), dtype=np.int64)
        for t in range(n_mc):
            parts = [rng.choice(pools[s], size=need[s], replace=True)
                     for s in range(25) if need[s] > 0 and len(pools[s]) > 0]
            d = np.concatenate(parts)
            if len(d) < n:          # stratum with no pool: top up from bg
                d = np.concatenate([d, rng.choice(bi, size=n - len(d))])
            draws[t] = d
        sp_obs = pos_all[fi]
        wb = pos_all // WINSZ
        wi_b = np.bincount(pos_all[bi] // WINSZ, minlength=nwin).astype(float)
        exp = n * wi_b / wi_b.sum()
        ok = exp > 0

        def chi2stat(idx):
            c = np.bincount(wb[idx], minlength=nwin).astype(float)
            return (((c - exp) ** 2) / np.maximum(exp, 1e-9))[ok].sum()
        obs_chi = chi2stat(fi)
        mc_chi = np.array([chi2stat(draws[t]) for t in range(n_mc)])
        rows.append(dict(site_set=name, test="window50kb_overdispersion",
                         n_sites=n, observed=round(float(obs_chi), 1),
                         null_mean=round(float(mc_chi.mean()), 1),
                         p_matched=float((1 + (mc_chi >= obs_chi).sum()) / (n_mc + 1)),
                         note="GC+depth matched MC, one-sided"))

        masks = {"core": (pos_all >= CORE_LO) & (pos_all <= CORE_HI)}
        for fk, iv in feats.items():
            masks[f"overlap:{fk}"] = in_intervals(pos_all, iv, L)
        for lab, m in masks.items():
            obs = int(m[fi].sum())
            mc = m[draws].sum(axis=1)
            pg = (1 + (mc >= obs).sum()) / (n_mc + 1)
            pl = (1 + (mc <= obs).sum()) / (n_mc + 1)
            rows.append(dict(site_set=name, test=lab, n_sites=n, observed=obs,
                             null_mean=round(float(mc.mean()), 2),
                             p_matched=float(min(1.0, 2 * min(pg, pl))),
                             note=f"two-sided MC; unmatched expectation "
                                  f"{n * m[bi].mean():.2f}"))
        dmc = np.array([ks_2samp(pos_all[draws[t]], pos_all[bi],
                                 method="asymp").statistic for t in range(200)])
        dobs = ks_2samp(sp_obs, pos_all[bi], method="asymp").statistic
        rows.append(dict(site_set=name, test="KS_coordinate", n_sites=n,
                         observed=round(float(dobs), 4),
                         null_mean=round(float(dmc.mean()), 4),
                         p_matched=float((1 + (dmc >= dobs).sum()) / 201),
                         note="GC+depth matched MC (200 draws), one-sided"))
        log(f"  matched null done: {name}")
    sm = pd.DataFrame(rows)
    sm["q_BH"] = bh(sm.p_matched.to_numpy())
    sm.to_csv(f"{TAB}/spatial_tests_gcmatched_261104.tsv", sep="\t", index=False)
    log(sm[sm.site_set.str.startswith("unassigned")].to_string(index=False))
    return sm


# ------------------------------------------------------------ stage 4 -----
def exact_spearman_null():
    r9 = np.arange(1, 10)
    P = np.array(list(permutations(r9)), dtype=np.int8)
    rho = 1 - 6 * ((P - r9) ** 2).sum(axis=1) / (9 * 80)
    return np.sort(np.abs(rho)), np.sort(np.unique(rho))


def spearman_rows(A, B):
    ra = np.apply_along_axis(rankdata, 1, A)
    rb = np.apply_along_axis(rankdata, 1, B)
    ra = ra - ra.mean(axis=1, keepdims=True)
    rb = rb - rb.mean(axis=1, keepdims=True)
    den = np.sqrt((ra ** 2).sum(axis=1) * (rb ** 2).sum(axis=1))
    return np.where(den > 0, (ra * rb).sum(axis=1) / np.maximum(den, 1e-12),
                    np.nan)


def stage_expression(S, cd, ctx):
    """114_pergene_corr_261102 procedure applied to minor-motif site sets."""
    L, pos_all = S["L"], S["pos_all"]
    Nv, Nm, Nf = S["Nv"], S["Nm"], S["Nf"]
    abs_null, rho_grid = exact_spearman_null()

    def exact_p(r):
        return ((len(abs_null)
                 - np.searchsorted(abs_null, np.abs(r) - 1e-12)) / len(abs_null))

    fc = pd.read_csv(FC, sep="\t", comment="#")
    fc.columns = list(fc.columns[:6]) + LIBS
    assert len(fc) == 8275, len(fc)
    rpk = fc.set_index("Geneid")[LIBS].div(fc.set_index("Geneid")["Length"], axis=0)
    tpm = rpk / rpk.sum() * 1e6
    g = fc[["Geneid", "Chr", "Start", "End", "Strand"]].copy()
    g["Chr"] = g.Chr.str.split(";").str[0]
    g["Start"] = g.Start.astype(str).str.split(";").str[0].astype(int)
    g["End"] = g.End.astype(str).str.split(";").str[-1].astype(int)
    g["Strand"] = g.Strand.astype(str).str[0]
    tss = pd.read_csv(TSSCSV)
    assert len(tss) == 3570, len(tss)

    # Jeong-TSS promoter window per gene, exactly as 114_
    cg = g[g.Chr == CHROM].copy()
    cg["gstart"] = np.where(cg.Strand == "+", cg.Start, cg.End)
    W = {}
    for tid, tp, st in zip(tss.tss_id, tss.tss_position, tss.strand):
        sub = cg[cg.Strand == st]
        dd = (sub.gstart - tp) if st == "+" else (tp - sub.gstart)
        ok = (dd >= 0) & (dd <= 500)
        if not ok.any():
            continue
        gid = sub.Geneid[ok].iloc[dd[ok].argmin()]
        a, b = (tp - 300, tp - 1) if st == "+" else (tp + 1, tp + 300)
        ab = int(tss.loc[tss.tss_id == tid, "abundance"].iloc[0])
        cur = W.get(gid)
        if cur is None or ab > cur[2]:
            W[gid] = (max(a, 1), min(b, L), ab)
    W = {k: (v[0], v[1]) for k, v in W.items()}

    order_pos = np.argsort(pos_all)
    pos_sorted = pos_all[order_pos]

    def gene_matrix(sel_mask):
        idx = order_pos[sel_mask[order_pos]]
        pp = pos_all[idx]
        nv, nm, nf = Nv[:, idx], Nm[:, idx], Nf[:, idx]
        gid, lev, dep, npos = [], [], [], []
        for kk, (a, b) in W.items():
            i, j = np.searchsorted(pp, [a - 1, b])
            if j <= i:
                continue
            dd = nv[:, i:j].sum(axis=1)
            ff = nf[:, i:j].sum(axis=1)
            mm = nm[:, i:j].sum(axis=1)
            gid.append(kk)
            dep.append(dd + ff)
            lev.append(np.where(dd > 0, mm / np.maximum(dd, 1), np.nan))
            npos.append(j - i)
        return (np.array(gid), np.array(lev), np.array(dep), np.array(npos))

    # candidate site sets: the unassigned sets, plus the 10 strongest motifs,
    # each as ALL its callable occurrences (the methylatable set of that motif)
    SETS = []
    bgr = S["sets"]["unassigned6mA_repro"]["bg"]
    SETS.append(("unassigned6mA_repro_sites",
                 S["sets"]["unassigned6mA_repro"]["fg"] & bgr))
    for _, r in cd[cd.site_set == "unassigned6mA_repro"].sort_values("p").head(10).iterrows():
        k, o = len(r.kmer), int(r.offset)
        cols = np.arange(WIN + o, WIN + o + k)
        code = np.array([B2I[b] for b in r.kmer], dtype=np.int8)
        has = (ctx[:, cols] == code[None, :]).all(axis=1)
        SETS.append((f"motif:{r.kmer}@{int(r.mod_index_in_kmer)}_all_occurrences",
                     bgr & has))

    rng = np.random.default_rng(SEED)
    res, per = [], []
    for label, mask in SETS:
        gid, lev, dep, npos = gene_matrix(mask)
        if len(gid) == 0:
            res.append(dict(set=label, n_positions=int(mask.sum()),
                            genes_with_position=0, n_genes=0,
                            note="no promoter contains a position"))
            continue
        ok = (dep >= MIN_DEPTH_EXPR).all(axis=1)
        Xa, gidok = lev[ok], gid[ok]
        T = tpm.reindex(gidok)[LIBS].to_numpy()
        keep = (np.isfinite(Xa).all(axis=1) & np.isfinite(T).all(axis=1)
                & (T.mean(axis=1) >= MIN_TPM))
        n_const = int((keep & (np.ptp(Xa, axis=1) == 0)).sum())
        keep &= np.ptp(Xa, axis=1) > 0      # Spearman undefined without variance
        gg, X, T = gidok[keep], Xa[keep], T[keep]
        n = len(gg)
        if n < 5:
            res.append(dict(set=label, n_positions=int(mask.sum()),
                            genes_with_position=int(len(gid)),
                            genes_depth_ok=int(ok.sum()), n_genes=n,
                            genes_constant_methylation=n_const,
                            note="n < 5: no test, no power"))
            continue
        order = np.argsort(g.set_index("Geneid").loc[gg, "Start"].to_numpy())
        gg, X, T = gg[order], X[order], T[order]
        rho = spearman_rows(X, T)
        p = np.array([exact_p(r_) for r_ in rho])
        q = bh(p)
        shifts = rng.choice(np.arange(1, n), size=min(N_SHIFTS, n - 1),
                            replace=False)
        nr = []
        for s_ in shifts:
            nr.append(spearman_rows(X, np.roll(T, s_, axis=0)))
        nr = np.concatenate(nr)
        ks = ks_2samp(rho, nr)
        thr = 0.05 / n
        crit = min([r_ for r_ in rho_grid[rho_grid > 0] if exact_p(r_) <= thr])
        res.append(dict(set=label, n_positions=int(mask.sum()),
                        genes_with_position=int(len(gid)),
                        genes_depth_ok=int(ok.sum()), n_genes=n,
                        genes_constant_methylation=n_const,
                        obs_mean_rho=round(float(rho.mean()), 4),
                        obs_sd_rho=round(float(rho.std()), 4),
                        null_mean_rho=round(float(nr.mean()), 4),
                        ks_D=round(float(ks.statistic), 4),
                        ks_p=float(ks.pvalue),
                        n_fdr5_observed=int((q < 0.05).sum()),
                        crit_abs_rho_fdr5=round(float(crit), 4),
                        min_attainable_two_sided_p=5.51e-06,
                        n_shifts=len(shifts)))
        for a_, b_, c_, d_ in zip(gg, rho, p, q):
            per.append(dict(set=label, gene=a_, rho=round(float(b_), 4),
                            p=float(c_), q_BH=float(d_)))
        log(f"  expression: {label} n={n}")
    pd.DataFrame(res).to_csv(f"{TAB}/expression_association_261104.tsv",
                             sep="\t", index=False)
    pd.DataFrame(per).to_csv(f"{TAB}/expression_pergene_261104.tsv",
                             sep="\t", index=False)
    log(pd.DataFrame(res).to_string(index=False))
    return pd.DataFrame(res)


def stage_match_decomposition(S, n_mc=2000):
    """Is the matched core excess driven by local GC or by depth?"""
    genome, L, pos_all, Nv = S["genome"], S["L"], S["pos_all"], S["Nv"]
    gcv = local_gc(genome, pos_all)
    dep = Nv.mean(axis=0)
    rng = np.random.default_rng(SEED)
    core = (pos_all >= CORE_LO) & (pos_all <= CORE_HI)
    rows = []
    for name in ("unassigned6mA_repro", "unassigned6mA_any"):
        fg = S["sets"][name]["fg"] & S["sets"][name]["bg"]
        bg = S["sets"][name]["bg"]
        bi, fi = np.flatnonzero(bg), np.flatnonzero(fg)
        n = len(fi)
        obs = int(core[fi].sum())
        for lab, var in (("none", None), ("local_GC_only", gcv),
                         ("mean_depth_only", dep), ("GC_x_depth", None)):
            if lab == "none":
                exp = n * core[bi].mean()
                mc = rng.binomial(n, core[bi].mean(), size=n_mc)
            else:
                if lab == "GC_x_depth":
                    q1 = np.searchsorted(np.quantile(gcv[bi], [.2, .4, .6, .8]), gcv)
                    q2 = np.searchsorted(np.quantile(dep[bi], [.2, .4, .6, .8]), dep)
                    strat = q1 * 5 + q2
                    ns = 25
                else:
                    strat = np.searchsorted(
                        np.quantile(var[bi], [.1, .2, .3, .4, .5, .6, .7, .8, .9]), var)
                    ns = 10
                need = np.bincount(strat[fi], minlength=ns)
                pools = {s: bi[strat[bi] == s] for s in range(ns)}
                exp = sum(need[s] * core[pools[s]].mean()
                          for s in range(ns) if need[s] and len(pools[s]))
                mc = np.array([
                    sum(int(core[rng.choice(pools[s], size=need[s])].sum())
                        for s in range(ns) if need[s] and len(pools[s]))
                    for _ in range(n_mc)])
            pg = (1 + (mc >= obs).sum()) / (n_mc + 1)
            pl = (1 + (mc <= obs).sum()) / (n_mc + 1)
            rows.append(dict(site_set=name, matched_on=lab, n_sites=n,
                             observed_core=obs,
                             expected_core=round(float(exp), 1),
                             obs_pct=round(100 * obs / n, 2),
                             exp_pct=round(100 * float(exp) / n, 2),
                             p=float(min(1.0, 2 * min(pg, pl)))))
            log(f"  decomposition {name} / {lab}: {obs} vs {exp:.1f}")
    d = pd.DataFrame(rows)
    d.to_csv(f"{TAB}/core_arm_matching_decomposition_261104.tsv", sep="\t",
             index=False)
    log(d.to_string(index=False))
    return d


if __name__ == "__main__":
    S = main()
    df, cd, ctx = stage_kmer(S)
    sp = stage_spatial(S, cd, ctx)
    sm = stage_spatial_matched(S, cd, ctx)
    dc = stage_match_decomposition(S)
    ex = stage_expression(S, cd, ctx)
