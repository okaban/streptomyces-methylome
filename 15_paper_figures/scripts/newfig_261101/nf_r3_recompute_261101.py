#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NF-R3 (local methylation) — recompute, author definition of 2026-11-01.

    python3 nf_r3_recompute_261101.py

Author decision (2026-11-01) fixed, and this script implements exactly:
    promoter = -300 .. TSS   (upstream of the TSS, TSS itself excluded)
    5' UTR   = TSS .. +50    (downstream of the TSS, TSS itself included)
    CDS      = inside an annotated CDS and not promoter/5'UTR
    other    = none of the above
    TSS population = the 3,570 experimental dRNA-seq TSSs of Jeong et al. 2016.

Why this is not just a re-run of rebuild_suppfig18a_261016.py
-------------------------------------------------------------
Three differences, all of which change numbers:

1. TSS population. 18a selects category in {P,S} AND gene_id non-null = 3,030
   TSSs. The author fixed the population at the full measured set, 3,570.
2. TSS strand. 18a ignores the `strand` column of jeong2016_all_tss.csv and
   treats "upstream" as "lower coordinate" for every TSS. 1,732 of 3,570 TSSs
   are on the minus strand, where upstream is the HIGHER coordinate. This
   script takes the signed distance along the TSS's own strand.
3. Sign. 18a's own comment reads `negative = upstream of the TSS`, but its
   `nearest` is -(p - tss_prev), i.e. negative where the position lies
   DOWNSTREAM of the nearest TSS, so its `prom` class is TSS..+300 and its
   `utr5` class is -50..TSS - the two windows are swapped relative to the
   declaration. Both defects are reported in c04_canon_corrections_261101.md.

Coordinate convention: `corrected` (both coordinates 1-based), the convention
that settled PROM.6mA.RATIO = 0.9105. BED start is 0-based, so the query is
start + 1.

Census rule (unchanged, canonical): per timepoint pool the 3 libraries of
pileup_canonical_261015 (via the npz cache), require pooled true depth
Nvalid_cov + Nfail >= 30, call a site at pooled Nmod/Nvalid >= 0.5.
01_integration/high_confidence_sites_weighted.csv is NOT used anywhere.

Part B (TSS-proximal level ratio) is NOT a census statistic: it is level-scale,
all-read (`pileup_nofilt_261015`, Nfail == 0), local-GC quintile matched,
culture unit n = 9, and is computed by importing the functions that produced
SuppTable22b unchanged. The 6mA focal rung is asserted against SuppTable22b.

Outputs (tables/):
    nf_r3_region_261101.tsv     region composition, 3 timepoints x 2 marks
    nf_r3_tssratio_261101.tsv   near/distal level ratio, both marks, 2 windows
"""
import importlib.util
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(HERE, "tables")
CACHE = ("/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/"
         "recompute_261015/pileup_cache_261015")
TSSCSV = ("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
          "18_tss_analyses/jeong2016_all_tss.csv")
GFF = ("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
       "39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff")
MOD = ("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
       "48_TSS_methylation_gradient/scripts/recompute_TSS_set_controls_261021.py")
SUPP22B = ("/Users/okaban/bioinfo/rna-seq/15_paper_figures/tables/supplementary/"
           "SuppTable22b_TSS_set_controls_261021.tsv")

_spec = importlib.util.spec_from_file_location("tssmod", MOD)
M = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M)

TIMEPOINTS = {"T1": ("1-1", "1-2", "1-3"), "T2": ("2-1", "2-3", "2-4"),
              "T3": ("3-2", "3-3", "3-4")}
CODE = {"6mA": 0, "4mC": 1}          # npz `code`: 0 = 'a', 1 = '21839', 2 = 'm'
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "K": "[GT]", "H": "[ACT]"}

# declared modified positions (canon; 0-based offsets into the motif)
SERIES = [
    ("TGGCCGGC", "4mC", [4], "TGGCCGGC C5"),
    ("AAGCCCG",  "4mC", [4], "AAGCCCG C5"),
    ("AAGCCCG",  "6mA", [0, 1], "AAGCCCG A1/A2"),
]
CLASSES = ["promoter", "utr5", "cds", "other"]


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def motif_positions(genome, motif, offsets):
    """Strand-oriented 0-based positions of the modified bases."""
    L = len(genome)
    pat = "".join(IUPAC[c] for c in motif)
    out = set()
    for m in re.finditer("(?=" + pat + ")", genome):
        for o in offsets:
            out.add((m.start() + o, True))
    for m in re.finditer("(?=" + pat + ")", rc(genome)):
        for o in offsets:
            out.add((L - 1 - (m.start() + o), False))
    p = np.array(sorted(out), dtype=object)
    return np.array([x[0] for x in p], int), np.array([x[1] for x in p], bool)


def pooled(libs, code):
    """Pooled Nvalid/Nmod/true-depth per (start, strand) for one timepoint."""
    frames = []
    for lib in libs:
        z = np.load(os.path.join(CACHE, f"{lib}.npz"))
        m = (z["code"] == code) & z["refok"]
        frames.append(pd.DataFrame(
            {"start": z["start"][m], "plus": z["plus"][m],
             "Nvalid": z["Nvalid"][m], "Nmod": z["Nmod"][m],
             "Nfail": z["Nfail"][m]}).set_index(["start", "plus"]))
    uni = frames[0].index
    for f in frames[1:]:
        uni = uni.union(f.index)
    s = sum(f.reindex(uni).fillna(0) for f in frames)
    s["dep"] = s.Nvalid + s.Nfail
    s["lvl"] = s.Nmod / s.Nvalid.replace(0, np.nan)
    return s


def signed_tss_distance(pos1, tss_pos, tss_strand):
    """Signed distance to the nearest TSS, measured along that TSS's strand.

    Negative = upstream of the TSS (the promoter side), positive = downstream.
    `pos1` is 1-based, as is `tss_pos` (coord_convention = corrected).
    """
    out = np.full(len(pos1), np.iinfo(np.int64).max // 4, dtype=np.int64)
    for strand, sign in (("+", 1), ("-", -1)):
        t = np.sort(tss_pos[tss_strand == strand].astype(np.int64))
        i = np.searchsorted(t, pos1)
        cand = []
        for j in (i - 1, i):
            k = np.clip(j, 0, len(t) - 1)
            cand.append(sign * (pos1 - t[k]))
        cand = np.vstack(cand)
        best = cand[np.argmin(np.abs(cand), axis=0), np.arange(len(pos1))]
        take = np.abs(best) < np.abs(out)
        out[take] = best[take]
    return out


def classify(pos0, d, cds_mask):
    """promoter / utr5 / cds / other, in that precedence."""
    prom = (d < 0) & (d >= -300)
    utr5 = (d >= 0) & (d <= 50)
    cds = cds_mask[pos0] & ~prom & ~utr5
    other = ~prom & ~utr5 & ~cds
    return {"promoter": prom, "utr5": utr5, "cds": cds, "other": other}


def part_a(genome, tss, cds_mask):
    rows = []
    store = {}
    for tp, libs in TIMEPOINTS.items():
        for code_name, code in CODE.items():
            store[(tp, code_name)] = pooled(libs, code)
    for tp in TIMEPOINTS:
        for motif, mark, offs, label in SERIES:
            pos0, plus = motif_positions(genome, motif, offs)
            s = store[(tp, mark)]
            idx = pd.MultiIndex.from_arrays([pos0, plus], names=["start", "plus"])
            d = s.reindex(idx)
            keep = (d["dep"] >= 30).to_numpy(dtype=bool)
            p = pos0[keep]
            called = (d["lvl"].to_numpy()[keep] >= 0.5)
            dist = signed_tss_distance(p + 1, tss["tss_position"].to_numpy(),
                                       tss["strand"].to_numpy())
            cls = classify(p, dist, cds_mask)
            for c in CLASSES:
                rows.append(dict(
                    timepoint=tp, mark=mark, series=label, region=c,
                    n_modified=int(cls[c][called].sum()),
                    n_measurable=int(cls[c].sum()),
                    n_modified_total=int(called.sum()),
                    n_measurable_total=int(len(p)),
                    pct_of_modified=round(100 * cls[c][called].sum()
                                          / max(called.sum(), 1), 2),
                    pct_of_measurable=round(100 * cls[c].sum()
                                            / max(len(p), 1), 2)))
    return pd.DataFrame(rows)


def part_b(genome):
    """Near/distal level ratio, all reads, GC-quintile matched, both marks."""
    cGC, cN, L = M.cum_gc(genome)
    j = pd.read_csv(TSSCSV)
    tss_all = np.array(sorted(j["tss_position"].astype(int)))
    assert len(np.unique(tss_all)) == 3570, len(np.unique(tss_all))

    tgg = set(int(m.start() + 4) for m in re.finditer("(?=TGGCCGGC)", genome))
    tgg |= set(L - 1 - (m.start() + 4)
               for m in re.finditer("(?=TGGCCGGC)", rc(genome)))

    universes = {
        "6mA, AAGCCCG A1/A2": (["a,AAGCCCG,0", "a,AAGCCCG,1"], None),
        "4mC, TGGCCGGC C5":   (["21839,GCCGGC,2"], tgg),
        "4mC, AAGCCCG C5":    (["21839,AAGCCCG,4"], None),
        # not a recognition sequence - carried only to identify which universe
        # the ledger row AUDIT.G04.4mC.MATCHED = 0.9472 was computed on
        "4mC, all GCCGGC C3": (["21839,GCCGGC,2"], None),
    }
    rows = []
    for label, (names, subset) in universes.items():
        pos0, strand0, Nv, Nm = M.load_universe(names)
        if subset is not None:
            k = np.array([p in subset for p in pos0])
            pos0, strand0, Nv, Nm = pos0[k], strand0[k], Nv[:, k], Nm[:, k]
        level = np.where(Nv > 0, Nm / np.maximum(Nv, 1), np.nan)
        for win in (100, 300):
            gc = M.local_gc(pos0, cGC, cN, L, win)
            near = M.min_distance(pos0 + 1, tss_all) <= win
            ratio, lo, hi, p, n = M.summarise(
                M.standardised_ratio(level, Nv, Nm, near, gc, nq=5))
            rows.append(dict(
                universe=label, universe_n=int(len(pos0)), window_bp=win,
                tss_set="jeong_measured", tss_n=3570, read_set="all reads",
                gc_matching="GC-matched", gc_quantiles=5,
                coord_convention="corrected",
                n_near=int(near.sum()), n_distal=int((~near).sum()),
                mean_level_near=round(float(np.nanmean(level[:, near])), 6),
                mean_level_distal=round(float(np.nanmean(level[:, ~near])), 6),
                ratio=round(ratio, 4), ci95_lo=round(lo, 4),
                ci95_hi=round(hi, 4), p_value=f"{p:.3g}", n_cultures=n))
    out = pd.DataFrame(rows)

    # the 6mA rungs must reproduce SuppTable22b, which the registry cites
    s = pd.read_csv(SUPP22B, sep="\t", dtype=str)
    s = s[s.tss_set != "tss_set"]
    for win, want_id in ((100, "PROM.6mA.RATIO"), (300, None)):
        w = s[(s.tss_set == "jeong_measured") & (s.universe == "both_A")
              & (s.window_bp == str(win)) & (s.coord_convention == "corrected")]
        assert len(w) == 1, (win, len(w))
        g = out[(out.universe == "6mA, AAGCCCG A1/A2")
                & (out.window_bp == win)].iloc[0]
        for col in ("ratio", "ci95_lo", "ci95_hi", "n_near", "n_distal"):
            assert float(g[col]) == float(w.iloc[0][col]), (
                win, col, g[col], w.iloc[0][col])
    return out


def main():
    os.makedirs(TAB, exist_ok=True)
    genome = M.read_reference(M.REF)
    tss = pd.read_csv(TSSCSV)
    assert len(tss) == 3570
    cds_mask = np.zeros(len(genome), bool)
    for line in open(GFF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[2] != "CDS":
            continue
        cds_mask[int(f[3]) - 1:int(f[4])] = True

    a = part_a(genome, tss, cds_mask)
    a.to_csv(os.path.join(TAB, "nf_r3_region_261101.tsv"), sep="\t", index=False)
    b = part_b(genome)
    b.to_csv(os.path.join(TAB, "nf_r3_tssratio_261101.tsv"), sep="\t", index=False)

    piv = (a[a.timepoint == "T1"]
           .pivot(index="series", columns="region", values="pct_of_modified"))
    print("== region composition, T1, % of modified sites ==")
    print(piv[CLASSES].to_string())
    print("\n== measurable-position background, T1, % ==")
    print(a[a.timepoint == "T1"].pivot(index="series", columns="region",
                                       values="pct_of_measurable")[CLASSES]
          .to_string())
    print("\n== totals per timepoint ==")
    print(a.groupby(["timepoint", "series"])[["n_modified_total",
                                              "n_measurable_total"]]
          .first().to_string())
    print("\n== near/distal level ratio ==")
    print(b[["universe", "universe_n", "window_bp", "n_near", "ratio",
             "ci95_lo", "ci95_hi", "p_value"]].to_string(index=False))


if __name__ == "__main__":
    main()
