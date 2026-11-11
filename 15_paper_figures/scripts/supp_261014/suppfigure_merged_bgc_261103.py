#!/usr/bin/env python3
"""Merged supplementary figure — biosynthetic gene clusters against the
methylation landscape, and why no cluster-level relationship is detectable.

Consolidation of two shipped supplementary figures (author approval 2026-11-03):

    old Supplementary Figure 6  "Predicted biosynthetic gene clusters are not
                                 aligned with the methylation landscape of
                                 either mark"
        generator: supp_261014/suppfigure13_bgc_by_mark.py
    old Supplementary Figure 7  "Why no cluster-level relationship between
                                 methylation and expression is detectable"
        generator: supp_261014/suppfigure14_bgc_relationships.py

The panels are re-laid on ONE canvas; neither original generator is called.
Every panel's computation is the code of its original generator, copied
unchanged, so no plotted number moves.  Panel letters are renumbered:

    old S6 a,b,c    ->  a,b,c   (chromosome-scale block, top)
    old S7 a,b,c,d  ->  d,e,f,g (cluster-level block, bottom)

Only the axis-label wording changed (units added in parentheses, FIG rule).

    SUPPFIG_OUT=<dir> PYTHONPATH=. python3 supp_261014/suppfigure_merged_bgc_261103.py
"""
import importlib
import json
import textwrap
import pathlib
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from _gccggc import DEPTH_MIN, FRAC_MIN, LIBS, genome          # noqa: E402
from _suppstyle import (BLUE, DATA_DIR, META_GREY, NAR_WIDTH_IN,  # noqa: E402
                        RED, apply_figure_style, figure_panel_letters, report,
                        save, source_note)

U = importlib.import_module("00_shared_utils")
EXPR, LOW, SIL = U.COL_BAR_DARK, U.COL_BAR_MID, U.COL_BAR_LIGHT
C_MARK = U.COL_4mC

FEATURECOUNTS = "/Users/okaban/bioinfo/rna-seq/03_quantification/data/featureCounts/featureCounts_M145_s2.txt"
TSS_CSV = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv"
CANON = pathlib.Path("/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/"
                     "recompute_261015/suppfig8_panelbcd_canonical_261022.tsv")
MIN_POS = 5
TPS = ("T1", "T2", "T3")

N = json.load(open(DATA_DIR / "numbers.json"))
bgc = pd.read_csv(DATA_DIR / "bgc_regions_29.csv")

BIN = 50_000
CORE = (1.50e6, 7.17e6)
ORIC = 4.270310e6
SUB = {"tgg": "pileup_GCCGGC", "c4": "pileup_AAGCCCG_4mC", "a01": "pileup_AAGCCCG_6mA"}


# ───────────── chromosome-scale block (old Supplementary Figure 6) ───────────
def motif_keys(g):
    return {
        "tgg": {(m.start() + 4, "+") for m in re.finditer("(?=TGGCCGGC)", g)}
               | {(m.start() + 3, "-") for m in re.finditer("(?=GCCGGCCA)", g)},
        "c4": {(m.start() + 4, "+") for m in re.finditer("(?=AAGCCCG)", g)}
              | {(m.start() + 2, "-") for m in re.finditer("(?=CGGGCTT)", g)},
        "a01": {(m.start() + i, "+") for m in re.finditer("(?=AAGCCCG)", g) for i in (0, 1)}
               | {(m.start() + 5 + i, "-") for m in re.finditer("(?=CGGGCTT)", g) for i in (0, 1)},
    }


def called(setname, tp, keys):
    acc = {}
    for lib in LIBS[tp]:
        d = pd.read_csv(DATA_DIR / SUB[setname] / f"{lib}.tsv", sep="\t", header=None,
                        names=["pos", "strand", "Nvalid", "Nmod", "Nfail"])
        for pos, st, nv, nm, nf in d.itertuples(index=False):
            k = (pos, st)
            if k not in keys[setname]:
                continue
            a = acc.setdefault(k, [0, 0, 0])
            a[0] += nv; a[1] += nm; a[2] += nf
    return [k[0] for k, v in acc.items()
            if v[0] + v[2] >= DEPTH_MIN and v[0] and v[1] / v[0] >= FRAC_MIN]


def collect():
    g = genome()
    keys = motif_keys(g)
    nb = int(np.ceil(len(g) / BIN))
    series = {}
    for tp in TPS:
        tgg = called("tgg", tp, keys)
        c4 = called("c4", tp, keys)
        a01 = called("a01", tp, keys)
        assert len(tgg) == N["TGG_called"][tp], f"{tp} TGGCCGGC {len(tgg)}"
        assert len(c4) == N["AAG_C4_called"][tp], f"{tp} AAGCCCG 4mC {len(c4)}"
        assert len(a01) == N["AAG_6mA_union"][tp], f"{tp} AAGCCCG 6mA {len(a01)}"
        series[tp] = {"m4": np.bincount(np.asarray(tgg + c4) // BIN, minlength=nb),
                      "m6": np.bincount(np.asarray(a01) // BIN, minlength=nb),
                      "n4": len(tgg) + len(c4), "n6": len(a01)}
    return len(g), nb, series


SHORT = {
    1: "leinamycin", 2: "isorenieratene", 3: "lobophorin CR4", 4: "coelichelin",
    5: "region 5 (unnamed)", 6: "informatipeptin", 7: "flaviolin", 8: "ectoine",
    9: "istamycin", 10: "desferrioxamine", 11: "CDA", 12: "actinorhodin",
    13: "albaflavenone", 14: "prejadomycin", 15: "spore pigment", 16: "kinamycin",
    17: "undecylprodigiosin", 18: "14-hydroxyisochainin", 19: "geosmin",
    20: "enduracidin", 21: "coelimycin P1", 22: "nogalamycin", 23: "SapB",
    24: "hopene", 25: "Bisenarsan", 26: "region 26 (unnamed)", 27: "germicidin",
    28: "dimethylallylindole", 29: "coelibactin",
}


def spread(xs, sep, lo, hi):
    x = np.array(xs, dtype=float)
    for i in range(1, len(x)):
        x[i] = max(x[i], x[i - 1] + sep)
    if x[-1] > hi:
        x[-1] = hi
        for i in range(len(x) - 2, -1, -1):
            x[i] = min(x[i], x[i + 1] - sep)
    x[0] = max(x[0], lo)
    for i in range(1, len(x)):
        x[i] = max(x[i], x[i - 1] + sep)
    return x


def shade(ax, glen):
    for x0, x1, col in [(0, CORE[0], "#f2f2f2"), (CORE[0], CORE[1], "#e6e6e6"),
                        (CORE[1], glen, "#f2f2f2")]:
        ax.add_patch(Rectangle((x0 / 1e6, 0), (x1 - x0) / 1e6, 1,
                               transform=ax.get_xaxis_transform(),
                               fc=col, ec="none", zorder=0))
    ax.axvline(ORIC / 1e6, color="#6b6b6b", lw=0.8, ls=(0, (4, 2)), zorder=1)


def track(ax, x, series, key, glen, ylab, colour):
    shade(ax, glen)
    stack = np.vstack([series[t][key] for t in TPS])
    lo, hi, mean = stack.min(0), stack.max(0), stack.mean(0)
    ax.fill_between(x, lo, hi, step="mid", color=colour, alpha=0.30, lw=0,
                    zorder=2, label="T1-T3 range")
    ax.step(x, mean, where="mid", lw=0.85, color=colour, zorder=3,
            label="mean of T1, T2, T3")
    ax.set_xlim(0, glen / 1e6)
    ax.set_ylim(0, hi.max() * 1.30)
    ax.set_ylabel(ylab)
    ax.legend(frameon=False, loc="upper left", ncol=2, fontsize=6.0,
              handlelength=1.3, columnspacing=1.2, borderpad=0.0,
              handletextpad=0.5, labelspacing=0.2)


# ───────────── cluster-level block (old Supplementary Figure 7) ──────────────
G = genome()
tgg = pd.DataFrame(
    sorted([(m.start() + 4, "+") for m in re.finditer("(?=TGGCCGGC)", G)]
           + [(m.start() + 3, "-") for m in re.finditer("(?=GCCGGCCA)", G)]),
    columns=["pos", "strand"])

fc = pd.read_csv(FEATURECOUNTS, sep="\t", comment="#")
bam_cols = [c for c in fc.columns if c.endswith(".bam")]
LIB_RE = re.compile(r"M145_(\d)_(\d)\.")
smap = {c: LIB_RE.search(c).group(1) + "-" + LIB_RE.search(c).group(2) for c in bam_cols}
fc = fc.rename(columns=smap)
GENES = fc[["Geneid", "Start", "End", "Strand"]].copy()
GENES["Start"] = GENES.Start.astype(str).str.split(";").str[0].astype(int)
GENES["End"] = GENES.End.astype(str).str.split(";").str[-1].astype(int)

mb = pd.read_csv(DATA_DIR / "MB02_BGC_crosswalk_supplementary.csv")
CAN = pd.read_csv(CANON, sep="\t")
assert set(CAN.panel) == {"b", "c", "d"}, sorted(set(CAN.panel))

tss = pd.read_csv(TSS_CSV)
ps = tss[tss.category.isin(["P", "S"]) & tss.gene_id.notna()]
prom = pd.DataFrame(
    [(int(t.tss_position) - 300, int(t.tss_position) - 1, t.gene_id)
     if t.strand in ("+", 1, "1") else
     (int(t.tss_position) + 1, int(t.tss_position) + 300, t.gene_id)
     for t in ps.itertuples()],
    columns=["s", "e", "gene"]).drop_duplicates()

bgc_genes = set()
for t in mb.itertuples():
    bgc_genes |= set(GENES[(GENES.Start >= t.start) & (GENES.End <= t.end)].Geneid)

tggpos = np.sort(tgg.pos.values)
prom = prom.assign(
    is_bgc=prom.gene.isin(bgc_genes),
    has_site=[int(np.searchsorted(tggpos, e, "right") - np.searchsorted(tggpos, s, "left")) > 0
              for s, e in zip(prom.s, prom.e)])
summ = prom.groupby("is_bgc").agg(n_promoters=("gene", "size"), with_site=("has_site", "sum"))
summ["pct"] = 100 * summ.with_site / summ.n_promoters
assert (int(summ.loc[True, "with_site"]), int(summ.loc[True, "n_promoters"])) == (18, 310), summ
assert (int(summ.loc[False, "with_site"]), int(summ.loc[False, "n_promoters"])) == (127, 2720), summ

panel_b = (CAN[CAN.panel == "b"]
           .pivot(index="timepoint", columns="locus",
                  values=["level_pooled_pct", "n_positions", "Nmod", "Nvalid"]))
PROM_KEY, REST_KEY = "promoter (-300..-1)", "rest of cluster"
for tp in TPS:
    for k in (PROM_KEY, REST_KEY):
        num = panel_b[("Nmod", k)][tp]
        den = panel_b[("Nvalid", k)][tp]
        assert abs(100 * num / den - panel_b[("level_pooled_pct", k)][tp]) < 0.02, (tp, k)

L4 = (CAN[CAN.panel == "c"]
      .merge(mb[["region", "expression_class", "coreTPM_T1", "coreTPM_T2", "coreTPM_T3"]],
             on="region")
      .rename(columns={"expression_class": "cls"}))
L4["tpm"] = [getattr(r, f"coreTPM_{r.timepoint}") for r in L4.itertuples()]
L4 = L4.dropna(subset=["level_pooled_pct"])
L4 = L4[L4.n_positions >= MIN_POS]
r4, p4 = spearmanr(np.log2(L4.tpm + 1), L4.level_pooled_pct)

DCLASSES = ["TGGCCGGC 4mC", "AAGCCCG C5 4mC", "AAGCCCG 6mA"]
dpan = CAN[(CAN.panel == "d") & CAN.mark_class.isin(DCLASSES)]
dpan = dpan.pivot(index=["mark_class", "timepoint"], columns="locus",
                  values=["level_pooled", "n_positions", "Nmod", "Nvalid"])
dpan = dpan.assign(ratio=dpan[("level_pooled", "interior")]
                   / dpan[("level_pooled", "flank")]).reset_index()
dpan.columns = ["_".join(c).rstrip("_") if isinstance(c, tuple) else c for c in dpan.columns]
dpan = dpan.sort_values(["mark_class", "timepoint"],
                        key=lambda s: s.map({m: i for i, m in enumerate(DCLASSES)})
                        if s.name == "mark_class" else s).reset_index(drop=True)

apply_figure_style(frame="boxed", sizes=(7.0, 6.2, 6.0))
H = 8.82  # 2026-11-09: 図下の 7 行の散文を凡例へ移し、0.78 in を詰めた


def build():
    glen, nb, series = collect()
    x = (np.arange(nb) * BIN + BIN / 2) / 1e6

    fig = plt.figure(figsize=(NAR_WIDTH_IN, H))
    gs_top = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 2.75], hspace=0.26,
                              left=0.118, right=0.985, top=0.9782, bottom=0.3872)
    gs_bot = fig.add_gridspec(2, 2, hspace=0.80, wspace=0.62,
                              left=0.105, right=0.985, top=0.3262, bottom=0.0564)

    # (a) 4mC track, (b) 6mA track
    ax_a = fig.add_subplot(gs_top[0])
    track(ax_a, x, series, "m4", glen, "4mC sites\n(per 50 kb)", RED)
    ax_a.set_xticklabels([])
    ax_b = fig.add_subplot(gs_top[1])
    track(ax_b, x, series, "m6", glen, "6mA sites\n(per 50 kb)", BLUE)
    ax_b.set_xticklabels([])

    # (c) cluster positions and expression class
    ax_c = fig.add_subplot(gs_top[2])
    shade(ax_c, glen)
    rows = [("expressed", EXPR, 2), ("low", LOW, 1), ("silent", SIL, 0)]
    colour_of = {st: c for st, c, _ in rows}
    y_of = {st: y for st, _, y in rows}
    for st, col, y in rows:
        for r in bgc[bgc.status == st].itertuples():
            ax_c.add_patch(Rectangle((r.start / 1e6, y - 0.22),
                                     max((r.end - r.start) / 1e6, 0.030), 0.44,
                                     fc=col, ec="none", zorder=4))
    ax_c.set_yticks([y for _, _, y in rows])
    ax_c.set_yticklabels([f"{st} ({len(bgc[bgc.status == st])})" for st, _, _ in rows])
    for lab, (_, col, _) in zip(ax_c.get_yticklabels(), rows):
        lab.set_color(col)
    ax_c.tick_params(axis="y", length=0)

    # 2026-11-07 補足図点検 4: 29 クラスタすべてに縦書きのラベルを 2 段で付けていたため
    # 名前が密集し、low / silent の淡色ラベルはほぼ読めず、引き出し線が互いに交差し、
    # ラベル上端がパネル (b) の領域に食い込んでいた。名前を出すのは expressed の
    # 18 件だけに絞り、1 段に収める。low (8) と silent (3) は行ラベルの件数と
    # Supplementary Table で参照でき、位置は矩形がそのまま示している。
    d = (bgc[bgc.status == "expressed"].sort_values("start").reset_index(drop=True))
    mid = ((d.start + d.end) / 2 / 1e6).to_numpy()
    lx = spread(mid, 0.300, 0.10, glen / 1e6 - 0.10)
    BASE = 2.62
    for xl, xm, rn, st in zip(lx, mid, d.region_number, d.status):
        col = colour_of[st]
        y0 = y_of[st] + 0.24
        ax_c.plot([xm, xm, xl], [y0, BASE - 0.26, BASE], lw=0.4, color=col,
                  alpha=0.8, zorder=3, solid_joinstyle="miter")
        ax_c.text(xl, BASE + 0.06, SHORT[rn], rotation=90, ha="center",
                  va="bottom", fontsize=6.0, color=col)
    ax_c.set_xlim(0, glen / 1e6)
    ax_c.set_ylim(-0.55, 8.6)
    for sp in ("left", "right", "top"):
        ax_c.spines[sp].set_visible(False)
    ax_c.set_xlabel("Chromosome position (Mb)")
    for xt, lab, col in [(0.75, "left arm", "#777777"), (3.25, "core", "#3a3a3a"),
                         (7.92, "right arm", "#777777")]:
        ax_c.text(xt, -0.44, lab, ha="center", va="bottom", fontsize=6, color=col)

    # (d) promoters carrying a TGGCCGGC position
    ax_d = fig.add_subplot(gs_bot[0, 0])
    cats = ["genes outside\nclusters", "genes inside\nclusters"]
    vals = [summ.loc[False, "pct"], summ.loc[True, "pct"]]
    cnts = [(int(summ.loc[False, "with_site"]), int(summ.loc[False, "n_promoters"])),
            (int(summ.loc[True, "with_site"]), int(summ.loc[True, "n_promoters"]))]
    ax_d.bar(cats, vals, width=0.5, color=C_MARK)
    for i, (v, (w_, n_)) in enumerate(zip(vals, cnts)):
        ax_d.text(i, v + 0.3, f"{w_}/{n_}\n({v:.1f} %)", ha="center", fontsize=6.0,
                  linespacing=1.4)
    ax_d.set_ylabel("Promoters carrying a\nTGGCCGGC position (%)")
    ax_d.set_ylim(0, 9)

    # (e) cluster promoter vs rest of cluster
    ax_e = fig.add_subplot(gs_bot[0, 1])
    xp = np.arange(3)
    w = 0.34
    n_prom = {int(panel_b[("n_positions", PROM_KEY)][t]) for t in TPS}
    n_rest = {int(panel_b[("n_positions", REST_KEY)][t]) for t in TPS}
    assert len(n_prom) == 1 and len(n_rest) == 1, (n_prom, n_rest)
    ax_e.bar(xp - w / 2, [panel_b[("level_pooled_pct", PROM_KEY)][t] for t in TPS], w,
             color=C_MARK, label=f"promoter \u2212300..\u22121 ({n_prom.pop()} positions)")
    ax_e.bar(xp + w / 2, [panel_b[("level_pooled_pct", REST_KEY)][t] for t in TPS], w,
             color=C_MARK, hatch="////", edgecolor="white", linewidth=0.0,
             label=f"rest of cluster ({n_rest.pop()} positions)")
    ax_e.set_xticks(xp)
    ax_e.set_xticklabels(["T1  12 h", "T2  24 h", "T3  50 h"])
    ax_e.set_ylabel("TGGCCGGC 4mC level (%)")
    ax_e.set_ylim(0, 100)
    ax_e.set_yticks([0, 25, 50, 75, 100])
    ax_e.legend(frameon=False, fontsize=6.0, loc="upper left", handlelength=1.4,
                bbox_to_anchor=(-0.02, -0.19), labelspacing=0.30, handleheight=1.0)

    # (f) 4mC level against cluster expression
    ax_f = fig.add_subplot(gs_bot[1, 0])
    CC = {"silent": (U.COL_BAR_LIGHT, "^", "none"),
          "low": (U.COL_BAR_MID, "s", U.COL_BAR_MID),
          "expressed": (U.COL_BAR_DARK, "o", U.COL_BAR_DARK)}
    for cl in ("silent", "low", "expressed"):
        g_ = L4[L4.cls == cl]
        edge, mark, face = CC[cl]
        ax_f.plot(np.log2(g_.tpm + 1), g_.level_pooled_pct, ls="none", marker=mark,
                  ms=3.6, mfc=face, mec=edge, mew=0.9, alpha=0.9, label=cl)
    ax_f.set_xlabel("Core-gene expression, log$_2$(TPM + 1)")
    ax_f.set_ylabel("TGGCCGGC 4mC level (%)")
    lo, hi = L4.level_pooled_pct.min(), L4.level_pooled_pct.max()
    # Headroom for the key, which sits above the points rather than on them; the
    # plotted values are unchanged, only the view window.
    ax_f.set_ylim(lo - 0.06 * (hi - lo), hi + 0.34 * (hi - lo))
    ax_f.legend(frameon=False, fontsize=6.0, loc="upper center", ncol=3,
                handlelength=0.9, columnspacing=1.1, borderaxespad=0.10)
    ax_f.margins(x=0.06)

    # (g) cluster interior / flank ratio
    ax_g = fig.add_subplot(gs_bot[1, 1])
    MARKS = {"TGGCCGGC 4mC": ("o", U.COL_4mC), "AAGCCCG C5 4mC": ("s", U.COL_4mC),
             "AAGCCCG 6mA": ("D", U.COL_6mA)}
    for i, t in enumerate(dpan.itertuples()):
        mk, col = MARKS[t.mark_class]
        ax_g.plot(t.ratio, i, ls="none", marker=mk, ms=3.6, color=col, zorder=3)
    ax_g.axvline(1.0, color=META_GREY, lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax_g.set_yticks(range(len(dpan)))
    ax_g.set_yticklabels([f"{t.mark_class}  {t.timepoint}" for t in dpan.itertuples()],
                         fontsize=6.0)
    ax_g.set_ylim(-0.6, len(dpan) - 0.4)
    ax_g.set_xlabel("Methylation level, cluster interior / flank (ratio)")
    rlo, rhi = dpan.ratio.min(), dpan.ratio.max()
    ax_g.set_xlim(min(rlo, 1.0) - 0.012, max(rhi, 1.0) + 0.012)
    ax_g.set_xticks([0.95, 0.98, 1.00, 1.02])

    figure_panel_letters(fig, [ax_a, ax_b, ax_c], "abc", dx=0.082, dy=0.0045)
    figure_panel_letters(fig, [ax_d, ax_e, ax_f, ax_g], "defg", dx=0.090, dy=0.0040)
    return fig, series


if __name__ == "__main__":
    fig, series = build()
    out = save(fig, "SuppFigure_merged_BGC_261103.png")
    for t in series:
        print(f"  {t}: 4mC {series[t]['n4']:,} sites, 6mA {series[t]['n6']:,} sites")
    report(out)
