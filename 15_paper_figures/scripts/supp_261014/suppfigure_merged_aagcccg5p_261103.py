#!/usr/bin/env python3
"""Merged supplementary figure — the 5' extension of AAGCCCG and what determines
whether the adenine is modified.

Consolidation of two shipped supplementary figures (author approval 2026-11-03):

    old Supplementary Figure 12  "The AAGCCCG recognition sequence extends 5' of
                                  the heptamer"          -- all-read calls
        generator: recompute_261015/aagcccg261016/make_fig_261016.py
    old Supplementary Figure 13  "The 5' context determines whether the AAGCCCG
                                  adenine is modified"   -- filtered site calls
        generator: supp_261014/suppfigure14_determinants_261020.py

The panels are re-laid on ONE canvas; neither original generator is called.  The
computation of every panel is the code of its original generator, copied
unchanged, so no plotted number moves.  Panel letters are renumbered:

    old S12 a,b,c  ->  a,b,c   (top row, all-read series)
    old S13 a,b,c  ->  d,e,f   (bottom row, filtered site calls)

Only the axis-label wording changed (units added in parentheses, FIG rule).

    SUPPFIG_OUT=<dir> PYTHONPATH=. python3 supp_261014/suppfigure_merged_aagcccg5p_261103.py
"""
import os
import pathlib
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from _suppstyle import (BLUE, GREY, META_GREY, NAR_WIDTH_IN, RED,  # noqa: E402
                        apply_figure_style, figure_panel_letters, report, save,
                        source_note)

RNA = pathlib.Path("/Users/okaban/bioinfo/rna-seq")
AAG = RNA / "15_paper_figures/scripts/recompute_261015/aagcccg261016"
S = RNA / "15_paper_figures/scripts/recompute_261015/pileup_cache_261015"
GEN = RNA / "11_epigenome_integration/data/NC_003888.3.fna"

# ───────────────────────── top row: all-read tables (old S12) ────────────────
LAD = pd.read_csv(AAG / "ladder_allread_261016.tsv", sep="\t")
PR = pd.read_csv(AAG / "H1_positional_profile_allread_261016.tsv", sep="\t")
CTRL = pd.read_csv(AAG / "minus3_specificity_controls_261016.tsv", sep="\t")
C4 = pd.read_csv(AAG / "fiveprime_contrast_261016.tsv", sep="\t")

TPS = ["T1", "T2", "T3"]
TPC = ["#2b2b2b", "#767676", "#b4b4b4"]      # FIG23: timepoints are neutral greys

ORDER = ["AAGCCCG (all instances)", "HAAGCCCG", "G-N-HAAGCCCG  (-3=G)",
         "C-N-G-N-HAAGCCCG  (-3=G,-5=C)"]
CHANS = ["6mA A0 | HAAGCCCG target", "6mA A0 | GAAGCCCG non-target",
         "4mC C4 | HAAGCCCG instances"]
# Storage keys in the .tsv files keep the old zero-based names (A0, C4); the
# DISPLAYED names follow the manuscript rule "1-based from the 5' end", so the
# adenine is A1 and the cytosine is C5 (D047, 2026-10-24).

# ──────────────────── bottom row: pileup scans (old S13) ─────────────────────
seq = "".join(l.strip() for l in GEN.read_text().split("\n")[1:]).upper()
COMP = str.maketrans("ACGT", "TGCA")
LIB = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"],
       "T3": ["3-2", "3-3", "3-4"]}
POOL = {}


def pooled(tp, code):
    k = (tp, code)
    if k in POOL:
        return POOL[k]
    acc = {}
    for lib in LIB[tp]:
        d = np.load(S / f"{lib}.npz")
        m = (d["code"] == code) & (d["refok"] == 1)
        st, pl, nv, nm, nf = d["start"][m], d["plus"][m], d["Nvalid"][m], d["Nmod"][m], d["Nfail"][m]
        key = st.astype(np.int64) * 2 + pl.astype(np.int64)
        for k2, v, mo, fa in zip(key, nv, nm, nf):
            a = acc.get(k2)
            if a is None:
                acc[k2] = [v, mo, fa]
            else:
                a[0] += v; a[1] += mo; a[2] += fa
    POOL[k] = acc
    return acc


def instances(motif, mod_offset):
    out = []
    for m in re.finditer(f"(?={motif})", seq):
        out.append((m.start() + mod_offset, 1, m.start()))
    rcm = motif.translate(COMP)[::-1]
    L = len(motif)
    for m in re.finditer(f"(?={rcm})", seq):
        out.append((m.start() + (L - 1 - mod_offset), 0, m.start()))
    return out


def flank_base(start, plus, L, off):
    if plus:
        i = start + (off if off < 0 else L - 1 + off)
        return seq[i] if 0 <= i < len(seq) else None
    i = start + (L - 1 - off if off < 0 else -off)
    if not (0 <= i < len(seq)):
        return None
    return seq[i].translate(COMP)


is_H = lambda st, plus, L: flank_base(st, plus, L, -1) in ("A", "C", "T")


def scan_cond(motif, mod_offset, code, label, cond=None, positions=range(-8, 9)):
    L = len(motif)
    rows = []
    for tp in ("T1", "T2", "T3"):
        acc = pooled(tp, code)
        dat = []
        for pos, plus, st in instances(motif, mod_offset):
            if cond and not cond(st, plus, L):
                continue
            a = acc.get(pos * 2 + plus)
            if a is None or (a[0] + a[2]) < 30:
                continue
            s = 1 if a[0] and a[1] / a[0] >= 0.5 else 0
            dat.append((st, plus, s))
        for off in positions:
            if off == 0:
                continue
            tab = {}
            for st2, plus2, s2 in dat:
                b = flank_base(st2, plus2, L, off)
                if b in "ACGT":
                    tab.setdefault(b, [0, 0])[s2] += 1
            nm = sum(v[1] for v in tab.values())
            nu = sum(v[0] for v in tab.values())
            for b, (u, mo) in sorted(tab.items()):
                if min(mo + u, nm + nu - mo - u) < 5:
                    continue
                orr, p = fisher_exact([[mo, nm - mo], [u, nu - u]])
                rows.append(dict(motif=label, timepoint=tp, position=off, base=b,
                                 n_mod_with=mo, n_mod_total=nm, n_unmod_with=u,
                                 n_unmod_total=nu,
                                 odds_ratio=round(float(orr), 4) if np.isfinite(orr) else np.nan,
                                 p=p))
    df = pd.DataFrame(rows)
    for tp2 in df.timepoint.unique():
        m2 = df.timepoint == tp2
        df.loc[m2, "q_BH"] = multipletests(df.loc[m2, "p"], method="fdr_bh")[1]
        df.loc[m2, "n_tests"] = int(m2.sum())
    sig = df[df.q_BH < 0.05].groupby(["position", "base"]).timepoint.nunique()
    df["survives_3tp"] = df.apply(lambda r: sig.get((r.position, r.base), 0) == 3, axis=1)
    return df


d6H = scan_cond("AAGCCCG", 0, 0, "HAAGCCCG 6mA(A1)", cond=is_H)
d4A = scan_cond("AAGCCCG", 4, 1, "AAGCCCG 4mC(pos5)", cond=is_H)

ctl_rows = []
for tp in ("T1", "T2", "T3"):
    acc = pooled(tp, 1)
    rec = []
    for pos, plus, st in instances("TGGCCGGC", 4):
        a = acc.get(pos * 2 + plus)
        if a is None or (a[0] + a[2]) < 30:
            continue
        s = 1 if a[0] and a[1] / a[0] >= 0.5 else 0
        up = "".join(flank_base(st, plus, 8, o) or "N" for o in range(-10, 0))
        dn = "".join(flank_base(st, plus, 8, o) or "N" for o in range(1, 11))
        f8 = "".join(flank_base(st, plus, 8, o) or "N"
                     for o in list(range(-8, 0)) + list(range(1, 9)))
        rec.append(dict(status=s, gc8=(f8.count("G") + f8.count("C")) / 16,
                        grun=int(bool(re.search(r"G{3,}", up + dn))),
                        failfrac=a[2] / (a[0] + a[2])))
    d = pd.DataFrame(rec)
    mo, un = d[d.status == 1], d[d.status == 0]
    ctl_rows.append(dict(timepoint=tp, n_mod=len(mo), n_unmod=len(un),
                         gc8_mod=mo.gc8.mean(), gc8_unmod=un.gc8.mean(),
                         grun_mod=mo.grun.mean(), grun_unmod=un.grun.mean(),
                         fail_mod=mo.failfrac.mean(), fail_unmod=un.failfrac.mean()))
ctl = pd.DataFrame(ctl_rows)

LADDER = HERE / "data" / "occupancy_ladder_261020.tsv"
lad = pd.read_csv(LADDER, sep="\t")
assert len(lad) == 15 and set(lad.timepoint) == {"T1", "T2", "T3"}, lad.shape
lad2 = lad.copy()
lad2["tier"] = lad2.definition.map({"AAGCCCG (all)": 1, "HAAGCCCG": 2, "+ -3=G": 3,
                                    "+ -5=C & -3=G": 4,
                                    "(ref) same-instance 4mC C5": "ref"})
lad2 = lad2[["timepoint", "tier", "definition", "n", "called", "occupancy_pct"]]

# ───────────────────────────────── figure ────────────────────────────────────
apply_figure_style(frame="boxed", sizes=(7.0, 6.2, 6.0))


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 5.80))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.60, 1.0, 1.0], hspace=0.62,
                          wspace=0.52, left=0.135, right=0.990, top=0.960,
                          bottom=0.1650)
    ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2]); ax_d = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1]); ax_f = fig.add_subplot(gs[1, 2])

    # (a) all-read extension ladder -----------------------------------------
    ax = ax_a
    n = [int(LAD[(LAD.definition == o) & (LAD.tp == "T1")].n_sites.iloc[0]) for o in ORDER]
    lab = [f"AAGCCCG  {n[0]}", f"HAAGCCCG  {n[1]}",
           f"+ G\u208b\u2083  {n[2]}", f"+ C\u208b\u2085 G\u208b\u2083  {n[3]}"]
    y = np.arange(len(ORDER))[::-1]
    w = 0.26
    for i, tp in enumerate(TPS):
        v = [LAD[(LAD.definition == o) & (LAD.tp == tp)].occupancy_pct.iloc[0] for o in ORDER]
        ax.barh(y - (i - 1) * w, v, w, color=TPC[i], label=tp, edgecolor="none")
    c4 = C4.C4_H_occ.mean()
    ax.axvline(c4, color=RED, lw=1.0, ls="--", label="4mC at C5 (same sites)")
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=6.0)
    ax.set_xlabel("6mA occupancy at A1 (%)\n[site count at right of label]")
    ax.set_xlim(0, 100)
    ax.legend(frameon=False, loc="lower right", fontsize=6.0, handlelength=1.0)

    # (b) positional profile, all-read ---------------------------------------
    ax = ax_b
    for i, tp in enumerate(TPS):
        s = PR[PR.tp == tp].sort_values("offset")
        ax.plot(s.offset, -np.log10(s.q.clip(lower=1e-20)), "o-", ms=2.2, lw=0.8,
                color=TPC[i], label=tp)
    ax.axvspan(-0.5, 6.5, color="#eeeeee", zorder=0)
    ax.text(3, 2.5, "AAGCCCG", ha="center", fontsize=6.0, color="#555555")
    ax.axhline(-np.log10(0.05), color=GREY, lw=0.7, ls=":", label="BH 5%")
    ax.plot([-5, -3], [-0.9, -0.9], ls="none", marker="^", ms=3.2, color="#333333",
            clip_on=False, label="passes registered scan")
    ax.set_xlabel("Offset from A1 (bp)")
    ax.set_ylabel("\u2212log\u2081\u2080 q (best base/offset)")
    ax.legend(frameon=False, fontsize=6.0, handlelength=1.0, loc="upper right",
              ncol=1, labelspacing=0.25, borderpad=0.1)

    # (c) -3 specificity control ---------------------------------------------
    ax = ax_c
    clab = ["6mA, target\nHAAGCCCG", "6mA, non-target\nGAAGCCCG", "4mC at C5,\nsame instances"]
    ccol = [BLUE, GREY, RED]
    bases = list("ACGT")
    xb = np.arange(4)
    for i, (ch, cc) in enumerate(zip(CHANS, ccol)):
        s = CTRL[(CTRL.channel == ch) & (CTRL.tp == "T1")].set_index("base_at_minus3")
        ax.plot(xb, [s.loc[b, "occupancy_pct"] for b in bases], "o-", ms=3, lw=1.1,
                color=cc, label=clab[i])
    ax.set_xticks(xb)
    ax.set_xticklabels(bases)
    ax.set_xlabel("Base at offset \u22123 (nt)")
    ax.set_ylabel("Occupancy (%), T1")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=6.0, handlelength=1.0, loc="center left")

    # (d) filtered occupancy ladder ------------------------------------------
    ax = ax_d
    tiers = [1, 2, 3, 4]
    ref = lad2[lad2.tier == "ref"].occupancy_pct
    for tp, mk in zip(("T1", "T2", "T3"), ("o", "s", "^")):
        g = lad2[(lad2.timepoint == tp) & (lad2.tier != "ref")].sort_values("tier")
        ax.plot(tiers, g.occupancy_pct, marker=mk, color=BLUE, lw=1.2, ms=3.6, label=tp)
    ax.axhspan(ref.min(), ref.max(), color=RED, alpha=0.18, lw=0)
    ax.text(1.05, ref.mean() + 3.0, "4mC, same instances", color=RED, fontsize=6.0,
            va="bottom", ha="left")
    ax.set_xticks(tiers)
    ax.set_xticklabels([f"{l}\nn={n_}" for l, n_ in zip(
        ["AAGCCCG", "HAAGCCCG", "+ \u22123 G", "+ \u22125 C"],
        lad2[(lad2.timepoint == "T1") & (lad2.tier != "ref")].sort_values("tier").n)],
        fontsize=6.0)
    ax.set_ylim(0, 100)
    ax.set_xlim(0.7, 4.4)
    ax.set_ylabel("6mA occupancy (%)")
    ax.legend(frameon=False, fontsize=6, loc="lower right", handlelength=1.2)

    # (e) positional scan, filtered ------------------------------------------
    ax = ax_e
    for d_, col, nm, mk in ((d6H, BLUE, "6mA (position 1)", "o"),
                            (d4A, RED, "4mC (position 5)", "s")):
        t1 = d_[d_.timepoint == "T1"]
        best = t1.loc[t1.groupby("position").q_BH.idxmin()].sort_values("position")
        ax.plot(best.position, best.odds_ratio, marker=mk, ms=3.0, lw=1.0, color=col,
                label=nm, alpha=0.85)
        sp = set(d_[d_.survives_3tp].position)
        for r in best[best.position.isin(sp)].itertuples():
            ax.plot(r.position, r.odds_ratio, marker=mk, ms=7, mfc="none", mec=col,
                    mew=1.4, lw=0)
            ax.annotate(r.base, (r.position, r.odds_ratio), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=6.5, color=col, weight="bold")
    ax.axhline(1, color=META_GREY, lw=0.8, ls=":")
    ax.set_yscale("log")
    ax.set_ylim(0.16, 7.5)
    ax.set_yticks([0.2, 1.0, 5.0])
    ax.set_yticklabels(["0.2", "1", "5"])
    ax.minorticks_off()
    ax.set_xlabel("Position relative to the motif (bp)")
    ax.set_ylabel("Odds ratio, most significant base")
    ax.legend(frameon=False, fontsize=6.0, loc="lower left", handlelength=1.1,
              bbox_to_anchor=(0.0, 0.0), ncol=1, labelspacing=0.3, borderpad=0.2)

    # (f) TGGCCGGC detection control -----------------------------------------
    ax = ax_f
    met = [("gc8_mod", "gc8_unmod"), ("grun_mod", "grun_unmod"), ("fail_mod", "fail_unmod")]
    x = np.arange(3)
    w = 0.36
    mv = [ctl[a].mean() for a, _ in met]
    uv = [ctl[b].mean() for _, b in met]
    ax.bar(x - w / 2, mv, w, color=RED, label="TGGCCGGC methylated\n1,886 positions")
    ax.bar(x + w / 2, uv, w, color=META_GREY, label="TGGCCGGC unmethylated\n205 positions")
    ax.set_xticks(x)
    ax.set_xticklabels(["\u00b18 bp GC", "G-run,\n\u00b110 bp", "call-failure\nfraction"],
                       fontsize=6.0)
    ax.set_ylim(0, 1.30)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    ax.set_ylabel("TGGCCGGC positions (fraction)")
    ax.legend(frameon=False, fontsize=6, loc="upper center", ncol=1, handlelength=1.2,
              labelspacing=0.25, borderpad=0.1)

    axes = [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f]
    figure_panel_letters(fig, axes, "abcdef", dx=0.058, dy=0.010)
    source_note(fig,
                "Panels a-c: all-read calls from the --no-filtering pileup, nine libraries, positions with \u2265 30 reads (aagcccg261016/*_261016.tsv).\n"
                "Panels d-f: pooled site calls \u2014 three libraries per timepoint, pooled true depth (Nvalid + Nfail) \u2265 30, reference base agreement,\n"
                "pooled Nmod/Nvalid \u2265 0.5. Panel e: Fisher tests at offsets \u22128..+8, Benjamini-Hochberg within each timepoint; a ring marks a base\n"
                "that is significant at all three timepoints.",
                y=0.008)
    print("panel f values (methylated, unmethylated):",
          {k: (round(mv[i], 2), round(uv[i], 2))
           for i, k in enumerate(["+-8 bp GC", "G-run within +-10 bp", "call-failure fraction"])})
    return fig


if __name__ == "__main__":
    report(save(build(), "SuppFigure_merged_AAGCCCG5p_261103.png",
                allow_over_data=("AAGCCCG", "4mC at C5 (same sites)",
                                 "4mC at C5,\nsame instances", "T1", "T2", "T3",
                                 "BH 5%", "passes registered scan",
                                 "4mC, same instances")))
