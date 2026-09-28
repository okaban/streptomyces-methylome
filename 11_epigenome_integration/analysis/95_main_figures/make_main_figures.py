#!/usr/bin/env python3
"""95 — the producing script for main Figures 1, 2 and 4.

Why this exists: these three figures were rebuilt on 2026-09-21 in throwaway
cells, so the PNGs shipped with NO producing script, ad-hoc colours, no font
ladder and widths of 182-196 mm. `Writing/check_figures.py` now fails on exactly
that. This script regenerates them from canonical inputs through the approved
shared style module (palette v2, font ladder, 174 mm NAR clamp).

Content changes adopted from the figure-design literature (see
Writing/FIGURE_DESIGN_RULES_260922.md):
  * Fig 1a is a SuperPlot (Lord et al. 2020, JCB): every site is drawn coloured
    by the biological replicate it came from, with the three replicate medians
    overlaid. The previous version pooled three replicates into one cloud.
  * Fig 1b and Fig 4a carry Wilson 95% CIs (Midway 2020, principle 5). Fig 4a
    compares 27/62 against 435/989 - without intervals two equal-height bars
    imply equal precision across a 16-fold difference in denominator.

Motif assignment and the 0-based `position` convention follow
90_per_timepoint_census_audit/per_timepoint_census.py exactly.

Usage:  python make_main_figures.py [--skip-pileup]
"""
import argparse, collections, csv, importlib.util, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
A = HERE.parent
REPO = A.parent.parent
TAB = HERE / "tables"; TAB.mkdir(exist_ok=True)
FIGDIR = HERE / "figures"; FIGDIR.mkdir(exist_ok=True)
SHIP = Path("/Users/okaban/obsidian/Research/rna-seq/Writing/fig_images")

REFFA = Path("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa")
PILEUP = Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup")
CANON = A / "01_integration/high_confidence_sites_weighted.csv"
CENSUS = A / "90_per_timepoint_census_audit/tables"
RESTRUCT = A / "94_figure_restructure_260921/tables"

T1_REPS = ["1-1", "1-2", "1-3"]
CODE_4mC, CODE_5mC = "21839", "m"
CORE = (1_500_000, 7_170_000)

REF = "".join(l.strip() for l in open(REFFA) if not l.startswith(">")).upper()
GENOME_LEN = len(REF)


def shared():
    p = REPO / "15_paper_figures/scripts/00_shared_utils.py"
    spec = importlib.util.spec_from_file_location("shared_utils", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["shared_utils"] = mod
    spec.loader.exec_module(mod)
    return mod


def motif(p, s, mod):
    """Identical to 90_/per_timepoint_census.py — 0-based position, strand-aware."""
    if mod == "4mC":
        if (REF[p - 2:p + 4] == "GCCGGC") if s == "+" else (REF[p - 3:p + 3] == "GCCGGC"):
            return "GCCGGC"
        if (REF[p - 4:p + 3] == "AAGCCCG") if s == "+" else (REF[p - 2:p + 5] == "CGGGCTT"):
            return "AAGCCCG"
    else:
        if s == "+" and (REF[p:p + 7] == "AAGCCCG" or REF[p - 1:p + 6] == "AAGCCCG"):
            return "AAGCCCG"
        if s == "-" and (REF[p - 6:p + 1] == "CGGGCTT" or REF[p - 5:p + 2] == "CGGGCTT"):
            return "AAGCCCG"
    return "other"


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def site_sets():
    """(timepoint, mod, motif) -> set of (pos, strand), from the canonical file."""
    S = collections.defaultdict(set)
    for r in csv.DictReader(open(CANON)):
        mod = "4mC" if "4mC" in r["mod_type"] else ("6mA" if "6mA" in r["mod_type"] else r["mod_type"])
        p, s = int(r["position"]), r["strand"]
        S[(r["timepoint"], mod, motif(p, s, mod))].add((p, s))
    return S


def per_replicate_table(sites):
    """4mC and 5mC percentage at each canonical T1 GCCGGC site, per replicate."""
    want = set(sites)
    rows = []
    for rep in T1_REPS:
        got = {}
        with open(PILEUP / f"{rep}_pileup.bed") as fh:
            for line in fh:
                c = line.split("\t", 11)
                if len(c) < 11 or c[3] not in (CODE_4mC, CODE_5mC):
                    continue
                key = (int(c[1]), c[5])
                if key in want:
                    got[(key, c[3])] = (float(c[10]), int(c[4]))
        n = 0
        for key in want:
            f4, f5 = got.get((key, CODE_4mC)), got.get((key, CODE_5mC))
            if f4 is None and f5 is None:
                continue
            rows.append(dict(replicate=rep, position=key[0], strand=key[1],
                             pct_4mC=(f4 or (np.nan, 0))[0], cov_4mC=(f4 or (np.nan, 0))[1],
                             pct_5mC=(f5 or (np.nan, 0))[0], cov_5mC=(f5 or (np.nan, 0))[1]))
            n += 1
        print(f"  {rep}: {n} of {len(want)} canonical sites covered", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "t1_per_replicate_sites.tsv", sep="\t", index=False)
    return df



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-pileup", action="store_true")
    args = ap.parse_args()
    su = shared()
    su.apply_unified_style() if hasattr(su, "apply_unified_style") else su.apply_style()

    S = site_sets()
    t1 = sorted(S[("T1", "4mC", "GCCGGC")])
    print(f"canonical T1 GCCGGC 4mC sites: {len(t1)}")

    cache = TAB / "t1_per_replicate_sites.tsv"
    if args.skip_pileup and cache.exists():
        rep = pd.read_csv(cache, sep="\t")
    else:
        print("reading T1 pileups (a few minutes)...")
        rep = per_replicate_table(t1)

    # ---- Figure 1 -----------------------------------------------------------
    fig1 = plt.figure(figsize=(su.mm_to_inch(174), su.mm_to_inch(76)))
    a1 = fig1.add_axes([0.095, 0.40, 0.145, 0.42])
    a2 = fig1.add_axes([0.400, 0.40, 0.170, 0.42])
    a3 = fig1.add_axes([0.750, 0.40, 0.165, 0.42])

    # (a) SuperPlot: sites coloured by replicate, replicate medians overlaid
    rng = np.random.default_rng(42)
    rep_cols = [su.COL_BAR_LIGHT, su.COL_BAR_MID, su.COL_BAR_DARK]
    med_rows = []
    for ch, xc, lab in ((("pct_4mC"), 0, "4mC"), (("pct_5mC"), 1, "5mC")):
        for i, r in enumerate(T1_REPS):
            v = rep.loc[rep.replicate == r, ch].dropna().values
            if len(v) == 0:
                continue
            x = xc + (i - 1) * 0.22 + rng.normal(0, 0.035, len(v))
            a1.scatter(x, v, s=1.2, c=rep_cols[i], alpha=0.35, linewidths=0, zorder=2)
            med_rows.append(dict(channel=ch, replicate=r, n=len(v), median=float(np.median(v))))
            a1.plot([xc + (i - 1) * 0.22], [np.median(v)], marker="o", ms=4.5,
                    mfc=rep_cols[i], mec=su.COL_DARK, mew=0.6, zorder=4)
    pd.DataFrame(med_rows).to_csv(TAB / "fig1a_superplot_medians.tsv", sep="\t", index=False)
    a1.set_xticks([0, 1]); a1.set_xticklabels(["4mC", "5mC"])
    a1.set_ylabel("modification frequency (%)")
    a1.set_ylim(-4, 104); a1.set_xlim(-0.55, 1.55)
    a1.set_title("The GCCGGC mark is 4mC,\nnot 5mC", loc="left")
    # The three large circles are the replicate medians; named in the caption
    # rather than in a legend box, so no text sits inside the data region.

    # (b) core fraction per timepoint with Wilson CI, against the motif null
    core_rows = []
    for tp in ("T1", "T2", "T3"):
        st = S[(tp, "4mC", "GCCGGC")]
        k = sum(1 for p, _ in st if CORE[0] <= p <= CORE[1])
        lo, hi = wilson(k, len(st))
        core_rows.append(dict(timepoint=tp, n=len(st), in_core=k, frac=k / len(st), lo=lo, hi=hi))
    inst = [m.start() for m in __import__("re").finditer("GCCGGC", REF)]
    null = sum(1 for p in inst if CORE[0] <= p <= CORE[1]) / len(inst)
    nlo, nhi = wilson(sum(1 for p in inst if CORE[0] <= p <= CORE[1]), len(inst))
    cf = pd.DataFrame(core_rows); cf.to_csv(TAB / "fig1b_core_fraction_wilson.tsv", sep="\t", index=False)
    x = np.arange(3)
    a2.bar(x, cf.frac * 100, color=su.COL_CORE, width=0.6, zorder=2)
    a2.errorbar(x, cf.frac * 100, yerr=[(cf.frac - cf.lo) * 100, (cf.hi - cf.frac) * 100],
                fmt="none", ecolor=su.COL_DARK, elinewidth=0.8, capsize=2, zorder=3)
    a2.axhspan(nlo * 100, nhi * 100, color=su.COL_GRAY, alpha=0.35, zorder=1)
    a2.axhline(null * 100, color=su.COL_GRAY, ls="--", lw=0.8, zorder=1)
    # Direct-label the null line in the margin the bars do not occupy, rather than
    # over the bars themselves (the previous placement overlapped them).
    a2.set_xlim(-0.55, 3.25)
    a2.text(2.62, null * 100, f"motif-matched null\n{null*100:.0f}% (95% CI)", ha="left",
            va="center", color=su.COL_GRAY, fontsize=su.FONT["annot"])
    a2.set_xticks(x); a2.set_xticklabels(su.TP_LABELS_NL)
    a2.spines["right"].set_visible(False)
    a2.set_ylabel("GCCGGC 4mC in core (%)")
    a2.set_ylim(0, 100)
    a2.set_title("Core-concentrated at\nevery timepoint", loc="left")

    # (c) promoter occupancy vs expression change, split by compartment
    d = pd.read_csv(RESTRUCT / "occupancy_vs_lfc.tsv", sep="\t")
    for rg, col, lab in ((0, su.COL_CORE, "core"), (1, su.COL_ARM, "arm")):
        sub = d[d.region == rg]
        a3.scatter(sub.occ2k, sub.LFC_T2vsT1, s=3, c=col, alpha=0.55, linewidths=0,
                   label=f"{lab} (n={len(sub)})")
    a3.axhline(0, color=su.COL_DARK, lw=0.7)
    raw = spearmanr(d.occ2k, d.LFC_T2vsT1)[0]
    a3.set_xlim(-20, 520)
    a3.set_xlabel("promoter GCCGGC\noccupancy (±2 kb, T1)")
    a3.set_ylabel("log$_2$ FC (T2 vs T1)")
    a3.set_title(("Coupling is positional:\n"
                  f"r {raw:.2f} → −0.09 controlled").replace("-", "−"), loc="left")
    # key below the panel: inside the axes it would sit on the point cloud
    a3.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.46), ncol=2,
              handletextpad=0.3, columnspacing=1.4, markerscale=2.2,
              fontsize=su.FONT["legend"])
    for ax, L in ((a1, "a"), (a2, "b"), (a3, "c")):
        ax.text(-0.30, 1.20, L, transform=ax.transAxes, fontweight="bold",
                fontsize=su.FONT["panel_letter"], va="bottom", ha="left")
    su.assert_no_text_collisions(fig1, "Figure1")
    su.save_figure(fig1, FIGDIR / "Figure1", formats=("png", "pdf"), width_class="full")
    plt.close(fig1)
    print("Figure1 written")

    # ---- Figure 2 -----------------------------------------------------------
    fig2 = plt.figure(figsize=(su.mm_to_inch(174), su.mm_to_inch(86)))
    b1 = fig2.add_axes([0.140, 0.68, 0.800, 0.22])
    b2 = fig2.add_axes([0.140, 0.16, 0.320, 0.26])
    b3 = fig2.add_axes([0.700, 0.16, 0.150, 0.26])

    # (a) where the sites are, timepoint by timepoint — ordered lightness ramp
    BIN = 100_000
    edges = np.arange(0, GENOME_LEN + BIN, BIN)
    for i, tp in enumerate(("T1", "T2", "T3")):
        pos = np.array([p for p, _ in S[(tp, "4mC", "GCCGGC")]])
        h, _ = np.histogram(pos, bins=edges)
        b1.plot(edges[:-1] / 1e6, h, lw=1.1, color=su.TIME_RAMP[i], ls=su.TIME_STYLE[i],
                label=su.TP_LABELS[i], zorder=2 + i)
    b1.axvspan(CORE[0] / 1e6, CORE[1] / 1e6, color=su.COL_GRAY, alpha=0.18, zorder=1)
    b1.set_xlim(0, GENOME_LEN / 1e6)
    b1.set_xticks(np.arange(0, 9, 2))
    b1.set_xlabel("chromosomal position (Mb)")
    b1.set_ylabel(f"GCCGGC 4mC sites\nper {BIN//1000} kb")
    b1.set_title("The site census is largely stable; the shaded band is the core",
                 loc="left")
    # key above the track: three lines fill the panel, so nothing can sit inside it
    b1.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=3,
              handletextpad=0.4, columnspacing=1.6, fontsize=su.FONT["legend"])

    # (b) persistence patterns, split by compartment
    pp = pd.read_csv(RESTRUCT / "gccggc_persistence_patterns.tsv", sep="\t")
    pp["pattern"] = pp["pattern"].astype(str).str.zfill(3)
    order = (pp.groupby("pattern")["n"].sum().sort_values(ascending=False).index.tolist())
    core_n = [int(pp[(pp.pattern == k) & (pp.region == "core")].n.sum()) for k in order]
    arm_n = [int(pp[(pp.pattern == k) & (pp.region == "arm")].n.sum()) for k in order]
    xs = np.arange(len(order))
    b2.bar(xs, core_n, color=su.COL_CORE, width=0.66, label="core", zorder=2)
    b2.bar(xs, arm_n, bottom=core_n, color=su.COL_ARM, width=0.66, label="arm", zorder=2)
    b2.set_xticks(xs); b2.set_xticklabels(order)
    # 2026-09-26: the raw pattern keys ("111", "110") were unreadable. Spell out
    # which timepoints the group was detected at.
    LAB = {"111": "all three", "110": "T1+T2", "011": "T2+T3", "101": "T1+T3",
           "100": "T1 only", "010": "T2 only", "001": "T3 only"}
    b2.set_xticks(xs)
    b2.set_xticklabels([LAB.get(k, k) for k in order], rotation=30, ha="right")
    b2.set_xlabel("timepoints at which the site was detected")
    b2.set_ylabel("distinct GCCGGC 4mC sites")
    b2.set_title(f"{core_n[0] + arm_n[0]:,} of {sum(core_n) + sum(arm_n):,} sites are seen\n"
                 f"at all three timepoints", loc="left")
    b2.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2,
              handletextpad=0.4, columnspacing=1.2, fontsize=su.FONT["legend"])

    # (c) the 62 Exposed promoters keep their mark
    ex = pd.read_csv(CENSUS / "exposed62_per_timepoint.tsv", sep="\t")
    b3.bar(np.arange(3), ex.within_293bp, color=su.COL_EXPOSED, width=0.6, zorder=2)
    b3.axhline(62, color=su.COL_DARK, ls=":", lw=0.8, zorder=1)
    b3.set_xticks(np.arange(3)); b3.set_xticklabels(su.TP_LABELS_NL)
    b3.set_ylabel("Exposed promoters with a\nGCCGGC 4mC site ≤293 bp")
    b3.set_ylim(0, 72)
    b3.set_title("Retained, not erased", loc="left")
    # Panel letters in FIGURE coordinates: an axes-relative offset collides with
    # a tall rotated y-label, whose extent depends on the label text.
    for ax, L in ((b1, "a"), (b2, "b"), (b3, "c")):
        bb = ax.get_position()
        lx = 0.022 if ax is not b3 else 0.605
        fig2.text(lx, bb.y1 + 0.025, L, fontweight="bold",
                  fontsize=su.FONT["panel_letter"], va="bottom", ha="left")
    su.assert_no_text_collisions(fig2, "Figure2")
    su.save_figure(fig2, FIGDIR / "Figure2_merged", formats=("png", "pdf"), width_class="full")
    plt.close(fig2)
    print("Figure2_merged written")

    # ---- Figure 4 (ships as Figure5.png) ------------------------------------
    # all_genes_features.tsv holds only 11 Exposed genes; the manuscript's universe
    # is the unified table (62 Exposed of 1,055 rows). Denominators here are genes
    # that HAVE an expression estimate (59 and 960) — the class sizes 62 and 989
    # quoted in the text include genes that cannot enter the numerator.
    reg = pd.read_csv(
        A / "52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv", sep="\t")
    E = reg[reg.is_exposed == 1]; Sh = reg[reg.is_exposed == 0]
    rows = []
    for nm, sub in (("Exposed", E), ("Shielded", Sh)):
        v = sub.LFC_T2vsT1.dropna()
        k = int((v.abs() >= 1).sum())
        lo, hi = wilson(k, len(v))
        rows.append(dict(cls=nm, n=len(v), n_de=k, frac=k / len(v), lo=lo, hi=hi,
                         median_abs_lfc=float(v.abs().median())))
    de = pd.DataFrame(rows); de.to_csv(TAB / "fig4a_de_fraction_wilson.tsv", sep="\t", index=False)

    fig4 = plt.figure(figsize=(su.mm_to_inch(174), su.mm_to_inch(72)))
    c1 = fig4.add_axes([0.085, 0.38, 0.135, 0.44])
    c2 = fig4.add_axes([0.400, 0.38, 0.185, 0.44])
    c3 = fig4.add_axes([0.760, 0.38, 0.185, 0.44])
    cols = [su.COL_EXPOSED, su.COL_SHIELDED]
    c1.bar(np.arange(2), de.frac * 100, color=cols, width=0.6, zorder=2)
    c1.errorbar(np.arange(2), de.frac * 100,
                yerr=[(de.frac - de.lo) * 100, (de.hi - de.frac) * 100],
                fmt="none", ecolor=su.COL_DARK, elinewidth=0.8, capsize=2.5, zorder=3)
    c1.set_xticks(np.arange(2))
    c1.set_xticklabels([f"Exposed\n({de.n_de[0]}/{de.n[0]})", f"Shielded\n({de.n_de[1]}/{de.n[1]})"])
    # criterion and denominator stated in the caption, not inside the panel
    c1.set_ylabel("responding genes (%)")
    c1.set_ylim(0, 72)
    c1.set_title("Marked promoters are no more\nlikely to respond", loc="left")

    for nm, sub, col in (("Exposed", E, su.COL_EXPOSED), ("Shielded", Sh, su.COL_SHIELDED)):
        v = np.sort(sub.LFC_T2vsT1.dropna().abs().values)
        c2.plot(v, np.arange(1, len(v) + 1) / len(v), color=col, lw=1.2, label=nm)
    c2.set_xlim(0, 6); c2.set_xlabel("|log$_2$ FC| (T2 vs T1)")
    c2.set_ylabel("cumulative fraction of genes")
    c2.set_title("and their response sizes\ncoincide", loc="left")
    c2.legend(frameon=False, loc="lower right", handletextpad=0.5,
              fontsize=su.FONT["legend"])

    fo = pd.read_csv(RESTRUCT / "fig4_correlation_forest.tsv", sep="\t")
    yy = np.arange(len(fo))[::-1]
    c3.errorbar(fo.r, yy, xerr=[fo.r - fo.lo, fo.hi - fo.r], fmt="o", ms=3.4,
                color=su.COL_BAR_DARK, ecolor=su.COL_BAR_DARK, elinewidth=0.9, capsize=2)
    c3.axvline(0, color=su.COL_DARK, lw=0.7)
    c3.set_yticks(yy)
    c3.set_yticklabels([f"{l} (n={n})" for l, n in zip(fo.label.str.replace(chr(10), " "), fo.n)])
    c3.set_xlabel("Spearman r, occupancy vs log$_2$ FC")
    c3.set_title("The coupling is carried by\nchromosomal position", loc="left")
    for ax, L, lx in ((c1, "a", 0.020), (c2, "b", 0.315), (c3, "c", 0.672)):
        bb = ax.get_position()
        fig4.text(lx, bb.y1 + 0.035, L, fontweight="bold",
                  fontsize=su.FONT["panel_letter"], va="bottom", ha="left")
    su.assert_no_text_collisions(fig4, "Figure4")
    su.save_figure(fig4, FIGDIR / "Figure5", formats=("png", "pdf"), width_class="full")
    plt.close(fig4)
    print("Figure5 (manuscript Figure 4) written")


if __name__ == "__main__":
    sys.exit(main())
