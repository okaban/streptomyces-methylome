#!/usr/bin/env python3
"""Supplementary Figure 14 — why no cluster-level methylation/expression
relationship is detectable.

Rebuild of the 212.7 mm original at 174 mm column width (2026-10-14).  Panel
values are re-derived from the deposited pileups, the Jeong 2016 experimental
TSS set and the 29-region BGC crosswalk, and asserted against the published
legend before anything is drawn.

2026-10-22 — panels b, c and d recomputed under the CANONICAL measurability rule.
Until now they used MIN_VALID = 10 (a per-library Nvalid floor), which is the
withdrawn depth filter: it carries a +15.5 pp upward bias because requiring
Nvalid alone to clear a floor selects positions whose calls were mostly valid,
i.e. mostly modified.  The canonical rule instead pools the three libraries of a
timepoint and asks for true depth (Nvalid + Nfail) >= 30, reference-base
agreement, and — where a SITE rather than a LEVEL is counted — pooled
Nmod/Nvalid >= 0.5.  Values now come from

    recompute_261015/recompute_suppfig8_bcd_261022.py
    -> recompute_261015/suppfig8_panelbcd_canonical_261022.tsv

which keeps the numerator and denominator of every plotted quantity.  The three
panels' conclusions are unchanged (promoter and cluster interior do not differ
from their comparators); the position counts rise because pooling admits
positions no single library could measure.  Panel a is unchanged.

    PYTHONPATH=. python3 suppfigure14_bgc_relationships.py
"""
import pathlib
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

from _gccggc import genome
from _suppstyle import (DATA_DIR, META_GREY, NAR_WIDTH_IN, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

FEATURECOUNTS = "/Users/okaban/bioinfo/rna-seq/03_quantification/data/featureCounts/featureCounts_M145_s2.txt"
TSS_CSV = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv"
CANON = pathlib.Path("/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/"
                     "recompute_261015/suppfig8_panelbcd_canonical_261022.tsv")
MIN_POS = 3             # panel c: regions with fewer measurable positions are not plotted
TPS = ("T1", "T2", "T3")

# Canonical mark palette (00_shared_utils, author decision 2026-09-26): 4mC is the
# muted red, 6mA the muted blue, genome-wide, in every main and supplementary
# figure.  Panels a-c all plot the same TGGCCGGC 4mC quantity, so they thread the
# 4mC colour; promoter vs rest of cluster is separated by hatch, not by hue.
import importlib  # noqa: E402
sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
U = importlib.import_module("00_shared_utils")
C_MARK = U.COL_4mC

G = genome()
tgg = pd.DataFrame(
    sorted([(m.start() + 4, "+") for m in re.finditer("(?=TGGCCGGC)", G)]
           + [(m.start() + 3, "-") for m in re.finditer("(?=GCCGGCCA)", G)]),
    columns=["pos", "strand"])

aag_a = sorted({m.start() + i for m in re.finditer("(?=AAGCCCG)", G) for i in (0, 1)}
               | {m.start() + 5 + i for m in re.finditer("(?=CGGGCTT)", G) for i in (0, 1)})
aag_a_set = set(aag_a)

fc = pd.read_csv(FEATURECOUNTS, sep="\t", comment="#")
bam_cols = [c for c in fc.columns if c.endswith(".bam")]
LIB_RE = re.compile(r"M145_(\d)_(\d)\.")
smap = {}
for c in bam_cols:
    m = LIB_RE.search(c)
    smap[c] = m.group(1) + "-" + m.group(2)
fc = fc.rename(columns=smap)
GENES = fc[["Geneid", "Start", "End", "Strand"]].copy()
GENES["Start"] = GENES.Start.astype(str).str.split(";").str[0].astype(int)
GENES["End"] = GENES.End.astype(str).str.split(";").str[-1].astype(int)

mb = pd.read_csv(DATA_DIR / "MB02_BGC_crosswalk_supplementary.csv")

# Canonical recompute (see module docstring).  Numerator/denominator are kept in
# the table so every plotted quantity can be read back as Nmod/Nvalid.
CAN = pd.read_csv(CANON, sep="\t")
assert set(CAN.panel) == {"b", "c", "d"}, sorted(set(CAN.panel))

# ---- promoter inventory (panel a) ---------------------------------------
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

# ---- cluster promoter vs rest of cluster (panel b) ----------------------
# Pooled TGGCCGGC 4mC level inside the 29 regions, split on the same promoter
# window panel a uses.  level = Nmod / Nvalid summed over measurable positions.
panel_b = (CAN[CAN.panel == "b"]
           .pivot(index="timepoint", columns="locus",
                  values=["level_pooled_pct", "n_positions", "Nmod", "Nvalid"]))
PROM_KEY, REST_KEY = "promoter (-300..-1)", "rest of cluster"
for tp in ("T1", "T2", "T3"):
    for k in (PROM_KEY, REST_KEY):
        num = panel_b[("Nmod", k)][tp]
        den = panel_b[("Nvalid", k)][tp]
        assert abs(100 * num / den - panel_b[("level_pooled_pct", k)][tp]) < 0.02, (tp, k)

# ---- cluster expression vs 4mC level (panel c) --------------------------
L4 = (CAN[CAN.panel == "c"]
      .merge(mb[["region", "expression_class", "coreTPM_T1", "coreTPM_T2", "coreTPM_T3"]],
             on="region")
      .rename(columns={"expression_class": "cls"}))
L4["tpm"] = [getattr(r, f"coreTPM_{r.timepoint}") for r in L4.itertuples()]
L4 = L4.dropna(subset=["level_pooled_pct"])
L4 = L4[L4.n_positions >= MIN_POS]
r4, p4 = spearmanr(np.log2(L4.tpm + 1), L4.level_pooled_pct)

# ---- cluster interior vs flank (panel d) --------------------------------
# One quantity per row: pooled methylation level of a motif's measurable
# positions, cluster interior divided by the 50 kb flanks.  The withdrawn table
# also carried an "unassigned 6mA" row; that class has no motif denominator, so
# its comparable quantity is a called-site density, not a level, and it is
# reported in the source note instead of being plotted on this axis.
DCLASSES = ["TGGCCGGC 4mC", "AAGCCCG C5 4mC", "AAGCCCG 6mA"]
dpan = CAN[(CAN.panel == "d") & CAN.mark_class.isin(DCLASSES)]
dpan = dpan.pivot(index=["mark_class", "timepoint"], columns="locus",
                  values=["level_pooled", "n_positions", "Nmod", "Nvalid"])
dpan = dpan.assign(ratio=dpan[("level_pooled", "interior")]
                   / dpan[("level_pooled", "flank")]).reset_index()
dpan.columns = ["_".join(c).rstrip("_") if isinstance(c, tuple) else c
                for c in dpan.columns]
dpan = dpan.sort_values(["mark_class", "timepoint"],
                        key=lambda s: s.map({m: i for i, m in enumerate(DCLASSES)})
                        if s.name == "mark_class" else s).reset_index(drop=True)

apply_figure_style(frame="open", sizes=(7.0, 6.2, 6.0))


def build():
    fig, axes = plt.subplots(2, 2, figsize=(NAR_WIDTH_IN, 4.6))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    cats = ["genes outside\nclusters", "genes inside\nclusters"]
    vals = [summ.loc[False, "pct"], summ.loc[True, "pct"]]
    cnts = [(int(summ.loc[False, "with_site"]), int(summ.loc[False, "n_promoters"])),
            (int(summ.loc[True, "with_site"]), int(summ.loc[True, "n_promoters"]))]
    ax_a.bar(cats, vals, width=0.5, color=C_MARK)
    for i, (v, (w, n)) in enumerate(zip(vals, cnts)):
        ax_a.text(i, v + 0.3, f"{w}/{n}\n({v:.1f} %)", ha="center", fontsize=6.0,
                  linespacing=1.4)
    ax_a.set_ylabel("Promoters carrying a\nTGGCCGGC position (%)")
    ax_a.set_ylim(0, 9)
    ax_a.set_title("",
                   loc="left", fontsize=6.6)

    # Panel b: zero baseline.  The withdrawn version started the axis at 60 %,
    # which made a 0.4 pp difference between the two bars look like a third of
    # the bar height (FIG: truncated baseline).
    xp = np.arange(3)
    w = 0.34
    # The measurable position count is the same at all three timepoints (19 in the
    # promoter window, 261 elsewhere in the clusters), so it is stated once in the
    # key rather than repeated over six bars.
    n_prom = {int(panel_b[("n_positions", PROM_KEY)][t]) for t in TPS}
    n_rest = {int(panel_b[("n_positions", REST_KEY)][t]) for t in TPS}
    assert len(n_prom) == 1 and len(n_rest) == 1, (n_prom, n_rest)
    ax_b.bar(xp - w / 2, [panel_b[("level_pooled_pct", PROM_KEY)][t] for t in TPS], w,
             color=C_MARK, label=f"promoter −300..−1 ({n_prom.pop()} positions)")
    ax_b.bar(xp + w / 2, [panel_b[("level_pooled_pct", REST_KEY)][t] for t in TPS], w,
             color=C_MARK, hatch="////", edgecolor="white", linewidth=0.0,
             label=f"rest of cluster ({n_rest.pop()} positions)")
    ax_b.set_xticks(xp)
    ax_b.set_xticklabels(["T1  12 h", "T2  24 h", "T3  50 h"])
    ax_b.set_ylabel("TGGCCGGC 4mC level (%)")
    ax_b.set_ylim(0, 100)
    ax_b.set_yticks([0, 25, 50, 75, 100])
    ax_b.set_title("", loc="left", fontsize=6.6)
    # Zero baseline leaves no headroom inside the panel, so the key sits in the
    # whitespace under the axes rather than on top of the bars.
    ax_b.legend(frameon=False, fontsize=6.0, loc="upper left", handlelength=1.4,
                bbox_to_anchor=(-0.02, -0.17), labelspacing=0.30, handleheight=1.0)

    # Panel c: expression class is ORDINAL, so it is encoded by a single-hue
    # lightness ramp AND by marker shape — never by a red/green pair, which the
    # withdrawn version used and which is indistinguishable in deuteranopia.
    CC = {"silent": (U.COL_BAR_LIGHT, "^", "none"),
          "low": (U.COL_BAR_MID, "s", U.COL_BAR_MID),
          "expressed": (U.COL_BAR_DARK, "o", U.COL_BAR_DARK)}
    for cl in ("silent", "low", "expressed"):
        g = L4[L4.cls == cl]
        edge, mark, face = CC[cl]
        ax_c.plot(np.log2(g.tpm + 1), g.level_pooled_pct, ls="none", marker=mark,
                  ms=3.6, mfc=face, mec=edge, mew=0.9, alpha=0.9, label=cl)
    ax_c.set_xlabel("Core-gene expression, log$_2$(TPM + 1)")
    ax_c.set_ylabel("TGGCCGGC 4mC level (%)")
    ax_c.set_title("", loc="left", fontsize=6.6)
    ax_c.legend(frameon=False, fontsize=6.0, loc="lower right", handlelength=0.9,
                borderaxespad=0.15, labelspacing=0.28)
    # Pad past the extreme marks so no glyph is cut by a spine.
    lo, hi = L4.level_pooled_pct.min(), L4.level_pooled_pct.max()
    ax_c.set_ylim(lo - 0.06 * (hi - lo), hi + 0.06 * (hi - lo))
    ax_c.margins(x=0.06)

    # Panel d: one quantity, three motif classes x three timepoints.  The row
    # label names the EIGHT-base recognition sequence TGGCCGGC; the withdrawn
    # version labelled this quantity "GCCGGC", the six-base core.
    MARKS = {"TGGCCGGC 4mC": ("o", U.COL_4mC), "AAGCCCG C5 4mC": ("s", U.COL_4mC),
             "AAGCCCG 6mA": ("D", U.COL_6mA)}
    for i, t in enumerate(dpan.itertuples()):
        mk, col = MARKS[t.mark_class]
        ax_d.plot(t.ratio, i, ls="none", marker=mk, ms=3.6, color=col, zorder=3)
    ax_d.axvline(1.0, color=META_GREY, lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax_d.set_yticks(range(len(dpan)))
    ax_d.set_yticklabels([f"{t.mark_class}  {t.timepoint}" for t in dpan.itertuples()],
                         fontsize=6.0)
    ax_d.set_ylim(-0.6, len(dpan) - 0.4)
    ax_d.set_xlabel("Methylation level, cluster interior / flank (ratio)")
    rlo, rhi = dpan.ratio.min(), dpan.ratio.max()
    ax_d.set_xlim(min(rlo, 1.0) - 0.012, max(rhi, 1.0) + 0.012)
    ax_d.set_xticks([0.95, 0.98, 1.00, 1.02])
    ax_d.set_title("", loc="left", fontsize=6.6)

    fig.subplots_adjust(left=0.105, right=0.985, top=0.925, bottom=0.300,
                        hspace=0.80, wspace=0.62)
    figure_panel_letters(fig, [ax_a, ax_b, ax_c, ax_d], "abcd",
                         dx=[0.090, 0.090, 0.090, 0.090], dy=0.026)
    source_note(fig,
                "Promoter window −300..−1 bp from an experimental primary or secondary TSS (Jeong et al. 2016); panel a uses the 3,030 of "
                "3,570 rows that are\ncategory P or S with a gene assignment. Panels b–d pool the three libraries of a timepoint: a position is "
                "measurable at pooled true depth\n(Nvalid + Nfail) ≥ 30 with the expected reference base, and the level is the pooled Nmod/Nvalid. "
                "29 antiSMASH regions; flank = 50 kb either\nside of a region. Panel c plots the "
                f"{len(L4)} region–timepoint pairs with ≥ {MIN_POS} measurable positions (Spearman ρ = {r4:+.2f}".replace("-", "\u2212") + f", p = {p4:.2f}). "
                "Unassigned 6mA has no\nmotif denominator, so it is not on panel d's axis; its called-site density, interior over flank, "
                "is 0.99, 0.82 and 0.91 at T1–T3.",
                y=0.010)
    return fig


if __name__ == "__main__":
    # 出荷名は SuppFigure8.png（補足図 14 -> 8 の改番、2026-10-24 に出力名を追随）。
    # SuppFigure14.png は現在「アデニン側の決定因子」図
    # （suppfigure14_determinants_261020.py）が占めている。
    report(save(build(), "SuppFigure8.png"))
