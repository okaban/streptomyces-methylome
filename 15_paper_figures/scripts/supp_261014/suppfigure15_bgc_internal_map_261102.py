# 2026-11-04 著者決定: 区別は軸ラベル・パネル内表記・形が担っているので、色で重ねて符号化しない（中立グレー）。
#!/usr/bin/env python3
"""Supplementary Figure 15 — what sets the number of marks inside a predicted cluster.

Co-author comment C07 (Asamizu, 2026-10-30): "a figure putting the methylation
positions and the transcription mapping on the BGC map is important".

A first version of this figure drew the modified positions as ticks on 29
length-normalised region rows. The author rejected it on 2026-11-02 for the right
reason: a reader cannot read uniformity off 382 ticks spread over 29 rows, and the
only legible content (per-region counts, expression class) duplicated the table.
This version plots the quantity a reader can actually judge.

Per antiSMASH region: the number of modified positions against the number of
recognition-sequence positions the region contains, one panel per mark, with the
core-gene expression class as the marker. Two things are legible at a glance:

  - the count is set by the region's recognition-sequence content. 4mC lies on a
    line through the origin (Spearman rho = 0.982, Pearson r = 0.990), with slope
    0.876 modified positions per recognition-sequence position - the genome-wide
    occupancy, reproduced region by region;
  - the expression classes are interleaved along that line. Residuals from the
    proportional fit do not differ by class (Kruskal-Wallis p = 0.165 for 4mC,
    0.709 for 6mA), so a cluster being transcribed does not move it off the line.

The 6mA panel is the same relation at a tenth of the counts (slope 0.160, range
0-7 per region), so its scatter is Poisson noise and is shown as the honest limit
of what 29 regions can resolve rather than as a weaker effect.

Region GC does not predict the count once the recognition sequences are counted
(Spearman rho = -0.114, p = 0.555 for 4mC): GC acts through motif density, and
conditioning on motif density leaves nothing.

Site rule: three libraries pooled per timepoint, pooled true depth
(Nvalid + Nfail) >= 30, pooled Nmod/Nvalid >= 0.5. T1 is drawn; the genome-wide
totals are asserted against the canonical 2,909 4mC and 416 6mA.

Environment: epi-rna.
"""
import importlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/newfig_261101")
U = importlib.import_module("00_shared_utils")
U.apply_unified_style()
matplotlib.rcParams["savefig.bbox"] = None
nf3 = importlib.import_module("nf_r3_recompute_261101")

BGC = ("/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014/"
       "data/bgc_regions_29.csv")
# 2026-11-04: 補足図 15→12 の再付番に追従。この生成器が書くのは補足図 12。
OUT = "/Users/okaban/bioinfo/rna-seq/15_paper_figures/figures_261013/SuppFigure12.png"
TAB = ("/Users/okaban/bioinfo/rna-seq/15_paper_figures/tables/supplementary/"
       "SuppTableS5d_bgc_internal_positions_261102.tsv")
TP = "T1"
CANON = {"4mC": 2909, "6mA": 416}


def called(genome, pool, motif, offsets, mark):
    pos0, plus = nf3.motif_positions(genome, motif, offsets)
    d = pool[mark].reindex(
        pd.MultiIndex.from_arrays([pos0, plus], names=["start", "plus"]))
    keep = (d["dep"].to_numpy() >= 30) & (d["lvl"].to_numpy() >= 0.5)
    return pos0[np.nan_to_num(keep, nan=False).astype(bool)]


def main():
    genome = nf3.M.read_reference(nf3.M.REF)
    libs = nf3.TIMEPOINTS[TP]
    pool = {m: nf3.pooled(libs, nf3.CODE[m]) for m in ("4mC", "6mA")}
    s4 = np.sort(np.concatenate([
        called(genome, pool, "TGGCCGGC", [4], "4mC"),
        called(genome, pool, "AAGCCCG", [4], "4mC")]))
    s6 = np.sort(called(genome, pool, "AAGCCCG", [0, 1], "6mA"))
    assert len(s4) == CANON["4mC"], len(s4)
    assert len(s6) == CANON["6mA"], len(s6)
    sites = {"4mC": s4, "6mA": s6}

    b = pd.read_csv(BGC).sort_values("start").reset_index(drop=True)
    occ = {"4mC": np.concatenate([nf3.motif_positions(genome, "TGGCCGGC", [4])[0],
                                  nf3.motif_positions(genome, "AAGCCCG", [4])[0]]),
           "6mA": nf3.motif_positions(genome, "AAGCCCG", [0, 1])[0]}
    rows = []
    for r in b.itertuples():
        seg = genome[r.start - 1:r.end]
        rec = dict(region_number=r.region_number, known_cluster=r.known_cluster,
                   status=r.status, length_kb=r.length_kb, timepoint=TP,
                   gc_percent=round(100 * (seg.count("G") + seg.count("C")) / len(seg), 2))
        for mark in ("4mC", "6mA"):
            o = occ[mark]
            s_ = sites[mark]
            inr = s_[(s_ >= r.start) & (s_ <= r.end)]
            rec[f"occurrences_{mark}"] = int(((o >= r.start) & (o <= r.end)).sum())
            rec[f"modified_{mark}"] = len(inr)
            rec[f"frac_positions_{mark}"] = ";".join(
                f"{(x - r.start) / (r.end - r.start):.4f}" for x in inr)
        rows.append(rec)
    cm = pd.DataFrame(rows)
    cm.to_csv(TAB, sep="\t", index=False)

    stats_out = {}
    for mark in ("4mC", "6mA"):
        o = cm[f"occurrences_{mark}"].to_numpy(float)
        n = cm[f"modified_{mark}"].to_numpy(float)
        rate = n.sum() / o.sum()
        resid = n - rate * o
        sp = stats.spearmanr(o, n)
        pe = stats.pearsonr(o, n)
        kw = stats.kruskal(*[resid[(cm.status == st).to_numpy()]
                             for st in ("expressed", "low", "silent")])
        gcr = stats.spearmanr(cm.gc_percent.to_numpy(float), n)
        stats_out[mark] = dict(rate=rate, rho=sp.statistic, rho_p=sp.pvalue,
                               r=pe.statistic, kw_p=kw.pvalue, gc_rho=gcr.statistic,
                               gc_p=gcr.pvalue)
        print(f"{mark}: slope={rate:.4f} rho={sp.statistic:.3f} (p={sp.pvalue:.3g}) "
              f"r={pe.statistic:.3f} residual-by-class KW p={kw.pvalue:.3g} "
              f"GC rho={gcr.statistic:.3f} (p={gcr.pvalue:.3g})")

    marker = {"expressed": ("o", U.COL_BAR_DARK), "low": ("s", U.COL_BAR_MID),
              "silent": ("^", U.COL_BAR_LIGHT)}
    fig = plt.figure(figsize=(6.85, 2.95), dpi=300)
    axes = [fig.add_axes([0.095, 0.215, 0.365, 0.700]),
            fig.add_axes([0.600, 0.215, 0.365, 0.700])]
    for ax, mark, col in zip(axes, ("4mC", "6mA"), ("#4a4a4a", "#4a4a4a")):
        o = cm[f"occurrences_{mark}"].to_numpy(float)
        n = cm[f"modified_{mark}"].to_numpy(float)
        st_ = stats_out[mark]
        xx = np.array([0.0, o.max() * 1.10])
        ax.plot(xx, st_["rate"] * xx, lw=0.9, ls=(0, (3, 1.8)), color=col, zorder=2)
        for name, (mk, gcol) in marker.items():
            m = (cm.status == name).to_numpy()
            ax.plot(o[m], n[m], mk, ms=3.6, mfc=gcol, mec="white", mew=0.5,
                    zorder=4, label=f"{name} ({int(m.sum())})")
        ax.text(0.04, 0.96,
                f"{mark}\n\u03c1 = {st_['rho']:.3f}\nresiduals by class\n"
                f"p = {st_['kw_p']:.2f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=6.3,
                color=col, linespacing=1.5)
        ax.set_xlim(0, o.max() * 1.10)
        ax.set_ylim(0, max(n.max(), st_["rate"] * o.max()) * 1.14)
        ax.set_xlabel("Recognition-sequence positions in the region (count)", fontsize=7.0)
        ax.set_ylabel(f"Modified {mark} positions (count)", fontsize=7.0)
        U.style_axes(ax, frame=True)
    axes[0].legend(loc="lower right", fontsize=6.0, frameon=False,
                   handletextpad=0.4, borderpad=0.1, labelspacing=0.25)
    U.assert_no_text_collisions(fig, "SuppFigure15")
    # 2026-11-04: 2 パネル図なので投稿規定どおり (a) (b) を振る。
    for ax, lab in zip(axes, ("a", "b")):
        bb = ax.get_position()
        fig.text(bb.x0 - 0.072, bb.y1 + 0.028, lab, fontsize=9,
                 fontweight="bold", ha="left", va="top")
    fig.savefig(OUT, dpi=300, bbox_inches=None, pad_inches=0.02)
    plt.close(fig)
    from PIL import Image
    print("wrote", OUT, Image.open(OUT).size)


if __name__ == "__main__":
    main()
