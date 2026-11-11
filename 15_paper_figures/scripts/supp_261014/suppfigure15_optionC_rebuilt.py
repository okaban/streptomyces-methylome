# 2026-11-05 D188: 補足図 12→9 の削減。出荷名 SuppFigure6.png -> SuppFigure5.png。
#!/usr/bin/env python3
"""Supplementary Figure 9 — rebuilt to lead with the molecule-level distribution.

Renumbered 15 -> 9; the output filename was corrected on 2026-10-22. Until then this
script still wrote SuppFigure15.png, which Supplementary Figure 15 (the duplex
schematic relocated from Figure 3a) now occupies — re-running it would have
overwritten that figure.

Adopted 2026-10-15 (author decision) in place of the earlier three-panel version in
suppfigure15_single_molecule.py, which is retained for reference but no longer ships.

Panel a asks the question in the reader's own terms: if the population splits into
fully modified and fully unmodified molecules (epialleles), the number of modified
sites per molecule piles up at the two ends; if methylation is placed independently
it follows a binomial.  Panels b-d quantify the residual coupling and bound it with
two controls.

All annotation is kept out of the data area: series are named in frameless legends
or on the category axis, never as text laid over bars or lines.
    python suppfigure15_optionC_rebuilt.py
"""
import json
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from _suppstyle import (BLUE, RED,  # noqa: E402
                        DATA_DIR, META_GREY, NAR_WIDTH_IN, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

H = pd.read_csv(DATA_DIR / "h26" / "H26.csv")
PCd = pd.read_csv(DATA_DIR / "h26" / "PC26.csv")
D26 = pd.read_csv(DATA_DIR / "h26" / "D26.csv")
PM = pd.read_csv(DATA_DIR / "epi261015" / "per_molecule_counts_HAAGCCCG_n4_T1_261015.csv")
EPI = pd.read_csv(DATA_DIR / "epi261015" / "epi_excess_HAAGCCCG_261015.csv")
NUM = json.load(open(DATA_DIR / "epi261015" / "epi_numbers_n4.json"))

six = H[H["mark"].str.contains("6mA")]
BIST = six.bistable_pred.mean()
WITHIN = PCd.excess.mean()
BETWEEN = six.excess_conc.mean()

C_OBS, C_STO, C_EPI = "#2b2b2b", "#6b6b6b", "#a8a8a8"  # FIG23: モデル区分は中立グレー
C_BTW = "#d4d4d4"
# データキー（CSV の mark 列）と表示ラベルは別物。キーを呼称規則に合わせて
# 書き換えると照合が空になるので、表示だけ「hexamer」を付ける。
# 2026-10-22: was 4mC navy / 6mA orange -- the canonical pair with the hues
# swapped. Bound to the one palette decision on record (00_shared_utils).
MC = {"GCCGGC 4mC": RED, "AAGCCCG 6mA": BLUE}
# 2026-10-22: 表示名を 8 塩基に直す。H26.csv の occupancy はライブラリあたり 0.69 前後で、
# これは TGGCCGGC 8 塩基の占有率（プール 90 %）の水準。六塩基 GCCGGC のシトシン全体では
# 1.8 % にしかならないので、この量は 8 塩基のものである。キーは照合用なので据え置く。
DISPLAY = {"GCCGGC 4mC": "TGGCCGGC 4mC", "AAGCCCG 6mA": "AAGCCCG 6mA"}

apply_figure_style(frame="boxed", sizes=(7.0, 6.2, 6.0))


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 2.95))  # 2026-11-09: 図下の散文を削り 0.35 in 詰めた
    gs = fig.add_gridspec(1, 4, width_ratios=[1.46, 1.20, 1.00, 1.10], wspace=0.95,
                          left=0.062, right=0.985, top=0.821, bottom=0.178)

    # (a) 1 分子が読めた 4 箇所のうち何箇所が 6mA か — 観測・独立配置・エピアレル
    ax_a = fig.add_subplot(gs[0, 0])
    tot = PM["observed"].sum()
    xs = PM["n_methylated"].values
    ax_a.bar(xs, 100 * PM.epiallele_prediction / tot, width=0.88, color="none",
             edgecolor=C_EPI, lw=0.9, hatch="////", zorder=2)
    ax_a.bar(xs, 100 * PM.observed / tot, width=0.56, color=C_OBS, zorder=3)
    ax_a.plot(xs, 100 * PM.stochastic_binomial / tot, color=C_STO, lw=1.3,
              marker="o", ms=3.2, zorder=4)
    ax_a.set_xticks(xs)
    ax_a.set_xlabel("Modified sites on one molecule\n(of the four AAGCCCG sites it covers)")
    ax_a.set_ylabel("Molecules (%)")
    # 2026-11-09 著者コメント [12]: 凡例が柱に重なっていた。上に凡例帯を確保する。
    ax_a.set_ylim(0, 88)
    ax_a.set_xlim(-0.62, 4.62)
    ax_a.set_title("", loc="left")
    ax_a.legend(handles=[Patch(facecolor=C_OBS, edgecolor="none", label="observed"),
                         Line2D([0], [0], color=C_STO, lw=1.3, marker="o", ms=3.0,
                                label="independent placement"),
                         Patch(facecolor="none", edgecolor=C_EPI, lw=0.9, hatch="////",
                               label="epiallele prediction")],
                frameon=False, fontsize=6.0, loc="upper left", handlelength=1.3,
                borderpad=0.0, labelspacing=0.26, handletextpad=0.5,
                bbox_to_anchor=(-0.02, 1.02))

    # (b) 残る連関を 2 つの参照点に対して置く — 水平棒にして軸側にラベル
    ax_b = fig.add_subplot(gs[0, 1])
    rows = [("between\nmolecules", BETWEEN, C_BTW),
            ("within one\ninstance", WITHIN, C_OBS),
            ("epiallele\nprediction", BIST, C_EPI)]
    for i, (lab, v, col) in enumerate(rows):
        ax_b.barh(i, v, height=0.52, color=col, zorder=3)
        ax_b.text(v + 0.014, i, f"{v:.3f}", va="center", ha="left", fontsize=6.2,
                  color="#4a4a4a", fontweight="bold")
    ax_b.axvline(0, color=META_GREY, lw=0.8)
    ax_b.set_yticks(range(len(rows)))
    ax_b.set_yticklabels([r[0] for r in rows], fontsize=6.0)
    ax_b.set_ylim(-0.60, len(rows) - 0.40)
    ax_b.set_xlim(0, 0.70)
    ax_b.set_xticks([0, 0.2, 0.4, 0.6])
    ax_b.set_xlabel("Excess concordance between\ntwo sites on one molecule")
    ax_b.set_title("", loc="left")

    # (c) 距離対照 — processive でない
    ax_c = fig.add_subplot(gs[0, 2])
    for mk in ("GCCGGC 4mC", "AAGCCCG 6mA"):
        dm = D26[(D26["mark"] == mk) & (D26.dist_max <= 15000)]
        g = dm.groupby("dist_max")["excess"]
        ax_c.errorbar(g.mean().index, g.mean().values, yerr=g.std().values,
                      color=MC[mk], marker="o", ms=3.0, lw=1.1, capsize=1.8,
                      elinewidth=0.8, label=DISPLAY[mk])
    ax_c.axhline(0, color=META_GREY, lw=0.8)
    ax_c.set_xscale("log")
    # 凡例を置く余白を上端に作る。目盛りは 0.04 までで打ち切り、
    # 0.04 より上はデータも目盛りラベルも無い帯にする。
    ax_c.set_ylim(0, 0.072)
    ax_c.set_yticks([0, 0.01, 0.02, 0.03, 0.04])
    ax_c.set_xlabel("Distance between the two sites (bp)\n200 bp to 15 kb")
    ax_c.set_ylabel("Excess concordance")
    ax_c.set_title("", loc="left")
    ax_c.legend(frameon=False, fontsize=6.0, loc="upper left", handlelength=0.9,
                borderpad=0.0, labelspacing=0.25, handletextpad=0.4)

    # (d) 配列一致対照 — 群名と中央値は凡例に置き、データ域にはテキストを置かない
    ax_d = fig.add_subplot(gs[0, 3])
    rng = np.random.default_rng(7)
    handles = []
    for i, (key, col, lab) in enumerate((("non-motif control", C_EPI, "control"),
                                         ("AAGCCCG instances (H context)", C_OBS, "5\u2032 not G"))):
        d = EPI[EPI.set == key]["excess"]
        ax_d.scatter(i + rng.uniform(-0.16, 0.16, len(d)), d, s=2.2, alpha=0.30,
                     color=col, linewidths=0)
        m = d.median()
        ax_d.hlines(m, i - 0.27, i + 0.27, color=col, lw=1.6, zorder=3)
        handles.append(Line2D([0], [0], color=col, lw=1.6, label=f"{lab}  {m:.3f}"))
    ax_d.set_xticks([])
    ax_d.set_xlim(-0.62, 1.62)
    ax_d.set_ylim(-0.03, 0.72)
    ax_d.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax_d.set_ylabel("Excess within-molecule structure\n(\u0394H above permutation null)")
    ax_d.set_title("", loc="left")
    ax_d.legend(handles=handles, frameon=False, fontsize=6.0, loc="upper center",
                handlelength=0.9, borderpad=0.0, labelspacing=0.30, handletextpad=0.4)

    figure_panel_letters(fig, [ax_a, ax_b, ax_c, ax_d], "abcd",
                         dx=[0.055, 0.068, 0.062, 0.070], dy=0.022)
    return fig


if __name__ == "__main__":
    report(save(build(), "SuppFigure5.png"))
