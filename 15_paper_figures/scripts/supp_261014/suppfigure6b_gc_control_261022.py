# 2026-11-05 D188: 補足図 12→9 の削減。出荷名 SuppFigure4.png -> SUPERSEDED_SuppFigure4_removed_D188.png。
# 2026-11-04 著者決定: 区別は軸ラベル・パネル内表記・形が担っているので、色で重ねて符号化しない（中立グレー）。
#!/usr/bin/env python3
"""補足図 6 — TSS 近傍のモチーフ枯渇に対する GC 含量の統制。

2026-10-22: この図は出荷 PNG のみが存在し figure_manifest.tsv の script 欄が MISSING
だった。さらに画像内に「GC content explains most of the apparent depletion — but not
all of GCCGGC」という断定的な全体タイトルが焼き込まれており FIG01（図画像に図全体の
タイトルを焼き込まない）と FIG17（解釈は凡例に置く）に違反していた。

生成元を復元するにあたり、描く量は recompute_261015/recompute_suppfigure12.py の
出力（suppfigure12_recomputed.tsv）から読む。同スクリプトは旧 PNG が主張していた
量を実際のゲノムと Jeong 2016 TSS 表から再計算している。2026-10-22 に照合した結果:

  一致（本スクリプト冒頭の assert で毎回確認する）
      観測数        2,539 / 305 / 142,109          旧図と完全一致
      GC 整合期待数 2,883.8 / 322.0 / 146,618.7    旧図と完全一致
      比            0.8804 / 0.9471 / 0.9692       小数 4 桁まで一致
      GC 非整合比   0.7525 / 1.0898 / 0.8675       小数 4 桁まで一致
      TSS 窓数 3,030、窓内 GC 69.54 % 対 全体 72.12 %  旧図と完全一致

  一致しない（assert していない。再計算値を描く）
      95% CI の端点が旧図の表示と第 3-4 位で最大 0.0014 ずれる
      （六塩基 0.8465-0.9154 対 旧 0.8478-0.9160、AAGCCCG 0.8438-1.0596 対
        0.8424-1.0601、三塩基 0.9642-0.9743 対 0.9646-0.9739）。旧図の CI が
      どの近似で引かれたかは生成元が無いため確認できない。本図は
      recompute_suppfigure12.py の Byar ポアソン近似をそのまま描き、
      凡例にも近似名を書く。点推定が動いていないので結論は変わらない。

したがって本図は旧図と同じ点推定を、タイトルなしで、CI は再計算値で描き直したもの。

配色は 00_shared_utils の正準パレットに合わせる（GCCGGC = 4mC の赤、AAGCCCG = 青、
三塩基の対照は灰）。GC 統制の有無は色ではなく記号で区別する（FIG06）。

Run: PYTHONPATH=. python3 suppfigure6b_gc_control_261022.py
"""
import importlib
import textwrap
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from _suppstyle import (META_GREY, NAR_WIDTH_IN, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
U = importlib.import_module("00_shared_utils")

TSV = pathlib.Path("/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/"
                   "recompute_261015/suppfigure12_recomputed.tsv")

# 表示名（2 行）と正準色。GCCGGC は六塩基コアであって 4mC 標的の 8 塩基 TGGCCGGC
# ではない。ここで数えているのは配列 availability の分母なので六塩基で正しい。
ROWS = [
    ("GCCGGC hexamer", "GCCGGC\n(hexamer)", "#2b2b2b"),
    ("AAGCCCG both strands", "AAGCCCG\n(both strands)", "#6b6b6b"),
    ("GGC or GCC trinucleotide", "GGC or GCC\n(trinucleotide)", U.COL_GRAY),
]

d = pd.read_csv(TSV, sep="\t").set_index("quantity")
# 再計算値が旧図の主張と一致していることを描く前に確かめる（一致しなければ落とす）。
for key, _, _ in ROWS:
    r = d.loc[key]
    assert abs(r.ratio - r.claimed_ratio) < 5e-4, (key, r.ratio, r.claimed_ratio)
    assert abs(r.ratio_unmatched - r.claimed_unmatched) < 5e-4, key
    assert int(r.observed) == int(r.claimed_observed), key

apply_figure_style(frame="boxed", sizes=(7.5, 6.5, 6.5))


def build():
    fig, ax = plt.subplots(figsize=(NAR_WIDTH_IN, 2.75))
    ys = np.arange(len(ROWS))[::-1]

    for y, (key, _, col) in zip(ys, ROWS):
        r = d.loc[key]
        # GC 統制あり: 95% CI 付きの塗り菱形。統制なし: 中抜き円。
        ax.plot([r.ci95_low, r.ci95_high], [y, y], color=col, lw=1.6,
                solid_capstyle="butt", zorder=3)
        ax.plot([r.ratio_unmatched], [y + 0.26], marker="o", ms=5.0, mfc="white",
                mec=col, mew=1.1, ls="none", zorder=4)
        ax.plot([r.ratio], [y], marker="D", ms=5.4, color=col, ls="none", zorder=5)
        ax.annotate("", xy=(r.ratio + 0.004, y + 0.055), xytext=(r.ratio_unmatched, y + 0.225),
                    arrowprops=dict(arrowstyle="-|>", color=META_GREY, lw=0.6,
                                    shrinkA=3.0, shrinkB=3.0), zorder=2)

    ax.axvline(1.0, color=META_GREY, lw=0.8, ls=(0, (4, 2)), zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for _, lab, _ in ROWS], linespacing=1.3)
    ax.set_ylim(-0.55, len(ROWS) - 0.45)
    ax.set_xlim(0.70, 1.27)
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
    ax.set_xlabel("Observed / expected motif count in TSS \u00b1300 bp (ratio)")

    ax.legend(handles=[
        Line2D([0], [0], marker="o", ls="none", mfc="white", mec=META_GREY, mew=1.1,
               ms=5.0, label="no GC control"),
        Line2D([0], [0], marker="D", ls="none", color=META_GREY, ms=5.4,
               label="GC-matched (95% CI)")],
        frameon=False, loc="lower right", handletextpad=0.5, borderaxespad=0.3,
        labelspacing=0.45)

    fig.subplots_adjust(left=0.165, right=0.985, top=0.960, bottom=0.400)
    gcw = d.loc["GC content, TSS +/-300 bp (%)"]
    source_note(fig, textwrap.fill(
        f"Expected counts from 601 bp windows stratified into GC quintiles; intervals are the "
        f"Byar Poisson approximation on the observed count. "
        f"{int(d.loc['TSS windows used'].observed):,} experimental primary or secondary TSS "
        f"(Jeong et al. 2016); GC content in those windows is {gcw.observed:.2f} % against "
        f"{gcw.expected_gc_matched:.2f} % genome-wide. GCCGGC here is the six-base core, used as "
        f"the sequence-availability denominator; the 4mC recognition sequence is the eight-base "
        f"TGGCCGGC. Recomputed by recompute_261015/recompute_suppfigure12.py.", 128),
        y=0.030)
    return fig


if __name__ == "__main__":
    report(save(build(), "SUPERSEDED_SuppFigure4_removed_D188.png"))
