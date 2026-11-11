# 2026-11-05 D188: 補足図 12→9 の削減。出荷名 SuppFigure8.png -> SuppFigure7.png。
# -*- coding: utf-8 -*-
"""Supplementary Figure 17 — candidate adenine MTase transcript levels against 6mA occupancy.

生成元が figure_manifest.tsv に無い孤児図だったため、2026-10-15 に生成元を復元して
スクリプト化した。復元は artifact 92979b18 の lineage から行ったが、復元コードの
focal_palette は灰色系を返し現物の青／黄褐／赤と一致しなかったため、配色は
00_shared_utils の定数に置き換えてある（FIG の「パレット定数を import する」規則にも適合）。

2026-10-15 の変更点:
  - 灰色破線の系列は「6mA occupancy at AAGCCCG」と注記されながら
    HAAGCCCG 565 分母の値 (45.84, 47.79, 46.02) を描いていた。注記と不一致。
    さらに T3 の 46.02 は正準値 44.96 ですらない旧値だった。
  - 著者決定 DECIDED.6mA.DENOMINATOR.261015 に従い AAGCCCG 1,334 分母の
    (24.21, 26.24, 25.11) に差し替え、注記が実際の量と一致するようにした。

TPM は featureCounts_M145_s2.txt から length 正規化した値（平均カウントではない）。
幾何（figsize / bbox_inches）は現物 1170x768 px を再現するため据え置き。

Run: PYTHONPATH=. python3 suppfigure17_mtase_expression.py
"""
import importlib
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')

OUT_DIR = os.environ.get(
    'SUPPFIG_OUT', '/Users/okaban/obsidian/Research/rna-seq/Writing/fig_images')

# axes.titlesize は既定だと font.size x 1.2 = 9.6 pt になり、表題がどの現実的な幅でも
# キャンバスをはみ出す（実測: 幅 3.90-4.20 in のいずれでも x1 > 1.0）。現物が収まって
# いたのは bbox_inches="tight" が後からキャンバスを広げていたためで、幾何を固定する
# 以上は明示的に下げる必要がある。
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 7, 'axes.titlesize': 8.5,
                     'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     # 2026-10-28 D121: 補足図も四辺の枠に揃える。ただしこの図は
                     # 第 2 軸 (ax2) を持ち、右スパインは C_OCC で着色して
                     # 第 2 軸の尺度を符号化している（L133）。右を枠線にすると
                     # その符号化が壊れるので、右だけは枠に含めず着色のまま残す。
                     'axes.spines.top': True, 'axes.spines.right': False})

# 候補酵素の転写量（TPM）。featureCounts から実行時に算出する。
# 2026-10-15 まで TPM はこのファイルにべた書きされており、再実行しても打ち込んだ値が
# 再現するだけでデータとの一致を検証できなかった。実測と小数 2 桁まで一致することは
# 確認済みだが、検証可能性のために読み込みに変更した。
FC = ('/Users/okaban/bioinfo/rna-seq/03_quantification/data/featureCounts/'
      'featureCounts_M145_s2.txt')
TP_LIBS = {'T1': ['1-1', '1-2', '1-3'], 'T2': ['2-1', '2-3', '2-4'],
           'T3': ['3-2', '3-3', '3-4']}


def _tpm_table():
    import re
    import pandas as pd
    fc = pd.read_csv(FC, sep='\t', comment='#').set_index('Geneid')
    ren = {c: re.search(r'M145_(\d_\d)', c).group(1).replace('_', '-')
           for c in fc.columns if 'M145_' in c}
    fc = fc.rename(columns=ren)
    libs = list(ren.values())
    assert len(libs) == 9, f'ライブラリ {len(libs)} 本（9 本を期待）'
    rpk = fc[libs].div(fc['Length'], axis=0)
    return rpk.div(rpk.sum(), axis=1) * 1e6, libs


TPM_TBL, LIBS = _tpm_table()


def _series(gene):
    return tuple(float(TPM_TBL.loc[gene, TP_LIBS[t]].mean()) for t in ('T1', 'T2', 'T3'))


# 生存候補 3 件（図 4c の漏斗を通過）
# 2026-11-09 (H04): NPPY は原稿にも出荷する補足表にも根拠が無く、本文は DPPY/HPPF で
# 記述している。裏の取れない注記を図から外し、HPPF だけ残す（本文 Discussion に記載あり）。
NAMES = ["SC_RS35335", "SC_RS28835", "SC_RS17645 (HPPF)"]
GENES = ["SC_RS35335", "SC_RS28835", "SC_RS17645"]
TPM = [_series(g) for g in GENES]
# 2026-10-22 — 色の帰属を直す。この図は候補遺伝子 3 本に COL_6mA（青）と COL_4mC（赤）
# を割り当て、実際の 6mA 占有率の系列を灰色（COL_BAR_MID）で描いていた。他の全図では
# 青 = 6mA、赤 = 4mC なので、ここだけ「赤 = SC_RS17645 という遺伝子」になっていて、
# 図をまたいで同じ実体を追えない。panels_261013/panel_4e.py が同じ構図で既に採って
# いる規約（転写量は中立色、占有率はマーク色）に合わせる:
#     6mA 占有率 -> COL_6mA（青。この図で唯一のメチル化マーク量）
#     候補遺伝子 -> 中立色 3 色（遺伝子名は図中に直接書いてあるので色は識別に不要）
COLS = [U.COL_BAR_DARK, U.COL_GRAY, U.COL_BAR_MID]  # FIG23: 候補遺伝子 3 本は中立グレー段
C_OCC = U.COL_6mA

# 348 aa のサイズ床だけで除外された 2 件。化学的な根拠ではなく長さで落ちているため、
# 除外が妥当かどうかは転写されているかで判断できる（図 4c の漏斗の最終段に対応）。
EXCL_NAMES = ["SC_RS36625 (DPPY, 251 aa)", "SC_RS19670 (DPPY, 248 aa)"]
EXCL_GENES = ["SC_RS36625", "SC_RS19670"]
EXCL_TPM = [_series(g) for g in EXCL_GENES]

# 全遺伝子 TPM 中央値（9 ライブラリの平均）。発現しているかの基準線。
MEDIAN_TPM = float(TPM_TBL[LIBS].median().mean())

# OCC.6mA.AAGCCCG.T1/T2/T3 — 分母は AAGCCCG 1,334 インスタンス（位置 1）
OCC = (24.21, 26.24, 25.11)

x = np.arange(3)
# 2026-11-07 補足図点検 1: 二軸（左 TPM / 右 6mA 占有率）なのに、どの系列が
# どちらの軸かを示す凡例が無く、実線・破線の規約も図中に書かれていなかった。
# さらに系列名が線の右端に直接置かれて破線に重なり、右軸は目盛だけで軸線が無かった。
# 二軸をやめ、上下 2 段（a 転写量 / b 6mA 占有率）に分ける。x は共有するので
# 時点どうしの対応は縦に読め、凡例で軸を説明する必要も無くなる。
fig = plt.figure(figsize=(4.55, 3.20), dpi=300)
ax = fig.add_axes([0.135, 0.500, 0.580, 0.430])     # (a) transcript level
ax2 = fig.add_axes([0.135, 0.165, 0.580, 0.275])    # (b) 6mA occupancy

# 2026-11-07: 全遺伝子中央値の基準線は、系列線やラベルと重なるわりに読み取れる
# ことが少ないので図から外した。値（MEDIAN_TPM）は凡例で述べる。
for name, v, c in zip(NAMES, TPM, COLS):
    ax.plot(x, v, "o-", color=c, lw=1.3, ms=4.2, zorder=3)
    ax.annotate(name, (2, v[2]), xytext=(4, 0), textcoords="offset points",
                color=c, fontsize=6, va="center")
# 348 aa の床だけで落ちた 2 件は経過がほぼ重なるので、まとめて 1 つの注記にする。
# 2026-11-07 著者指示（作図規約）: 色だけで区別が付く場合は線種を変えない。
# 除外 2 件は COL_BAR_LIGHT、候補 3 件は DARK / GRAY / MID で明度が分かれるので、
# 破線をやめて全系列を実線にする。
for v in EXCL_TPM:
    ax.plot(x, v, "o-", color=U.COL_BAR_LIGHT, lw=0.9, ms=3.4,
            mfc=U.COL_BAR_LIGHT, mec="white", mew=0.5, zorder=2)
ax.annotate("excluded by the 348 aa floor\n(SC_RS36625, SC_RS19670)",
            (2.05, 16.5), xytext=(4, 0), textcoords="offset points",
            color=U.COL_BAR_LIGHT, fontsize=5.6, va="center", linespacing=1.4)
ax.set_xticks(x); ax.set_xticklabels([]); ax.set_xlim(-0.3, 4.6)
ax.set_ylabel("Transcript level (TPM)"); ax.set_ylim(0, 46)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

ax2.plot(x, OCC, "s-", color=C_OCC, lw=1.1, ms=3.6, zorder=2)
ax2.set_xticks(x); ax2.set_xticklabels(["T1", "T2", "T3"]); ax2.set_xlim(-0.3, 4.6)
ax2.set_ylabel("6mA occupancy (%)", color=C_OCC)
ax2.set_ylim(0, 40); ax2.set_yticks([0, 20, 40])
ax2.tick_params(axis="y", colors=C_OCC)
ax2.spines["left"].set_color(C_OCC)
for s in ("top", "right"):
    ax2.spines[s].set_visible(False)
ax2.annotate("AAGCCCG, n = 1,334 instances", (2, OCC[2]), xytext=(4, 0),
             textcoords="offset points", color=C_OCC, fontsize=6, va="center")

for a_, lab in ((ax, "a"), (ax2, "b")):
    box = a_.get_position()
    fig.text(max(0.004, box.x0 - 0.115), min(0.950, box.y1 + 0.006), lab, fontsize=8.5,
             fontweight="bold", ha="left", va="bottom", color="black")

U.assert_no_text_collisions(fig, "SuppFigure7")
out = os.path.join(OUT_DIR, "SuppFigure10.png")
fig.savefig(out, dpi=300)
print("saved:", out, plt.imread(out).shape, "median TPM", round(MEDIAN_TPM, 2))
