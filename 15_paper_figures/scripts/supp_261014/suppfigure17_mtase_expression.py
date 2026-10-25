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
                     'axes.spines.top': False, 'axes.spines.right': False})

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
NAMES = ["SC_RS35335 (NPPY)", "SC_RS28835 (NPPY)", "SC_RS17645 (HPPF)"]
GENES = ["SC_RS35335", "SC_RS28835", "SC_RS17645"]
TPM = [_series(g) for g in GENES]
# 2026-10-22 — 色の帰属を直す。この図は候補遺伝子 3 本に COL_6mA（青）と COL_4mC（赤）
# を割り当て、実際の 6mA 占有率の系列を灰色（COL_BAR_MID）で描いていた。他の全図では
# 青 = 6mA、赤 = 4mC なので、ここだけ「赤 = SC_RS17645 という遺伝子」になっていて、
# 図をまたいで同じ実体を追えない。panels_261013/panel_4e.py が同じ構図で既に採って
# いる規約（転写量は中立色、占有率はマーク色）に合わせる:
#     6mA 占有率 -> COL_6mA（青。この図で唯一のメチル化マーク量）
#     候補遺伝子 -> 中立色 3 色（遺伝子名は図中に直接書いてあるので色は識別に不要）
COLS = [U.COL_BAR_DARK, U.COL_ORANGE, U.COL_BAR_MID]
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
# 幾何は固定する。現物は bbox_inches="tight" で保存されており、出力幅が表題の長さに
# 依存していた（1170 px = 99 mm）。同じ幅を figsize で明示し tight は使わない。
fig, ax = plt.subplots(figsize=(4.55, 2.70), dpi=300)
ax2 = ax.twinx()
# 2026-10-22: この基準線は axhline で軸いっぱいに引かれていたため、x > 2 の
# ラベル帯（データは T3 = x=2 で終わり、その右は系列名を置くためだけの余白）を
# 横断し、"SC_RS28835 (NPPY)" の上を通って文字を読めなくしていた。SC_RS28835 の
# T3 は 9.1 TPM で全遺伝子中央値 9.1 TPM とほぼ同値なので、ラベルを上下に逃がす
# 余地がない。線の描画範囲をデータのある区間に限定して解決する。
X_DATA_END = 2.05
ax.plot([-0.3, X_DATA_END], [MEDIAN_TPM, MEDIAN_TPM], color=U.COL_GRAY, lw=0.7,
        ls=(0, (3, 2)), zorder=1, solid_capstyle="butt", clip_on=False)
# 注記は線の終端の真上に置く（x=1.0-2.05, y>9.1 の帯は空いている）。
ax.annotate(f"all-gene median {MEDIAN_TPM:.1f} TPM", (X_DATA_END, MEDIAN_TPM),
            xytext=(0, 3), textcoords="offset points", color=U.COL_GRAY,
            fontsize=6.0, ha="right", va="bottom")
for name, v, c in zip(NAMES, TPM, COLS):
    ax.plot(x, v, "o-", color=c, lw=1.3, ms=4.2, zorder=3)
    ax.annotate(name, (2, v[2]), xytext=(4, 0), textcoords="offset points",
                color=c, fontsize=6, va="center")
# 除外 2 件は経過がほぼ重なるため、個別注記ではなく 1 つの凡例にまとめる。
for v in EXCL_TPM:
    ax.plot(x, v, "o--", color=U.COL_GRAY, lw=0.9, ms=3.4, mfc="white",
            mew=0.9, zorder=2)
ax.annotate("SC_RS36625, SC_RS19670\n(DPPY, 251 and 248 aa:\nexcluded by the 348 aa floor)",
            (2, EXCL_TPM[1][2]), xytext=(5, -19), textcoords="offset points",
            color=U.COL_GRAY, fontsize=6, va="top", linespacing=1.25)
ax2.plot(x, OCC, "s--", color=C_OCC, lw=1.1, ms=3.6, zorder=2)

ax.set_xticks(x); ax.set_xticklabels(["T1", "T2", "T3"]); ax.set_xlim(-0.3, 4.6)
ax.set_ylabel("Transcript level (TPM)"); ax.set_ylim(0, 46)
ax2.set_ylabel("6mA occupancy (%)", color=C_OCC)
ax2.set_ylim(0, 40); ax2.set_yticks([0, 10, 20, 30, 40])
ax2.tick_params(axis="y", colors=C_OCC)
ax2.spines["right"].set_color(C_OCC)
ax2.spines["top"].set_visible(False)
ax2.annotate("6mA occupancy at AAGCCCG\n(position 1, n = 1,334)", (0.05, OCC[0]),
             xytext=(0, -12), textcoords="offset points", color=C_OCC,
             fontsize=6, ha="left", va="top")
ax.set_title("", loc="left")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.subplots_adjust(left=0.115, right=0.865, top=0.905, bottom=0.255)

# 出荷名は SuppFigure11.png（補足図 17 -> 11 の改番、2026-10-24 に出力名を追随）。
U.assert_no_text_collisions(fig, "SuppFigure11")
out = os.path.join(OUT_DIR, "SuppFigure11.png")
fig.savefig(out, dpi=300)
print("saved:", out, plt.imread(out).shape)
