# -*- coding: utf-8 -*-
"""Supplementary Figure 7 — region composition of modified sites, by mark and timepoint.

生成元が無い図だったため 2026-10-15 にスクリプト化した。旧版は画像として出荷されており、
焼き込まれた部位数が T3 の旧値（1,877 / 977 / 457）のまま、凡例の正準値（1,866 / 973 / 453）と
食い違っていた。この再実装でその不整合が消える。

区分規則と入力の来歴は data/region261015/README.md。要点だけ再掲すると、promoter は
**遺伝子体の外**かつ最近接 TSS の上流 1–500 bp で、旧版の「遺伝子体に関係なく promoter 優先」
とは百分率が比較できない。補足図 18 パネル a は旧版規則なので割合を直接並べてはいけない。

Run: PYTHONPATH=. python3 suppfigure7_region_composition.py
"""
import textwrap

import matplotlib.pyplot as plt
import pandas as pd

from _suppstyle import (DATA_DIR, GREY, NAR_WIDTH_IN, TAN,
                        apply_figure_style, figure_panel_letters, report, save,
                        source_note)

import importlib  # noqa: E402
import sys as _sys  # noqa: E402
_sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
_U = importlib.import_module("00_shared_utils")

apply_figure_style(frame="open", sizes=(7.0, 6.2, 6.0))

CSV = DATA_DIR / 'region261015' / 'region_composition_by_mark_261024.csv'
ORDER = ['promoter', "5'UTR", 'CDS-internal', 'intergenic']
# mark 値は入力表の照合キー。表示名と分けて持つ（表のラベルを変えるとデータが黙って消える）。
# 2026-10-22: 表示名を 8 塩基に直す。この行の n_sites は 1,886 / 1,892 / 1,866 で
# numbers.json の TGG_called（TGGCCGGC 8 塩基の正準センサス）と完全に一致する。
# 六塩基 GCCGGC はこの図の量ではない（六塩基は availability の分母で 32,102 箇所）。
# CSV の mark 列は照合キーなので変えない。
# 2026-10-24: 著者コメント #7 により 4 パネル 2x2 に変更し、AAGCCCG の両マーク同時
# （同一インスタンスで C5 の 4mC と A1/A2 の 6mA がともに呼ばれたもの）を 4 枚目に足す。
# 入力表は region261015/build_region_composition_261024.py が生成する。同スクリプトは
# 旧 9 行を小数第 2 位まで完全再現することを assert しているので、4 行目は他の 3 パネルと
# 同じ区分規則・同じ部位判定規則で出ている。
MARKS = [('GCCGGC hexamer 4mC', 'TGGCCGGC 4mC'),
         ('AAGCCCG 4mC', 'AAGCCCG 4mC'),
         ('AAGCCCG 6mA', 'AAGCCCG 6mA'),
         ('AAGCCCG 4mC + 6mA', 'AAGCCCG 4mC + 6mA, same instance')]
TPS = ['T1', 'T2', 'T3']
# 4 区分は順序尺度ではないので濃淡ではなく別色にする（積み上げ棒では色が唯一の識別子）。
# 2026-10-22: この 4 色は RED / BLUE / TAN / GREY で、先頭 2 色が 4mC と 6mA の
# マーク色そのものだった。パネル見出しが "TGGCCGGC 4mC" / "AAGCCCG 6mA" である図の中で
# 赤が「promoter」、青が「5'UTR」も指してしまい、色が実体を横断して追えなくなる。
# ゲノム区分はマークではないので、マーク色と重ならない 4 色に移す。
# 2026-10-24（著者コメント #7「色の違いが判別しづらい。リッチなカラーで統一」）:
# 彩度を上げ、明度差を広げる。面積の 9 割を占める CDS-internal だけを低彩度の淡色に
# 落とし、1-6 % しかない 3 区分に彩度を寄せる。小さい区切りほど目立つ配色にしないと、
# 色相を変えても判別は上がらない。マーク色 COL_4mC/#8A3D37・COL_6mA/#68849D とは
# 色相でも明度でも離してある。
REGION_VIOLET = '#4C2A85'   # promoter
REGION_TEAL = '#0F7B6C'     # 5'UTR
REGION_SAND = '#E8C39E'     # CDS-internal（面積が大きいので低彩度）
REGION_SLATE = '#46525C'    # intergenic
COLS = [REGION_VIOLET, REGION_TEAL, REGION_SAND, REGION_SLATE]

df = pd.read_csv(CSV)
fig, axes2d = plt.subplots(2, 2, figsize=(NAR_WIDTH_IN, 4.95), dpi=300)
axes = axes2d.ravel()
fig.subplots_adjust(left=0.085, right=0.988, top=0.918, bottom=0.270,
                    wspace=0.16, hspace=0.52)

for ax, (key, label) in zip(axes, MARKS):
    sub = df[df['mark'] == key].set_index('timepoint').loc[TPS]
    assert len(sub) == 3, f"入力表に mark='{key}' の 3 時点が揃っていない"
    bottom = [0.0] * 3
    for region, col in zip(ORDER, COLS):
        vals = sub[region].to_numpy()
        ax.bar(range(3), vals, bottom=bottom, width=0.62, color=col,
               linewidth=0, zorder=3, label=region if ax is axes[0] else None)
        bottom = [b + v for b, v in zip(bottom, vals)]
    for i, tp in enumerate(TPS):
        ax.text(i, 102.5, f"n = {int(sub.loc[tp, 'n_sites']):,}", ha='center',
                va='bottom', fontsize=6.0, color=GREY)
    ax.set_xticks(range(3)); ax.set_xticklabels(TPS)
    ax.set_ylim(0, 112); ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xlim(-0.6, 2.6)
    ax.set_title(label, loc='left', pad=14)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
for ax in (axes[0], axes[2]):
    ax.set_ylabel('% of modified sites')
for ax in (axes[1], axes[3]):
    ax.set_yticklabels([])

# 凡例の handle は axes[0] に付けたラベルから取り、描画は図の下端に置く。
# axes[n].legend() を引数なしで呼ぶとそのパネル自身のラベルを探して空になる。
_h, _l = axes[0].get_legend_handles_labels()
assert len(_h) == len(ORDER), f"凡例の handle が {len(_h)} 個しかない"
fig.legend(_h, _l, loc='lower center', bbox_to_anchor=(0.5, 0.168), ncol=4,
           frameon=False, fontsize=6.0, handlelength=1.1,
           columnspacing=1.5, handletextpad=0.45)

figure_panel_letters(fig, list(axes), 'abcd',
                     dx=[0.070, 0.030, 0.070, 0.030], dy=0.020)
source_note(fig, textwrap.fill(
            'Canonical per-timepoint census (three libraries pooled, pre-pooling true depth '
            'Nvalid_cov + Nfail >= 30, pooled Nmod/Nvalid >= 0.5). Regions assigned in priority '
            'order: promoter (outside any gene body, 1-500 bp upstream of the nearest of 3,030 '
            'experimental primary or secondary TSSs), 5\u2032UTR (0 to +100 bp of a TSS), '
            'CDS-internal, intergenic. This differs from the earlier rule, which assigned promoter '
            'first irrespective of gene body, so percentages are not comparable between versions. '
            'Panel d counts AAGCCCG instances at which both marks are called in the same timepoint '
            '(4mC at C5 and 6mA at A1 or A2 of that instance); regions are assigned at the C5 '
            'position, the same basis as panels b and c.',
            148), y=0.020)
# 出荷名は SuppFigure4.png（補足図 7 -> 4 の改番、2026-10-24 に出力名を追随）。
# 拡張子を省くと拡張子なしのファイルが出るので .png まで書く。
report(save(fig, 'SuppFigure4.png'))
