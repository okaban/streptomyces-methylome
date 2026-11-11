# 2026-11-04 著者決定: 区別は軸ラベル・パネル内表記・形が担っているので、色で重ねて符号化しない（中立グレー）。
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NF-R4 — how the methyltransferase candidates were narrowed, and what that
cost the three candidates named by Lutz et al. 2019.

    python3 nf_r4_figure_261101.py

Writes tables/nf_r4_funnel_261101.tsv and the PNG.

The figure exists to answer two open findings:
  V19 - the manuscript eliminates SCO3510 and SCO6885, which Lutz et al. 2019
        proposed as N4-cytosine MTases, without showing why. Both fall at the
        same stage, the 348 aa size floor, and the figure puts them there by
        name.
  V21 - the expression axis must be TPM. `MTase_candidate_elimination_261008.csv`
        carries expression_T1/T2/T3 columns that are DESeq2-normalised counts,
        not TPM (SC_RS17645: 245.11/51.74/145.14 there, against TPM
        9.95/1.56/5.47 which is what canon records). This script reads TPM from
        113_mtase_completeness_261024 and asserts the canon value before
        drawing, so the defective columns cannot reach a panel.

Stage numbering follows 113_mtase_completeness_261024 (S2 = 21, S3 = 12,
S4 = 5, S5 = 3). The entry set is the 129 annotated SAM-dependent MTases
(registry MTASE.CATALOG.TOTAL). The hatched extension on S2 and S3 is the one
protein an annotation-independent Pfam-A scan adds, SC_RS19860, which the
product-name entry misses; it changes no later stage.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("SUPPFIG_OUT", HERE)
SUPP = "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014"
sys.path.insert(0, SUPP)
sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from _suppstyle import (META_GREY, NAR_WIDTH_IN, apply_figure_style,  # noqa: E402
                        figure_panel_letters, report, save)
import importlib  # noqa: E402
U = importlib.import_module("00_shared_utils")

A = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
T107 = f"{A}/107_MTase_candidate_triage/tables/MTase_candidate_elimination_261008.csv"
T113 = f"{A}/113_mtase_completeness_261024/tables/mtase_completeness_261024.tsv"
TAB = os.path.join(HERE, "tables")
NAME = "NF-R4_MTase_candidates_261101.png"

SIZE_FLOOR = 348          # registry MTASE.SIZE_FLOOR.AA
LUTZ = {"SC_RS17645": "SCO3104", "SC_RS19670": "SCO3510",
        "SC_RS36625": "SCO6885"}
EXTRA = "SC_RS19860"      # entrant found only by the annotation-independent scan
STAGES = [
    ("Entry", "Annotated MTases", 1, 0),
    ("S1", "Translation present", 2, 0),
    ("S2", "Motif IV present", 3, 1),
    ("S3", "Motif IV in catalytic domain", 4, 1),
    ("S4", "DNA amino-MTase family", 5, 0),
    ("S5", f"Length \u2265 {SIZE_FLOOR} aa", 6, 0),
]
C_PASS = U.COL_BAR_DARK
C_FAIL = U.COL_BAR_LIGHT


def load():
    t107 = pd.read_csv(T107)
    t113 = pd.read_csv(T113, sep="\t")
    assert len(t107) == 129, len(t107)
    tp = t113.set_index("locus_tag")[["TPM_T1", "TPM_T2", "TPM_T3",
                                      "length_aa"]]
    # V21 guard: the TPM source must reproduce the canonical SC_RS17645 series
    got = tp.loc["SC_RS17645", ["TPM_T1", "TPM_T2", "TPM_T3"]].round(2).tolist()
    assert got == [9.95, 1.56, 5.47], got
    return t107, tp


def funnel_table(t107):
    rows = []
    for sid, label, thr, extra in STAGES:
        n = int((t107.surviving_stage >= thr).sum())
        rows.append(dict(stage=sid, criterion=label.replace("\n", " "),
                         n_annotation_based=n,
                         n_annotation_independent=n + extra))
    return pd.DataFrame(rows)


def panel_a(ax, fun):
    y = np.arange(len(fun))[::-1]
    ax.barh(y, fun.n_annotation_based, height=0.62, color=C_PASS, zorder=3)
    for yy, n, nn in zip(y, fun.n_annotation_based,
                         fun.n_annotation_independent):
        txt = f"{n}" if n == nn else f"{n} ({nn})"
        ax.text(nn + 4.0, yy, txt, va="center", ha="left", fontsize=6.2,
                color=U.COL_BAR_DARK)
    ax.set_yticks(y)
    ax.set_yticklabels([s[1] for s in STAGES], fontsize=5.8)
    ax.set_xlim(0, 172)
    ax.set_ylim(-0.6, len(fun) - 0.4)
    ax.set_xlabel("Proteins remaining (count)")
    ax.margins(y=0.04)


def six(t107, tp):
    """The five proteins that reach S4, plus the annotation-independent one."""
    s4 = t107[t107.surviving_stage >= 5][["locus_tag", "SCO", "length"]]
    rows = [dict(locus_tag=r.locus_tag, sco=r.SCO, length=int(r.length))
            for r in s4.itertuples()]
    rows.append(dict(locus_tag=EXTRA, sco="",
                     length=int(tp.loc[EXTRA, "length_aa"])))
    for r in rows:
        r["passes"] = r["length"] >= SIZE_FLOOR and r["locus_tag"] != EXTRA
        r["lutz"] = r["locus_tag"] in LUTZ
        for tpt in ("T1", "T2", "T3"):
            r[tpt] = float(tp.loc[r["locus_tag"], f"TPM_{tpt}"])
    return pd.DataFrame(rows).sort_values("length", ascending=False)


def label_of(r):
    """2 行の水平ラベル（2026-11-07 C04: 回転ラベルをやめて読みやすくする）。"""
    second = f"({r.sco})" if r.sco else ""
    if r.lutz:
        second = (second + " \u2020").strip()
    return f"{r.locus_tag}\n{second}" if second else r.locus_tag


def panel_b(ax, d):
    """Length of the six candidates against the size floor.

    2026-11-07 共著者コメント C04「書くなら、縦軸と横軸を入れ替えては？あと、
    このグラフでどこに注目したらいいのか不明」。蛋白質を x、長さを y に入れ替え、
    下段の発現パネルと x を共有させて、床を割る 3 本が縦に揃うようにした。
    """
    x = np.arange(len(d))
    col = [C_PASS if p else C_FAIL for p in d.passes]
    ax.vlines(x, 0, d.length, color=col, lw=1.0, zorder=2)
    ax.scatter(x, d.length, s=22, c=col, zorder=3)
    ax.axhline(SIZE_FLOOR, color="#6b6b6b", lw=0.9, ls=(0, (3, 2)), zorder=1)
    # 2026-11-09 [4]: 下限線の説明は凡例文へ（図中に文字を置かない）。
    ax.set_xticks(x)
    ax.set_xticklabels([])
    ax.set_xlim(-0.6, len(d) - 0.4)
    ax.set_ylim(0, d.length.max() * 1.12)
    ax.set_ylabel("Length (aa)")


TP_STYLE = [("T1", "T1, 12 h", "o", U.COL_BAR_LIGHT),
            ("T2", "T2, 24 h", "^", U.COL_BAR_MID),
            ("T3", "T3, 50 h", "s", U.COL_BAR_DARK)]


def panel_c(ax, d):
    # 2026-11-09 著者コメント [4]: 記号が重なって増減が読めない。時点 3 本の
    # 灰色棒に替え、遺伝子ごとに T1→T3 の向きが見えるようにする。
    x = np.arange(len(d))
    BW = 0.26
    FLOOR = 0.05          # 対数軸の棒は底から描けないので下限を置く
    for k, (col, lab, mk, c) in enumerate(TP_STYLE):
        v = np.maximum(d[col].to_numpy(dtype=float), FLOOR)
        ax.bar(x + (k - 1) * BW, v - FLOOR, width=BW * 0.86, bottom=FLOOR,
               color=U.TIME_RAMP[k], linewidth=0, label=lab, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.05, 900)
    ax.set_yticks([0.1, 1, 10, 100])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.set_ylabel("Expression (TPM)")
    ax.set_xticks(x)
    ax.set_xticklabels([label_of(r) for r in d.itertuples()], fontsize=5.8)
    ax.set_xlim(-0.6, len(d) - 0.4)
    ax.legend(frameon=False, fontsize=5.6, loc="upper left",
              bbox_to_anchor=(0.005, 1.02), ncol=3, handlelength=0.8,
              handletextpad=0.35, borderpad=0.0, columnspacing=0.9)


def main():
    apply_figure_style()
    t107, tp = load()
    fun = funnel_table(t107)
    os.makedirs(TAB, exist_ok=True)
    d = six(t107, tp)
    fun.to_csv(f"{TAB}/nf_r4_funnel_261101.tsv", sep="\t", index=False)
    d.to_csv(f"{TAB}/nf_r4_candidates_261101.tsv", sep="\t", index=False)
    print(fun.to_string(index=False))
    print()
    print(d.to_string(index=False))

    # 2026-11-07 共著者コメント C06「これ、グラフで書くものではなく表とかで
    # まとめるものでは？本来、あるいは、本文で書いてもいいレベル。とにかく fig
    # グラフでまとめる必要はないかと」。旧パネル (a)（絞り込み funnel）は図から外し、
    # 段階ごとの残存数は本文 R4 の 1 文と tables/nf_r4_funnel_261101.tsv に残す。
    # funnel_table() は引き続き実行して表を書き出す（図に描かないだけ）。
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 1.95))
    ax_b = fig.add_axes([0.140, 0.4638, 0.845, 0.4381])
    ax_c = fig.add_axes([0.140, 0.1304, 0.845, 0.2877])
    axes = [ax_b, ax_c]
    panel_b(ax_b, d)
    panel_c(ax_c, d)
    figure_panel_letters(fig, axes, ["a", "b"], dx=(0.128, 0.128), dy=0.012)
    p = save(fig, NAME)
    report(p)


if __name__ == "__main__":
    main()
