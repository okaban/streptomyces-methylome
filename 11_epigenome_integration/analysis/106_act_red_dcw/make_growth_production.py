#!/usr/bin/env python3
"""Supplementary Figure 25 — growth and antibiotic production at the three sampled timepoints.

Source of the values: the author's own measurements, pasted from the Notion lab
notebook pages "M145培養液におけるAct/Redの定量" (ID 021) and "M145株の乾燥重量測定"
(ID 022) on 2026-09-25, and written to tables/act_red_dcw_raw.tsv. The author's
exported charts (notion_source/*.svg) are kept for reference but are NOT shipped:
none of them carries x-axis tick labels, so the reader cannot tell which bar is
which timepoint. This script redraws them in the project's grey BAR_RAMP, which
matches the palette of those exports (#293039 is common to both).

Reproduces the author's reported statistics exactly: one-way ANOVA
F(2,9) = 199.04 (dry cell weight) and 192.79 (actinorhodin).
"""
import csv, importlib.util, os, sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

H = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "su", os.path.join(H, "..", "..", "..", "15_paper_figures", "scripts", "00_shared_utils.py"))
su = importlib.util.module_from_spec(spec); sys.modules["su"] = su; spec.loader.exec_module(su)

rows = list(csv.DictReader(open(os.path.join(H, "tables", "act_red_dcw_raw.tsv"), encoding="utf-8"),
                           delimiter="\t"))
RAW = {}
for r in rows:
    RAW.setdefault(r["measurement"], {}).setdefault(r["timepoint"], []).append(float(r["value"]))

PANELS = [("dry_cell_weight_g_L", "Dry cell weight (g/L)"),
          ("actinorhodin_mg_L", "Actinorhodin (mg/L)"),
          ("undecylprodigiosin_mg_L", "Undecylprodigiosin (mg/L)")]
TP = ("T1", "T2", "T3")


def stars(p):
    for thr, s in ((0.05, "ns"), (0.01, "*"), (0.001, "**"), (1e-4, "***")):
        if p >= thr:
            return s
    return "****"


def main():
    su.apply_unified_style()
    fig = plt.figure(figsize=(su.mm_to_inch(174), su.mm_to_inch(66)))
    axes = [fig.add_axes([0.085 + 0.315 * i, 0.22, 0.205, 0.56]) for i in range(3)]
    stat_rows = []
    for ax, (key, lab) in zip(axes, PANELS):
        g = [np.array(RAW[key][t]) for t in TP]
        F, p = stats.f_oneway(*g)
        tk = stats.tukey_hsd(*g)
        m = [v.mean() for v in g]
        sd = [v.std(ddof=1) for v in g]
        x = np.arange(3)
        ax.bar(x, m, yerr=sd, capsize=2.5, width=0.62, color=su.BAR_RAMP, edgecolor="none",
               error_kw=dict(lw=0.8, ecolor=su.COL_DARK))
        for xi, v in zip(x, g):
            ax.scatter(np.full(len(v), xi) + np.linspace(-0.14, 0.14, len(v)), v, s=8,
                       facecolor="white", edgecolor=su.COL_DARK, linewidth=0.6, zorder=3)
        top = max(np.array(m) + np.array(sd))
        ax.set_ylim(0, top * 1.42)
        for k, (i, j) in enumerate([(0, 1), (1, 2), (0, 2)]):
            y = top * (1.08 + 0.11 * k)
            ax.plot([i, i, j, j], [y, y * 1.015, y * 1.015, y], lw=0.7, color=su.COL_DARK)
            ax.text((i + j) / 2, y * 1.025, stars(tk.pvalue[i, j]), ha="center", va="bottom",
                    fontsize=su.FONT["annot"])
            stat_rows.append(dict(measurement=key, comparison=f"{TP[i]}_vs_{TP[j]}",
                                  p_tukey=float(f"{tk.pvalue[i, j]:.3g}"),
                                  anova_F=round(float(F), 2), anova_p=float(f"{p:.3g}")))
        ax.set_xticks(x)
        ax.set_xticklabels(["12 h", "24 h", "50 h"])
        ax.set_ylabel(lab)
        ax.set_xlabel("time after inoculation")
    for ax, L in zip(axes, "abc"):
        ax.text(-0.30, 1.10, L, transform=ax.transAxes, fontsize=su.FONT["panel_letter"],
                fontweight="bold", va="bottom", ha="left")
    su.assert_no_text_collisions(fig, "GrowthProduction")
    su.save_figure(fig, os.path.join(H, "figures", "growth_production"),
                   formats=("png", "pdf"), width_class="full")
    out = os.path.join(H, "tables", "act_red_dcw_stats.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(stat_rows[0]))
        w.writeheader(); w.writerows(stat_rows)
    for r in stat_rows[::3]:
        print(r["measurement"], "ANOVA F =", r["anova_F"], "p =", r["anova_p"])


if __name__ == "__main__":
    sys.exit(main())
