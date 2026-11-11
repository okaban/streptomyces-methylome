#!/usr/bin/env python3
"""Real figures for the P1 reviewer-robustness re-analysis bundle.
Reads the TSVs written by the P1_* scripts; muted Okabe-Ito palette."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/analysis/78_reviewer_robustness")
T = os.path.join(BASE, "tables"); F = os.path.join(BASE, "figures")
# 2026-11-03 FIG23（著者決定 2026-11-02）: 意味を持つのは 4mC 赤・6mA 青・デュアル紫だけ。
# この図の区分（near/far、観測/帰無、感度曲線）はマークではないので中立グレー段に置き換える。
# 名前は呼び出し側を変えないために残す。
OK = dict(blue="#444444", vermillion="#2b2b2b", green="#c8c8c8",
          grey="#b4b4b4", purple="#8c8c8c", orange="#d4d4d4")
# 4mC / 6mA のマーク色は原稿共通の定数から取る（00_shared_utils、著者決定 2026-09-26）。
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("_su", os.path.expanduser(
    "~/bioinfo/rna-seq/15_paper_figures/scripts/00_shared_utils.py"))
_su = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_su)
MARK_4mC, MARK_6mA = _su.COL_4mC, _su.COL_6mA
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


# ---- Fig A: cross-talk specificity + probability QC ------------------------
# 2026-09-20: read the C4 tables (4mC is at C4; the C3/C5 QC panel measured
# positions that carry ~0 modification). The near/far test itself is unchanged
# because it scans offsets around each 6mA call, not fixed motif positions.
rates = pd.read_csv(os.path.join(T, "P1_1_crosstalk_rates_C4.tsv"), sep="\t")
probs = pd.read_csv(os.path.join(T, "P1_1_prob_distributions_C4.tsv"), sep="\t")

# 2026-10-22 — この図（Supplementary Figure 5）だけ、出荷されている PNG が
# スクリプトの出力そのものではなかった。
#
#   スクリプト出力 : 3000 x 960 px = 254.0 mm 幅（本誌の 174 mm 上限を超える）
#   出荷 PNG       : 2040 x 653 px = 172.7 mm 幅 = 上の 0.68 倍の縮小ラスタ
#
# 誰かが手で 0.68 倍に縮めて fig_images に置いたもので、figure_manifest.tsv の
# sha256 はさらにそのどちらとも違う第 3 の値だった。縮小はフォントも 0.68 倍に
# するので、出荷 PNG の実効文字サイズは本文 6.1 pt・注記 5.4 pt・凡例 4.8 pt で、
# FIG09 の 6 pt 床を凡例と注記が割っていた（スクリプトの数値リテラルだけを見る
# 監査では見えない種類の違反）。
#
# 対処: 最初から NAR 段幅 6.8 in で描き、その寸法で 6 pt 以上になる字送りを与える。
# 以後ラスタ縮小は不要で、スクリプトの出力がそのまま出荷物になる。
NAR_W_IN = 6.8                      # 2040 px at 300 dpi = 172.7 mm
_rc_saved = {k: plt.rcParams[k] for k in
             ("font.size", "axes.titlesize", "axes.labelsize",
              "xtick.labelsize", "ytick.labelsize", "legend.fontsize")}
plt.rcParams.update({"font.size": 7.0, "axes.titlesize": 7.0, "axes.labelsize": 7.0,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "legend.fontsize": 6.0})
fig, ax = plt.subplots(1, 3, figsize=(NAR_W_IN, 3.2 * NAR_W_IN / 10))
aag = rates[rates.group == "AAGCCCG_6mA"].iloc[0]
non = rates[rates.group == "nonAAGCCCG_6mA_generic"].iloc[0]
# panel a: near vs far co-call rate
x = np.arange(2); w = 0.35
ax[0].bar(x - w/2, [aag.near_rate, non.near_rate], w, label="near (±3,4 bp)", color=OK["vermillion"])
ax[0].bar(x + w/2, [aag.far_rate, non.far_rate], w, label="far (20–50 bp)", color=OK["grey"])
ax[0].set_yscale("log"); ax[0].set_xticks(x)
ax[0].set_xticklabels(["AAGCCCG\n6mA", "non-AAGCCCG\n6mA (generic)"])
ax[0].set_ylabel("4mC co-call rate / offset (log)")
ax[0].set_title("4mC near a 6mA call", loc="center", fontsize=7.0)
ax[0].legend(fontsize=6.0, frameon=False, loc="upper right")
# fold annotations sit just above each bar-pair. Headroom is opened so the
# title band above stays clear (R2-11: title/letter must not touch the label).
ax[0].set_ylim(top=ax[0].get_ylim()[1]*3)
ax[0].text(0 - w/2, aag.near_rate*1.35, f"{aag.near_over_far:.0f}×", ha="center", fontsize=6.5)
ax[0].text(1 - w/2, non.near_rate*1.35, f"{non.near_over_far:.2f}×", ha="center", fontsize=6.5)
# panel b: motif specificity
spec = aag.near_rate / non.near_rate
ax[1].bar([0, 1], [non.near_rate, aag.near_rate], color=[OK["grey"], OK["vermillion"]])
ax[1].set_yscale("log"); ax[1].set_xticks([0, 1])
ax[1].set_xticklabels(["generic\n6mA", "AAGCCCG\n6mA"])
ax[1].set_ylabel("Near 4mC co-call rate (log)")
ax[1].set_title("near vs far co-call rate", loc="center", fontsize=7.0)  # 値は凡例へ（主張タイトルの削除、2026-10-15）
# panel c: probability QC
# 2026-10-22: この 2 本は 6mA = Okabe-Ito blue、4mC = vermillion で、原稿全体で
# 五つ目の 4mC/6mA 配色になっていた。00_shared_utils の著者決定（4mC は赤、
# 6mA は青）に合わせる。Okabe-Ito 自体は CVD 安全だが、図をまたいで同じマークが
# 同じ色であることのほうが読者には効く。形（o / s）も併用して色だけに頼らない。
ax[2].plot(probs.prob_bin_center, probs.AAGCCCG_A_6mA/probs.AAGCCCG_A_6mA.sum(),
           "-o", ms=3, color=MARK_6mA, label="A₀/A₁ (6mA)")
ax[2].plot(probs.prob_bin_center, probs.AAGCCCG_C_4mC/probs.AAGCCCG_C_4mC.sum(),
           "-s", ms=3, color=MARK_4mC, label="C₄ (4mC)")
ax[2].set_xlabel("Per-read modification probability")
ax[2].set_ylabel("Fraction of calls")
ax[2].set_title("call confidence at AAGCCCG", loc="center", fontsize=7.0)
ax[2].legend(fontsize=6.0, frameon=False)
# uniform panel letters in the outer margin (R2-11: letter must not collide
# with the in-plot fold annotation). qc helper anchors a/b/c identically.
try:
    for _ax, _L in zip(ax, "abc"):
        qc_add_panel_label(_ax, _L)
except NameError:
    for _ax, _L in zip(ax, "abc"):
        _ax.text(-0.12, 1.06, _L, transform=_ax.transAxes,
                 fontsize=9, fontweight="bold", va="bottom", ha="right")
fig.tight_layout(); fig.savefig(os.path.join(F, "P1_1_crosstalk.png"), dpi=300); plt.close(fig)
plt.rcParams.update(_rc_saved)      # 以降の 3 図は従来どおりの字送りで描く

# ---- Fig B: 5mC vs 4mC (P1-2) ---------------------------------------------
fig, ax = plt.subplots(figsize=(3.4, 3.2))
ax.bar([0, 1], [0.0014, 0.0226], color=[OK["grey"], OK["vermillion"]])
ax.axhline(0.01, ls="--", color="k", lw=0.8)
ax.text(1.4, 0.0105, "0.01% claim", fontsize=7, va="bottom", ha="right")
ax.set_xticks([0, 1]); ax.set_xticklabels(["5mC (m)", "4mC (21839)"])
ax.set_ylabel("Genome-wide modified fraction (%)\ncytosine calls, cov≥10")
ax.set_title("Cytosine modification is m4C, not 5mC", loc="left", fontweight="bold", fontsize=9)
for i, v in enumerate([0.0014, 0.0226]):
    ax.text(i, v + 0.0008, f"{v:.4f}%", ha="center", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(F, "P1_2_5mC_vs_4mC.png"), dpi=300); plt.close(fig)

# ---- Fig C: spatial autocorrelation (P1-3) --------------------------------
# rebuild a quick null for the figure (cheap)
from scipy import stats
A = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/analysis")
dose = pd.read_csv(os.path.join(A, "41_GCCGGC_dose_response/tables/GCCGGC_gene_site_counts.tsv"), sep="\t")
trans = pd.read_csv(os.path.join(A, "42_GCCGGC_temporal_derepression/tables/gene_methylation_transitions.tsv"), sep="\t")[["gene_id", "LFC_T2vsT1"]]
d = dose.merge(trans, left_on="locus_tag", right_on="gene_id").dropna(subset=["LFC_T2vsT1", "gene_start", "dose_group"]).sort_values("gene_start").reset_index(drop=True)
code = d.dose_group.astype(str).map({"0":0,"1":1,"2":2,"3":3,"4+":4}).values
lfc = d.LFC_T2vsT1.values
obs, _ = stats.spearmanr(code, lfc)
rng = np.random.default_rng(42)
null = np.array([stats.spearmanr(code, np.roll(lfc, rng.integers(1, len(lfc))))[0] for _ in range(3000)])
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.2))
ax[0].hist(null, bins=40, color=OK["grey"], alpha=0.8)
ax[0].axvline(obs, color=OK["vermillion"], lw=2, label=f"observed ρ = {obs:+.3f}\n(block-perm p = 0.22)")
ax[0].axvline(-np.percentile(np.abs(null),95), color="k", ls=":", lw=0.8)
ax[0].axvline(np.percentile(np.abs(null),95), color="k", ls=":", lw=0.8, label="null 95%")
ax[0].set_xlabel("Spearman ρ (dose vs LFC)"); ax[0].set_ylabel("Permutations")
# 2026-09-13 (FIG-05): neutral titles, no editorial labels ("→ artefact",
# "SURVIVES"); panel letters a-d run across the two stacked rows (S18).
ax[0].set_title("a  Genome-wide dose vs log$_2$FC (n = 7,496)\n    spatial-block permutation null",
                loc="left", fontweight="bold", fontsize=9)
ax[0].set_ylim(top=ax[0].get_ylim()[1] * 1.28)  # headroom for the legend
ax[0].legend(fontsize=7, frameon=False, loc="upper center")
# panel b: regulatory survives
reg = pd.read_csv(os.path.join(T, "P1_3b_permissive_spatial.tsv"), sep="\t").iloc[0]
ax[1].bar([0, 1], [reg.p_partial_naive, reg.block_perm_p_partial], color=[OK["grey"], OK["green"]])
ax[1].axhline(0.05, ls="--", color="k", lw=0.8); ax[1].set_yscale("log")
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["naive p", "block-perm p"])
ax[1].set_ylabel("p (partial Spearman | region)")
ax[1].set_title(f"b  Regulatory genes (n = {int(reg.n):,})\n    partial Spearman ρ = {reg.rho_partial:+.2f}\n    (TSS–GCCGGC distance vs log$_2$FC | region)",
                loc="left", fontweight="bold", fontsize=8)
ax[1].text(0.5, 0.95, f"naive p = {reg.p_partial_naive:.1e}\nblock-perm p = {reg.block_perm_p_partial:.4f}",
           transform=ax[1].transAxes, fontsize=7.5, ha="center", va="top")
fig.tight_layout(); fig.savefig(os.path.join(F, "P1_3_spatial_autocorr.png"), dpi=300); plt.close(fig)

# ---- Fig D: equivalence + threshold sensitivity (P1-4/5) ------------------
eq = pd.read_csv(os.path.join(T, "P1_4_equivalence.tsv"), sep="\t")
sens = pd.read_csv(os.path.join(T, "P1_5_threshold_sensitivity.tsv"), sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.2))
row = eq[eq.contrast == "|LFC| T2vsT1"].iloc[0]
ax[0].axvspan(-0.147, 0.147, color=OK["green"], alpha=0.15, label="negligible band")
ax[0].errorbar([row.cliff_delta], [0], xerr=[[row.cliff_delta-row.cliff_lo],[row.cliff_hi-row.cliff_delta]],
               fmt="o", color=OK["vermillion"], capsize=4)
ax[0].axvline(0, color="k", lw=0.6)
ax[0].set_yticks([]); ax[0].set_xlim(-0.5, 0.5)
ax[0].set_xlabel("Cliff's δ (Exposed − Shielded |LFC|)")
ax[0].set_title(f"c  Permissive equivalence (TOST)\n    Cliff's δ = {row.cliff_delta:+.2f}, 95% CI [{row.cliff_lo:+.2f}, {row.cliff_hi:+.2f}]",
                loc="left", fontweight="bold", fontsize=9)
ax[0].legend(fontsize=6.0, frameon=False, loc="upper right")
ax[1].plot(sens.window_bp, sens.cliff_delta, "-o", color=OK["blue"], ms=4)
ax[1].axhspan(-0.147, 0.147, color=OK["green"], alpha=0.15)
ax[1].axvline(293, ls="--", color=OK["vermillion"], lw=1, label="293 bp")
ax[1].set_xlabel("Exposed window (bp)"); ax[1].set_ylabel("Cliff's δ (|LFC|)")
ax[1].set_ylim(-0.4, 0.4)
ax[1].set_title("d  Cliff's δ across Exposed-window\n    definitions (100–500 bp)", loc="left", fontweight="bold", fontsize=9)
ax[1].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(F, "P1_4_5_equivalence.png"), dpi=300); plt.close(fig)

# ---- S18 slot: rows C (a,b) and D (c,d) stacked into one PNG -------------
from PIL import Image
_top = Image.open(os.path.join(F, "P1_3_spatial_autocorr.png"))
_bot = Image.open(os.path.join(F, "P1_4_5_equivalence.png"))
_W = max(_top.width, _bot.width); _gap = 40
_st = Image.new("RGB", (_W, _top.height + _bot.height + _gap), "white")
_st.paste(_top, (0, 0)); _st.paste(_bot, (0, _top.height + _gap))
_st.save(os.path.join(F, "P1_3_4_5_stacked_S18.png"), dpi=(300, 300))

print("Figures ->", F)
for fn in ["P1_1_crosstalk.png","P1_2_5mC_vs_4mC.png","P1_3_spatial_autocorr.png","P1_4_5_equivalence.png"]:
    p = os.path.join(F, fn); print("  ", fn, "OK" if os.path.exists(p) else "MISSING")
