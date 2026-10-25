"""Fair copy-number control figure (reviewer request):
(a) chromosome-wide T1 coverage (copy-number proxy) overlaid with GCCGGC
    methylated-site density -> shows the oriC coverage gradient transparently;
(b) coverage-stratified %methylated of GCCGGC motifs, core vs arm
    -> core > arm at MATCHED coverage = enrichment is not a copy-number artefact;
(c) per-site methylation frequency, core vs arm -> intensity is uniform.
"""
import os, subprocess
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis"
PILE = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015/1-1_pileup.bed"
OUT = f"{B}/78_reviewer_robustness"
SITES = f"{B}/37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv"
strat = pd.read_csv(f"{OUT}/tables/C1_coverage_stratified.tsv", sep="\t")
OK = dict(verm="#D55E00", blue="#0072B2", grey="#999999", green="#009E73")
CORE_LO, CORE_HI = 1_500_000, 7_170_000; WIN = 50_000; ORIC = 4_270_000

# --- per-50kb mean coverage (T1, code 21839) via awk ---
awk = r'''$4=="21839"{w=int($2/%d); s[w]+=$10; c[w]++} END{for(w in s) print w"\t"s[w]/c[w]}''' % WIN
res = subprocess.run(["awk","-F","\t",awk,PILE], capture_output=True, text=True)
cov = pd.DataFrame([l.split("\t") for l in res.stdout.strip().split("\n")],
                   columns=["win","mean_cov"]).astype({"win":int,"mean_cov":float}).sort_values("win")
cov["mb"] = cov.win*WIN/1e6

# --- per-50kb GCCGGC methylated-site density (T1) ---
st = pd.read_csv(SITES, sep="\t"); st = st[st.timepoint=="T1"].copy()
st["win"] = (st.position//WIN).astype(int)
dens = st.groupby("win").size().rename("n_meth").reset_index(); dens["mb"]=dens.win*WIN/1e6

fig = plt.figure(figsize=(11,7))
gs = fig.add_gridspec(2,2, height_ratios=[1,1], hspace=0.42, wspace=0.28)

# (a) chromosome tracks
axa = fig.add_subplot(gs[0,:])
axa.axvspan(0, CORE_LO/1e6, color=OK["grey"], alpha=0.12)
axa.axvspan(CORE_HI/1e6, 8.667, color=OK["grey"], alpha=0.12, label="arm")
axa.fill_between(cov.mb, cov.mean_cov, color=OK["grey"], alpha=0.55, step="mid", label="coverage (copy-number proxy)")
axa.axvline(ORIC/1e6, color="k", ls=":", lw=1); axa.text(ORIC/1e6, axa.get_ylim()[1]*0.92, " oriC", fontsize=8)
axa.set_ylabel("mean coverage (×)", color="#555")
axa.set_xlabel("chromosomal position (Mb)")
axa.set_title("a  T1 coverage (copy-number) vs GCCGGC methylation density", loc="left", fontweight="bold", fontsize=10)
axb = axa.twinx()
axb.bar(dens.mb, dens.n_meth, width=WIN/1e6, color=OK["verm"], alpha=0.85, label="GCCGGC m4C sites")
axb.set_ylabel("GCCGGC m4C sites / 50 kb", color=OK["verm"])
h1,l1 = axa.get_legend_handles_labels(); h2,l2 = axb.get_legend_handles_labels()
axa.legend(h1+h2, l1+l2, fontsize=7, loc="upper left", framealpha=0.9)

# (b) coverage-stratified %meth core vs arm
axc = fig.add_subplot(gs[1,0])
x = np.arange(len(strat)); w=0.38
axc.bar(x-w/2, strat.core_pct_meth, w, color=OK["verm"], label="core")
axc.bar(x+w/2, strat.arm_pct_meth, w, color=OK["blue"], label="arm")
axc.set_xticks(x); axc.set_xticklabels(strat.cov_bin.str.replace("-inf","+"), fontsize=8)
axc.set_xlabel("coverage band (×)"); axc.set_ylabel("% of GCCGGC motifs methylated")
axc.set_title("b  matched-coverage: core > arm\n(enrichment is not copy-number)", loc="left", fontweight="bold", fontsize=10)
axc.legend(fontsize=8, frameon=False)

# (c) per-site frequency core vs arm (T1)
axd = fig.add_subplot(gs[1,1])
core_f = st[st.region=="core"].frequency.values; arm_f = st[st.region=="arm"].frequency.values
parts = axd.violinplot([core_f, arm_f], showmedians=True)
for i,b in enumerate(parts['bodies']): b.set_facecolor([OK["verm"],OK["blue"]][i]); b.set_alpha(0.6)
axd.set_xticks([1,2]); axd.set_xticklabels([f"core\n(n={len(core_f)})", f"arm\n(n={len(arm_f)})"])
axd.set_ylabel("per-site m4C frequency (%)")
axd.set_title(f"c  methylation intensity uniform\n(median {np.median(core_f):.0f}% vs {np.median(arm_f):.0f}%)", loc="left", fontweight="bold", fontsize=10)

fig.savefig(f"{OUT}/figures/C1_copynumber_control.png", dpi=200, bbox_inches="tight")
print("saved figures/C1_copynumber_control.png")
print(f"coverage windows: {len(cov)}, methylated-site windows: {len(dens)}")
