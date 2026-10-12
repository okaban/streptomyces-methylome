"""Reviewer-requested figures, reframe-consistent (2026-06-15).

Fig2c/d : Shielded(989)/Exposed(62) nearest-distance + expression-variability (62/989 set).
Pie     : methylation occupancy by genomic feature (GCCGGC vs AAGCCCG; promoter/5'UTR/CDS/intergenic)
          -> addresses 'occupancy pie charts' + 'intergenic' reviewer comments.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import mannwhitneyu

B = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis")
OUT = B/"77_reviewer_figures/figures"; OUT.mkdir(parents=True, exist_ok=True)
W = 293
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

# ---- data: regulatory genes with non-circular Exposed label ----
reg = pd.read_csv(B/"52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv", sep="\t").dropna(subset=['tss']).copy()
reg['tss'] = reg['tss'].astype(int)
g = pd.read_csv(B/"37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv", sep="\t")
def nd(tss, pos):
    if len(pos)==0: return np.nan
    i=np.clip(np.searchsorted(pos,tss),1,len(pos)-1); return min(abs(tss-pos[i-1]),abs(tss-pos[i]))
posT1 = np.sort(g[g.timepoint=='T1'].position.values)
reg['dT1'] = reg['tss'].apply(lambda t: nd(t,posT1))
reg['Exposed'] = reg['dT1'] <= W
E = reg[reg.Exposed]; S = reg[~reg.Exposed]

# ===== Figure 2c/d =====
fig, ax = plt.subplots(1, 2, figsize=(8.4, 4.0))
# (c) nearest-distance violin (log10)
sd = np.log10(S['dT1'].clip(lower=1)); ed = np.log10(E['dT1'].clip(lower=1))
parts = ax[0].violinplot([sd, ed], showmedians=True)
ax[0].axhline(np.log10(W), color="red", ls="--", lw=1, label=f"293 bp boundary")
ax[0].set_xticks([1,2]); ax[0].set_xticklabels([f"Shielded\n(n={len(S)})", f"Exposed\n(n={len(E)})"])
ax[0].set_ylabel("log₁₀ nearest GCCGGC 4mC distance (bp)")
u,p = mannwhitneyu(S['dT1'].dropna(), E['dT1'].dropna())
ax[0].set_title(f"(c) Promoter-proximal methylation\nMWU p={p:.1e}"); ax[0].legend(frameon=False, fontsize=8)
# (d) expression variability |LFC|
ev = [S['LFC_T2vsT1'].abs().dropna(), E['LFC_T2vsT1'].abs().dropna(),
      S['LFC_T3vsT1'].abs().dropna(), E['LFC_T3vsT1'].abs().dropna()]
bp = ax[1].boxplot(ev, positions=[1,1.7,3,3.7], widths=0.55, showfliers=False, patch_artist=True)
for i,patch in enumerate(bp['boxes']): patch.set_facecolor(["#9AA7B0","#A64B44"][i%2])
ax[1].set_xticks([1.35,3.35]); ax[1].set_xticklabels(["T2 vs T1","T3 vs T1"])
ax[1].set_ylabel("|log₂ fold-change|")
ax[1].set_title("(d) Expression variability (T2/T3 vs T1)")
ax[1].annotate("n.s. (Cliff's δ=−0.01,\nequivalent → permissive)", xy=(2.2,ax[1].get_ylim()[1]*0.9),
               ha="center", fontsize=8, color="#666")
fig.tight_layout(); fig.savefig(OUT/"Figure2cd_distance_variability.pdf"); fig.savefig(OUT/"Figure2cd_distance_variability.png", dpi=150); plt.close(fig)
print(f"Fig2c/d: Exposed={len(E)} Shielded={len(S)} ; distance MWU p={p:.2e}")

# ===== Occupancy pie charts (by genomic feature) =====
# unified region definition: experimental TSS, 300 bp promoter, precedence
# promoter > 5'-UTR > CDS > intergenic, painted per base (mutually exclusive)
cls = pd.read_csv(B/"62_GO_KEGG_enrichment/tables/F1_classified_meth_sites_unified.tsv", sep="\t")
order = ['promoter','5UTR_approx','CDS_internal','intergenic']
lbls = ['Promoter','5′UTR','CDS','Intergenic']
colors = ['#A64B44','#C0803A','#3A6B8C','#3E7256']  # unified calm qualitative
fig, ax = plt.subplots(1, 2, figsize=(9, 4.2))
MOT_MOD = {'GCCGGC': '4mC', 'AAGCCCG': '4mC/6mA'}  # GCCGGC = 4mC; AAGCCCG = dual
for i,mot in enumerate(['GCCGGC','AAGCCCG']):
    sub = cls[cls.motif==mot]
    counts = [ (sub.category==c).sum() for c in order ]
    tot = sum(counts)
    wedges,_,autotexts = ax[i].pie(counts, labels=lbls, colors=colors,
                            autopct=lambda p:f"{p:.0f}%",
                            startangle=90, textprops={'fontsize':9})
    # percent labels sit inside the (saturated/dark) wedges → white for contrast
    for t in autotexts:
        t.set_color('white'); t.set_fontweight('bold')
    ax[i].set_title(f"{mot} {MOT_MOD[mot]} sites (n={tot})")
fig.suptitle("Genomic-feature distribution of methylation sites", y=1.02, fontsize=11)
fig.tight_layout()
fig.savefig(OUT/"SuppFig_occupancy_pie.pdf", bbox_inches="tight")
fig.savefig(OUT/"SuppFig_occupancy_pie.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("pie: GCCGGC", {c:(cls[(cls.motif=='GCCGGC')].category==c).sum() for c in order})
print("pie: AAGCCCG", {c:(cls[(cls.motif=='AAGCCCG')].category==c).sum() for c in order})
print("saved:", [p.name for p in sorted(OUT.glob('*.pdf'))])
