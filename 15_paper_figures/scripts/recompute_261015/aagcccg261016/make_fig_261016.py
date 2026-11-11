"""AAGCCCG 6mA occupancy: the recognition sequence extends 5' of HAAGCCCG.
Rendered with the house supplementary style (_suppstyle). SUPPFIG_OUT is set by the caller
so nothing is written under Writing/ while the manuscript build is live."""
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014')
from _suppstyle import apply_figure_style, figure_panel_letters, source_note, save, NAR_WIDTH_IN, RED, BLUE, GREY, TAN

D = os.path.dirname(os.path.abspath(__file__))
LAD  = pd.read_csv(f'{D}/ladder_allread_261016.tsv', sep='\t')
PR   = pd.read_csv(f'{D}/H1_positional_profile_allread_261016.tsv', sep='\t')
CTRL = pd.read_csv(f'{D}/minus3_specificity_controls_261016.tsv', sep='\t')
C4   = pd.read_csv(f'{D}/fiveprime_contrast_261016.tsv', sep='\t')

apply_figure_style()
fig, axes = plt.subplots(1, 3, figsize=(NAR_WIDTH_IN, 2.35))
TPS = ['T1', 'T2', 'T3']; TPC = ['#2b2b2b', '#767676', '#b4b4b4']  # FIG23: 時点は中立グレー段

# (a) extension ladder (horizontal: the definitions need room)
ax = axes[0]
order = ['AAGCCCG (all instances)', 'HAAGCCCG', 'G-N-HAAGCCCG  (-3=G)', 'C-N-G-N-HAAGCCCG  (-3=G,-5=C)']
n = [int(LAD[(LAD.definition == o) & (LAD.tp == 'T1')].n_sites.iloc[0]) for o in order]
lab = [f'AAGCCCG  {n[0]}', f'HAAGCCCG  {n[1]}', f'+ G\u208b\u2083  {n[2]}', f'+ C\u208b\u2085 G\u208b\u2083  {n[3]}']
y = np.arange(len(order))[::-1]; w = 0.26
for i, tp in enumerate(TPS):
    v = [LAD[(LAD.definition == o) & (LAD.tp == tp)].occupancy_pct.iloc[0] for o in order]
    ax.barh(y - (i - 1) * w, v, w, color=TPC[i], label=tp, edgecolor='none')
# 2026-10-24, D047 (consistency audit CA02): displayed position names follow the
# manuscript's declared rule, "1-based from the 5' end of each sequence", so for
# AAGCCCG the modified cytosine is C5 and the two adenines are A1 and A2 — the
# same labels Supplementary Figure 1e prints. This figure had called them C4 and
# A0, which collided with Supplementary Figure 1e's C4, a different and almost
# unmodified cytosine. Column names and data keys below (C4_H_occ, the `chans`
# strings) are storage keys in the .tsv files and are deliberately NOT renamed.
c4 = C4.C4_H_occ.mean()
ax.axvline(c4, color=RED, lw=1.0, ls='--', label='4mC at C5 (same sites)')
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=6.0)
ax.set_xlabel('6mA occupancy at A1 (%)  [site count at right of label]'); ax.set_xlim(0, 100)
ax.legend(frameon=False, loc='lower right', fontsize=6.0, handlelength=1.0)
# (b) positional profile
ax = axes[1]
for i, tp in enumerate(TPS):
    s = PR[PR.tp == tp].sort_values('offset')
    ax.plot(s.offset, -np.log10(s.q.clip(lower=1e-20)), 'o-', ms=2.2, lw=0.8, color=TPC[i], label=tp)
ax.axvspan(-0.5, 6.5, color='#eeeeee', zorder=0)
ax.text(3, 2.5, 'AAGCCCG', ha='center', fontsize=6.0, color='#555555')
ax.axhline(-np.log10(0.05), color=GREY, lw=0.7, ls=':', label='BH 5%')
ax.plot([-5, -3], [-0.9, -0.9], ls='none', marker='^', ms=3.2, color='#333333',
        clip_on=False, label='passes registered scan')
ax.set_xlabel('Offset from A1 (bp)')
ax.set_ylabel('\u2212log\u2081\u2080 q (best base/offset)')
ax.legend(frameon=False, fontsize=6.0, handlelength=1.0, loc='upper right', ncol=2, columnspacing=0.8)

# (c) -3 specificity control
ax = axes[2]
chans = ['6mA A0 | HAAGCCCG target', '6mA A0 | GAAGCCCG non-target', '4mC C4 | HAAGCCCG instances']
clab  = ['6mA, target\nHAAGCCCG', '6mA, non-target\nGAAGCCCG', '4mC at C5,\nsame instances']
ccol  = [BLUE, GREY, RED]
bases = list('ACGT'); xb = np.arange(4)
for i, (ch, cc) in enumerate(zip(chans, ccol)):
    s = CTRL[(CTRL.channel == ch) & (CTRL.tp == 'T1')].set_index('base_at_minus3')
    ax.plot(xb, [s.loc[b, 'occupancy_pct'] for b in bases], 'o-', ms=3, lw=1.1, color=cc, label=clab[i])
ax.set_xticks(xb); ax.set_xticklabels(bases)
ax.set_xlabel('Base at offset \u22123'); ax.set_ylabel('Occupancy (%), T1')
ax.set_ylim(0, 100); ax.legend(frameon=False, fontsize=6.0, handlelength=1.0, loc='center left')

fig.subplots_adjust(left=0.125, right=0.985, top=0.93, bottom=0.235, wspace=0.60)
figure_panel_letters(fig, axes, ['a', 'b', 'c'], dx=0.058)
source_note(fig, 'all-read calls, 9 libraries, sites with \u226530 reads; aagcccg261016/*_261016.tsv', fontsize=6.0)
save(fig, 'SuppFig_AAGCCCG_recognition_extension_261016', allow_over_data=('AAGCCCG', '4mC at C5 (same sites)', '4mC at C5,\nsame instances', 'T1', 'T2', 'T3', 'BH 5%', 'passes registered scan'))
