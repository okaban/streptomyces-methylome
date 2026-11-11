import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Global state mutations
plt.rcParams.update({
    'font.size': 8,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
})

META_GREY = '#6B7783'

def apply_figure_style(sizes=(8, 7, 6)):
    plt.rcParams.update({
        'font.size': sizes[0],
        'axes.labelsize': sizes[1],
        'xtick.labelsize': sizes[2],
        'ytick.labelsize': sizes[2],
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

def focal_palette(names, focal=None, focal_color="#c0392b", other="muted"):
    muted_colors = ['#7f8c8d', '#95a5a6', '#bdc3c7']
    result = {}
    muted_idx = 0
    for n in names:
        if n == focal:
            result[n] = focal_color
        else:
            result[n] = muted_colors[muted_idx % len(muted_colors)]
            muted_idx += 1
    return result

def set_frame(ax, style="open"):
    if style == "open":
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    elif style == "lb":
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

apply_figure_style(sizes=(8, 7, 6))
names = ["SC_RS35335 (NPPY)", "SC_RS28835 (NPPY)", "SC_RS17645 (HPPF)"]
vals = [(38.45, 39.41, 26.27), (14.04, 7.53, 9.07), (9.95, 1.56, 5.47)]
occ = (45.84, 47.79, 46.02)
x = np.arange(3)
pal = focal_palette(names, focal=names[2], focal_color="#c0392b", other="muted")
cols = dict(zip(names, pal)) if isinstance(pal, list) else pal

fig, ax = plt.subplots(figsize=(3.7, 2.6), dpi=300)
ax2 = ax.twinx()
for g, v in zip(names, vals):
    ax.plot(x, v, "o-", color=cols[g], lw=1.3, ms=4.2, zorder=3)
    ax.annotate(g, (2, v[2]), xytext=(4, 0), textcoords="offset points",
                color=cols[g], fontsize=6, va="center")
ax2.plot(x, occ, "s--", color=META_GREY, lw=1.1, ms=3.6, zorder=2)
ax.set_xticks(x)
ax.set_xticklabels(["T1", "T2", "T3"])
ax.set_xlim(-0.3, 4.6)
ax.set_ylabel("Transcript level (TPM)")
ax.set_ylim(0, 46)
ax2.set_ylabel("6mA occupancy (%)", color=META_GREY)
ax2.set_ylim(0, 62)
ax2.tick_params(axis="y", colors=META_GREY)
ax2.spines["right"].set_color(META_GREY)
ax2.spines["top"].set_visible(False)
ax2.annotate("6mA occupancy at AAGCCCG", (0.05, occ[0]), xytext=(0, -12),
             textcoords="offset points", color=META_GREY, fontsize=6,
             ha="left", va="top")
ax.set_title("Candidates are transcribed; occupancy does not track them", loc="left")
set_frame(ax, "open")
fig.tight_layout()
out = "SuppFigure17.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
print("saved:", out)