"""
Shared utilities for paper figure generation.
S. coelicolor A3(2) M145 Nanopore methylome analysis.
"""

import sys
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO

# ── Base paths ──────────────────────────────────────────────────────────────
# Session-agnostic resolution: prefer first existing path
def _resolve_base():
    import os
    candidates = [
        '/sessions/busy-eloquent-bell/mnt/rna-seq',
        '/sessions/zen-modest-davinci/mnt/rna-seq',
    ]
    for c in candidates:
        if os.path.isdir(c):
            return Path(c)
    # fallback: use the directory two levels above this file (scripts/ → 15_paper_figures/ → rna-seq/)
    return Path(__file__).resolve().parents[2]

BASE = _resolve_base()
EPIGENOME = BASE / '11_epigenome_integration' / 'analysis'
METHYL = BASE
FIG_DIR = BASE / '15_paper_figures' / 'figures' / 'main'
FIG_SUP_DIR = BASE / '15_paper_figures' / 'figures' / 'supplementary'
TABLE_SUP_DIR = BASE / '15_paper_figures' / 'tables' / 'supplementary'

# ── Colors — UNIFIED PALETTE "Rich & Calm" (approved 2026-07) ────────────────
# One hex per semantic role across ALL main + supplementary figures.
#   4mC (4mC) = rich muted red, 6mA = rich muted blue, dual = purple.
#   Bars use a grayscale ramp keyed to SuppFig15 (author's hand-made reference).
#   Motif sequence-logo colours (MOTIF_COLORS) are intentionally left colourful.
# v2 (2026-09-22): hues kept, RELATIVE LUMINANCE re-spaced. Measured on v1:
#   4mC vs 6mA  greyscale contrast 1.01  -> identical in black-and-white print
#   6mA vs dual deuteranope distance 0.036 -> identical to a red-green CVD reader,
#                and these two sit side by side in the KEGG figure.
# v1 put the five categorical hues at almost the same luminance (0.118-0.143), so
# they were separated by chroma alone -- the two failure modes Crameri, Shephard &
# Heron (2020, Nat Commun 11:5444) name. v2 measures: 4mC-6mA 1.92 / 0.335,
# 6mA-dual 1.83 / 0.393, 4mC-dual 3.50 / 0.699, core-arm 2.22 / 0.350.
# NOTE the one pair luminance cannot fix: COL_4mC vs COL_CORE (red vs green) stays
# at 1.21 / 0.053 -- never use them as opposing categories in the SAME panel.
COL_4mC = '#8A3D37'   # deep muted red   — 4mC modification mark      (L 0.09)
COL_6mA = '#68849D'   # mid muted blue   — 6mA modification mark      (L 0.22)
COL_BOTH = '#68849D'  # mid muted blue   — AAGCCCG, the dual 4mC+6mA motif (L 0.22)
#   2026-09-26 author decision: AAGCCCG is BLUE and GCCGGC is RED across every figure.
#   The lavender used since 09-22 read as grey on screen. 6mA-as-a-mark shares this blue,
#   which is safe because AAGCCCG is the only 6mA-carrying motif shown anywhere.
# superseded v1, kept so old outputs can be identified: 4mC #A64B44, 6mA #3A6B8C,
# dual #6E5495, core #3E7256, arm #C0803A
COL_DUAL = COL_BOTH   # alias

# 2026-09-26 author decision: the DUAL-TARGETED gene set is purple — red (GCCGGC)
# plus blue (AAGCCCG) mixing to purple is the mnemonic. It replaces the neutral
# dark grey used since 09-22. Hue sits between the two mark colours; luminance is
# spaced away from both so the three columns stay separable in greyscale.
COL_DUAL_SET = '#6E4B8E'   # medium purple — genes near BOTH motifs (a gene set)

# Grayscale bar ramp (SuppFig15 reference: bar #293039, grid #EAEAEA)
COL_BAR_DARK = '#293039'   # T3 / late / single-category bars
COL_BAR_MID = '#6B7783'    # T2 / mid
COL_BAR_LIGHT = '#AAB3BB'  # T1 / early
COL_GRID = '#EAEAEA'       # gridlines
BAR_RAMP = [COL_BAR_LIGHT, COL_BAR_MID, COL_BAR_DARK]  # T1→T2→T3

# 2026-09-26: the three timepoint TRACKS in Figure 2a were drawn with BAR_RAMP,
# three greys that are hard to tell apart as overlaid lines. Bars can carry a grey
# ramp (they are separated in space); lines cannot. TIME_RAMP is a single-hue
# sequential blue with wide luminance spacing, paired with distinct linestyles so
# the tracks survive greyscale printing too. Blue is free here: this panel plots
# GCCGGC only, so it cannot be confused with the 6mA mark colour.
TIME_RAMP = ['#9FC3E0', '#4C86B8', '#1F3E5A']   # T1 → T2 → T3 (lines)
TIME_STYLE = ['-', '--', '-']

COL_GRAY = '#9AA7B0'  # neutral grey (Shielded / unassigned)
COL_DARK = '#22282E'  # axis / text near-black
COL_GREEN = '#3A6B51'   # deep muted green — activation / core
COL_ORANGE = '#C69062'  # calm amber — repression / arm

# Gatekeeper / categorical accents (muted, deduplicated)
COL_ACTIVATION = '#3A6B51'    # deep muted green — activation bloc / core
COL_REPRESSION = '#C69062'    # calm amber — repression bloc / arm
COL_CORE = '#3A6B51'          # chromosomal core (= activation)
COL_ARM = '#C69062'           # chromosomal arm (= repression)
COL_EXPOSED = '#895981'       # muted plum — exposed TFs (distinct from dual)
COL_SHIELDED = '#9AA7B0'      # neutral grey — shielded TFs
COL_ARTIFACT = '#C69062'      # amber — Simpson's artifact (= repression tone)
COL_SIGNAL = '#68849D'        # blue — genuine signal (= 6mA tone)

# ── Timepoint labels ────────────────────────────────────────────────────
TP_LABELS = ['T1 (12 h)', 'T2 (24 h)', 'T3 (50 h)']
TP_LABELS_NL = ['T1\n(12 h)', 'T2\n(24 h)', 'T3\n(50 h)']

MOTIF_COLORS = {
    'CCGG': '#E53935',
    'GGCCGG': '#C62828',
    'TGGCCGGC': '#B71C1C',
    'AAGCCCG': '#7E57C2',
    'GATC': '#1565C0',
    'CCGKCA': '#0D47A1',
    'unassigned': '#B0BEC5',
    'other': '#78909C',
}

# ── Style ───────────────────────────────────────────────────────────────────
STYLE = {
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'pdf.fonttype': 42,      # TrueType (editable in Illustrator)
    'ps.fonttype': 42,
    'svg.fonttype': 'none',  # Text as text, not paths
}


def apply_style():
    """Apply publication-quality matplotlib style."""
    plt.rcParams.update(STYLE)
    warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')


# ── Unified role-based font ladder (NAR full-width 6.85") ─────────────────────
# One size per role across ALL figures, so panels read consistently at 174 mm.
FONT = {
    'panel_title': 8,    # ax.set_title
    'axis_label': 7.5,   # ax.set_xlabel / set_ylabel
    'tick': 7,           # tick labels
    'legend': 6.5,       # legend entries
    'annot': 6.5,        # in-panel annotations / stat text
    'panel_letter': 14,  # a b c d (bold, top-left)
}

# rcParams flavour of the ladder — apply per script for the compressed panels.
STYLE_UNIFIED = {
    'font.family': 'Arial',
    'font.size': FONT['annot'],
    'axes.titlesize': FONT['panel_title'],
    'axes.labelsize': FONT['axis_label'],
    'xtick.labelsize': FONT['tick'],
    'ytick.labelsize': FONT['tick'],
    'legend.fontsize': FONT['legend'],
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.edgecolor': COL_DARK,
    'axes.labelcolor': COL_DARK,
    'text.color': COL_DARK,
    'xtick.color': COL_DARK,
    'ytick.color': COL_DARK,
    'grid.color': COL_GRID,
    'grid.linewidth': 0.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
}


def apply_unified_style():
    """Apply the unified role-based font ladder + palette-consistent axis colours.

    Use in place of apply_style() for figures being brought to the unified
    2026-07 style (rich/calm palette, grayscale bars, lowercase panel letters).
    Font sizes are sized for NAR full-width (174 mm) compressed panels.
    """
    plt.rcParams.update(STYLE_UNIFIED)
    warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')


def style_axes(ax, grid=False, grid_axis='both'):
    """Normalise a single Axes to the unified look: hide top/right spines,
    colour the remaining spines/ticks near-black, optional light gridlines."""
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_edgecolor(COL_DARK)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=COL_DARK, labelcolor=COL_DARK)
    if grid:
        ax.grid(True, axis=grid_axis, color=COL_GRID, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
    return ax


def add_panel_label(ax, label, dx=-26, dy=12, fontsize=13, lower=True,
                    x=None, y=None):
    """Add a bold panel label (a, b, c, ...) at a UNIFORM visual position.

    The label is anchored to the axes' top-left CORNER and offset by a constant
    number of typographic points (dx, dy), so every panel across every figure
    gets an identically-sized, identically-placed letter regardless of how wide
    that panel's y-axis label/tick text happens to be. This is the fix for the
    previous per-panel `x`/`y` axes-fraction overrides, which made the letters
    drift and appear different sizes. NAR convention is lowercase (matches the
    figure legends' "(a) ... (b) ..."), so `lower=True` by default.
    (x, y kwargs are accepted for backward-compatibility but ignored.)
    """
    lab = label.lower() if lower else label.upper()
    # Anchor in the figure's OUTER margin at the panel's top-left, so the letter
    # clears the y-axis label/ticks no matter how wide they are, and lands in the
    # same visual spot for every panel. We use the axes' tightbbox (includes the
    # y-label) left edge, then place the letter a small constant gap to its left
    # and at the axes' top. Falls back to an axes-fraction offset if the renderer
    # can't supply a tightbbox yet.
    fig = ax.figure
    try:
        fig.canvas.draw()  # ensure a renderer exists for tightbbox
        rend = fig.canvas.get_renderer()
        tb = ax.get_tightbbox(rend).transformed(fig.transFigure.inverted())
        ax_pos = ax.get_position()
        xf = max(0.002, tb.x0 - 0.006)      # just left of the widest left-side text
        # Sit at the axes top, but never below the top of the y-label text, and
        # add a small constant lift so the letter clears a tall wrapped y-label.
        yf = max(ax_pos.y1, tb.y1) + 0.012
        fig.text(xf, yf, lab, fontsize=fontsize, fontweight='bold',
                 va='bottom', ha='left')
    except Exception:
        ax.annotate(lab, xy=(0, 1), xycoords='axes fraction',
                    xytext=(dx, dy), textcoords='offset points',
                    fontsize=fontsize, fontweight='bold', va='bottom', ha='left',
                    annotation_clip=False)


# ── NAR figure width limits (mm → inch) ──────────────────────────────────────
# Nucleic Acids Research column widths: single 86 mm, intermediate 120 mm,
# full 174 mm. A figure MUST NOT exceed its width class or the journal
# down-scales it, crushing fonts. save_figure() below enforces the ceiling.
NAR_WIDTH_MM = {'single': 86, 'intermediate': 120, 'full': 174}


def save_figure(fig, path, formats=('pdf', 'svg'), width_class='full', clamp=True):
    """Save figure at a NAR-compliant physical width (no silent tight re-expansion).

    IMPORTANT — why this does NOT use bbox_inches='tight':
    matplotlib's `tight` re-crops the canvas to the CONTENT bounding box and so
    IGNORES the declared figsize. A figure declared at 174 mm whose content spans
    more is re-expanded past the NAR full-width limit (that is how Figure1/2/8
    came out at 7.1-7.4"). To honour the declared/clamped width we save with the
    figsize as-is and clear rcParams['savefig.bbox'] for the duration of the save
    (passing bbox_inches=None alone is a no-op — print_figure defers None to the
    rcParam, which apply_style() sets to 'tight').

    width_class : 'single' | 'intermediate' | 'full' — NAR column ceiling.
    clamp       : if True and the figure is wider than the ceiling, scale the
                  whole figure down (aspect preserved) so width == ceiling.
                  Figures already within the ceiling are left untouched.

    If content genuinely overruns the clamped width (axis-external labels, wide
    legends, suptitles placed outside the axes), it will clip — that is the
    signal to tighten the layout in-script, not to let tight hide it. For a
    guaranteed-no-clip PNG fallback see save_figure_raster().
    """
    ceiling = mm_to_inch(NAR_WIDTH_MM[width_class])
    w, h = fig.get_size_inches()
    if clamp and w > ceiling + 1e-6:
        fig.set_size_inches(ceiling, h * ceiling / w)
        w, h = fig.get_size_inches()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    prev_bbox = plt.rcParams['savefig.bbox']
    plt.rcParams['savefig.bbox'] = None
    try:
        for fmt in formats:
            out = path.with_suffix(f'.{fmt}')
            fig.savefig(out, format=fmt, dpi=300, bbox_inches=None, pad_inches=0)
            print(f'  Saved: {out}  (width={w:.2f}" <= {ceiling:.2f}" [{width_class}])')
    finally:
        plt.rcParams['savefig.bbox'] = prev_bbox
    plt.close(fig)


def save_figure_raster(fig, path, width_class='full', dpi=300):
    """Guaranteed-no-clip PNG save at NAR width (raster-fit fallback).

    Renders with bbox_inches='tight' (full content preserved, no clipping), then
    if the rendered width exceeds the NAR ceiling downscales the RASTER to exactly
    the ceiling width (Lanczos, aspect preserved; never upscales). Use for figures
    whose layout places content outside the axes box and cannot be quickly
    reflowed. PNG only — for vector output use save_figure() and fix the layout.
    """
    from PIL import Image
    target_px = round(mm_to_inch(NAR_WIDTH_MM[width_class]) * dpi)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = path.with_suffix('.png')
    tmp = str(out) + '.__tight.png'
    prev_bbox = plt.rcParams['savefig.bbox']
    plt.rcParams['savefig.bbox'] = 'tight'
    try:
        fig.savefig(tmp, dpi=dpi, bbox_inches='tight')
    finally:
        plt.rcParams['savefig.bbox'] = prev_bbox
    im = Image.open(tmp)
    w, h = im.size
    scaled = w > target_px
    if scaled:
        im = im.resize((target_px, round(h * target_px / w)), Image.LANCZOS)
    im.convert('RGB').save(out, dpi=(dpi, dpi))
    Path(tmp).unlink(missing_ok=True)
    print(f"  Saved: {out}  (width={im.size[0] / dpi:.2f}\" <= "
          f"{mm_to_inch(NAR_WIDTH_MM[width_class]):.2f}\" [{width_class}], "
          f"{'downscaled' if scaled else 'native'})")
    plt.close(fig)


# ── Data loaders ────────────────────────────────────────────────────────────

def load_methylation_census():
    """Load 4mC and 6mA census data.

    Returns:
        df_4mc, df_6ma: DataFrames with columns:
            chrom, position, strand, mod_type, timepoint, frequency,
            sequence, center_base, at_AAGCCCG, final_motif
    """
    path_4mc = EPIGENOME / '23_expanded_motif_search' / '4mC_final_census.csv'
    path_6ma = EPIGENOME / '23_expanded_motif_search' / '6mA_final_census.csv'
    df_4mc = pd.read_csv(path_4mc)
    df_6ma = pd.read_csv(path_6ma)
    return df_4mc, df_6ma


def load_methylation_hc_all():
    """Load all HC methylation sites with timepoint info (including overlaps).

    Returns:
        DataFrame with columns:
            chrom, position, strand, mod_type, timepoint,
            n_reps, total_coverage, weighted_mod_freq, unweighted_mod_freq
    """
    path = EPIGENOME / '01_integration' / 'high_confidence_sites_weighted.csv'
    df = pd.read_csv(path)
    print(f'  HC all sites: {len(df)} rows '
          f'(4mC: {(df.mod_type=="4mC").sum()}, 6mA: {(df.mod_type=="6mA").sum()})')
    return df


def load_methylation_unique_positions():
    """Load unique methylation positions (collapsed across timepoints).

    Returns:
        df_unique: DataFrame with columns: position, mod_type, final_motif, strand
    """
    df_4mc, df_6ma = load_methylation_census()
    df_all = pd.concat([df_4mc, df_6ma], ignore_index=True)
    df_unique = (df_all.drop_duplicates(subset=['position', 'mod_type'])
                 [['chrom', 'position', 'strand', 'mod_type', 'final_motif']]
                 .sort_values('position')
                 .reset_index(drop=True))
    return df_unique


def load_tss_jeong2016():
    """Load TSS data, filtered to Jeong2016 dRNA-seq experimental TSS only.

    Returns:
        DataFrame with 2,703 genes (Jeong2016 experimental TSS)
    """
    path = EPIGENOME / '18_tss_analyses' / 'comprehensive_tss_table.csv'
    df = pd.read_csv(path)
    df_j = df[df['tss_source'] == 'Jeong2016_dRNA-seq'].copy()
    print(f'  Jeong2016 TSS loaded: {len(df_j)} genes')
    return df_j


def load_reference_genome():
    """Load reference genome as SeqRecord.

    Returns:
        SeqRecord (NC_003888.3, 8,667,507 bp)
    """
    path = METHYL / 'data' / 'ref.fa'
    record = SeqIO.read(path, 'fasta')
    print(f'  Genome: {record.id}, {len(record.seq):,} bp')
    return record


def load_reference_gbk():
    """Load reference GenBank file for gene annotations.

    Returns:
        SeqRecord with features
    """
    path = METHYL / 'data' / 'ref.gbk'
    record = SeqIO.read(path, 'genbank')
    print(f'  GenBank: {record.id}, {len(record.features)} features')
    return record


def parse_meme_pwm(meme_file):
    """Parse MEME meme.txt file and extract position weight matrices.

    Args:
        meme_file: Path to meme.txt

    Returns:
        dict: {motif_name: {'pwm': pd.DataFrame, 'nsites': int, 'evalue': str, 'width': int}}
    """
    meme_file = Path(meme_file)
    text = meme_file.read_text()

    motifs = {}
    # Find all motif blocks
    pattern = (r'Motif (\S+) MEME-(\d+) position-specific probability matrix\n'
               r'-+\n'
               r'letter-probability matrix: alength= 4 w= (\d+) nsites= (\d+) E= (\S+)\s*\n'
               r'((?:\s*[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*\n)+)')
    for m in re.finditer(pattern, text):
        name = m.group(1)
        rank = int(m.group(2))
        width = int(m.group(3))
        nsites = int(m.group(4))
        evalue = m.group(5)
        rows = []
        for line in m.group(6).strip().split('\n'):
            vals = [float(x) for x in line.split()]
            rows.append(vals)
        pwm = pd.DataFrame(rows, columns=['A', 'C', 'G', 'T'])
        motifs[f'MEME-{rank}'] = {
            'name': name,
            'pwm': pwm,
            'nsites': nsites,
            'evalue': evalue,
            'width': width,
        }

    print(f'  Parsed {len(motifs)} motifs from {meme_file.name}')
    return motifs


def load_tf_binding_sites():
    """Load TF binding site master table.

    Returns:
        DataFrame with BS_source-based tier classification added.
    """
    path = (BASE / '13_TF_binding-site' / 'analysis' /
            '01_master_TF_list_260206_v1' /
            'master_TF_binding_sites_M145_ARCHIVED_260224.tsv')
    df = pd.read_csv(path, sep='\t')

    # Tier classification
    tier1_sources = ['RegPrecise', 'ZorroAranda2022_Curated_Strong', 'FIMO_curated_motif']
    df['tier'] = df['BS_source'].apply(lambda x: 1 if x in tier1_sources else 2)

    n1 = (df['tier'] == 1).sum()
    n2 = (df['tier'] == 2).sum()
    print(f'  TF BS loaded: {len(df)} total (Tier 1: {n1}, Tier 2: {n2})')
    return df


def load_mtase_expression():
    """Load MTase genes with expression data.

    Returns:
        DataFrame with 22 MTase genes and DESeq2 expression data.
    """
    path = EPIGENOME / '11_rm_system_identification' / 'mtase_genes_with_expression.csv'
    df = pd.read_csv(path)
    print(f'  MTase genes: {len(df)}')
    return df


def load_rebase_conservation():
    """Load REBASE motif conservation matrix.

    Returns:
        DataFrame with conservation rates per motif.
    """
    path = EPIGENOME / '21_genuswide_motif_conservation' / 'rebase_motif_conservation_matrix.csv'
    df = pd.read_csv(path)
    return df


def load_5mc_4mc_data():
    """Load 5mC vs 4mC comparison data at CCGG sites.

    Returns:
        df_detailed, df_summary: detailed per-sample and genome-wide summary
    """
    base = EPIGENOME / '24_5mC_vs_4mC_CCGG' / 'tables'
    df_det = pd.read_csv(base / 'CCGG_4mC_vs_5mC_detailed.csv')
    df_sum = pd.read_csv(base / 'genome_wide_4mC_vs_5mC_summary.csv')
    return df_det, df_sum


def load_motif_summary():
    """Load comprehensive motif summary.

    Returns:
        DataFrame with motif, mod, sites, m145_oe, rebase, etc.
    """
    path = EPIGENOME / '23_expanded_motif_search' / 'comprehensive_motif_summary.csv'
    df = pd.read_csv(path)
    return df


def load_integrated_expression():
    """Load integrated methylation-expression table.

    Returns:
        DataFrame with 8,083 genes.
    """
    path = (EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
            '01_integration' / 'integrated_methyl_expression.csv')
    df = pd.read_csv(path)
    return df


# ── Gatekeeper model data loaders ──────────────────────────────────────────

def load_spatial_profile():
    """Load TSS methylation gradient spatial profile data (analysis/48).

    Returns:
        DataFrame with columns: site_type, timepoint, gene_category,
        bin_center_bp, density_per_kb_per_gene, ci_low, ci_high, n_genes, n_sites
    """
    path = EPIGENOME / '48_TSS_methylation_gradient' / 'tables' / 'spatial_profile_data.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Spatial profile: {len(df)} rows')
    return df


def load_all_genes_features():
    """Load shielded/exposed classification for all regulatory genes (analysis/52).

    Returns:
        DataFrame with 1,017 regulatory genes and features:
        locus_tag, gene_name, old_locus_tag, product, tf_family, start, end, strand,
        tss, region, is_exposed, baseMean, nearest_methyl_distance,
        n_methyl_sites_2kb, gene_length, LFC_T2vsT1, LFC_T3vsT1, n_FIMO_hits
    """
    path = EPIGENOME / '52_shielded_exposed_boundary' / 'tables' / 'all_genes_features.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  All genes features: {len(df)} genes ({df["is_exposed"].sum()} exposed)')
    return df


def load_exposed_regulators():
    """Load 62 exposed TF full annotation (analysis/51).

    Returns:
        DataFrame with 62 exposed TFs and comprehensive annotation.
    """
    path = EPIGENOME / '51_exposed_regulators_characteristics' / 'tables' / 'exposed_regulators_full_table.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Exposed regulators: {len(df)} TFs')
    return df


def load_coexpression_matrix():
    """Load 62×62 co-expression matrix (analysis/55).

    Returns:
        DataFrame (62×62) with locus_tags as both index and columns.
    """
    path = EPIGENOME / '55_exposed_regulatory_module' / 'tables' / 'coexpression_matrix.tsv'
    df = pd.read_csv(path, sep='\t', index_col=0)
    print(f'  Co-expression matrix: {df.shape}')
    return df


def load_module_summary():
    """Load regulatory module summary (analysis/55).

    Returns:
        DataFrame with module IDs, members, evidence scores.
    """
    path = EPIGENOME / '55_exposed_regulatory_module' / 'tables' / 'module_summary.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Modules: {len(df)}')
    return df


def load_temporal_classification():
    """Load temporal dynamics classification for 62 exposed TFs (analysis/57).

    Returns:
        DataFrame with bloc (activation/repression/unassigned), phase_ratio,
        temporal_class, T1/T2/T3 z-scores, LFC values.
    """
    path = EPIGENOME / '57_temporal_dynamics_exposed_TF' / 'tables' / 'temporal_classification.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Temporal classification: {len(df)} TFs '
          f'(act={sum(df["bloc"]=="activation")}, rep={sum(df["bloc"]=="repression")})')
    return df


def load_bloc_family_enrichment():
    """Load TF family enrichment by bloc (analysis/58).

    Returns:
        DataFrame with tf_family, activation_count, repression_count,
        odds_ratio, p_value, enriched_in, functional_category, p_adj.
    """
    path = EPIGENOME / '58_exposed_TF_functional_prediction' / 'tables' / 'bloc_family_enrichment.tsv'
    df = pd.read_csv(path, sep='\t')
    df = df.dropna(subset=['tf_family'])
    print(f'  Bloc family enrichment: {len(df)} families')
    return df


def load_bloc_comparison():
    """Load bloc-level summary statistics (analysis/57).

    Returns:
        DataFrame with activation/repression aggregate metrics.
    """
    path = EPIGENOME / '57_temporal_dynamics_exposed_TF' / 'tables' / 'bloc_comparison.tsv'
    df = pd.read_csv(path, sep='\t')
    return df


def load_geographic_redistribution():
    """Load GCCGGC geographic redistribution data (analysis/37).

    Returns:
        df_summary: core/arm distribution by timepoint
        df_sites: individual GCCGGC site positions by timepoint
    """
    # 2026-09-21 (BLOCKER-0): the 37_ table is position-DEDUPLICATED across
    # timepoints (first-appearance: 1,289/407/21 = sites NEW at each TP). Read the
    # canonical per-timepoint file instead (1,289/1,595/1,073) via the shared
    # loader in 90_per_timepoint_census_audit/canonical_sites.py.
    import importlib.util as _ilu
    _spec = _ilu.spec_from_file_location(
        'canonical_sites', EPIGENOME / '90_per_timepoint_census_audit' / 'canonical_sites.py')
    _cs = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_cs)
    df_sites = _cs.gccggc_by_timepoint()
    df_summary = (df_sites.groupby('timepoint')
                  .agg(n_sites=('position', 'size'),
                       core_fraction=('region', lambda r: (r == 'core').mean()))
                  .reset_index())
    print(f'  GCCGGC sites (canonical per-timepoint): T1={len(df_sites[df_sites.timepoint=="T1"])}, '
          f'T2={len(df_sites[df_sites.timepoint=="T2"])}, '
          f'T3={len(df_sites[df_sites.timepoint=="T3"])}')
    return df_summary, df_sites


def load_gene_methylation_transitions(motif='GCCGGC'):
    """Load gene-level methylation transition data (analysis/42 or 44).

    Args:
        motif: 'GCCGGC' (analysis/42) or 'AAGCCCG' (analysis/44)

    Returns:
        df_transitions: per-gene methylation transition and expression data
        df_geographic: geographic stratified test results
    """
    if motif == 'GCCGGC':
        base = EPIGENOME / '42_GCCGGC_temporal_derepression' / 'tables'
    elif motif == 'AAGCCCG':
        base = EPIGENOME / '44_AAGCCCG_temporal_causality' / 'tables'
    else:
        raise ValueError(f'Unknown motif: {motif}')
    df_trans = pd.read_csv(base / 'gene_methylation_transitions.tsv', sep='\t')
    df_geo = pd.read_csv(base / 'geographic_stratified_tests.tsv', sep='\t')
    print(f'  {motif} transitions: {len(df_trans)} genes')
    return df_trans, df_geo


def load_roc_analysis(analysis='boundary'):
    """Load ROC analysis results (analysis/52 or 53).

    Args:
        analysis: 'boundary' (analysis/52) or 'sequence' (analysis/53)

    Returns:
        DataFrame with feature, AUC, CI, threshold, etc.
    """
    if analysis == 'boundary':
        path = EPIGENOME / '52_shielded_exposed_boundary' / 'tables' / 'ROC_analysis.tsv'
    elif analysis == 'sequence':
        path = EPIGENOME / '53_TSS_sequence_determinants' / 'tables' / 'ROC_analysis.tsv'
    else:
        raise ValueError(f'Unknown analysis: {analysis}')
    df = pd.read_csv(path, sep='\t')
    print(f'  ROC ({analysis}): {len(df)} features')
    return df


def load_expression_quintile():
    """Load expression quintile analysis (analysis/52).

    Returns:
        DataFrame with quintile, n_genes, n_exposed, frac_exposed, etc.
    """
    path = EPIGENOME / '52_shielded_exposed_boundary' / 'tables' / 'expression_quintile.tsv'
    df = pd.read_csv(path, sep='\t')
    return df


def load_tfbs_spatial_profile():
    """Load TF BS methylation spatial profile (analysis/49).

    Returns:
        DataFrame with methylation density around TF binding sites.
    """
    path = EPIGENOME / '49_TF_BS_methylation_protection' / 'tables' / 'spatial_profile_data.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  TFBS spatial profile: {len(df)} rows')
    return df


def load_tcs_pairs():
    """Load TCS pair analysis (analysis/55).

    Returns:
        DataFrame with 7 TCS pairs showing exposed/shielded asymmetry.
    """
    path = EPIGENOME / '55_exposed_regulatory_module' / 'tables' / 'TCS_pairs_analysis.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  TCS pairs: {len(df)}')
    return df


def load_neighborhood_results():
    """Load neighborhood effect analysis (analysis/56).

    Returns:
        df_permutation: permutation test results
        df_stats: statistical test results
    """
    base = EPIGENOME / '56_exposed_TF_neighborhood' / 'tables'
    df_perm = pd.read_csv(base / 'permutation_results.tsv', sep='\t')
    df_stats = pd.read_csv(base / 'statistical_tests.tsv', sep='\t')
    return df_perm, df_stats


def load_exposed_methylation_status():
    """Load AAGCCCG/GCCGGC methylation status for 62 exposed TFs (analysis/58_causal).

    Returns:
        DataFrame with per-TF AAGCCCG and GCCGGC methylation counts.
    """
    path = EPIGENOME / '58_AAGCCCG_exposed_TF_causal' / 'tables' / 'exposed_TF_methylation_status.tsv'
    df = pd.read_csv(path, sep='\t')
    print(f'  Exposed methylation status: {len(df)} TFs')
    return df


# ── Utility functions ───────────────────────────────────────────────────────

def mm_to_inch(mm):
    """Convert millimeters to inches."""
    return mm / 25.4


def format_pvalue(p):
    """Format p-value for display."""
    if pd.isna(p):
        return 'N/A'
    if p < 0.001:
        return f'{p:.1e}'
    if p < 0.01:
        return f'{p:.3f}'
    if p < 0.05:
        return f'{p:.2f}'
    return f'{p:.2f}'


if __name__ == '__main__':
    apply_style()
    print('=== Testing shared utilities ===')
    print()

    print('Loading methylation census...')
    df_4mc, df_6ma = load_methylation_census()
    print(f'  4mC: {len(df_4mc)} rows, 6mA: {len(df_6ma)} rows')

    print('Loading unique positions...')
    df_uniq = load_methylation_unique_positions()
    print(f'  Unique: {len(df_uniq)} positions')

    print('Loading Jeong2016 TSS...')
    df_tss = load_tss_jeong2016()

    print('Loading reference genome...')
    ref = load_reference_genome()

    print('Loading MEME 4mC...')
    m4 = parse_meme_pwm(EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
                         '07_motif_analysis' / 'meme_4mC' / 'meme.txt')
    for k, v in m4.items():
        print(f'  {k}: {v["name"]} (w={v["width"]}, n={v["nsites"]}, E={v["evalue"]})')

    print('Loading MEME 6mA...')
    m6 = parse_meme_pwm(EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
                         '07_motif_analysis' / 'meme_6mA' / 'meme.txt')
    for k, v in m6.items():
        print(f'  {k}: {v["name"]} (w={v["width"]}, n={v["nsites"]}, E={v["evalue"]})')

    print()
    print('=== All utilities OK ===')

# ── Text-collision gate ───────────────────────────────────────────────────────
# Rule (2026-09-22): text may sit INSIDE the axes only when it identifies an
# adjacent mark (direct labelling, which Wong's "Axes, ticks and grids" column
# frames as keeping navigational elements unobtrusive so the DATA keeps visual
# priority, and which the CPT visual-communication tutorial prefers over a
# legend box). Narration — n, statistics, null levels — goes outside the data
# region: title, axis label, margin, or caption. Nothing may overlap data marks
# or other text, and nothing may fall outside the canvas.
# Every figure script must call assert_no_text_collisions() before saving;
# Writing/check_figures.py fails a script that does not.

def _renderer(fig):
    # Draw twice: after a save the first draw can still report stale extents for
    # axis labels, which made the canvas-overflow test fire on labels that were
    # comfortably inside. Two draws settle the layout deterministically.
    fig.canvas.draw()
    fig.canvas.draw()
    if hasattr(fig.canvas, "get_renderer"):
        return fig.canvas.get_renderer()
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    return FigureCanvasAgg(fig).get_renderer()


def text_collisions(fig, data_pad=1.0):
    """Report (overlapping text pairs, text outside canvas, text over data marks)."""
    import matplotlib as mpl
    r = _renderer(fig)
    ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    # Structural text — titles, axis labels, tick labels — is positioned by
    # matplotlib relative to its own axes. Two of them abutting is a layout
    # artefact, not a defect; only collisions that involve free-floating text
    # (annotations, legend entries) or cross different axes are reported.
    def _struct(ax):
        # loc="left"/"right" titles are SEPARATE artists (_left_title/_right_title),
        # not ax.title. Omitting them made the checker report a left-aligned panel
        # title as an annotation sitting on the data.
        out = {ax.title, ax.xaxis.label, ax.yaxis.label}
        for attr in ("_left_title", "_right_title"):
            t = getattr(ax, attr, None)
            if t is not None:
                out.add(t)
        return ticks.get(ax, set()) | out

    structural = set()
    own = {}
    for ax in fig.axes:
        st = _struct(ax)
        structural |= st
        for t in st:
            own[t] = ax
    # Legend entries are packed by matplotlib's offsetbox; their reported extents
    # can nominally intersect even when the drawn rows are clearly separated, so
    # legend-internal pairs are excluded. Legend text that covers DATA or leaves
    # the canvas is still reported below.
    legend_texts = set()
    for ax in list(fig.axes) + [fig]:
        lg = getattr(ax, "legend_", None)
        if lg is not None:
            legend_texts |= set(lg.get_texts())
    pairs = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            if not ba.overlaps(bb):
                continue
            # Report a pair only when at least one member is free-floating text
            # (an annotation or value label the script placed itself). Titles,
            # axis labels and tick labels are positioned by matplotlib, which
            # does not draw them on top of each other; their reported extents
            # can nominally intersect without any visible collision.
            if a in structural and b in structural:
                continue
            if a in legend_texts and b in legend_texts:
                continue
            pairs.append((a.get_text()[:28], b.get_text()[:28]))
    outside = [t.get_text()[:28] for t, b in texts
               if b.x0 < -1 or b.y0 < -1 or b.x1 > fig.bbox.x1 + 1 or b.y1 > fig.bbox.y1 + 1]
    over_data = []
    for ax in fig.axes:
        # Test against the actual drawn positions, not each artist's bounding box:
        # a scatter cloud's bbox is as wide as the axes, so a bbox test silently
        # passed text sitting on top of points.
        pts = []
        for c in ax.collections:
            try:
                off = c.get_offsets()
                if off is not None and len(off):
                    pts.extend(ax.transData.transform(off))
            except Exception:
                pass
        for l in ax.lines:
            try:
                xy = l.get_xydata()
                if xy is not None and len(xy):
                    pts.extend(ax.transData.transform(xy))
            except Exception:
                pass
        boxes = []
        for m in ax.patches:
            try:
                bb = m.get_window_extent(r)
                if bb.width > 0 and bb.height > 0 and bb.width < ax.bbox.width * 0.98:
                    boxes.append(bb)
            except Exception:
                pass
        for t, bt in texts:
            if t in structural:          # never an annotation-over-data finding
                continue
            hit_patch = any(bt.overlaps(bb) for bb in boxes)
            hit_point = any(bt.x0 - 1 <= px <= bt.x1 + 1 and bt.y0 - 1 <= py <= bt.y1 + 1
                            for px, py in pts)
            if hit_patch or hit_point:
                over_data.append(t.get_text()[:28])
    return pairs, outside, sorted(set(over_data))


def assert_no_text_collisions(fig, name="figure", allow_over_data=()):
    """Raise unless the figure is free of text collisions. allow_over_data lists
    label texts that are deliberately placed on a mark (e.g. a value inside a bar,
    which must then meet the contrast floor)."""
    pairs, outside, over = text_collisions(fig)
    over = [o for o in over if o not in allow_over_data]
    msg = []
    if pairs:
        msg.append(f"text overlaps text: {pairs}")
    if outside:
        msg.append(f"text outside the canvas: {outside}")
    if over:
        msg.append(f"text over data marks: {over}")
    if msg:
        raise AssertionError(f"[{name}] " + " | ".join(msg))
    print(f"  [{name}] text-collision check: clean")
