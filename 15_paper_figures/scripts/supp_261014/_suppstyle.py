"""Shared style for the 174 mm-wide supplementary figure rebuilds (2026-10-14).

The supplementary PDF scales every image to the A4 text column, so a figure
wider than the column has its text shrunk below legibility.  These rebuilds fix
the rendered width at NAR_WIDTH_IN and save with ``bbox_inches=None`` so the
PNG width is deterministic (6.80 in x 300 dpi = 2040 px = 172.7 mm).

Font floor is 6 pt at the rendered size (FIG09); panel letters are placed in
FIGURE coordinates after composition (FIG03); no whole-figure title is baked
into the image (FIG01).
"""
import os
from pathlib import Path

import matplotlib as mpl

DATA_DIR = Path(__file__).resolve().parent / "data"
# SUPPFIG_OUT lets a dry run render somewhere else before the shipped PNG is replaced.
OUT_DIR = Path(os.environ.get(
    "SUPPFIG_OUT", "/Users/okaban/obsidian/Research/rna-seq/Writing/fig_images"))

NAR_WIDTH_MM = 174.0
NAR_WIDTH_IN = 6.85          # 174.0 mm exactly; 4110 px at 600 dpi
DPI = int(os.environ.get('FIG_DPI', 300))

# 2026-10-22: the mark colours are now IMPORTED from 00_shared_utils rather than
# restated here.  Until today the shipped set carried four different 4mC/6mA
# colour pairs across its 18 figures — this module's #9e4a4a/#45708f, the main
# figures' #8A3D37/#68849D, SuppFigure 2/9's #1b4965/#c0552b (which also had the
# hues SWAPPED, 4mC blue and 6mA orange), and SuppFigure 14's #1b6ca8/#b5651d —
# so a reader could not follow one mark from one figure to the next.  There is
# one author decision on record (00_shared_utils, 2026-09-26: 4mC/GCCGGC is the
# muted red, 6mA/AAGCCCG the muted blue, with the v2 luminance spacing measured
# for greyscale print and deuteranopia), and it governs every figure.  The names
# RED and BLUE are kept so the call sites do not have to change.
import importlib.util as _ilu  # noqa: E402

_SU = Path(__file__).resolve().parent.parent / "00_shared_utils.py"
_spec = _ilu.spec_from_file_location("_shared_utils_palette", _SU)
_U = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_U)

RED = _U.COL_4mC     # 4mC  — #8A3D37 deep muted red
BLUE = _U.COL_6mA    # 6mA  — #68849D mid muted blue
TAN = "#d9a05b"
GREY = "#9a9a9a"
LIGHT = "#e8e8e8"
META_GREY = "#888888"

FONT_FLOOR = 6.0


def apply_figure_style(*, frame="boxed", sizes=(7.5, 6.5, 6.0), grid=False):
    """Publication rc for a column-width figure.  sizes = (base, secondary, tick).

    frame の既定は 2026-10-28 に "open" から "boxed" に変えた（著者決定 D121:
    補足図も本文図と同じく四辺の枠に揃える）。"open" を明示すれば従来の
    上右なし体裁に戻せるが、投稿物の中で作図規約を混在させないこと。
    """
    import glob
    import os
    import sys

    import matplotlib.font_manager as fm

    try:
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = frame == "boxed"
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 2.4, "ytick.major.size": 2.4,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.pad": 2.0, "ytick.major.pad": 2.0,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": DPI,
        "savefig.bbox": None,          # deterministic width — do not use "tight"
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.0,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def figure_panel_letters(fig, axes, letters, dx=0.040, dy=0.012, fontsize=8.5):
    """Place panel letters in FIGURE coordinates, after the layout is final."""
    dxs = dx if isinstance(dx, (list, tuple)) else [dx] * len(axes)
    dys = dy if isinstance(dy, (list, tuple)) else [dy] * len(axes)
    for ax, lab, dx, dy in zip(axes, letters, dxs, dys):
        box = ax.get_position()
        fig.text(max(0.002, box.x0 - dx), min(0.999, box.y1 + dy), lab,
                 fontsize=fontsize, fontweight="bold", ha="left", va="bottom",
                 color="black")


def source_note(fig, text, y=0.012, x=0.012, fontsize=6.0):
    fig.text(x, y, text, fontsize=fontsize, color="0.42", ha="left", va="bottom",
             linespacing=1.45)


def _collision_checker():
    """Load assert_no_text_collisions from the main-figure utils.

    The main figures (panels_261013/*.py) have called this guard since 2026-09-22;
    the supplementary pipeline was built separately on 2026-10-14 and did not
    inherit it, which let text-over-data ship in supplementary figures.  Importing
    the same function keeps one definition of the rule for both pipelines.
    """
    import importlib.util
    src = Path(__file__).resolve().parent.parent / "00_shared_utils.py"
    spec = importlib.util.spec_from_file_location("_shared_utils_for_supp", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.assert_no_text_collisions


def _prune_offview_ticks(fig):
    """Drop ticks that lie outside their axis view limits.

    matplotlib's auto locator routinely places one tick beyond the data range
    (a tick at 9.0 when the limit is 8.667). It draws no tick MARK there, but it
    still creates the LABEL and renders it past the axes edge -- on a figure whose
    axes reach within ~15 px of the canvas, that label leaves the canvas entirely.
    Removing the off-view ticks is visually a no-op apart from the spilled label.
    """
    for ax in fig.axes:
        for axis, getlim in ((ax.xaxis, ax.get_xlim), (ax.yaxis, ax.get_ylim)):
            lo, hi = sorted(getlim())
            ticks = [t for t in axis.get_ticklocs() if lo - 1e-9 <= t <= hi + 1e-9]
            if len(ticks) != len(axis.get_ticklocs()):
                axis.set_ticks(ticks)


def save(fig, name, allow_over_data=()):
    """Render to OUT_DIR after the text-collision guard passes.

    FIGURE_DESIGN_RULES_260922.md line 140 forbids text overlapping other text,
    text over data marks, and text outside the canvas.  The check runs before the
    file is written, so a violating figure is never produced.  allow_over_data
    lists label strings deliberately placed on a mark (e.g. a value inside a bar),
    which must then meet the contrast floor.
    """
    _prune_offview_ticks(fig)
    _collision_checker()(fig, name=name, allow_over_data=allow_over_data)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    # 2026-11-09: 余白を 3 mm 以下に落とすため tight 描画→174 mm 幅に正規化する
    # 経路（本文図と同じ save_figure_raster）に揃えた。幅は常に 174 mm 相当。
    _U.save_figure_raster(fig, path, width_class="full", dpi=DPI)
    return path


def report(path):
    from PIL import Image
    im = Image.open(path)
    w, h = im.size
    print(f"{path.name}: {w}x{h} px, {w / DPI * 25.4:.1f} x {h / DPI * 25.4:.1f} mm")
