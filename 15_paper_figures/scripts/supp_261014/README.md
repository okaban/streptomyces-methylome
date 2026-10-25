# Supplementary figure rebuilds at 174 mm column width (2026-10-14)

`build_pdf_weasy.py` styles images with `img { max-width:100%; height:auto; }`, so the
supplementary PDF scales every figure to the A4 text column. Six supplementary figures
were far wider than the column, which shrank their text below legibility (7 pt became
3.2–5.9 pt). These scripts rebuild them at a fixed 172.7 mm (2,040 px at 300 dpi), inside
the 174 mm cap enforced by FIG08 in `Writing/check_figure_conventions.py`.

| figure | script | old width | new width |
|---|---|---|---|
| SuppFigure3.png  | `suppfigure3_motif_positions.py`            | 387.1 → see note | 172.7 mm |
| SuppFigure4.png  | `suppfigure4_core_share.py`                 | 240.1 mm | 172.7 mm |
| SuppFigure6.png  | `suppfigure6_motif_exclusive_asymmetry.py`  | 387.1 mm | 172.7 mm |
| SuppFigure13.png | `suppfigure13_bgc_vs_density.py`            | 298.0 mm | 172.7 mm |
| SuppFigure14.png | `suppfigure14_bgc_relationships.py`         | 212.7 mm | 172.7 mm |
| SuppFigure15.png | `suppfigure15_single_molecule.py`           | 258.0 mm | 172.7 mm |

(SuppFigure3 was 321.6 mm.)

Run any script from this directory:

    cd 15_paper_figures/scripts/supp_261014
    PYTHONPATH=. python suppfigure6_motif_exclusive_asymmetry.py

Each writes directly to `Writing/fig_images/<same filename>`. Set `SUPPFIG_OUT=<dir>` to
render elsewhere for a dry run. `_suppstyle.py` fixes the width and saves with
`bbox_inches=None` so the PNG width is deterministic — do not switch back to
`bbox_inches="tight"`, which makes the rendered width depend on the text.

## Provenance of the originals

None of the six had a producing script on disk; `figure_manifest.tsv` recorded MISSING
for five. The originals were produced in kernel cells and saved as artifacts under
earlier supplementary numbering; each was identified by matching the SHA-256 of the
shipped PNG against the artifact store and recovering the cell source from its lineage:

| shipped as | original artifact |
|---|---|
| SuppFigure3.png  | `SuppFigure10.png` |
| SuppFigure4.png  | `SuppFigure11_corrected.png` |
| SuppFigure6.png  | `SuppFigure13.png` |
| SuppFigure13.png | `SuppFigure26.png` |
| SuppFigure14.png | `SuppFig_BGC_relationships.png` |
| SuppFigure15.png | `Fig_H26_single_molecule.png` |

The recovered cells read intermediate tables from session scratch directories
(`work/c4m`, `work/ext`, `/tmp/pr2`) that no longer exist. Those tables were re-derived
from the deposited modkit pileups and the nine modification-tagged BAMs and are staged
under `data/` so the scripts run standalone; see the README in each data subdirectory
for the exact extraction command.

## Numerical identity

Every re-derived quantity is asserted against the canonical census
(`data/numbers.json`) or the published legend before the figure is drawn:
TGGCCGGC 1,886 / 1,892 / 1,877; AAGCCCG 4mC 1,023 / 977 / 977; AAGCCCG 6mA union
416 / 450 / 457; 16,043 eligible GCCGGC duplexes; core shares 63.7 / 63.6 / 63.4,
69.9 / 68.9 / 68.7 and 71.2 / 71.3 / 70.2 against nulls 63.45 % and 69.72 %;
promoter inventory 18/310 vs 127/2720; ρ = −0.090, p = 0.53, n = 51; single-molecule
0.470 / 0.171 / 0.016 and 3.0 % / 3.3 %. No value was changed.

## Convention changes made while rebuilding

* whole-figure titles baked into SuppFigure3, 6 and 13 were removed (FIG01); the same
  sentence is already the first line of each legend.
* `(a)`/`(b)` prefixes inside SuppFigure4 panel titles were replaced by panel letters
  placed in figure coordinates after composition (FIG02, FIG03).
* SuppFigure14 panel a was two-coloured for a single series; it is now single-colour
  (FIG07).
* the smallest rendered font is 6.0 pt everywhere (FIG09); the originals used 5.0–5.8 pt
  for tick labels and source notes.
* "GCCGGC 4mC sites" was written as "TGGCCGGC 4mC sites" in SuppFigure4 and
  SuppFigure13 and as "TGGCCGGC 4mC level" in SuppFigure14, matching those figures' own
  legends and the recognition-sequence convention. SuppFigure6 keeps "GCCGGC", which is
  correct there: that panel is about the hexamer duplex.
