#!/usr/bin/env python3
"""Figure 4b ladder — all four rungs on ONE coordinate convention (`corrected`).

    python3 recompute_F4b_ladder_261022.py

Why this exists
---------------
The Figure 4b ladder shows the 6mA TSS-proximal level ratio under the 2x2 of
(read set) x (local-GC matching). On 2026-10-22 the focal rung (all reads,
GC-matched) was corrected from 0.9212 to 0.9105 when the TSS-distance coordinate
convention was fixed, but the other three rungs stayed frozen constants inside
`15_paper_figures/scripts/panels_261013/build_F4c_promoter_input.py`, whose own
docstring recorded that they were on the `ledger` convention. One figure then
carried two conventions that differ by one base, so the rungs were not
comparable. This script recomputes all four on `corrected`.

It is a driver, not a new implementation: every quantity is computed by the
functions of `recompute_TSS_set_controls_261021.py`, which produced
SuppTable22b, imported from that file unchanged. The only thing added here is
the second read set (see below), because SuppTable22b sweeps window, universe
and convention on the all-read pileup only.

Metric definition (identical to SuppTable22b)
---------------------------------------------
level(position, library) = Nmod / Nvalid

Within each stratum q of local GC (window +/-100 bp, the TSS window), with
direct standardisation of the distal arm to the near arm's GC distribution:

    ratio(l) = SUM_q n_near(q) * mean_near(q,l)  /  SUM_q n_near(q) * mean_distal(q,l)

ratio   = exp(mean_l log ratio(l)), 95 % CI from t(0.975, 8), p from a two-sided
          one-sample t test of log ratio == 0.  Unit of replication = culture
          (9 = 3 timepoints x 3 biological replicates).
GC matching : nq = 5 (quintiles) for "GC-matched"; nq = 1 (one stratum, i.e. the
          plain near/distal mean ratio) for "unmatched".

The two read sets
-----------------
all reads = `pileup_nofilt_261015`, modkit 0.6.4 `--no-filtering` (Nfail == 0:
          every read contributes a call).
filtered  = `pileup_canonical_261015`, the default modkit confidence filter
          (Nfail > 0: low-confidence calls are discarded from Nvalid). This is
          the pileup whose filter inflates level statistics by +15.5 pp
          (canon: F7). It is RAW bedMethyl - its name column carries the
          modification code ('a'), not a motif annotation - so the motif
          positions are taken from the all-read universe and the filtered rows
          are selected at exactly those (start, strand) keys with code 'a'.
          Both read sets therefore span the same motif positions; only which
          reads enter Nvalid/Nmod differs, which is the definition of the rung.

NOT the canonical census filter: no methylation threshold and no depth floor is
applied to either read set. Thresholding a LEVEL-scale statistic on pooled
Nmod/Nvalid >= 0.5 would select positions by the quantity being measured.

Fixed parameters (the combination the manuscript quotes)
--------------------------------------------------------
tss_set = jeong_measured (3,570 experimental dRNA-seq TSSs, Jeong et al. 2016)
universe = both_A (both adenines of AAGCCCG) - this is the universe that
          reproduces the focal 0.9105; pos1_only/corrected gives 0.9130
window_bp = 100, coord_convention = corrected (both coordinates 1-based)

The focal rung is asserted against SuppTable22b. A mismatch stops the script
rather than writing a table.
"""
import importlib.util
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MODPATH = os.path.join(HERE, "recompute_TSS_set_controls_261021.py")
FILTERED = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015"
SUPP22B = ("/Users/okaban/bioinfo/rna-seq/15_paper_figures/tables/supplementary/"
           "SuppTable22b_TSS_set_controls_261021.tsv")
OUT = os.path.join(os.path.dirname(HERE), "tables", "F4b_ladder_corrected_261022.tsv")

_spec = importlib.util.spec_from_file_location("recompute_tss_261021", MODPATH)
M = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M)                      # main() is __main__-guarded

NAMES = ["a,AAGCCCG,0", "a,AAGCCCG,1"]           # universe = both_A
WINDOW = 100
TSS_SET = "jeong_measured"
UNIVERSE = "both_A"
CONV = "corrected"


def load_filtered(pos0, strand0):
    """Nvalid/Nmod from the confidence-filtered pileup at the given motif keys.

    Returns (Nvalid, Nmod, n_missing_cells): arrays shaped (9, len(pos0)), and
    the number of (library, position) cells with no 6mA row in the filtered
    pileup. A missing cell leaves Nvalid = 0, which the level computation turns
    into NaN, exactly as load_universe() does for a zero-coverage row.
    """
    keyfile = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    for s, st in zip(pos0, strand0):
        keyfile.write(f"{s}\t{st}\n")
    keyfile.close()
    idx = {(int(s), st): i for i, (s, st) in enumerate(zip(pos0, strand0))}
    Nvalid = np.zeros((len(M.LIBS), len(pos0)))
    Nmod = np.zeros((len(M.LIBS), len(pos0)))
    missing = 0
    for li, lib in enumerate(M.LIBS):
        bed = os.path.join(FILTERED, f"{lib}_pileup.bed")
        awk = ('NR==FNR{k[$1"\t"$2]=1; next} '
               '$4=="a" && (($2"\t"$6) in k){print $2"\t"$6"\t"$10"\t"$12}')
        res = subprocess.run(["awk", "-F", "\t", awk, keyfile.name, bed],
                             capture_output=True, text=True, check=True)
        seen = set()
        for line in res.stdout.splitlines():
            s, st, nv, nm = line.split("\t")
            i = idx.get((int(s), st))
            if i is not None:
                Nvalid[li, i], Nmod[li, i] = float(nv), float(nm)
                seen.add(i)
        missing += len(pos0) - len(seen)
    os.unlink(keyfile.name)
    return Nvalid, Nmod, missing


def rungs(conv=CONV):
    seq = M.read_reference(M.REF)
    cGC, cN, L = M.cum_gc(seq)

    pos0, strand0, Nvalid_all, Nmod_all = M.load_universe(NAMES)
    Nvalid_flt, Nmod_flt, missing = load_filtered(pos0, strand0)
    sys.stderr.write(f"[universe] {len(pos0)} motif positions present in all 9 "
                     f"all-read libraries; filtered pileup missing "
                     f"{missing} / {9 * len(pos0)} cells\n")

    j = pd.read_csv(f"{M.TSSDIR}/jeong2016_all_tss.csv")
    tss = np.array(sorted(j["tss_position"].astype(int)))
    assert len(np.unique(tss)) == 3570, len(np.unique(tss))

    gc = M.local_gc(pos0, cGC, cN, L, WINDOW)
    query = pos0 if conv == "ledger" else pos0 + 1      # as in the 261021 module
    near = M.min_distance(query, tss) <= WINDOW

    rows = []
    for read_set, (Nv, Nm) in (("all reads", (Nvalid_all, Nmod_all)),
                               ("filtered", (Nvalid_flt, Nmod_flt))):
        level = np.where(Nv > 0, Nm / np.maximum(Nv, 1), np.nan)
        for match, nq in (("GC-matched", 5), ("unmatched", 1)):
            ratio, lo, hi, p, n = M.summarise(
                M.standardised_ratio(level, Nv, Nm, near, gc, nq=nq))
            rows.append(dict(
                variant=f"{read_set}, {match}",
                read_set=read_set, gc_matching=match, gc_quantiles=nq,
                tss_set=TSS_SET, tss_n=int(len(np.unique(tss))),
                universe=UNIVERSE, universe_n=int(len(pos0)),
                window_bp=WINDOW, coord_convention=conv,
                n_near=int(near.sum()), n_distal=int((~near).sum()),
                mean_level_near=round(float(np.nanmean(level[:, near])), 6),
                mean_level_distal=round(float(np.nanmean(level[:, ~near])), 6),
                reads_mod_near=int(Nm[:, near].sum()),
                reads_valid_near=int(Nv[:, near].sum()),
                reads_mod_distal=int(Nm[:, ~near].sum()),
                reads_valid_distal=int(Nv[:, ~near].sum()),
                ratio=round(ratio, 4), ci95_lo=round(lo, 4), ci95_hi=round(hi, 4),
                p_value=f"{p:.3g}", n_cultures=n))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    # `--conv ledger` is a diagnostic only: it answers whether the three frozen
    # constants the panel carried (0.8052 / 0.8945 / 0.8475) are reproducible
    # under the convention their docstring claimed. It never writes the table.
    conv = "ledger" if "--conv=ledger" in sys.argv else CONV
    df = rungs(conv)
    if conv != CONV:
        print(df[["variant", "ratio", "ci95_lo", "ci95_hi", "p_value",
                  "n_near"]].to_string(index=False))
        print(f"\n[diagnostic only, coord_convention={conv}: nothing written]")
        raise SystemExit(0)

    # the focal rung must reproduce SuppTable22b, which the registry and the
    # manuscript both cite; otherwise stop without writing anything
    s = pd.read_csv(SUPP22B, sep="\t", dtype=str)
    key = s[(s.tss_set == TSS_SET) & (s.universe == UNIVERSE)
            & (s.window_bp == str(WINDOW)) & (s.coord_convention == CONV)]
    assert len(key) == 1, f"{len(key)} rows match the focal key in SuppTable22b"
    want = key.iloc[0]
    got = df[df.variant == "all reads, GC-matched"].iloc[0]
    for col in ("ratio", "ci95_lo", "ci95_hi", "n_near", "n_distal", "n_cultures"):
        assert float(got[col]) == float(want[col]), (
            f"focal rung {col}: recomputed {got[col]} vs SuppTable22b {want[col]}")
    assert abs(float(got.p_value) - float(want.p_value)) < 1e-9, (
        f"focal rung p: recomputed {got.p_value} vs SuppTable22b {want.p_value}")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    print(df[["variant", "ratio", "ci95_lo", "ci95_hi", "p_value",
              "n_near", "n_distal", "mean_level_near",
              "mean_level_distal"]].to_string(index=False))
    print(f"\n[write] {OUT}")
