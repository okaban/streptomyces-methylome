#!/usr/bin/env python3
"""Canonical per-timepoint site tables for figure scripts (2026-09-21, BLOCKER-0).

Every manuscript figure that shows per-timepoint GCCGGC / AAGCCCG site sets must
read the CANONICAL high-confidence file

    01_integration/high_confidence_sites_weighted.csv   (0-based positions,
    depth >= 10 per replicate, 3 replicates, weighted modification frequency >= 50 %)

and NOT the position-deduplicated tables 07_/23_/37_ (first-appearance: a site is
kept only at the FIRST timepoint it is called, so "T2" there means "new at T2";
see 90_per_timepoint_census_audit/README.md).

Motif assignment is byte-identical to 90_per_timepoint_census_audit/per_timepoint_census.py
(which produced tables/per_timepoint_census.tsv: GCCGGC 4mC 1,289/1,595/1,073;
AAGCCCG 6mA 418/451/443; AAGCCCG 4mC(C4) 698/851/574).

Public helpers
--------------
load_canonical()              -> DataFrame: chrom, position, strand, mod (4mC|6mA), timepoint,
                                 frequency (= weighted_mod_freq), motif (GCCGGC|AAGCCCG|other), region
gccggc_by_timepoint()         -> GCCGGC 4mC rows, column layout compatible with the retired
                                 37_defense_island_GCCGGC/tables/GCCGGC_sites_by_timepoint.tsv
                                 (chrom, position, strand, timepoint, frequency, final_motif, region)
aagcccg_pooled_by_timepoint() -> AAGCCCG rows with 4mC (C4) and 6mA (A0/A1) pooled
                                 (one row per (position, strand, mod, timepoint); 1,116/1,302/1,017)
aagcccg_by_timepoint(mod)     -> AAGCCCG rows for one modification ('4mC' -> C4; '6mA' -> A0/A1)
"""
from pathlib import Path
import pandas as pd

ANALYSIS = Path(__file__).resolve().parent.parent
CANONICAL_CSV = ANALYSIS / "01_integration" / "high_confidence_sites_weighted.csv"
REF_FA = Path("/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa")
# Core boundary: use the value the figure scripts use (77_reviewer_figures ARM_RIGHT),
# not a re-rounded 7_170_000 — a second hardcoded boundary in the loader meant to
# unify them defeats the purpose. Measured 2026-09-23: zero sites in
# high_confidence_sites_weighted.csv fall in the 2,493 bp between the two values,
# so no published number changes; this is consistency, not a correction.
CORE_LO, CORE_HI = 1_500_000, 7_167_507      # manuscript core 1.5-7.17 Mb (Methods)

_REF = None


def reference():
    global _REF
    if _REF is None:
        _REF = "".join(l.strip() for l in open(REF_FA) if not l.startswith(">")).upper()
    return _REF


def motif(p, s, mod, ref=None):
    """Identical to 90_/per_timepoint_census.py::motif (0-based position p)."""
    ref = ref if ref is not None else reference()
    if mod == "4mC":
        if (ref[p-2:p+4] == "GCCGGC") if s == "+" else (ref[p-3:p+3] == "GCCGGC"):
            return "GCCGGC"
        if (ref[p-4:p+3] == "AAGCCCG") if s == "+" else (ref[p-2:p+5] == "CGGGCTT"):
            return "AAGCCCG"
    else:
        if s == "+" and (ref[p:p+7] == "AAGCCCG" or ref[p-1:p+6] == "AAGCCCG"):
            return "AAGCCCG"
        if s == "-" and (ref[p-6:p+1] == "CGGGCTT" or ref[p-5:p+2] == "CGGGCTT"):
            return "AAGCCCG"
    return "other"


def load_canonical():
    hc = pd.read_csv(CANONICAL_CSV)
    hc["position"] = hc["position"].astype(int)
    hc["mod"] = hc["mod_type"].astype(str)
    ref = reference()
    hc["motif"] = [motif(p, s, m, ref) for p, s, m in zip(hc["position"], hc["strand"], hc["mod"])]
    hc["frequency"] = hc["weighted_mod_freq"].astype(float)
    hc["region"] = ["core" if CORE_LO <= p <= CORE_HI else "arm" for p in hc["position"]]
    return hc[["chrom", "position", "strand", "mod", "timepoint", "frequency", "motif", "region",
               "n_reps", "total_coverage"]]


def gccggc_by_timepoint():
    hc = load_canonical()
    g = hc[(hc["mod"] == "4mC") & (hc["motif"] == "GCCGGC")].copy()
    g["final_motif"] = "GCCGGC"
    g = g[["chrom", "position", "strand", "timepoint", "frequency", "final_motif", "region"]]
    assert not g.duplicated(["position", "strand", "timepoint"]).any()
    return g.sort_values(["timepoint", "position"]).reset_index(drop=True)


def aagcccg_by_timepoint(mod):
    hc = load_canonical()
    a = hc[(hc["mod"] == mod) & (hc["motif"] == "AAGCCCG")].copy()
    a["final_motif"] = "AAGCCCG"
    return a[["chrom", "position", "strand", "mod", "timepoint", "frequency", "final_motif", "region"]] \
        .sort_values(["timepoint", "position"]).reset_index(drop=True)


def aagcccg_pooled_by_timepoint():
    a = pd.concat([aagcccg_by_timepoint("4mC"), aagcccg_by_timepoint("6mA")], ignore_index=True)
    assert not a.duplicated(["position", "strand", "mod", "timepoint"]).any()
    return a.sort_values(["timepoint", "position"]).reset_index(drop=True)


if __name__ == "__main__":
    g = gccggc_by_timepoint()
    print("GCCGGC 4mC:", g.groupby("timepoint").size().to_dict(),
          g.groupby("timepoint")["region"].apply(lambda r: round((r == "core").mean() * 100, 1)).to_dict())
    a = aagcccg_pooled_by_timepoint()
    print("AAGCCCG pooled:", a.groupby("timepoint").size().to_dict(),
          a.groupby("timepoint")["region"].apply(lambda r: round((r == "core").mean() * 100, 1)).to_dict())
