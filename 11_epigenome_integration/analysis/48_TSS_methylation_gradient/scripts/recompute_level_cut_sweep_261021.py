#!/usr/bin/env python3
"""
SuppTable29 - occupancy as a function of the pooled level cut.
Written 2026-10-21. Recomputes from the genome-wide pileup cache.

Context
-------
The headline occupancies (TGGCCGGC 4mC 90.11 %, AAGCCCG 6mA 24.29 %) and the
central contrast ("6mA roughly one third of 4mC") are all functions of ONE
threshold: a position is called a site when the pooled Nmod/Nvalid >= 0.5.
The sweep shown in the manuscript varies the per-read probability cutoff, which
is a different threshold. The level cut itself has never been swept. This
script sweeps it from 0.1 to 0.9.

Site-calling rule (verified against canonical_sites_genomewide_261015.tsv)
-------------------------------------------------------------------------
Per timepoint, the 3 libraries of that timepoint are POOLED by summing Nvalid,
Nmod and Nfail at each (position, strand, mod-code).

    eligible  : Nvalid_pooled + Nfail_pooled >= 30   (true depth)
                AND refok in at least one of the 3 libraries
    level     : Nmod_pooled / Nvalid_pooled
    site      : eligible AND level >= cut

At cut = 0.5 this reproduces canonical_sites_genomewide_261015.tsv exactly for
both motifs at all three timepoints (T1 1,886 / T2 1,892 / T3 1,866 TGGCCGGC
4mC; T1 323 / T2 350 / T3 335 AAGCCCG 6mA position 1), with zero positions
differing in either direction. Requiring refok in ALL three libraries instead
gives 1,884 / 1,892 / 1,861 and does NOT reproduce the canonical list.

Occupancy, stated as numerator / denominator
--------------------------------------------
    occupancy_eligible_pct = 100 * n_sites / n_eligible      <- canonical rule
    occupancy_allpos_pct   = 100 * n_sites / n_positions      <- all motif copies

The two headline numbers do not use the same denominator (see the script's
stdout summary), which is why 90.11 % and 24.29 % cannot both be reproduced
with one convention.

Modified-base positions (verified against the canonical site list)
------------------------------------------------------------------
    TGGCCGGC 4mC       : plus-strand match + 4 ; GCCGGCCA match + 3 on minus
    AAGCCCG  6mA pos 1 : plus-strand match + 0 ; CGGGCTT  match + 6 on minus
    HAAGCCCG 6mA pos 1 : as above, restricted to instances whose 5' flanking
                         base in motif orientation is not G
"""
import os, sys, numpy as np, pandas as pd

CACHE = "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/recompute_261015/pileup_cache_261015"
REF   = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna"
CANON = "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/recompute_261015/canonical_sites_genomewide_261015.tsv"
OUT   = os.environ.get("SUPP29_OUT",
        "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/48_TSS_methylation_gradient/tables")
TP    = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"], "T3": ["3-2", "3-3", "3-4"]}
CUTS  = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
MIN_TRUE_DEPTH = 30
CODE = {"6mA": 0, "4mC": 1, "5mC": 2}


def read_reference(path):
    parts = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def find_all(seq, pat):
    out, i = [], seq.find(pat)
    while i != -1:
        out.append(i)
        i = seq.find(pat, i + 1)
    return out


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def motif_positions(seq):
    """name -> (sorted [(pos0, is_plus)], mod name)"""
    m = {}
    p = [(x + 4, True) for x in find_all(seq, "TGGCCGGC")] + \
        [(x + 3, False) for x in find_all(seq, revcomp("TGGCCGGC"))]
    m["TGGCCGGC_4mC"] = (sorted(set(p)), "4mC")
    a = [(x + 0, True) for x in find_all(seq, "AAGCCCG")] + \
        [(x + 6, False) for x in find_all(seq, "CGGGCTT")]
    m["AAGCCCG_6mA_pos1"] = (sorted(set(a)), "6mA")
    h = [(x + 0, True) for x in find_all(seq, "AAGCCCG") if x >= 1 and seq[x - 1] in "ACT"] + \
        [(x + 6, False) for x in find_all(seq, "CGGGCTT") if x + 7 < len(seq) and seq[x + 7] in "AGT"]
    m["HAAGCCCG_6mA_pos1"] = (sorted(set(h)), "6mA")
    return m


def fetch(lib, msets):
    """Pull Nvalid/Nmod/Nfail/refok at the requested positions out of one cache file."""
    z = np.load(f"{CACHE}/{lib}.npz")
    key_all = z["start"].astype(np.int64) * 2 + z["plus"].astype(np.int64)
    by_code = {}
    for code in sorted({CODE[mod] for _, mod in msets.values()}):
        sel = z["code"] == code
        k = key_all[sel]
        order = np.argsort(k, kind="stable")
        k = k[order]
        assert np.all(np.diff(k) > 0), f"{lib}: duplicate (start,strand) within code {code}"
        by_code[code] = (k, z["Nvalid"][sel][order], z["Nmod"][sel][order],
                         z["Nfail"][sel][order], z["refok"][sel][order])
    out = {}
    for name, (pl, mod) in msets.items():
        q = np.array([p * 2 + int(b) for p, b in pl], dtype=np.int64)
        k, nv, nm, nf, ro = by_code[CODE[mod]]
        i = np.searchsorted(k, q)
        ic = np.clip(i, 0, len(k) - 1)
        hit = (i < len(k)) & (k[ic] == q)
        d = {x: np.zeros(len(q), dtype=np.int64) for x in ("Nvalid", "Nmod", "Nfail")}
        d["refok"] = np.zeros(len(q), dtype=bool)
        for x, src in (("Nvalid", nv), ("Nmod", nm), ("Nfail", nf)):
            d[x][hit] = src[ic][hit]
        d["refok"][hit] = ro[ic][hit]
        out[name] = d
    return out


def main():
    seq = read_reference(REF)
    msets = motif_positions(seq)
    for name, (pl, mod) in msets.items():
        sys.stderr.write(f"[motif] {name}: {len(pl)} positions ({mod})\n")

    per_lib = {lib: fetch(lib, msets) for libs in TP.values() for lib in libs}

    canon = pd.read_csv(CANON, sep="\t")
    rows, check = [], []
    for tp, libs in TP.items():
        for name, (pl, mod) in msets.items():
            nv = sum(per_lib[l][name]["Nvalid"] for l in libs)
            nm = sum(per_lib[l][name]["Nmod"] for l in libs)
            nf = sum(per_lib[l][name]["Nfail"] for l in libs)
            ok = np.logical_or.reduce([per_lib[l][name]["refok"] for l in libs])
            eligible = (nv + nf >= MIN_TRUE_DEPTH) & ok
            level = np.where(nv > 0, nm / np.maximum(nv, 1), np.nan)
            for cut in CUTS:
                site = eligible & (level >= cut)
                rows.append(dict(
                    motif=name, mod=mod, timepoint=tp, level_cut=cut,
                    n_positions=len(pl), n_eligible=int(eligible.sum()),
                    n_sites=int(site.sum()),
                    occupancy_eligible_pct=round(100 * site.sum() / eligible.sum(), 2),
                    occupancy_allpos_pct=round(100 * site.sum() / len(pl), 2),
                    pooled_reads_mod=int(nm[eligible].sum()),
                    pooled_reads_valid=int(nv[eligible].sum()),
                    median_level_eligible=round(float(np.nanmedian(level[eligible])), 4),
                    is_canonical_cut=(cut == 0.5)))
            # reproduction check of the canonical site list at cut 0.5
            P = np.array([p for p, _ in pl]); B = np.array([b for _, b in pl])
            mine = {(int(P[i]), bool(B[i])) for i in np.where(eligible & (level >= 0.5))[0]}
            c = {(p, s == "+") for p, s in zip(
                canon.loc[(canon.timepoint == tp) & (canon["mod"] == mod), "pos"],
                canon.loc[(canon.timepoint == tp) & (canon["mod"] == mod), "strand"])} & set(pl)
            check.append((tp, name, len(mine), len(c), len(mine - c), len(c - mine)))

    df = pd.DataFrame(rows)
    # per (timepoint, cut) contrast: is 4mC occupancy above 6mA occupancy?
    piv = df.pivot_table(index=["timepoint", "level_cut"], columns="motif",
                         values="occupancy_eligible_pct")
    cont = (piv["TGGCCGGC_4mC"] / piv["AAGCCCG_6mA_pos1"]).rename("ratio_4mC_over_6mA").reset_index()
    cont["ordering_4mC_gt_6mA"] = cont["ratio_4mC_over_6mA"] > 1
    cont["ratio_4mC_over_6mA"] = cont["ratio_4mC_over_6mA"].round(3)
    df = df.merge(cont, on=["timepoint", "level_cut"], how="left")
    df = df.sort_values(["motif", "timepoint", "level_cut"]).reset_index(drop=True)

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "SuppTable29_level_cut_sweep_261021.tsv")
    df.to_csv(path, sep="\t", index=False)

    sys.stderr.write("[check] cut 0.5 vs canonical_sites_genomewide_261015.tsv\n")
    for tp, name, nmine, nc, extra, missing in check:
        sys.stderr.write(f"         {tp} {name:20s} mine={nmine} canon={nc} "
                         f"extra={extra} missing={missing}\n")
    sys.stderr.write(f"[write] {path}  ({len(df)} rows)\n")
    return df


if __name__ == "__main__":
    d = main()
    print(d[d.is_canonical_cut][["motif", "timepoint", "n_sites", "n_eligible", "n_positions",
                                 "occupancy_eligible_pct", "occupancy_allpos_pct"]].to_string(index=False))
    print()
    print(d[["timepoint", "level_cut", "ratio_4mC_over_6mA", "ordering_4mC_gt_6mA"]]
          .drop_duplicates().sort_values(["timepoint", "level_cut"]).to_string(index=False))
