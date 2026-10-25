#!/usr/bin/env python3
"""
SuppTable22 - 6mA depletion near TSS, recomputed for four TSS sets.
Written 2026-10-21. Recomputes from raw pileup; reads no derived table.

Context
-------
The ledger (Writing/canon_rebuild_notes.md L29-36) records four TSS-set results
for AAGCCCG 6mA depletion near TSS, but none of the backing files it names
(6mA_TSS_depletion_fully_corrected.csv, 6mA_depletion_by_TSS_set.csv, ...)
exists anywhere on disk. This script recomputes the four rows from the all-read
pileup and the TSS tables, and prints the recorded values alongside for
comparison. It does not write to the manuscript, legends or value registry.

Metric definition (numerator / denominator stated explicitly)
-------------------------------------------------------------
level(position, library) = Nmod / Nvalid            (all-read pileup; Nfail == 0)

Within each local-GC quintile q of the position universe U:
    mean_near(q,l)   = mean over near positions in q of level(position,l)
    mean_distal(q,l) = mean over distal positions in q of level(position,l)

ratio(l) = SUM_q n_near(q) * mean_near(q,l)      <- numerator
           --------------------------------
           SUM_q n_near(q) * mean_distal(q,l)    <- denominator

i.e. direct standardisation of the distal arm to the near arm's local-GC
distribution, computed independently in each of the 9 cultures.

ratio   = exp( mean_l log ratio(l) )                      (geometric mean, n = 9)
95 % CI = exp( mean +/- t(0.975, 8) * SE(log ratio) )
p       = two-sided one-sample t-test of log ratio(l) == 0, df = 8
Unit of replication = culture (9 cultures = 3 timepoints x 3 biological reps).

Variants swept
--------------
window      : +/-100 bp and +/-300 bp around the TSS
universe    : both_A    = both adenines of AAGCCCG (motif positions 1 and 2)
              pos1_only = motif position 1 only (the canonical 6mA position)
coord_conv  : ledger    = distance between the 0-based BED start of the adenine
                         and the 1-based TSS coordinate (reproduces the ledger's
                         near-position counts exactly; it is an off-by-one, so
                         the effective window is [TSS-99, TSS+101] in 1-based
                         coordinates)
              corrected = both coordinates 1-based
"""
import os, sys, numpy as np, pandas as pd
from scipy import stats

NOFILT = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015"
REF    = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna"
TSSDIR = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses"
OUT    = os.environ.get("SUPP22_OUT",
         "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/48_TSS_methylation_gradient/tables")
LIBS   = ["1-1", "1-2", "1-3", "2-1", "2-3", "2-4", "3-2", "3-3", "3-4"]
CHROM  = "NC_003888.3"
BEDCOL = ["chrom", "start", "end", "name", "score", "strand", "ts", "te", "color",
          "Nvalid", "pct", "Nmod", "Ncanon", "Nother", "Ndel", "Nfail", "Ndiff", "Nnocall"]

# values recorded in Writing/canon_rebuild_notes.md L29-36, for comparison only
RECORDED = {
    "jeong_measured":       dict(tss_n=3570, n_near=236, ratio=0.9212, p=1.71e-06,
                                 ci=(0.9086, 0.9337)),
    "comprehensive_jeong":  dict(tss_n=2703, n_near=192, ratio=0.9152, p=8.35e-06, ci=(None, None)),
    "gff_gene_start_proxy": dict(tss_n=4977, n_near=256, ratio=1.0135, p=0.0409,  ci=(None, None)),
    "mixed":                dict(tss_n=7675, n_near=438, ratio=0.9684, p=0.00133, ci=(None, None)),
}


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


def load_universe(names):
    """Positions of the requested pileup rows present in ALL 9 libraries."""
    per = {}
    for lib in LIBS:
        d = pd.read_csv(f"{NOFILT}/{lib}_nofilt.bed", sep="\t", header=None, names=BEDCOL)
        assert int(d["Nfail"].sum()) == 0, f"{lib}: all-read pileup must have Nfail == 0"
        per[lib] = (d[d["name"].isin(names)]
                    .drop_duplicates(subset=["start", "strand"])
                    .sort_values("start").reset_index(drop=True))
    keys = sorted(set.intersection(*[set(zip(p["start"], p["strand"])) for p in per.values()]))
    idx = {k: i for i, k in enumerate(keys)}
    Nvalid = np.zeros((len(LIBS), len(keys)))
    Nmod = np.zeros((len(LIBS), len(keys)))
    for li, lib in enumerate(LIBS):
        p = per[lib]
        for s, st, nv, nm in zip(p["start"], p["strand"], p["Nvalid"], p["Nmod"]):
            i = idx.get((s, st))
            if i is not None:
                Nvalid[li, i], Nmod[li, i] = nv, nm
    return np.array([k[0] for k in keys]), np.array([k[1] for k in keys]), Nvalid, Nmod


def cum_gc(seq):
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    gc = ((a == ord("G")) | (a == ord("C"))).astype(np.int64)
    acgt = np.isin(a, [ord(c) for c in "ACGT"]).astype(np.int64)
    return np.concatenate([[0], np.cumsum(gc)]), np.concatenate([[0], np.cumsum(acgt)]), len(a)


def local_gc(pos0, cGC, cN, L, half):
    a = np.clip(pos0 - half, 0, L)
    b = np.clip(pos0 + half + 1, 0, L)
    return (cGC[b] - cGC[a]) / np.maximum(cN[b] - cN[a], 1)


def min_distance(query, tss):
    t = np.unique(np.asarray(tss))
    k = np.searchsorted(t, query)
    best = np.full(len(query), np.iinfo(np.int64).max)
    for off in (-1, 0):
        kk = np.clip(k + off, 0, len(t) - 1)
        best = np.minimum(best, np.abs(query - t[kk]))
    return best


def standardised_ratio(level, Nvalid, Nmod, near, gc, nq=5):
    q = pd.qcut(gc, nq, labels=False, duplicates="drop")
    ratios = []
    for li in range(level.shape[0]):
        num = den = 0.0
        for qq in np.unique(q):
            s = q == qq
            nn, dd = s & near, s & ~near
            if nn.sum() == 0 or dd.sum() == 0:
                continue
            w = nn.sum()
            num += w * np.nanmean(level[li, nn])
            den += w * np.nanmean(level[li, dd])
        ratios.append(num / den)
    return np.array(ratios)


def summarise(r):
    lg = np.log(r)
    n = len(lg)
    se = lg.std(ddof=1) / np.sqrt(n)
    tc = stats.t.ppf(0.975, n - 1)
    return (float(np.exp(lg.mean())), float(np.exp(lg.mean() - tc * se)),
            float(np.exp(lg.mean() + tc * se)),
            float(2 * stats.t.sf(abs(lg.mean() / se), n - 1)), n)


def main():
    seq = read_reference(REF)
    cGC, cN, L = cum_gc(seq)

    # --- cross-check the pileup's motif universe against the reference ------
    ref_plus = find_all(seq, "AAGCCCG")
    ref_minus = [p + 6 for p in find_all(seq, "CGGGCTT")]
    sys.stderr.write(f"[check] reference AAGCCCG instances: {len(ref_plus)} (+) "
                     f"{len(ref_minus)} (-) = {len(ref_plus) + len(ref_minus)}\n")

    universes = {
        "both_A":    load_universe(["a,AAGCCCG,0", "a,AAGCCCG,1"]),
        "pos1_only": load_universe(["a,AAGCCCG,0"]),
    }
    p1 = universes["pos1_only"][0]
    assert set(p1) <= set(ref_plus) | set(ref_minus), "pileup positions absent from reference motif set"

    j = pd.read_csv(f"{TSSDIR}/jeong2016_all_tss.csv")
    comp = pd.read_csv(f"{TSSDIR}/comprehensive_tss_table.csv")
    compc = comp[comp["chrom"] == CHROM]
    tss_sets = {
        "jeong_measured":       np.array(sorted(j["tss_position"].astype(int))),
        "comprehensive_jeong":  np.array(sorted(compc.loc[compc["tss_source"] == "Jeong2016_dRNA-seq", "tss"].astype(int))),
        "gff_gene_start_proxy": np.array(sorted(compc.loc[compc["tss_source"] == "GFF_annotation", "tss"].astype(int))),
        "mixed":                np.array(sorted(compc["tss"].astype(int))),
    }
    # the ledger reports mixed as 7,675; the table holds 7,680 rows that collapse
    # to 7,675 unique positions (5 TSS shared by a Jeong TSS and a GFF gene start)
    assert len(np.unique(tss_sets["mixed"])) == 7675, "mixed unique-position count changed"

    rows = []
    for uni_name, (pos0, strand0, Nvalid, Nmod) in universes.items():
        level = np.where(Nvalid > 0, Nmod / np.maximum(Nvalid, 1), np.nan)
        for W in (100, 300):
            gc = local_gc(pos0, cGC, cN, L, W)       # GC window matched to the TSS window
            for conv in ("ledger", "corrected"):
                query = pos0 if conv == "ledger" else pos0 + 1
                for set_name, tss in tss_sets.items():
                    d = min_distance(query, tss)
                    near = d <= W
                    if near.sum() == 0:
                        continue
                    ratio, lo, hi, p, n = summarise(
                        standardised_ratio(level, Nvalid, Nmod, near, gc))
                    rec = RECORDED[set_name]
                    rows.append(dict(
                        tss_set=set_name,
                        tss_n=len(np.unique(tss)),
                        universe=uni_name,
                        universe_n=len(pos0),
                        window_bp=W,
                        coord_convention=conv,
                        n_near=int(near.sum()),
                        n_distal=int((~near).sum()),
                        mean_level_near=round(float(np.nanmean(level[:, near])), 6),
                        mean_level_distal=round(float(np.nanmean(level[:, ~near])), 6),
                        reads_mod_near=int(Nmod[:, near].sum()),
                        reads_valid_near=int(Nvalid[:, near].sum()),
                        reads_mod_distal=int(Nmod[:, ~near].sum()),
                        reads_valid_distal=int(Nvalid[:, ~near].sum()),
                        ratio=round(ratio, 4),
                        ci95_lo=round(lo, 4),
                        ci95_hi=round(hi, 4),
                        p_value=f"{p:.3g}",
                        n_cultures=n,
                        recorded_ratio=rec["ratio"],
                        recorded_p=rec["p"],
                        recorded_n_near=rec["n_near"],
                        n_near_matches_record=int(near.sum()) == rec["n_near"],
                        ratio_matches_record_3dp=round(ratio, 3) == round(rec["ratio"], 3),
                    ))

    out = pd.DataFrame(rows).sort_values(
        ["universe", "window_bp", "coord_convention", "tss_set"]).reset_index(drop=True)
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "SuppTable22_TSS_set_controls_261021.tsv")
    out.to_csv(path, sep="\t", index=False)
    sys.stderr.write(f"[write] {path}  ({len(out)} rows)\n")
    return out


if __name__ == "__main__":
    df = main()
    key = df[(df.universe == "both_A") & (df.window_bp == 100) & (df.coord_convention == "ledger")]
    print(key[["tss_set", "tss_n", "n_near", "recorded_n_near", "ratio", "recorded_ratio",
               "ci95_lo", "ci95_hi", "p_value", "recorded_p"]].to_string(index=False))
