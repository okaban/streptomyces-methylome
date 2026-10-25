"""Pooled GCCGGC hexamer 4mC levels from the deposited pileups (see data/pileup_GCCGGC/README.md)."""
import numpy as np
import pandas as pd

from _suppstyle import DATA_DIR

GENOME = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna"
LIBS = {"T1": ["1-1", "1-2", "1-3"],
        "T2": ["2-1", "2-3", "2-4"],
        "T3": ["3-2", "3-3", "3-4"]}
PILEUP_DIR = DATA_DIR / "pileup_GCCGGC"
DEPTH_MIN = 30
FRAC_MIN = 0.5


def genome():
    seq = []
    with open(GENOME) as fh:
        fh.readline()
        for line in fh:
            seq.append(line.strip())
    return "".join(seq).upper()


def find_all(s, sub):
    out, i = [], s.find(sub)
    while i != -1:
        out.append(i)
        i = s.find(sub, i + 1)
    return out


def motif_positions(g):
    """Strand-oriented modified-cytosine positions."""
    gcc = find_all(g, "GCCGGC")
    duplex = [(i + 2, i + 3) for i in gcc]                       # (+ strand C, - strand C)
    tgg = [(i + 4, "+") for i in find_all(g, "TGGCCGGC")]
    tgg += [(i + 3, "-") for i in find_all(g, "GCCGGCCA")]
    return duplex, sorted(tgg)


def pooled(timepoint):
    """dict (pos, strand) -> (Nvalid, Nmod, Nfail) summed over the three libraries."""
    acc = {}
    for lib in LIBS[timepoint]:
        d = pd.read_csv(PILEUP_DIR / f"{lib}.tsv", sep="\t", header=None,
                        names=["pos", "strand", "Nvalid", "Nmod", "Nfail"])
        for pos, strand, nv, nm, nf in d.itertuples(index=False):
            key = (pos, strand)
            cur = acc.get(key)
            if cur is None:
                acc[key] = [nv, nm, nf]
            else:
                cur[0] += nv
                cur[1] += nm
                cur[2] += nf
    return acc


def level(acc, key):
    """Pooled modified fraction, or nan when the position is not measurable."""
    v = acc.get(key)
    if v is None or v[0] == 0:
        return np.nan
    return v[1] / v[0]


def eligible(acc, key):
    v = acc.get(key)
    return v is not None and (v[0] + v[2]) >= DEPTH_MIN


def called(acc, keys):
    return [k for k in keys
            if eligible(acc, k) and acc[k][0] > 0 and acc[k][1] / acc[k][0] >= FRAC_MIN]
