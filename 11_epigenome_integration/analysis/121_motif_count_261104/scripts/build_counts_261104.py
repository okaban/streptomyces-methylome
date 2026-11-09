#!/usr/bin/env python3
"""121 - per-gene recognition-sequence counts (possible vs actually modified).

Motif definitions and the site-calling rule are taken verbatim from
48_TSS_methylation_gradient/scripts/recompute_level_cut_sweep_261021.py
(pooled per timepoint: eligible = Nvalid+Nfail >= 30 AND refok in >=1 library;
 site = eligible AND Nmod/Nvalid >= 0.5).

Writes: tables/per_site_261104.tsv, tables/per_gene_counts_261104.tsv
Memory: one pileup npz at a time (~0.3 GB), nothing else large.
"""
import os, re, sys, resource
import numpy as np, pandas as pd

BASE = "/Users/okaban/bioinfo/rna-seq"
HERE = f"{BASE}/11_epigenome_integration/analysis/121_motif_count_261104"
CACHE = f"{BASE}/15_paper_figures/scripts/recompute_261015/pileup_cache_261015"
REF = f"{BASE}/11_epigenome_integration/data/NC_003888.3.fna"
ANN = f"{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
TSS_CSV = f"{BASE}/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv"
TP = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"], "T3": ["3-2", "3-3", "3-4"]}
MIN_TRUE_DEPTH, CUT = 30, 0.5
CODE = {"6mA": 0, "4mC": 1, "5mC": 2}
PROM = 300


def read_reference(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def find_all(seq, pat):
    out, i = [], seq.find(pat)
    while i != -1:
        out.append(i); i = seq.find(pat, i + 1)
    return out


def motif_positions(seq):
    """name -> ([(pos0, is_plus)], mod). Offsets verified in recompute_level_cut_sweep."""
    m = {}
    m["TGGCCGGC_4mC"] = (sorted(set([(x + 4, True) for x in find_all(seq, "TGGCCGGC")] +
                                    [(x + 3, False) for x in find_all(seq, "GCCGGCCA")])), "4mC")
    m["AAGCCCG_4mC"] = (sorted(set([(x + 4, True) for x in find_all(seq, "AAGCCCG")] +
                                   [(x + 2, False) for x in find_all(seq, "CGGGCTT")])), "4mC")
    m["AAGCCCG_6mA"] = (sorted(set([(x + 0, True) for x in find_all(seq, "AAGCCCG")] +
                                   [(x + 6, False) for x in find_all(seq, "CGGGCTT")])), "6mA")
    return m


def fetch(lib, msets):
    z = np.load(f"{CACHE}/{lib}.npz")
    key_all = z["start"].astype(np.int64) * 2 + z["plus"].astype(np.int64)
    by_code, out = {}, {}
    for code in sorted({CODE[mod] for _, mod in msets.values()}):
        sel = z["code"] == code
        k = key_all[sel]; order = np.argsort(k, kind="stable"); k = k[order]
        assert np.all(np.diff(k) > 0), f"{lib}: duplicate key in code {code}"
        by_code[code] = (k, z["Nvalid"][sel][order], z["Nmod"][sel][order],
                         z["Nfail"][sel][order], z["refok"][sel][order])
    for name, (pl, mod) in msets.items():
        q = np.array([p * 2 + int(b) for p, b in pl], dtype=np.int64)
        k, nv, nm, nf, ro = by_code[CODE[mod]]
        i = np.searchsorted(k, q); ic = np.clip(i, 0, len(k) - 1)
        hit = (i < len(k)) & (k[ic] == q)
        d = {x: np.zeros(len(q), np.int64) for x in ("Nvalid", "Nmod", "Nfail")}
        d["refok"] = np.zeros(len(q), bool)
        for x, src in (("Nvalid", nv), ("Nmod", nm), ("Nfail", nf)):
            d[x][hit] = src[ic][hit]
        d["refok"][hit] = ro[ic][hit]
        out[name] = d
    del z, key_all, by_code
    return out


def main():
    seq = read_reference(REF)
    sys.stderr.write(f"[ref] {len(seq)} bp\n")
    msets = motif_positions(seq)
    for n, (pl, mod) in msets.items():
        sys.stderr.write(f"[motif] {n}: {len(pl)} positions ({mod})\n")

    # ---- per-site modification calls, one library at a time -------------------
    acc = {tp: {n: {x: np.zeros(len(pl), np.int64) for x in ("Nvalid", "Nmod", "Nfail")}
                for n, (pl, _) in msets.items()} for tp in TP}
    for tp in TP:
        for n in msets:
            acc[tp][n]["refok"] = np.zeros(len(msets[n][0]), bool)
    for tp, libs in TP.items():
        for lib in libs:
            d = fetch(lib, msets)
            for n in msets:
                for x in ("Nvalid", "Nmod", "Nfail"):
                    acc[tp][n][x] += d[n][x]
                acc[tp][n]["refok"] |= d[n]["refok"]
            del d
        sys.stderr.write(f"[pooled] {tp}\n")

    rows = []
    for n, (pl, mod) in msets.items():
        pos = np.array([p for p, _ in pl]); plus = np.array([b for _, b in pl])
        r = pd.DataFrame({"motif": n, "mod": mod, "pos": pos,
                          "strand": np.where(plus, "+", "-")})
        for tp in TP:
            a = acc[tp][n]
            elig = (a["Nvalid"] + a["Nfail"] >= MIN_TRUE_DEPTH) & a["refok"]
            lev = np.where(a["Nvalid"] > 0, a["Nmod"] / np.maximum(a["Nvalid"], 1), np.nan)
            r[f"eligible_{tp}"] = elig
            r[f"level_{tp}"] = np.round(lev, 4)
            r[f"site_{tp}"] = elig & (lev >= CUT)
        rows.append(r)
    sites = pd.concat(rows, ignore_index=True)
    sites.to_csv(f"{HERE}/tables/per_site_261104.tsv", sep="\t", index=False)
    for n in msets:
        s = sites[sites.motif == n]
        sys.stderr.write(f"[call] {n}: n={len(s)} " + " ".join(
            f"{tp} elig={int(s['eligible_'+tp].sum())} site={int(s['site_'+tp].sum())}"
            for tp in TP) + "\n")
    del acc

    # ---- gene / promoter intervals -------------------------------------------
    ann = pd.read_csv(ANN, sep="\t")
    ann = ann[ann.contig == "NC_003888.3"].copy()
    assert len(ann) == 7863, len(ann)
    ann["b0"] = ann.start.astype(int) - 1       # 0-based inclusive
    ann["b1"] = ann.end.astype(int)             # 0-based exclusive
    ann["glen"] = ann.b1 - ann.b0

    tss = pd.read_csv(TSS_CSV)
    tss = tss[tss.gene_id.notna() & (tss.gene_id.astype(str).str.strip() != "")]
    # one TSS per gene: the most abundant
    tss = tss.sort_values("abundance", ascending=False).drop_duplicates("gene_id")
    tmap = dict(zip(tss.gene_id, tss.tss_position.astype(int)))
    ann["tss_exp"] = ann.gene_id.map(tmap)
    ann["has_tss"] = ann.tss_exp.notna()
    anchor = np.where(ann.has_tss, ann.tss_exp.fillna(0).astype(int),
                      np.where(ann.strand == "+", ann.b0, ann.b1 - 1))
    ann["anchor"] = anchor
    ann["p0"] = np.where(ann.strand == "+", ann.anchor - PROM, ann.anchor + 1)
    ann["p1"] = np.where(ann.strand == "+", ann.anchor, ann.anchor + 1 + PROM)
    ann["p0"] = ann.p0.clip(0, len(seq)); ann["p1"] = ann.p1.clip(0, len(seq))

    # ---- GC ------------------------------------------------------------------
    g = np.frombuffer(seq.encode(), dtype=np.uint8)
    isgc = ((g == ord("G")) | (g == ord("C"))).astype(np.int64)
    cum = np.concatenate([[0], np.cumsum(isgc)])

    def gc(a, b):
        a = np.asarray(a); b = np.asarray(b)
        n = np.maximum(b - a, 1)
        return (cum[b] - cum[a]) / n
    ann["gene_gc"] = gc(ann.b0.values, ann.b1.values)
    ann["prom_gc"] = gc(ann.p0.values, ann.p1.values)
    ann["prom_len"] = ann.p1.values - ann.p0.values

    # ---- counting ------------------------------------------------------------
    def count_in(intervals_a, intervals_b, positions):
        """# of positions inside each [a,b) interval. intervals may overlap; O(n log n)."""
        p = np.sort(np.asarray(positions))
        return (np.searchsorted(p, intervals_b) - np.searchsorted(p, intervals_a))

    out = ann[["gene_id", "old_locus_tag", "start", "end", "strand", "glen",
               "gene_gc", "prom_gc", "prom_len", "has_tss", "b0", "b1", "p0", "p1"]].copy()
    for n in msets:
        s = sites[sites.motif == n]
        out[f"n_{n}_gene"] = count_in(out.b0.values, out.b1.values, s.pos.values)
        out[f"n_{n}_prom"] = count_in(out.p0.values, out.p1.values, s.pos.values)
        for tp in TP:
            m = s.loc[s[f"site_{tp}"], "pos"].values
            out[f"nmod_{n}_gene_{tp}"] = count_in(out.b0.values, out.b1.values, m)
            out[f"nmod_{n}_prom_{tp}"] = count_in(out.p0.values, out.p1.values, m)
    # combined columns
    for reg in ("gene", "prom"):
        out[f"n_ANY_{reg}"] = sum(out[f"n_{n}_{reg}"] for n in msets)
        out[f"n_AAGCCCG_{reg}"] = out[f"n_AAGCCCG_4mC_{reg}"]  # one instance per copy
        for tp in TP:
            out[f"nmod_ANY_{reg}_{tp}"] = sum(out[f"nmod_{n}_{reg}_{tp}"] for n in msets)
    out.to_csv(f"{HERE}/tables/per_gene_counts_261104.tsv", sep="\t", index=False)
    sys.stderr.write(f"[write] per_gene_counts_261104.tsv  ({len(out)} genes)\n")
    sys.stderr.write(f"[genes with TSS] {int(out.has_tss.sum())}\n")
    sys.stderr.write("[peak RSS GB] %.2f\n" %
                     (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3))


if __name__ == "__main__":
    main()
