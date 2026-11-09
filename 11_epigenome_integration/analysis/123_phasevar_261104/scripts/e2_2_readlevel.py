#!/usr/bin/env python
"""E2-2: per-read homopolymer tract length at candidate phase-variation loci and at
matched control loci (same base, same reference length, outside MTase/BREX genes).
Control loci supply the nanopore length-error background."""
import os, sys, resource, re
import pandas as pd, numpy as np, pysam
from multiprocessing import Pool

BASE = "/Users/okaban/bioinfo/rna-seq"
OUT  = f"{BASE}/11_epigenome_integration/analysis/123_phasevar_261104"
BAMDIR = "/Users/okaban/bioinfo/methyl/260102_M145/analysis"
LIBS = ["1-1","1-2","1-3","2-1","2-3","2-4","3-2","3-3","3-4"]
CHROM = "NC_003888.3"
ANCHOR = 25      # bp of aligned flank required each side
MINMAPQ = 20
CAND_MINLEN = 7  # candidate tracts: homopolymers >= 7 bp in target genes
CTRL_PAD = 10000 # control loci must be >= this far from any target gene

def build_loci():
    rep = pd.read_csv(f"{OUT}/tables/e2_1_repeats_all_261104.tsv", sep="\t")
    hp  = rep[rep.repeat_class == "homopolymer"].copy()
    cand = hp[hp.repeat_len >= CAND_MINLEN].copy()
    # dedupe tracts reported twice (divergent gene pair / both strands)
    cand["key"] = cand.genome_start.astype(str) + "_" + cand.genome_end.astype(str)
    agg = (cand.groupby("key")
              .agg(genome_start=("genome_start","first"), genome_end=("genome_end","first"),
                   ref_len=("repeat_len","first"),
                   genes=("gene_id", lambda s: "|".join(sorted(set(s)))),
                   sco=("old_locus_tag", lambda s: "|".join(sorted(set(map(str,s))))),
                   gname=("gene_name", lambda s: "|".join(sorted(set(map(str,s))))),
                   region=("region", lambda s: "|".join(sorted(set(s)))))
              .reset_index(drop=True))
    agg["locus_class"] = "candidate"
    # genomic base of the tract (plus strand), read from reference
    gw = pd.read_csv(f"{OUT}/tables/e2_1_genomewide_homopolymers_261104.tsv", sep="\t")
    bmap = {(r.genome_start, r.genome_end): r.base for r in gw.itertuples()}
    agg["base"] = [bmap.get((s,e), "?") for s,e in zip(agg.genome_start, agg.genome_end)]

    # target gene intervals to exclude from controls
    rep_g = rep[["gene_start","gene_end"]].drop_duplicates()
    iv = [(int(a)-CTRL_PAD, int(b)+CTRL_PAD) for a,b in rep_g.itertuples(index=False)]
    def near_target(p):
        return any(a <= p <= b for a,b in iv)

    ctrls = []
    for L in sorted(agg.ref_len.unique()):
        for base in sorted(agg[agg.ref_len==L].base.unique()):
            pool = gw[(gw.len == L) & (gw.base == base)]
            for r in pool.itertuples():
                if near_target(r.genome_start): continue
                ctrls.append(dict(genome_start=r.genome_start, genome_end=r.genome_end,
                                  ref_len=L, genes="", sco="", gname="", region="",
                                  locus_class="control", base=base))
    ctl = pd.DataFrame(ctrls)
    loci = pd.concat([agg, ctl], ignore_index=True)
    loci["locus_id"] = [f"{c}_{b}{l}_{s}" for c,b,l,s in
        zip(loci.locus_class.str[:4], loci.base, loci.ref_len, loci.genome_start)]
    return loci

def tract_segment(read, la0, ra0):
    """Return the read bases strictly between reference positions la0 and ra0
    (0-based, both must be aligned as M/=/X). CIGAR walk with early exit."""
    rpos = read.reference_start
    qpos = 0
    lq = rq = None
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):                      # M / = / X
            if lq is None and rpos <= la0 < rpos + ln:
                lq = qpos + (la0 - rpos)
            if rpos <= ra0 < rpos + ln:
                rq = qpos + (ra0 - rpos)
                break
            rpos += ln; qpos += ln
        elif op in (1, 4):                       # I / S
            qpos += ln
        elif op in (2, 3):                       # D / N
            if rpos <= la0 < rpos + ln or rpos <= ra0 < rpos + ln:
                return None                      # anchor deleted in this read
            rpos += ln
        elif op == 5:                            # H
            pass
        if rpos > ra0 and rq is None:
            return None
    if lq is None or rq is None or rq <= lq:
        return None
    return read.query_sequence[lq + 1: rq].upper()

def scan_bam(args):
    lib, loci = args
    bam = pysam.AlignmentFile(f"{BAMDIR}/{lib}_mapped.bam", "rb")
    rows = []
    for r in loci.itertuples():
        ts, te, base = int(r.genome_start), int(r.genome_end), r.base
        la0, ra0 = ts - 2, te                    # 0-based flanking anchors
        n_span = 0
        for read in bam.fetch(CHROM, max(0, la0 - ANCHOR), ra0 + ANCHOR):
            if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
            if read.mapping_quality < MINMAPQ: continue
            if read.reference_start > la0 - ANCHOR or read.reference_end < ra0 + ANCHOR + 1:
                continue
            seg = tract_segment(read, la0, ra0)
            if seg is None: continue
            n_span += 1
            runs = [len(m.group(0)) for m in re.finditer(base + "+", seg)] or [0]
            rows.append((lib, r.locus_id, r.locus_class, base, int(r.ref_len),
                         len(seg), max(runs), int(seg == base * len(seg))))
        rows.append((lib, r.locus_id, "_NSPAN", base, int(r.ref_len), -1, -1, n_span))
    bam.close()
    return rows

def main():
    loci = pd.read_csv(f"{OUT}/tables/e2_2_loci_261104.tsv", sep="\t")
    sys.stderr.write(f"loci: {loci.locus_class.value_counts().to_dict()}\n")
    with Pool(9) as p:
        res = p.map(scan_bam, [(l, loci) for l in LIBS])
    recs = [x for sub in res for x in sub]
    df = pd.DataFrame(recs, columns=["lib","locus_id","locus_class","base","ref_len",
                                     "seg_len","max_run","pure"])
    nspan = df[df.locus_class == "_NSPAN"]
    df = df[df.locus_class != "_NSPAN"]
    df.to_csv(f"{OUT}/tables/e2_2_perread_lengths_261104.tsv.gz", sep="\t",
              index=False, compression="gzip")
    nspan.groupby(["locus_id"]).pure.sum().rename("n_spanning").reset_index().to_csv(
        f"{OUT}/tables/e2_2_locus_depth_261104.tsv", sep="\t", index=False)
    print("per-read measurements:", len(df))
    print(df.groupby("locus_class").size().to_string())
    print("candidate loci depth:")
    c = df[df.locus_class=="candidate"].groupby("locus_id").size()
    print(c.to_string())
    print("control loci: n=%d median depth=%.0f" % (
        df[df.locus_class=="control"].locus_id.nunique(),
        df[df.locus_class=="control"].groupby("locus_id").size().median()))
    print("peak_rss_gb", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**3,
          "children", resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024**3)

if __name__ == "__main__":
    main()
