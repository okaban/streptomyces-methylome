#!/usr/bin/env python
"""E2-1: identify homopolymer runs and short tandem repeats (SSR) in M145 MTase
catalogue genes + BREX locus genes (gene body and promoter), as candidate
slipped-strand phase-variation switches."""
import os, re, sys, resource
import pandas as pd

BASE = "/Users/okaban/bioinfo/rna-seq"
OUT  = f"{BASE}/11_epigenome_integration/analysis/123_phasevar_261104"
REF  = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"

HP_MIN   = 6     # homopolymer run minimum length reported
SSR_UNIT = range(2, 7)
SSR_MINCOPY = 3
SSR_MINLEN  = 9
PROM = 300

def load_genome(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"): continue
            seq.append(line.strip())
    return "".join(seq).upper()

COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(COMP)[::-1]

def find_homopolymers(seq, minlen=HP_MIN):
    out = []
    for m in re.finditer(r"(A+|C+|G+|T+)", seq):
        if m.end() - m.start() >= minlen:
            out.append((m.start(), m.end(), m.group(1)[0], m.end()-m.start()))
    return out

def find_ssr(seq):
    """perfect tandem repeats, unit 2..6, >=3 copies, >=9 bp; greedy longest, non-redundant"""
    hits = []
    for u in SSR_UNIT:
        pat = re.compile(r"((.{%d})\2{%d,})" % (u, SSR_MINCOPY-1))
        for m in pat.finditer(seq):
            full, unit = m.group(1), m.group(2)
            if len(full) < SSR_MINLEN: continue
            if len(set(unit)) == 1: continue      # pure homopolymer, handled above
            # reject units that are themselves a smaller repeat (e.g. ATAT as unit 4)
            red = False
            for d in range(1, u):
                if u % d == 0 and unit == unit[:d] * (u // d): red = True
            if red: continue
            hits.append((m.start(), m.end(), unit, u, len(full)//u, len(full)))
    # drop hits fully contained in a longer hit
    hits.sort(key=lambda h: (h[0], -(h[1]-h[0])))
    keep = []
    for h in hits:
        if any(h[0] >= k[0] and h[1] <= k[1] for k in keep): continue
        keep.append(h)
    return keep

def main():
    g = load_genome(REF)
    L = len(g)
    sys.stderr.write(f"genome {L} bp\n")

    ann = pd.read_csv(f"{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/"
                      "gene_annotation_basic.tsv", sep="\t", low_memory=False)
    ann = ann[ann.contig == "NC_003888.3"].copy()

    mt = pd.read_csv(f"{BASE}/11_epigenome_integration/analysis/33_GCCGGC_RM_identification/"
                     "tables/M145_MTase_complete_catalog.tsv", sep="\t")
    brex = pd.read_csv(f"{BASE}/17_brex_261103/tables/brex_locus_main_261103.tsv", sep="\t")

    targets = {}
    for gid in mt.locus_tag.astype(str):
        targets[gid] = "MTase_catalog"
    for gid in brex.gene_id.astype(str):
        targets[gid] = "BREX_locus" if gid not in targets else "MTase_catalog+BREX_locus"

    tss = pd.read_csv(f"{BASE}/11_epigenome_integration/analysis/18_tss_analyses/"
                      "jeong2016_all_tss.csv")
    tss = tss.dropna(subset=["gene_id"])
    tssmap = {}
    for gid, sub in tss.groupby("gene_id"):
        # strongest TSS for the gene
        r = sub.sort_values("abundance", ascending=False).iloc[0]
        tssmap[gid] = (int(r.tss_position), r.strand)

    sub = ann[ann.gene_id.isin(targets)].copy()
    missing = sorted(set(targets) - set(sub.gene_id))
    sys.stderr.write(f"targets={len(targets)} found_in_annotation={len(sub)} missing={len(missing)}\n")

    rows = []
    for _, r in sub.iterrows():
        gid, st, en, strand = r.gene_id, int(r.start), int(r.end), r.strand
        regions = []
        # gene body, in transcribed orientation
        regions.append(("gene_body", st, en))
        # promoter
        if gid in tssmap:
            t, _s = tssmap[gid]
            p = (t - PROM, t) if strand == "+" else (t, t + PROM)
            regions.append(("promoter_TSS", max(1, p[0]), min(L, p[1])))
        else:
            p = (st - PROM, st - 1) if strand == "+" else (en + 1, en + PROM)
            regions.append(("promoter_CDSupstream", max(1, p[0]), min(L, p[1])))

        for rtype, a, b in regions:
            if b <= a: continue
            fwd = g[a-1:b]
            s = fwd if strand == "+" else rc(fwd)
            for (i, j, base, ln) in find_homopolymers(s):
                gs, ge = (a + i, a + j - 1) if strand == "+" else (b - j + 1, b - i)
                rows.append(dict(gene_id=gid, old_locus_tag=r.old_locus_tag,
                    gene_name=r.gene_name, product=r.product, source=targets[gid],
                    gene_start=st, gene_end=en, strand=strand, region=rtype,
                    repeat_class="homopolymer", unit=base, unit_len=1,
                    copies=ln, repeat_len=ln,
                    genome_start=gs, genome_end=ge,
                    offset_in_region=i,
                    offset_in_gene=(i if rtype == "gene_body" else None),
                    frame_pos=(i % 3 if rtype == "gene_body" else None),
                    frameshift_capable=True))
            for (i, j, unit, u, cp, ln) in find_ssr(s):
                gs, ge = (a + i, a + j - 1) if strand == "+" else (b - j + 1, b - i)
                rows.append(dict(gene_id=gid, old_locus_tag=r.old_locus_tag,
                    gene_name=r.gene_name, product=r.product, source=targets[gid],
                    gene_start=st, gene_end=en, strand=strand, region=rtype,
                    repeat_class="SSR", unit=unit, unit_len=u,
                    copies=cp, repeat_len=ln,
                    genome_start=gs, genome_end=ge,
                    offset_in_region=i,
                    offset_in_gene=(i if rtype == "gene_body" else None),
                    frame_pos=(i % 3 if rtype == "gene_body" else None),
                    frameshift_capable=(u % 3 != 0)))

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/tables/e2_1_repeats_all_261104.tsv", sep="\t", index=False)

    # genome-wide homopolymer length background (for context + control pool in E2-2)
    hp = find_homopolymers(g, minlen=5)
    bg = pd.DataFrame(hp, columns=["start0", "end0", "base", "len"])
    bg["genome_start"] = bg.start0 + 1
    bg["genome_end"]   = bg.end0
    bg[["genome_start", "genome_end", "base", "len"]].to_csv(
        f"{OUT}/tables/e2_1_genomewide_homopolymers_261104.tsv", sep="\t", index=False)

    with open(f"{OUT}/logs/e2_1.log", "w") as fh:
        fh.write(f"genome_bp\t{L}\n")
        fh.write(f"target_genes\t{len(targets)}\n")
        fh.write(f"genes_in_annotation\t{len(sub)}\n")
        fh.write(f"missing_gene_ids\t{','.join(missing)}\n")
        fh.write(f"repeats_total\t{len(df)}\n")
        if len(df):
            for k, v in df.groupby(["region", "repeat_class"]).size().items():
                fh.write(f"count\t{k[0]}\t{k[1]}\t{v}\n")
            fh.write(f"homopolymer_len_max\t{df[df.repeat_class=='homopolymer'].repeat_len.max()}\n")
        fh.write(f"genomewide_hp_ge5\t{len(bg)}\n")
        for n in (6, 7, 8, 9, 10, 11, 12):
            fh.write(f"genomewide_hp_ge{n}\t{(bg.len >= n).sum()}\n")

    print("repeats", len(df))
    print(df.groupby(["region", "repeat_class"]).size().to_string())
    print("peak_rss_gb", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3)

main()
