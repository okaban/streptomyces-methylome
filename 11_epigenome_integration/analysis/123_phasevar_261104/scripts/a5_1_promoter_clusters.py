#!/usr/bin/env python
"""A5-1: promoters (-300..TSS) carrying clustered modifiable positions, i.e. a
pap/Dam-GATC-like switch architecture. Modifiable position = an instance of a
methylated recognition sequence (TGGCCGGC 4mC-C5; AAGCCCG 4mC-C5 / 6mA-A1),
counted on both strands."""
import re, sys, resource
import numpy as np, pandas as pd

BASE = "/Users/okaban/bioinfo/rna-seq"
OUT  = f"{BASE}/11_epigenome_integration/analysis/123_phasevar_261104"
REF  = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
PROM = 300
MOTIFS = {"TGGCCGGC": "4mC", "AAGCCCG": "4mC+6mA"}
NPERM  = 10000
rng = np.random.default_rng(20261104)

COMP = str.maketrans("ACGTN","TGCAN")
def rc(s): return s.translate(COMP)[::-1]

def load_genome(p):
    s=[]
    for ln in open(p):
        if ln[0]==">": continue
        s.append(ln.strip())
    return "".join(s).upper()

def count_motifs(seq):
    """instances on both strands; returns list of (offset, motif, strand)"""
    hits=[]
    for m, _ in MOTIFS.items():
        for mm in re.finditer("(?=%s)" % m, seq):
            hits.append((mm.start(), m, "+"))
        r = rc(m)
        if r != m:
            for mm in re.finditer("(?=%s)" % r, seq):
                hits.append((mm.start(), m, "-"))
    return sorted(hits)

def main():
    g = load_genome(REF); L=len(g)
    tss = pd.read_csv(f"{BASE}/11_epigenome_integration/analysis/18_tss_analyses/"
                      "jeong2016_all_tss.csv")
    sys.stderr.write(f"TSS rows {len(tss)}\n")

    rows=[]
    for r in tss.itertuples():
        t = int(r.tss_position)
        if r.strand == "+": a,b = t-PROM, t-1
        else:               a,b = t+1, t+PROM
        if a < 1 or b > L: continue
        fwd = g[a-1:b]
        s = fwd if r.strand=="+" else rc(fwd)
        gc = (s.count("G")+s.count("C"))/len(s)
        hits = count_motifs(s)
        pos = [h[0] for h in hits]
        minsp = min([pos[i+1]-pos[i] for i in range(len(pos)-1)], default=np.nan)
        rows.append(dict(tss_id=r.tss_id, gene_id=r.gene_id, sco_gene=r.sco_gene,
            tss_position=t, strand=r.strand, prom_start=a, prom_end=b, gc=gc,
            n_sites=len(hits),
            n_TGGCCGGC=sum(1 for h in hits if h[1]=="TGGCCGGC"),
            n_AAGCCCG=sum(1 for h in hits if h[1]=="AAGCCCG"),
            min_spacing=minsp,
            offsets=";".join(str(p) for p in pos)))
    pr = pd.DataFrame(rows)
    pr.to_csv(f"{OUT}/tables/a5_1_promoter_sites_261104.tsv", sep="\t", index=False)

    N = len(pr); total = int(pr.n_sites.sum())
    obs = pr.n_sites.value_counts().sort_index()

    # --- Null 1 (as specified): random placement of the observed total across promoters
    lam = total / N
    from scipy.stats import poisson
    exp1 = {k: N*poisson.pmf(k, lam) for k in range(0, int(pr.n_sites.max())+1)}

    # --- Null 2 (GC-matched 300-bp genomic windows): controls the real confounder
    starts = rng.integers(1, L-PROM, size=NPERM)
    nullrows=[]
    for a in starts:
        s = g[a-1:a-1+PROM]
        gc = (s.count("G")+s.count("C"))/PROM
        nullrows.append((gc, len(count_motifs(s))))
    nu = pd.DataFrame(nullrows, columns=["gc","n_sites"])
    bins = np.quantile(pr.gc, np.linspace(0,1,11)); bins[0]-=1e-9; bins[-1]+=1e-9
    pr["gcbin"] = pd.cut(pr.gc, bins, labels=False)
    nu["gcbin"] = pd.cut(nu.gc, bins, labels=False)

    # GC-matched expected count distribution and per-promoter tail p
    pvals=[]; expmatch=[]
    for i, r in pr.iterrows():
        pool = nu[nu.gcbin == r.gcbin].n_sites.values
        if len(pool) < 50: pool = nu.n_sites.values
        pvals.append(float((pool >= r.n_sites).mean()))
        expmatch.append(float(pool.mean()))
    pr["p_gcmatched"] = pvals; pr["exp_gcmatched"] = expmatch
    from statsmodels.stats.multitest import multipletests
    pr["q_gcmatched"] = multipletests(pr.p_gcmatched, method="fdr_bh")[1]
    pr.drop(columns=["gcbin"]).to_csv(
        f"{OUT}/tables/a5_1_promoter_sites_261104.tsv", sep="\t", index=False)

    clus = pr[pr.n_sites >= 2].copy()
    tight = clus[clus.min_spacing <= 50]
    sig = pr[(pr.n_sites >= 2) & (pr.q_gcmatched < 0.05)]
    clus.sort_values(["n_sites","min_spacing"], ascending=[False,True]).to_csv(
        f"{OUT}/tables/a5_1_clustered_promoters_261104.tsv", sep="\t", index=False)

    with open(f"{OUT}/tables/a5_1_summary_261104.tsv","w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"promoters_tested\t{N}\n")
        fh.write(f"total_modifiable_sites\t{total}\n")
        fh.write(f"mean_sites_per_promoter\t{lam:.4f}\n")
        fh.write(f"null2_gcmatched_mean_sites_per_window\t{nu.n_sites.mean():.4f}\n")
        for k in sorted(obs.index):
            fh.write(f"obs_n_promoters_with_{k}_sites\t{obs[k]}\n")
            fh.write(f"null1_poisson_expected_{k}_sites\t{exp1.get(k,0):.2f}\n")
            fh.write(f"null2_gcmatched_expected_{k}_sites\t{N*(nu.n_sites==k).mean():.2f}\n")
        fh.write(f"promoters_ge2_sites\t{len(clus)}\n")
        fh.write(f"null1_poisson_expected_ge2\t{N*(1-poisson.cdf(1,lam)):.2f}\n")
        fh.write(f"null2_gcmatched_expected_ge2\t{N*(nu.n_sites>=2).mean():.2f}\n")
        fh.write(f"promoters_ge2_sites_within_50bp\t{len(tight)}\n")
        fh.write(f"promoters_ge2_and_qlt0.05_gcmatched\t{len(sig)}\n")
        fh.write(f"promoters_ge3_sites\t{(pr.n_sites>=3).sum()}\n")
        fh.write(f"null2_gcmatched_expected_ge3\t{N*(nu.n_sites>=3).mean():.2f}\n")
        fh.write(f"max_sites_in_a_promoter\t{pr.n_sites.max()}\n")
        fh.write(f"min_q_gcmatched\t{pr.q_gcmatched.min():.4g}\n")
    print(open(f"{OUT}/tables/a5_1_summary_261104.tsv").read())
    print("peak_rss_gb", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**3)

main()
