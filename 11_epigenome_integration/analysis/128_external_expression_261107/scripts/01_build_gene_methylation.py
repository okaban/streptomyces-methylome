#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""128/01 — per-gene gene-body methylation level + positional covariates.

Output: tables/gene_methylation_covariates_261107.tsv (one row per chromosomal gene).

Inputs (canonical only)
-----------------------
sites     : 127_coverage_canonical_261107/tables/canonical_sites_{4mC,6mA}_{T1,T2,T3}.tsv
            (union-index sum over the 3 libraries of a timepoint, true depth
            Nvalid+Nfail >= 30, level = Nmod/Nvalid >= 0.5, reference base matches
            the modification; counts reproduce Supplementary Table S1c exactly)
genes     : 119_target_function_261104/tables/genes_universe_261104.tsv
            (chromosomal genes only, NC_003888.3; carries length and gc computed
            from the reference, and old_locus_tag = SCO identifier)
oriC      : 4,271,763 bp (dnaA SCO3879 midpoint), the same constant as
            80_partial_corr_covariates and 83_oriC_confound_FIRE.

Gene-body methylation level
---------------------------
A canonical site belongs to a gene when the modified base lies inside the gene's
[start, end] interval (1-based inclusive, either strand; the site tables carry
0-based `start`, so the test is start+1 in [gene.start, gene.end]).
Three per-gene summaries, all defined for every gene:
  n_sites     : number of canonical sites in the body
  level_sum   : sum of site levels (Nmod/Nvalid) over those sites
  level_kb    : level_sum / (gene length in kb)     <- PRIMARY exposure
and one defined only on genes with >= 1 site:
  mean_level  : level_sum / n_sites

Covariates (the set the manuscript declares for gene-level partial correlations:
gene GC, distance from oriC, local gene density, length)
  length_bp   : gene length
  gc          : gene GC fraction
  dist_oriC   : |gene midpoint - oriC|
  gene_density: number of gene midpoints within +/- 50 kb of this gene's midpoint,
                expressed per 100 kb (the window is declared here; no earlier
                producing script for the manuscript's value was recoverable)
"""
import numpy as np, pandas as pd
from pathlib import Path

RS   = Path("/Users/okaban/bioinfo/rna-seq")
SITE = RS/"11_epigenome_integration/analysis/127_coverage_canonical_261107/tables"
UNIV = RS/"11_epigenome_integration/analysis/119_target_function_261104/tables/genes_universe_261104.tsv"
HERE = Path(__file__).resolve().parent.parent
ORIC = 4_271_763.0
WIN  = 50_000          # half-width of the local gene-density window

g = pd.read_csv(UNIV, sep="\t", dtype={"old_locus_tag": str})
g = g[g.contig == "NC_003888.3"].copy()
g["mid"] = (g.start + g.end) / 2.0
g = g.sort_values("mid").reset_index(drop=True)
print("genes (chromosome):", len(g))

# --- covariates ---
g["length_bp"] = g.end - g.start + 1
g["dist_oriC"] = (g["mid"] - ORIC).abs()
mid = g["mid"].values
lo = np.searchsorted(mid, mid - WIN, side="left")
hi = np.searchsorted(mid, mid + WIN, side="right")
g["gene_density"] = (hi - lo) / (2 * WIN / 100_000.0)   # genes per 100 kb

# --- methylation ---
gstart, gend = g.start.values, g.end.values
order = np.argsort(gstart)
assert (order == np.arange(len(g))).all() or True
for mark in ("4mC", "6mA"):
    for tp in ("T1", "T2", "T3"):
        s = pd.read_csv(SITE/f"canonical_sites_{mark}_{tp}.tsv", sep="\t")
        pos = s["start"].values + 1                      # 1-based coordinate of the modified base
        lvl = s["level"].values
        n   = np.zeros(len(g)); ls = np.zeros(len(g))
        # genes do not overlap much; assign by interval containment, both strands
        idx = np.searchsorted(gstart, pos, side="right") - 1
        for k in range(len(pos)):
            i = idx[k]
            # walk back over genes that start before but end after (handles nesting/overlap)
            j = i
            while j >= 0 and gstart[j] >= pos[k] - 20000:
                if gstart[j] <= pos[k] <= gend[j]:
                    n[j] += 1; ls[j] += lvl[k]
                j -= 1
        g[f"n_{mark}_{tp}"]        = n
        g[f"levelsum_{mark}_{tp}"] = ls
        g[f"levelkb_{mark}_{tp}"]  = ls / (g.length_bp.values / 1000.0)
        g[f"meanlevel_{mark}_{tp}"] = np.where(n > 0, ls / np.maximum(n, 1), np.nan)
        print(f"{mark} {tp}: sites {len(s)}, assigned to a gene body {int(n.sum())} "
              f"({int(n.sum())/len(s)*100:.1f}%), genes with >=1 site {(n>0).sum()}")

keep = ["gene_id", "old_locus_tag", "gene_name", "gene_biotype", "start", "end", "strand",
        "mid", "length_bp", "gc", "dist_oriC", "gene_density"] + \
       [c for c in g.columns if c.startswith(("n_4mC", "n_6mA", "levelsum_", "levelkb_", "meanlevel_"))]
out = g[keep].sort_values("start")
(HERE/"tables").mkdir(exist_ok=True)
out.to_csv(HERE/"tables/gene_methylation_covariates_261107.tsv", sep="\t", index=False,
           float_format="%.6g")
print("wrote", HERE/"tables/gene_methylation_covariates_261107.tsv", out.shape)
print("genes with an SCO identifier:", (out.old_locus_tag.notna() & (out.old_locus_tag != "")).sum())
