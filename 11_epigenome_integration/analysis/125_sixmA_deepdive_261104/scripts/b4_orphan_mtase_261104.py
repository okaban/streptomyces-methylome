#!/usr/bin/env python3
"""B4 - orphan MTase vs R-M-paired MTase, and recognition-site avoidance.

Reuses the 2026-11-02 adjacency scan (Writing/mtase_nuclease_adjacency_261102.tsv);
nothing in that scan is recomputed. Adds (a) an orphan/paired split of the
129-gene annotation catalogue at a stated threshold, (b) the status of the
three adenine-MTase candidates that clear the 348 aa floor, and (c) a
recognition-site avoidance test that does not need an enzyme assignment.
"""
import re, numpy as np, pandas as pd
RNA="/Users/okaban/bioinfo/rna-seq"
AN=f"{RNA}/11_epigenome_integration/analysis"
OUT=f"{AN}/125_sixmA_deepdive_261104"; TAB=f"{OUT}/tables"
W="/Users/okaban/obsidian/Research/rna-seq/Writing"
ADJ=pd.read_csv(f"{W}/mtase_nuclease_adjacency_261102.tsv",sep="\t")
CAT=pd.read_csv(f"{AN}/33_GCCGGC_RM_identification/tables/M145_MTase_complete_catalog.tsv",sep="\t")
print("catalogue genes:",len(CAT),"  adjacency pairs:",len(ADJ))

# ---- (a) orphan / paired split of the 129-gene catalogue -------------------
# threshold: a partner gene with nuclease-domain evidence at ANY tier within
# +-3 genes AND <=10,000 bp edge-to-edge (the 261102 scan's own window).
with_nuc = set(ADJ.loc[ADJ.domain_evidence_tier.notna() & (ADJ.domain_evidence_tier!="none"),
                       "mtase_locus_tag"])
cat_lt=set(CAT.locus_tag)
paired = cat_lt & with_nuc
rows=[dict(threshold="+-3 genes and <=10,000 bp edge-to-edge, nuclease domain at any tier",
           catalogue_n=len(cat_lt),
           paired_with_nuclease=len(paired),
           orphan=len(cat_lt)-len(paired),
           orphan_pct=round(100*(len(cat_lt)-len(paired))/len(cat_lt),1))]
# stricter: immediate neighbour (offset +-1) only
imm=set(ADJ.loc[(ADJ.gene_offset.abs()==1)&ADJ.domain_evidence_tier.notna()
                &(ADJ.domain_evidence_tier!="none"),"mtase_locus_tag"])
rows.append(dict(threshold="immediate neighbour (offset +-1) only, nuclease domain at any tier",
    catalogue_n=len(cat_lt),paired_with_nuclease=len(cat_lt&imm),
    orphan=len(cat_lt)-len(cat_lt&imm),
    orphan_pct=round(100*(len(cat_lt)-len(cat_lt&imm))/len(cat_lt),1)))
# strictest: Pfam-GA tier
ga=set(ADJ.loc[ADJ.domain_evidence_tier.astype(str).str.contains("Pfam-GA"),"mtase_locus_tag"])
rows.append(dict(threshold="+-3 genes / <=10,000 bp, nuclease domain at Pfam gathering threshold",
    catalogue_n=len(cat_lt),paired_with_nuclease=len(cat_lt&ga),
    orphan=len(cat_lt)-len(cat_lt&ga),
    orphan_pct=round(100*(len(cat_lt)-len(cat_lt&ga))/len(cat_lt),1)))
O=pd.DataFrame(rows); O.to_csv(f"{TAB}/b4_orphan_split_261104.tsv",sep="\t",index=False)
print("\n=== orphan / R-M-paired split of the 129-gene catalogue ===")
print(O.to_string(index=False))

# ---- (b) the three candidates that clear the 348 aa floor ------------------
CAND={"SC_RS35335":("SCO6627","pglX, complete BREX-2 locus"),
      "SC_RS28835":("SCO5331","pglS, isolated pglX paralogue"),
      "SC_RS17645":("SCO3104","Lutz 2019 candidate")}
sub=ADJ[ADJ.mtase_locus_tag.isin(CAND)][
    ["mtase_locus_tag","mtase_old_locus_tag","mtase_aa","partner_locus_tag",
     "partner_old_locus_tag","partner_aa","nuclease_class","domain_evidence_tier",
     "gene_offset","intervening_bp","orientation"]]
sub.to_csv(f"{TAB}/b4_candidate_neighbourhoods_261104.tsv",sep="\t",index=False)
print("\n=== nuclease neighbourhoods of the 3 adenine-MTase candidates ===")
print(sub.to_string(index=False))
for lt,(sco,note) in CAND.items():
    n=len(sub[sub.mtase_locus_tag==lt])
    print(f"  {lt} ({sco}, {note}): {n} nuclease-evidence neighbour(s)")

# ---- (c) recognition-site avoidance (no enzyme assignment needed) ----------
REF="/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
g="".join(l.strip() for l in open(REF) if not l.startswith(">")).upper()
L=len(g); comp=str.maketrans("ACGT","TGCA"); grc=g.translate(comp)[::-1]
def count(m):
    return len(re.findall("(?="+m+")",g))+len(re.findall("(?="+m+")",grc))
def kcounts(k):
    from collections import Counter
    c=Counter()
    for s in (g,grc):
        for i in range(len(s)-k+1):
            w=s[i:i+k]
            if "N" not in w: c[w]+=1
    return c
KC={k:kcounts(k) for k in range(1,9)}
def karlin(m):
    """maximal-order Markov expectation (Karlin rho): E = (prod of the two
    (k-1)-mers) / (the central (k-2)-mer), counts from both strands."""
    k=len(m)
    num=KC[k-1][m[:-1]]*KC[k-1][m[1:]]
    den=KC[k-2][m[1:-1]]
    return num/den if den else np.nan
rows=[]
for m in ["AAGCCCG","TGGCCGGC","GCCGGC","GGCCGG","AAGCCC","AGCCCG"]:
    o=count(m); e=karlin(m)
    rows.append(dict(motif=m,k=len(m),observed=o,expected_maxorder_markov=round(e,1),
                     rho=round(o/e,3) if e==e else np.nan))
A=pd.DataFrame(rows); A.to_csv(f"{TAB}/b4_site_avoidance_261104.tsv",sep="\t",index=False)
print("\n=== recognition-site avoidance (Karlin rho = obs/exp, maximal-order Markov, both strands) ===")
print(A.to_string(index=False))
# genome-wide rho distribution for same-length words, for calibration
for k in (7,8):
    rs=[]
    for w,o in KC[k].items():
        num=KC[k-1][w[:-1]]*KC[k-1][w[1:]]; den=KC[k-2][w[1:-1]]
        if den and o>=50: rs.append(o/(num/den))
    rs=np.array(rs)
    tgt=[m for m in A.motif if len(m)==k]
    for m in tgt:
        r=float(A.loc[A.motif==m,"rho"])
        print(f"  {m}: rho={r:.3f}, percentile among {len(rs):,} {k}-mers (count>=50) = "
              f"{100*(rs<r).mean():.1f}  [5th/50th/95th pct of rho = "
              f"{np.percentile(rs,5):.2f}/{np.percentile(rs,50):.2f}/{np.percentile(rs,95):.2f}]")
