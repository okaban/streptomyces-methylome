#!/usr/bin/env python3
"""121 - robustness of the AAGCCCG count ~ expression association."""
import os, numpy as np, pandas as pd
from scipy import stats
BASE="/Users/okaban/bioinfo/rna-seq"
HERE=f"{BASE}/11_epigenome_integration/analysis/121_motif_count_261104"
ANN=f"{BASE}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
FC=f"{BASE}/03_quantification/data/featureCounts/featureCounts_M145_s2.txt"
def pspear(x,y,Z):
    X=np.column_stack([stats.rankdata(x),stats.rankdata(y)]+[stats.rankdata(Z[:,j]) for j in range(Z.shape[1])])
    C=np.corrcoef(X-X.mean(0),rowvar=False);P=np.linalg.inv(C)
    r=-P[0,1]/np.sqrt(P[0,0]*P[1,1]);df=len(x)-2-Z.shape[1]
    return float(r),float(2*stats.t.sf(abs(r*np.sqrt(df/(1-r*r))),df))
ann=pd.read_csv(ANN,sep='\t');ann=ann[ann.contig=='NC_003888.3']
g=pd.read_csv(f"{HERE}/tables/per_gene_counts_261104.tsv",sep='\t').merge(ann[['gene_id','gene_biotype','product']],on='gene_id')
fc=pd.read_csv(FC,sep='\t',comment='#');lib=[c for c in fc.columns if c.endswith('.bam')]
lab=[os.path.basename(c).split('.')[0].replace('M145_','') for c in lib]
rpk=fc[lib].to_numpy(float)/fc['Length'].to_numpy(float)[:,None];tpm=rpk/rpk.sum(0)*1e6
T=pd.DataFrame(tpm,columns=lab);T['gene_id']=fc.Geneid.values
T['tpm_T1']=T[[c for c in lab if c.startswith('1_')]].mean(1)
g=g.merge(T[['gene_id','tpm_T1']],on='gene_id')
rows=[]
def add(tag,d,pred):
    x=d[pred].to_numpy(float);y=np.log2(d.tpm_T1+1).to_numpy(float)
    Z=np.column_stack([np.log10(d.glen.to_numpy(float)),d.gene_gc.to_numpy(float)])
    r0=stats.spearmanr(x,y);r1=pspear(x,y,Z)
    rows.append(dict(subset=tag,predictor=pred,n=len(d),naive_rho=round(r0.statistic,4),
                     naive_p=r0.pvalue,adj_rho=round(r1[0],4),adj_p=r1[1]))
for pred in ['n_AAGCCCG_4mC_gene','n_TGGCCGGC_4mC_gene']:
    add('all_chromosomal_genes',g,pred)
    add('protein_coding_only',g[g.gene_biotype=='protein_coding'],pred)
    add('exclude_count_ge3',g[g[pred]<3],pred)
    add('expressed_TPM_ge1',g[g.tpm_T1>=1],pred)
    add('TSS_defined_promoter',g[g.has_tss],pred)
R=pd.DataFrame(rows);R.to_csv(f"{HERE}/tables/robustness_261104.tsv",sep='\t',index=False)
print(R.to_string(index=False))
hi=g[g.n_AAGCCCG_4mC_gene>=3][['gene_id','old_locus_tag','glen','n_AAGCCCG_4mC_gene','tpm_T1','product']]
hi.sort_values('tpm_T1',ascending=False).to_csv(f"{HERE}/tables/high_count_genes_261104.tsv",sep='\t',index=False)
