import numpy as np, pickle, time, itertools, pandas as pd
from multiprocessing import Pool
from collections import Counter
t0=time.time()
GEN="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna"
g="".join(l.strip() for l in open(GEN) if not l.startswith(">")).upper()
comp=str.maketrans("ACGT","TGCA"); rcg=g.translate(comp)[::-1]
MAP=np.full(256,-1,np.int64)
for i,b in enumerate("ACGT"): MAP[ord(b)]=i
GA=MAP[np.frombuffer(g.encode(),np.uint8)]
RA=MAP[np.frombuffer(rcg.encode(),np.uint8)]
S=pickle.load(open("sets.pkl","rb")); W=25
K=(3,4,5); SP=range(2,13); MINFG=15; MINGEN=30

def enum(seqs):
    C=Counter()
    for w in seqs:
        A=W
        for k1 in K:
            for k2 in K:
                for sp in SP:
                    pl=k1+sp+k2
                    for st in range(A-pl+1,A+1):
                        if st<0 or st+pl>len(w): continue
                        off=A-st
                        if k1<=off<k1+sp: continue
                        C[(k1,sp,k2,w[st:st+k1],w[st+k1+sp:st+pl],off)]+=1
    return C

def kmer_codes(a,k):
    n=len(a); out=np.zeros(n-k+1,np.int64); bad=np.zeros(n-k+1,bool)
    for i in range(k):
        seg=a[i:n-k+1+i]; out=out*4+np.maximum(seg,0); bad|=(seg<0)
    out[bad]=-1; return out

def code(s):
    v=0
    for ch in s: v=v*4+"ACGT".index(ch)
    return v

def census(cls):
    k1,sp,k2=cls
    tot=np.zeros(4**(k1+k2),np.int64)
    for a in (GA,RA):
        Lc=kmer_codes(a,k1); Rc=kmer_codes(a,k2)
        n=len(a)-(k1+sp+k2)+1
        l=Lc[:n]; r=Rc[k1+sp:k1+sp+n]
        m=(l>=0)&(r>=0)
        j=l[m]*(4**k2)+r[m]
        tot+=np.bincount(j,minlength=4**(k1+k2))
    return cls,tot

if __name__=="__main__":
    fgC=enum(S["fg"]); bgC=enum(S["bg"])
    print("fg patterns",len(fgC),"bg",len(bgC),"t",round(time.time()-t0,1),flush=True)
    cand=set((k1,sp,k2) for (k1,sp,k2,_,_,_),c in fgC.items() if c>=MINFG)
    cand|=set((k1,sp,k2) for (k1,sp,k2,_,_,_),c in bgC.items() if c>=MINFG)
    classes=sorted(cand)
    print("classes to census:",len(classes),flush=True)
    with Pool(8) as p:
        res=dict(p.map(census,classes))
    print("census done t",round(time.time()-t0,1),flush=True)
    rows=[]
    for label,C in (("fg",fgC),("bg",bgC)):
        for (k1,sp,k2,l,r,off),c in C.items():
            if c<MINFG or (k1,sp,k2) not in res: continue
            gc=int(res[(k1,sp,k2)][code(l)*(4**k2)+code(r)])
            if gc<MINGEN: continue
            rows.append(dict(set=label,pattern=f"{l}{'N'*sp}{r}",k1=k1,spacer=sp,k2=k2,
                             A_offset=off,n_sites=c,n_genome=gc,hit_rate=100*c/gc))
    t=pd.DataFrame(rows).sort_values("hit_rate",ascending=False)
    t.to_csv("bipartite_scan_T1_261015.tsv",sep="\t",index=False)
    print("FG TOP"); print(t[t.set=="fg"].head(15).to_string(index=False))
    print("BG TOP"); print(t[t.set=="bg"].head(8).to_string(index=False))
    print("total s",round(time.time()-t0,1))
