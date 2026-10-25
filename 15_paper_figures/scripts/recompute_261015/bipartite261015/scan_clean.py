"""補正後の未帰属 6mA 部位集合に対する二分割モチーフ探索＋連続 k-mer 走査。
判定はヒット率（そのモチーフのゲノム出現のうち高信頼 6mA である割合）。
ベンチマーク: TGGCCGGC 4mC = 90%（本物）、CGGCAACC = 9.0%（陰性）。
対照は GC 含量・鎖方向を整合させた同数のランダムアデニン。"""
import numpy as np, pandas as pd, re, sys, pathlib, time
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0,'supp_261014' if pathlib.Path('supp_261014').exists() else '../../supp_261014')
HERE=pathlib.Path(__file__).resolve().parent
GEN="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna"
g="".join(l.strip() for l in open(GEN) if not l.startswith(">")).upper()
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
rcg=rc(g); L=len(g)
MAP=np.full(256,-1,np.int64)
for i,b in enumerate("ACGT"): MAP[ord(b)]=i
GA=MAP[np.frombuffer(g.encode(),np.uint8)]; RA=MAP[np.frombuffer(rcg.encode(),np.uint8)]
sites=np.load(HERE.parent/"corrected_unassigned_T1_6mA.npy")
W=25
fg=[]
for st,plus in sites:
    if st-W<0 or st+W+1>L: continue
    w=g[st-W:st+W+1]
    fg.append(w if plus else rc(w))
assert all(w[W]=="A" for w in fg), "中心がアデニンでない窓がある"
rng=np.random.default_rng(0)
arr=np.frombuffer(g.encode(),np.uint8)
pool={True:np.flatnonzero(arr==ord("A")),False:np.flatnonzero(arr==ord("T"))}
bg=[]
for (st,plus),w in zip(sites,fg):
    tgt=(w.count("G")+w.count("C"))/len(w)
    for _ in range(400):
        p=int(rng.choice(pool[bool(plus)]))
        if p-W<0 or p+W+1>L: continue
        v=g[p-W:p+W+1]; v=v if plus else rc(v)
        if abs((v.count("G")+v.count("C"))/len(v)-tgt)<=0.02: bg.append(v); break
print("fg",len(fg),"bg",len(bg),flush=True)
K=(3,4,5); SP=range(2,13); MINFG=15; MINGEN=30
def enum_sp(seqs):
    C=Counter()
    for w in seqs:
        for k1 in K:
            for k2 in K:
                for sp in SP:
                    pl=k1+sp+k2
                    for s0 in range(W-pl+1,W+1):
                        if s0<0 or s0+pl>len(w): continue
                        off=W-s0
                        if k1<=off<k1+sp: continue
                        C[(k1,sp,k2,w[s0:s0+k1],w[s0+k1+sp:s0+pl],off)]+=1
    return C
def enum_contig(seqs,ks=range(4,9)):
    C=Counter()
    for w in seqs:
        for k in ks:
            for s0 in range(W-k+1,W+1):
                if s0<0 or s0+k>len(w): continue
                C[(k,w[s0:s0+k],W-s0)]+=1
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
def census_sp(cls):
    k1,sp,k2=cls; tot=np.zeros(4**(k1+k2),np.int64)
    for a in (GA,RA):
        Lc=kmer_codes(a,k1); Rc=kmer_codes(a,k2)
        n=len(a)-(k1+sp+k2)+1
        l=Lc[:n]; r=Rc[k1+sp:k1+sp+n]; m=(l>=0)&(r>=0)
        tot+=np.bincount(l[m]*(4**k2)+r[m],minlength=4**(k1+k2))
    return cls,tot
def census_k(k):
    tot=np.zeros(4**k,np.int64)
    for a in (GA,RA):
        c=kmer_codes(a,k); tot+=np.bincount(c[c>=0],minlength=4**k)
    return k,tot
if __name__=="__main__":
    t0=time.time()
    fgS,bgS=enum_sp(fg),enum_sp(bg); fgK,bgK=enum_contig(fg),enum_contig(bg)
    cls=sorted({(a,b,c) for (a,b,c,_,_,_),n in list(fgS.items())+list(bgS.items()) if n>=MINFG})
    with Pool(7) as p:
        resS=dict(p.map(census_sp,cls)); resK=dict(p.map(census_k,list(range(4,9))))
    rows=[]
    for lab,C in (("fg",fgS),("bg",bgS)):
        for (k1,sp,k2,l,r,off),n in C.items():
            if n<MINFG or (k1,sp,k2) not in resS: continue
            gn=int(resS[(k1,sp,k2)][code(l)*(4**k2)+code(r)])
            if gn<MINGEN: continue
            rows.append(dict(set=lab,kind="bipartite",pattern=f"{l}{'N'*sp}{r}",
                             A_offset=off,n_sites=n,n_genome=gn,hit_rate_pct=100*n/gn))
    for lab,C in (("fg",fgK),("bg",bgK)):
        for (k,s,off),n in C.items():
            if n<MINFG: continue
            gn=int(resK[k][code(s)])
            if gn<MINGEN: continue
            rows.append(dict(set=lab,kind="contiguous",pattern=s,A_offset=off,
                             n_sites=n,n_genome=gn,hit_rate_pct=100*n/gn))
    t=pd.DataFrame(rows).sort_values("hit_rate_pct",ascending=False)
    t.to_csv(HERE.parent/"motif_scan_clean_T1_261015.tsv",sep="\t",index=False)
    for kind in ("bipartite","contiguous"):
        for lab in ("fg","bg"):
            s=t[(t.kind==kind)&(t.set==lab)]
            print(f"--- {kind} {lab} (n={len(s)})"); print(s.head(8).to_string(index=False))
    print("s",round(time.time()-t0,1))
