import pandas as pd, numpy as np, re, time, pickle
t0=time.time()
GEN="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna"
g="".join(l.strip() for l in open(GEN) if not l.startswith(">")).upper()
L=len(g); comp=str.maketrans("ACGT","TGCA")
rc=lambda s: s.translate(comp)[::-1]
PD="/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015"
COLS=["start","code","strand","Nvalid","Nmod","Nfail"]
def load6(lib):
    d=pd.read_csv(f"{PD}/{lib}_pileup.bed",sep=r"\s+",header=None,usecols=[1,3,5,9,11,15],names=COLS)
    return d[d.code=="a"].groupby(["start","strand"])[["Nvalid","Nmod","Nfail"]].sum()
acc=[load6(l) for l in ("1-1","1-2","1-3")]; idx=acc[0].index
nv=sum(a["Nvalid"].reindex(idx).fillna(0) for a in acc)
nm=sum(a["Nmod"].reindex(idx).fillna(0) for a in acc)
dep=nv+sum(a["Nfail"].reindex(idx).fillna(0) for a in acc)
df=pd.DataFrame({"dep":dep,"lvl":nm/nv.replace(0,np.nan)}).reset_index()
df=df[(df.dep>=30)&(df.lvl>=0.9)]
aag=set()
for m in re.finditer("(?=AAGCCCG)",g): aag |= {m.start(),m.start()+1}
rcg=rc(g)
for m in re.finditer("(?=AAGCCCG)",rcg): aag |= {L-1-m.start(), L-1-(m.start()+1)}
un=df[~df.start.isin(aag)]
print("pooled T1 hi-conf:",len(df),"assigned:",int(df.start.isin(aag).sum()),"unassigned:",len(un))
W=25
fg=[]
for s,st in zip(un.start,un.strand):
    if s-W<0 or s+W+1>L: continue
    w=g[s-W:s+W+1]
    fg.append(rc(w) if st=="-" else w)
rng=np.random.default_rng(0)
arr=np.frombuffer(g.encode(),dtype=np.uint8)
pool={"+":np.flatnonzero(arr==ord("A")),"-":np.flatnonzero(arr==ord("T"))}
gcs=[(w.count("G")+w.count("C"))/len(w) for w in fg]
bg=[]
for st,target in zip(un.strand,gcs):
    for _ in range(400):
        p=int(rng.choice(pool[st]))
        if p-W<0 or p+W+1>L: continue
        w=g[p-W:p+W+1]
        if st=="-": w=rc(w)
        if abs((w.count("G")+w.count("C"))/len(w)-target)<=0.02: bg.append(w); break
print("fg",len(fg),"bg",len(bg),"prep s",round(time.time()-t0,1))
for name,S in (("unassigned_T1",fg),("background_T1",bg)):
    with open(name+".fa","w") as f:
        for i,w in enumerate(S): f.write(f">{name[0]}{i}\n{w}\n")
pickle.dump({"fg":fg,"bg":bg},open("sets.pkl","wb"))
# ---- timing: one genome count for a spaced pattern
t1=time.time()
pat="AACG"+"."*7+"GTGC"
n=len(re.findall("(?="+pat+")",g))+len(re.findall("(?="+pat+")",rcg))
print("one genome regex count:",n,"in",round(time.time()-t1,3),"s")
# ---- timing: foreground pattern enumeration for 50 windows
t2=time.time()
from collections import Counter
C=Counter()
for w in fg[:50]:
    A=W
    for k1 in (3,4,5):
        for k2 in (3,4,5):
            for sp in range(2,13):
                plen=k1+sp+k2
                for st in range(A-plen+1, A+1):
                    if st<0 or st+plen>len(w): continue
                    off=A-st
                    if not (off<k1 or off>=k1+sp):   # A must sit in a block
                        continue
                    C[(w[st:st+k1],sp,w[st+k1+sp:st+plen],off)]+=1
print("enum 50 windows:",len(C),"patterns in",round(time.time()-t2,2),"s")
