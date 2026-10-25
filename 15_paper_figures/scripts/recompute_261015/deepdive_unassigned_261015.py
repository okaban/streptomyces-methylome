"""未帰属 6mA の 2 つの有意所見（遺伝子間 2.1 倍・ホモポリマー 3.9 倍）の深掘り。

対照を「GC 整合ランダムアデニン」から「同じ深度条件を満たしながら呼ばれなかったアデニン」に
置き換える。後者が正しい分母で、カバレッジ・検出可能性の交絡を外から持ち込まない。
"""
import numpy as np, pandas as pd, pathlib, re, sys
from scipy import stats
import statsmodels.api as sm
sys.path.insert(0,'supp_261014'); from _gccggc import genome
HERE=pathlib.Path(__file__).resolve().parent; C=HERE/"pileup_cache_261015"
g=genome(); L=len(g)
TP={"T1":["1-1","1-2","1-3"],"T2":["2-1","2-3","2-4"],"T3":["3-2","3-3","3-4"]}
aag_A=set()
for m in re.finditer("AAGCCCG",g): aag_A|={(m.start(),True),(m.start()+1,True)}
for m in re.finditer("CGGGCTT",g): aag_A|={(m.start()+6,False),(m.start()+5,False)}
def pooled6(tp):
    fr=[]
    for lib in TP[tp]:
        z=np.load(C/f"{lib}.npz"); m=(z["code"]==0)
        fr.append(pd.DataFrame({"start":z["start"][m],"plus":z["plus"][m],"Nvalid":z["Nvalid"][m],
            "Nmod":z["Nmod"][m],"Nfail":z["Nfail"][m],"refok":z["refok"][m]}).set_index(["start","plus"]))
    uni=fr[0].index
    for f in fr[1:]: uni=uni.union(f.index)
    s=sum(f[["Nvalid","Nmod","Nfail"]].reindex(uni).fillna(0) for f in fr)
    ok=fr[0]["refok"].reindex(uni)
    for f in fr[1:]: ok=ok.fillna(f["refok"].reindex(uni))
    s["refok"]=ok.fillna(False).astype(bool); s=s.reset_index()
    s["depth"]=s.Nvalid+s.Nfail; s["lvl"]=s.Nmod/s.Nvalid.replace(0,np.nan)
    return s[s.refok & (s.depth>=30) & (s.Nvalid>0)]
e=pooled6("T1")
e["called"]=e.lvl>=0.5
e["in_motif"]=[(p,b) in aag_A for p,b in zip(e.start,e.plus)]
e=e[~e.in_motif]
print("適格アデニン(非モチーフ):",len(e)," うち呼ばれた:",int(e.called.sum()),flush=True)
# ---- 注釈
TSS=("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv")
t=pd.read_csv(TSS); sel=t[t.category.astype(str).str.upper().isin(["P","S"])&t.gene_id.notna()
      &(t.gene_id.astype(str).str.strip().str.len()>0)]
tss=np.sort(sel.tss_position.astype(int).values)
GFF=("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/"
     "data/GCF_000203835.1_ASM20383v1_genomic.gff")
cds=np.zeros(L,bool); rna=np.zeros(L,bool); mge=np.zeros(L,bool)
for line in open(GFF):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    if len(f)<9: continue
    s0,e0,ty,att=int(f[3])-1,int(f[4]),f[2],f[8]
    if ty=="CDS": cds[s0:e0]=True
    if ty in ("rRNA","tRNA","ncRNA","tmRNA"): rna[s0:e0]=True
    if ty=="CDS" and re.search("transposase|integrase|recombinase|IS[0-9]|mobile",att,re.I): mge[s0:e0]=True
pos=e.start.to_numpy()
i=np.searchsorted(tss,pos)
lo=np.where(i>0,pos-tss[np.clip(i-1,0,None)],10**9); hi=np.where(i<len(tss),tss[np.clip(i,None,len(tss)-1)]-pos,10**9)
e["d_tss"]=np.minimum(lo,hi); e["promoter"]=e.d_tss<=300
e["in_cds"]=cds[pos]; e["rna_gene"]=rna[pos]; e["mge"]=mge[pos]
e["intergenic"]=(~e.in_cds)&(~e.promoter)
gc=np.frombuffer(g.encode(),np.uint8)
isgc=((gc==ord("G"))|(gc==ord("C"))).astype(np.int32)
cum=np.concatenate([[0],np.cumsum(isgc)])
e["gc50"]=(cum[np.minimum(pos+26,L)]-cum[np.maximum(pos-25,0)])/(np.minimum(pos+26,L)-np.maximum(pos-25,0))
# ホモポリマー（±5 bp に 5 連続同一塩基）と、その塩基種
hp=np.zeros(len(e),bool); hpb=np.array([""]*len(e),dtype=object)
for idx,p in enumerate(pos):
    s=g[max(0,p-5):p+6]
    for b in "ACGT":
        if b*5 in s: hp[idx]=True; hpb[idx]=b; break
e["homopolymer"]=hp; e["hp_base"]=hpb
# 21-mer の一意性（マッパビリティ代理）
K=21; code=np.full(256,-1,np.int64)
for ii,b in enumerate("ACGT"): code[ord(b)]=ii
a=code[np.frombuffer(g.encode(),np.uint8)]
km=np.zeros(L-K+1,np.int64); bad=np.zeros(L-K+1,bool)
for ii in range(K):
    seg=a[ii:L-K+1+ii]; km=km*4+np.maximum(seg,0); bad|=(seg<0)
km[bad]=-1
canon=np.minimum(km,km)  # 片鎖で十分（反復検出が目的）
u,cnt=np.unique(canon,return_counts=True)
occ=cnt[np.searchsorted(u,canon)]
st=np.clip(pos-K//2,0,L-K)
e["kmer21_occ"]=occ[st]
e["repeat"]=e.kmer21_occ>1
# ---- 検定: 呼ばれた vs 呼ばれなかった（同じ適格集合内）
rows=[]
for feat in ["promoter","in_cds","intergenic","rna_gene","mge","homopolymer","repeat"]:
    a1=int(e[e.called][feat].sum()); b1=int(e.called.sum())-a1
    c1=int(e[~e.called][feat].sum()); d1=int((~e.called).sum())-c1
    orr,p=stats.fisher_exact([[a1,b1],[c1,d1]])
    rows.append(dict(feature=feat,called_pct=round(100*a1/(a1+b1),2),notcalled_pct=round(100*c1/(c1+d1),2),
                     odds_ratio=round(orr,3),p=f"{p:.3g}",n_called=a1))
for col in ["gc50","depth","Nvalid","kmer21_occ","d_tss"]:
    u2,p=stats.mannwhitneyu(e[e.called][col],e[~e.called][col])
    rows.append(dict(feature=col+" (median)",called_pct=round(float(e[e.called][col].median()),4),
                     notcalled_pct=round(float(e[~e.called][col].median()),4),odds_ratio=np.nan,
                     p=f"{p:.3g}",n_called=int(e.called.sum())))
r=pd.DataFrame(rows); print(r.to_string(index=False),flush=True)
# ---- 多変量: GC と深度で調整しても残るか
X=pd.DataFrame({"gc50":e.gc50,"log_depth":np.log10(e.depth),"intergenic":e.intergenic.astype(int),
                "homopolymer":e.homopolymer.astype(int),"promoter":e.promoter.astype(int),
                "repeat":e.repeat.astype(int)})
X=sm.add_constant(X)
m=sm.Logit(e.called.astype(int).values,X.values).fit(disp=0)
co=pd.DataFrame({"term":X.columns,"coef":m.params,"OR":np.exp(m.params),"p":m.pvalues})
print("\n--- ロジスティック回帰（called ~ GC + 深度 + 位置）"); print(co.to_string(index=False),flush=True)
# ---- ホモポリマー部位の質
hpr=[]
for lab,sub in (("homopolymer",e[e.called&e.homopolymer]),("non-homopolymer",e[e.called&~e.homopolymer])):
    hpr.append(dict(group=lab,n=len(sub),median_level=round(float(sub.lvl.median()),3),
                    median_Nvalid=float(sub.Nvalid.median()),median_depth=float(sub.depth.median()),
                    pct_in_cds=round(100*sub.in_cds.mean(),1)))
print("\n--- ホモポリマー部位の質"); print(pd.DataFrame(hpr).to_string(index=False))
print("\nホモポリマー塩基種（呼ばれた部位）:",e[e.called&e.homopolymer].hp_base.value_counts().to_dict())
print("ホモポリマー塩基種（呼ばれなかった適格）:",e[~e.called&e.homopolymer].hp_base.value_counts().to_dict())
# ---- 時点間再現性を層別
U={}
for tp in TP:
    s=pooled6(tp); U[tp]=set(zip(s[s.lvl>=0.5].start,s[s.lvl>=0.5].plus))-aag_A
cal=e[e.called]
for lab,mask in (("intergenic",cal.intergenic),("CDS",cal.in_cds),("homopolymer",cal.homopolymer),
                 ("repeat",cal.repeat),("all",pd.Series(True,index=cal.index))):
    k=set(zip(cal[mask].start,cal[mask].plus))
    if not k: continue
    print(f"再現性 {lab:12s} n={len(k):5d}  T2 でも未帰属 {100*len(k&U['T2'])/len(k):5.1f}%  3 時点 {100*len(k&U['T2']&U['T3'])/len(k):5.1f}%")
r.to_csv(HERE/"deepdive_unassigned_261015.tsv",sep="\t",index=False)
co.to_csv(HERE/"deepdive_logit_261015.tsv",sep="\t",index=False)
e[e.called].to_csv(HERE/"unassigned_T1_sites_annotated_261015.tsv",sep="\t",index=False)
