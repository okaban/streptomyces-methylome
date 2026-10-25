"""正典のモチーフ定義規則（F1 最大化）を ±8 塩基の窓で走らせ、AAGCCCG より長い定義が勝つかを見る。

宇宙 = T1 で深度条件を満たす参照アデニン（鎖指定）。陽性 = level>=0.5 で呼ばれたもの。
候補 = 呼ばれた部位の窓から観測された制約集合（連続 k-mer と、間隙を許す位置集合）。
精度 = その定義に合致する適格アデニンのうち呼ばれた割合、再現率 = 呼ばれた部位のうち合致する割合。
"""
import numpy as np, pandas as pd, re, sys, pathlib
from collections import Counter
sys.path.insert(0,'supp_261014'); from _gccggc import genome
HERE=pathlib.Path("recompute_261015"); C=HERE/"pileup_cache_261015"
g=genome(); L=len(g); comp=str.maketrans("ACGT","TGCA")
W=8
fr=[]
for lib in ("1-1","1-2","1-3"):
    z=np.load(C/f"{lib}.npz"); m=(z["code"]==0)&z["refok"]
    fr.append(pd.DataFrame({"start":z["start"][m],"plus":z["plus"][m],"Nvalid":z["Nvalid"][m],
        "Nmod":z["Nmod"][m],"Nfail":z["Nfail"][m]}).set_index(["start","plus"]))
uni=fr[0].index
for f in fr[1:]: uni=uni.union(f.index)
s=sum(f.reindex(uni).fillna(0) for f in fr)
s["dep"]=s.Nvalid+s.Nfail; s["lvl"]=s.Nmod/s.Nvalid.replace(0,np.nan)
e=s[(s.dep>=30)&(s.Nvalid>0)].reset_index()
e["called"]=(e.lvl>=0.5).values
pos=e.start.to_numpy(); plus=e.plus.to_numpy(); called=e.called.to_numpy()
print("適格アデニン",len(e)," 呼ばれた",int(called.sum()),flush=True)
# 窓行列（鎖方向、±W）
gb=np.frombuffer(g.encode(),np.uint8)
cmpl=np.zeros(256,np.uint8)
for a,b in zip(b"ACGT",b"TGCA"): cmpl[a]=b
off=np.arange(-W,W+1)
idx=np.where(plus[:,None], pos[:,None]+off[None,:], pos[:,None]-off[None,:])
ok=(idx>=0)&(idx<L); idx=np.clip(idx,0,L-1)
Wm=gb[idx]
Wm=np.where(plus[:,None], Wm, cmpl[Wm])
Wm[~ok]=ord("N")
del idx
# 候補生成: 呼ばれた部位の窓から、連続 k-mer（offset 0 を含む）を数える
MINSUP=20
cand=Counter()
Wc=Wm[called]
for k in range(3,13):
    for st in range(-W, W-k+2):
        if not (st<=0<st+k): continue
        cols=slice(st+W, st+W+k)
        sub=Wc[:,cols]
        if (sub==ord("N")).any(): continue
        keys=[bytes(r).decode() for r in sub]
        for key,n in Counter(keys).items():
            if n>=MINSUP: cand[(st,key)]+=n
# 明示的に評価する定義（間隙つき）
EXPLICIT={
 "AAGCCCG":          [(0,"A"),(1,"A"),(2,"G"),(3,"C"),(4,"C"),(5,"C"),(6,"G")],
 "HAAGCCCG":         [(-1,"ACT"),(0,"A"),(1,"A"),(2,"G"),(3,"C"),(4,"C"),(5,"C"),(6,"G")],
 "G-N-HAAGCCCG":     [(-3,"G"),(-1,"ACT"),(0,"A"),(1,"A"),(2,"G"),(3,"C"),(4,"C"),(5,"C"),(6,"G")],
 "C-N-G-N-HAAGCCCG": [(-5,"C"),(-3,"G"),(-1,"ACT"),(0,"A"),(1,"A"),(2,"G"),(3,"C"),(4,"C"),(5,"C"),(6,"G")],
}
def evaluate(constraints):
    m=np.ones(len(e),bool)
    for o,bases in constraints:
        col=Wm[:,o+W]
        mm=np.zeros(len(e),bool)
        for b in bases: mm |= (col==ord(b))
        m &= mm
    tp=int((m&called).sum()); fp=int((m&~called).sum()); fn=int((~m&called).sum())
    P=tp/(tp+fp) if tp+fp else 0.0; R=tp/(tp+fn) if tp+fn else 0.0
    F=2*P*R/(P+R) if P+R else 0.0
    return dict(n_match=tp+fp, tp=tp, precision_pct=round(100*P,2), recall_pct=round(100*R,2), f1=round(F,4))
rows=[]
for name,c in EXPLICIT.items():
    r=evaluate(c); r["definition"]=name; r["kind"]="explicit"; rows.append(r)
for (st,key),n in cand.items():
    if n<MINSUP: continue
    c=[(st+i, key[i]) for i in range(len(key))]
    r=evaluate(c); r["definition"]=("."*(st+W) if False else f"[{st}]{key}"); r["kind"]="contiguous"; rows.append(r)
t=pd.DataFrame(rows).sort_values("f1",ascending=False)
t.to_csv(HERE/"motif_definition_f1_261016.tsv",sep="\t",index=False)
print("\n--- 明示的定義"); print(t[t.kind=="explicit"][["definition","n_match","tp","precision_pct","recall_pct","f1"]].to_string(index=False))
print("\n--- F1 上位 12（連続 k-mer を含む全候補）")
print(t.head(12)[["definition","kind","n_match","tp","precision_pct","recall_pct","f1"]].to_string(index=False))
