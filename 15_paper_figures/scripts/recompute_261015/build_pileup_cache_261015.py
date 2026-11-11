"""9 ライブラリのパイルアップを 1 回だけ npz 化し、参照塩基一致フラグを付ける。
以後の再計算はこのキャッシュから読む。列: start, plus, code(0=a,1=4mC,2=5mC), Nvalid, Nmod, Nfail, refok
ライブラリ単位で並列（既定 7 プロセス、1 プロセスあたり実測 1.5 GB 程度）。"""
import numpy as np, pandas as pd, pathlib, sys, time, os
from multiprocessing import Pool
sys.path.insert(0,'supp_261014'); from _gccggc import genome
PD=pathlib.Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015")
OUT=pathlib.Path(__file__).resolve().parent/"pileup_cache_261015"
OUT.mkdir(exist_ok=True)
LIBS=["1-1","1-2","1-3","2-1","2-3","2-4","3-2","3-3","3-4"]
COLS=["start","code","strand","Nvalid","Nmod","Nfail"]
CMAP={"a":0,"21839":1,"m":2}
EXP={0:(ord("A"),ord("T")),1:(ord("C"),ord("G")),2:(ord("C"),ord("G"))}

def one(lib):
    t0=time.time(); f=PD/f"{lib}_pileup.bed"
    if not f.exists(): return f"MISSING {lib}"
    g=np.frombuffer(genome().encode(),dtype=np.uint8)
    d=pd.read_csv(f,sep="\t",header=None,usecols=[1,3,5,9,11,15],names=COLS,
                  dtype={"start":np.int32,"Nvalid":np.int32,"Nmod":np.int32,"Nfail":np.int32})
    code=d.code.map(CMAP).fillna(255).to_numpy(np.uint8)
    plus=(d.strand.to_numpy()=="+")
    base=g[d.start.to_numpy()]
    refok=np.zeros(len(d),bool)
    for c,(pb,mb) in EXP.items():
        m=(code==c); refok |= m & np.where(plus, base==pb, base==mb)
    np.savez_compressed(OUT/f"{lib}.npz", start=d.start.to_numpy(), plus=plus, code=code,
                        Nvalid=d.Nvalid.to_numpy(), Nmod=d.Nmod.to_numpy(),
                        Nfail=d.Nfail.to_numpy(), refok=refok)
    return f"{lib} rows={len(d)} refok={100*refok.mean():.1f}% {time.time()-t0:.0f}s"

if __name__=="__main__":
    n=int(os.environ.get("NPROC","7")); t0=time.time()
    with Pool(n) as p:
        for r in p.imap_unordered(one,LIBS): print(r,flush=True)
    print("cache done in",round(time.time()-t0,1),"s")
