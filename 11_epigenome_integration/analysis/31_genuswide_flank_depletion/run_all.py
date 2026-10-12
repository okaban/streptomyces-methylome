import os, sys, json, time, traceback
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scan_genome import scan

ROOT = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/ncbi_genomes_all_species/extracted/ncbi_dataset/data"
M145 = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna"

jobs = []
for d in sorted(os.listdir(ROOT)):
    sub = os.path.join(ROOT, d)
    if not os.path.isdir(sub):
        continue
    fns = [f for f in os.listdir(sub) if f.endswith("_genomic.fna") or f.endswith("_genomic.fna.gz")]
    if not fns:
        jobs.append((None, d)); continue
    jobs.append((os.path.join(sub, sorted(fns)[0]), d))
jobs.append((M145, "M145_NC_003888.3"))

def work(job):
    path, acc = job
    if path is None:
        return {"accession": acc, "error": "no genomic fasta in assembly directory"}
    try:
        return scan(path, acc)
    except Exception as e:
        return {"accession": acc, "error": f"{type(e).__name__}: {e}"}

if __name__ == "__main__":
    t0 = time.time()
    with Pool(11) as pool:
        res = []
        for i, r in enumerate(pool.imap_unordered(work, jobs, chunksize=4), 1):
            res.append(r)
            if i % 100 == 0:
                print(i, round(time.time() - t0, 1), flush=True)
    with open("raw_counts.json", "w") as fh:
        json.dump(res, fh)
    errs = [r for r in res if "error" in r]
    print("total", len(res), "errors", len(errs), "elapsed_s", round(time.time() - t0, 1))
    for e in errs[:20]:
        print(e)
