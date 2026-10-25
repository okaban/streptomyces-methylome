import pysam, sys
from collections import defaultdict
lib = sys.argv[1]
targ = defaultdict(set)   # pos -> {setname}
for ln in open('both.bed'):
    f = ln.split()
    targ[int(f[1])].add(f[3])
bam = pysam.AlignmentFile(f'/Users/okaban/bioinfo/methyl/260102_M145/analysis/{lib}_mapped.bam', 'rb')
WANT = {'tgg': 21839, 'a01': 'a'}
out = {s: open(f'pr2b/{s}_{lib}.tsv', 'w') for s in WANT}
for r in bam.fetch('NC_003888.3'):
    try:
        mb = r.modified_bases
    except Exception:
        continue
    if not mb:
        continue
    strand = '-' if r.is_reverse else '+'
    d = dict(r.get_aligned_pairs(matches_only=True))
    byq = defaultdict(dict)
    for (base, mstrand, code), vals in mb.items():
        for qpos, score in vals:
            byq[qpos][code] = score
    for qpos, codes in byq.items():
        rp = d.get(qpos)
        if rp is None:
            continue
        for s in targ.get(rp, ()):
            want = WANT[s]
            if want not in codes:
                continue
            tot = sum(codes.values())
            canon = 255 - tot
            p = codes[want]
            state = int(p > canon and all(p >= v for k, v in codes.items() if k != want))
            out[s].write(f"{r.query_name}\t{rp}\t{strand}\t{p/255:.4f}\t{state}\n")
for f in out.values():
    f.close()
