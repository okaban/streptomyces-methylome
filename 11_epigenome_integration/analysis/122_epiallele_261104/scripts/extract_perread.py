"""Per-read modification calls at the three motif classes.
Convention matched to 15_paper_figures/scripts/supp_261014/extract_perread_calls.py
(pysam `modified_bases`, argmax against canonical = 255 - sum(mod scores)).
"""
import pysam, sys, os
from collections import defaultdict
lib = sys.argv[1]; outdir = sys.argv[2]
targ = defaultdict(set)
for ln in open(os.path.join(outdir, 'targets.bed')):
    f = ln.split()
    targ[int(f[1])].add(f[3])
WANT = {'TGG_4mC': 21839, 'AAG_4mC': 21839, 'AAG_6mA_A0': 'a', 'AAG_6mA_A1': 'a'}
bam = pysam.AlignmentFile(f'/Users/okaban/bioinfo/methyl/260102_M145/analysis/{lib}_mapped.bam','rb')
out = {s: open(f'{outdir}/pr/{s}_{lib}.tsv','w') for s in WANT}
nr = 0
for r in bam.fetch('NC_003888.3'):
    try: mb = r.modified_bases
    except Exception: continue
    if not mb: continue
    nr += 1
    strand = '-' if r.is_reverse else '+'
    d = dict(r.get_aligned_pairs(matches_only=True))
    byq = defaultdict(dict)
    for (base, mstrand, code), vals in mb.items():
        for qpos, score in vals: byq[qpos][code] = score
    for qpos, codes in byq.items():
        rp = d.get(qpos)
        if rp is None: continue
        for s in targ.get(rp, ()):
            want = WANT[s]
            if want not in codes: continue
            tot = sum(codes.values()); canon = 255 - tot; p = codes[want]
            state = int(p > canon and all(p >= v for k, v in codes.items() if k != want))
            out[s].write(f"{r.query_name}\t{rp}\t{strand}\t{p/255:.4f}\t{state}\n")
for f in out.values(): f.close()
print(lib, 'reads_with_mods', nr, flush=True)
