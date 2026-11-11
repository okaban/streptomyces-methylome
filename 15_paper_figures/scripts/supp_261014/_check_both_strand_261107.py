"""2026-11-07 点検: Supplementary Figure 3b の右上象限（両鎖とも >= 50%）の点が
何件で、そのすべてが回文 TGGCCGGCCA かを、図と同じ入力・同じ規則で数え直す。
結果は _check_both_strand_261107.txt に書く。"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _gccggc import eligible, genome, level, motif_positions, pooled  # noqa: E402

g = genome()
duplex, _ = motif_positions(g)
acc = pooled("T1")
both, elig = [], 0
for p, m in duplex:
    if eligible(acc, (p, "+")) and eligible(acc, (m, "-")):
        elig += 1
        if level(acc, (p, "+")) >= 0.5 and level(acc, (m, "-")) >= 0.5:
            both.append((p, m))

ctx = {}
for p, _m in both:
    s = p - 2                      # GCCGGC hexamer start
    dec = g[s - 2:s + 8]           # the 10-mer TGGCCGGCCA frame
    ctx[dec] = ctx.get(dec, 0) + 1

out = Path(__file__).with_suffix(".txt")
with out.open("w") as fh:
    fh.write(f"eligible duplex hexamers (T1): {elig}\n")
    fh.write(f"both strands >= 50%: {len(both)}\n")
    fh.write(f"TGGCCGGCCA: {ctx.get('TGGCCGGCCA', 0)}\n")
    for k, v in sorted(ctx.items(), key=lambda kv: -kv[1]):
        fh.write(f"  {k}\t{v}\n")
print(out)
