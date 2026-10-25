#!/usr/bin/env python3
"""
P1-1  k-mer cross-talk artefact test for the AAGCCCG same-strand 4mC+6mA
      dual modification.

Reviewer concern (already flagged in the manuscript): because Nanopore calls
derive from overlapping k-mers, a strong modification call at one base may bias
the call at a base 3-4 bp away. The same-strand A0/A1 (6mA) <-> C4 (4mC)
co-occurrence sits exactly inside that window, so the giant OR could be a
calling artefact rather than genuine co-deposition.

DISCRIMINATING LOGIC
  Cross-talk is a *generic* property of the basecaller k-mer model: a strong
  6mA call should inflate nearby C calls regardless of sequence context.
  Genuine co-deposition is *motif-specific*: only AAGCCCG should show it.

  We therefore measure, per read, the rate at which a high-confidence 6mA call
  has a high-confidence 4mC call at a NEAR offset (+/-3, +/-4 bp; inside the
  k-mer window) vs a FAR offset (+/-20..50 bp; outside it), separately for
  6mA calls that ARE vs ARE NOT in an AAGCCCG motif (the latter dominated by
  GATC/Dam and other contexts = the generic-cross-talk baseline).

  near/far ratio in NON-AAGCCCG 6mA  -> magnitude of any generic short-range bleed
  AAGCCCG near rate vs NON-AAGCCCG near rate -> motif specificity (real vs artefact)

Also: probability-distribution QC at AAGCCCG A/C positions (graded = bleed-like;
bimodal high-confidence = genuine state calls).

Input (read_id-sorted, one row per read x position):
  65_per_read_comod/1-1_modonly.tsv  (T1 replicate 1)   cols: read_id, ref_position, call_prob, call_code, fail
Reference: /Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa
"""
# DATA_PATH_EXEMPT: 65_per_read_comod は 2026-09-21 に superseded 宣言されているが、その理由は methylation_site_sequences.csv の sequence 窓経由でオフセットを求めると 1 bp ずれる点にある。このスクリプトは位置集合を参照ゲノムから re.finditer で直接構築しており、その欠陥を継承しない。modonly.tsv からは ref_position / call_prob / call_code の生の列しか読まない。2026-10-15 に再計算で検証: AAGCCCG near/far 365.0、非モチーフ near/far 1.0915、モチーフ/非モチーフ 435.2（出所表 429.7 との 1.3% 差は参照塩基照合の有無で説明され、非モチーフ n も同率ずれる）。
import os, re, sys
from collections import defaultdict
import numpy as np

PROB = 0.75                      # match manuscript high-confidence threshold
REF = "/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa"
MODONLY = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/"
                             "analysis/65_per_read_comod/1-1_modonly.tsv")
OUT = os.path.expanduser("~/bioinfo/rna-seq/11_epigenome_integration/"
                         "analysis/78_reviewer_robustness/tables")
CALL_6mA, CALL_4mC = "a", "21839"
MOTIF = "AAGCCCG"
A_POS, C_POS = {0, 1}, {4}   # 2026-09-20: 4mC is at C4; C3/C5 carry ~0 (0.15%/0.05%)
NEAR = (-4, -3, 3, 4)
FAR = tuple(o for o in range(-50, 51) if 20 <= abs(o) <= 50)


def load_genome(fa):
    parts = []
    with open(fa) as f:
        for ln in f:
            if ln.startswith(">"):
                if parts:
                    break
            else:
                parts.append(ln.strip().upper())
    return "".join(parts)


def motif_sets(genome):
    """Return (set of AAGCCCG 6mA genomic positions, set of 4mC positions)."""
    a_set, c_set = set(), set()
    L = len(MOTIF)
    rc = MOTIF[::-1].translate(str.maketrans("ACGT", "TGCA"))
    for m in re.finditer(f"(?={MOTIF})", genome):
        s = m.start()
        for p in A_POS: a_set.add(s + p)
        for p in C_POS: c_set.add(s + p)
    for m in re.finditer(f"(?={rc})", genome):
        s = m.start()
        for p in A_POS: a_set.add(s + (L - 1 - p))
        for p in C_POS: c_set.add(s + (L - 1 - p))
    return a_set, c_set


def main():
    genome = load_genome(REF)
    aag_A, aag_C = motif_sets(genome)
    print(f"AAGCCCG genomic positions: 6mA(A0/A1)={len(aag_A)}  4mC(C4)={len(aag_C)}")

    # tallies
    cnt = defaultdict(int)          # category counters
    prob_hist_A = np.zeros(20)      # AAGCCCG A positions, all calls
    prob_hist_C = np.zeros(20)      # AAGCCCG C positions, all calls
    bins = np.linspace(0, 1, 21)

    def flush(read_6mA, read_4mC):
        """read_6mA, read_4mC: dict pos->prob (high-conf only) for one read."""
        if not read_6mA:
            return
        c4_pos = read_4mC  # dict
        for p in read_6mA:
            is_aag = p in aag_A
            grp = "AAG" if is_aag else "NON"
            cnt[f"{grp}_n6mA"] += 1
            near_hit = sum(1 for o in NEAR if (p + o) in c4_pos)
            far_hit = sum(1 for o in FAR if (p + o) in c4_pos)
            cnt[f"{grp}_near_hits"] += near_hit
            cnt[f"{grp}_far_hits"] += far_hit
            if near_hit:
                cnt[f"{grp}_near_anyread"] += 1

    n_lines = 0
    cur = None
    r6, r4 = {}, {}                 # high-conf calls for current read
    a6_all, c4_all = {}, {}         # ALL calls (for prob QC) at AAGCCCG pos
    with open(MODONLY) as f:
        f.readline()
        for line in f:
            n_lines += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            rid, pos, prob, code, fail = parts
            if code != CALL_6mA and code != CALL_4mC:
                continue
            pos = int(pos); prob = float(prob)
            if rid != cur:
                if cur is not None:
                    flush(r6, r4)
                cur = rid; r6, r4 = {}, {}
            # probability QC at AAGCCCG positions (all calls)
            if code == CALL_6mA and pos in aag_A:
                prob_hist_A[min(int(prob * 20), 19)] += 1
            if code == CALL_4mC and pos in aag_C:
                prob_hist_C[min(int(prob * 20), 19)] += 1
            # high-confidence filter for co-occurrence
            if prob >= PROB and fail != "true":
                if code == CALL_6mA:
                    r6[pos] = prob
                else:
                    r4[pos] = prob
        flush(r6, r4)

    print(f"Streamed {n_lines:,} mod rows.")

    # ---- rates -------------------------------------------------------------
    def rates(grp):
        n = cnt[f"{grp}_n6mA"]
        if n == 0:
            return dict(n=0)
        near_rate = cnt[f"{grp}_near_hits"] / (len(NEAR) * n)   # per offset position
        far_rate = cnt[f"{grp}_far_hits"] / (len(FAR) * n)
        any_near = cnt[f"{grp}_near_anyread"] / n
        return dict(n=n, near_rate=near_rate, far_rate=far_rate,
                    near_over_far=(near_rate / far_rate if far_rate else float('inf')),
                    frac_6mA_with_near_4mC=any_near)

    aag = rates("AAG"); non = rates("NON")
    print("\n" + "=" * 68)
    print("CROSS-TALK TEST  (per-offset 4mC co-call rate around 6mA calls)")
    print("=" * 68)
    for lab, r in [("AAGCCCG 6mA", aag), ("non-AAGCCCG 6mA (generic)", non)]:
        print(f"\n[{lab}]  n_6mA(high-conf)={r['n']:,}")
        if r['n']:
            print(f"  4mC co-call rate per NEAR offset (+/-3,4 bp): {r['near_rate']:.5f}")
            print(f"  4mC co-call rate per FAR  offset (20-50 bp): {r['far_rate']:.5f}")
            print(f"  near/far ratio (generic short-range bleed index): {r['near_over_far']:.2f}")
            print(f"  fraction of 6mA with >=1 near 4mC: {r['frac_6mA_with_near_4mC']:.4f}")
    if aag.get('near_rate') and non.get('near_rate'):
        print(f"\n  >>> AAGCCCG near-rate / non-AAGCCCG near-rate (motif specificity) = "
              f"{aag['near_rate'] / non['near_rate']:.1f}x")

    # ---- prob QC -----------------------------------------------------------
    centers = (bins[:-1] + bins[1:]) / 2
    def summ(h):
        tot = h.sum()
        if tot == 0: return (0, 0, 0)
        frac_hi = h[centers >= 0.9].sum() / tot       # near 1.0
        frac_mid = h[(centers >= 0.3) & (centers < 0.7)].sum() / tot
        return tot, frac_hi, frac_mid
    tA, hiA, midA = summ(prob_hist_A)
    tC, hiC, midC = summ(prob_hist_C)
    print("\n" + "=" * 68)
    print("PROBABILITY QC at AAGCCCG positions (bleed -> graded/mid; real -> hi)")
    print("=" * 68)
    print(f"  A0/A1 (6mA): n={int(tA):,}  frac prob>=0.9={hiA:.3f}  frac mid(0.3-0.7)={midA:.3f}")
    print(f"  C4 (4mC): n={int(tC):,}  frac prob>=0.9={hiC:.3f}  frac mid(0.3-0.7)={midC:.3f}")

    # ---- save --------------------------------------------------------------
    import pandas as pd
    pd.DataFrame([
        dict(group="AAGCCCG_6mA", **aag),
        dict(group="nonAAGCCCG_6mA_generic", **non),
    ]).to_csv(os.path.join(OUT, "P1_1_crosstalk_rates_C4.tsv"), sep="\t", index=False)
    pd.DataFrame(dict(prob_bin_center=centers,
                      AAGCCCG_A_6mA=prob_hist_A,
                      AAGCCCG_C_4mC=prob_hist_C)).to_csv(
        os.path.join(OUT, "P1_1_prob_distributions_C4.tsv"), sep="\t", index=False)
    print("\nDONE -> 78_reviewer_robustness/tables/P1_1_*.tsv")


if __name__ == "__main__":
    main()
