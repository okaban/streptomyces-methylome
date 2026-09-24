#!/usr/bin/env python3
"""Build the replacement Supplementary Table 9 (per-read AAGCCCG 4mC+6mA
co-modification, FULL read denominator, T1) from comod_full_denominator_T1.tsv.

Adds a 95% CI for the odds ratio by the log-OR normal approximation
(Woolf): SE = sqrt(1/a + 1/b + 1/c + 1/d); CI = exp(log OR +/- 1.96 SE).
a = n_both, b = n_A_only, c = n_C_only, d = n_neither.
P(6mA|4mC) = a/(a+c); P(6mA) = (a+b)/N.  Output: SuppTable9_comod_full_denominator_T1.{tsv,md}
"""
import math, sys, pandas as pd
from pathlib import Path
H = Path(__file__).resolve().parent
# usage: build_suppTable9.py [suffix]   suffix "" -> original (C offsets 3,5); "_C4" -> corrected offset-4 scan
SUF = sys.argv[1] if len(sys.argv) > 1 else ""
t = pd.read_csv(H / f"comod_full_denominator_T1{SUF}.tsv", sep="\t")
# The Source line is DERIVED from the input actually read, not hardcoded: the
# earlier hardcoded line survived every variant build, so the C4 and the
# orientation-fixed tables both claimed to come from the uncorrected scan.
SCRIPTS = {"": "11_epigenome_integration/analysis/79_comod_full_denominator/comod_full_denominator.py",
           "_C4": "11_epigenome_integration/analysis/79_comod_full_denominator/comod_full_denominator_C4.py",
           "_C4_orientation_fixed": "11_epigenome_integration/analysis/91_comod_stratified_and_hemi/"
                                    "comod_full_denominator_C4_orientation_fixed.py"}
assert SUF in SCRIPTS, f"unknown suffix {SUF!r}: add its producing script to SCRIPTS"
SRC_LINE = (f"Source: 11_epigenome_integration/analysis/79_comod_full_denominator/"
            f"comod_full_denominator_T1{SUF}.tsv (script {SCRIPTS[SUF]}; "
            f"table built by build_suppTable9.py).")
rows = []
for r in t.itertuples():
    a, b, c, d = r.n_both, r.n_A_only, r.n_C_only, r.n_neither
    OR = (a * d) / (b * c)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    lo, hi = math.exp(math.log(OR) - 1.96*se), math.exp(math.log(OR) + 1.96*se)
    rows.append(dict(threshold=r.threshold, n_both=a, n_A_only=b, n_C_only=c, n_neither=d, N=a+b+c+d,
                     OR=OR, CI95_lower=lo, CI95_upper=hi, p_fisher=r.p_value,
                     P_6mA_given_4mC=a/(a+c), P_6mA=(a+b)/(a+b+c+d), fold_vs_independence=r.fold_enrichment))
out = pd.DataFrame(rows)
assert abs(out.OR - t.OR).max() < 1e-3, "OR mismatch vs source"
out.to_csv(H / f"SuppTable9_comod_full_denominator_T1{SUF}.tsv", sep="\t", index=False, float_format="%.6g")
def fp(p):
    if p == 0: return "< 1e-300"
    m, e = f"{p:.1e}".split("e"); return f"{m} × 10^{int(e)}"
md = ["| Threshold | n_both | n_A_only (6mA only) | n_C_only (4mC only) | n_neither | OR | 95% CI | *p* (Fisher) | P(6mA \\| 4mC) | P(6mA) |",
      "|---|---|---|---|---|---|---|---|---|---|"]
for r in out.itertuples():
    md.append(f"| ≥ {r.threshold:.2f} | {r.n_both:,} | {r.n_A_only:,} | {r.n_C_only:,} | {r.n_neither:,} | {r.OR:.2f} | [{r.CI95_lower:.2f}, {r.CI95_upper:.2f}] | {fp(r.p_fisher)} | {r.P_6mA_given_4mC:.3f} | {r.P_6mA:.3f} |")
hdr = ("**Supplementary Table 9 | Per-read probability threshold-sensitivity sweep for the AAGCCCG 4mC+6mA co-modification (full read denominator, T1).** "
       f"Unit: (motif instance × read) pairs in which the read covers both an A (A0/A1) and a C ({'C4' if SUF.startswith('_C4') else 'C3/C5'}) position of the same AAGCCCG instance "
       f"(N = {int(out.N.iloc[0]):,} pairs; 1,334 instances; T1 replicates 1-1, 1-2, 1-3). A pair is called 6mA-positive (A) / 4mC-positive (C) if the maximum "
       "per-read modification probability over its A / C positions is ≥ threshold. " + ("C position = motif offset 4 (AAGCC*CG), the cytosine carrying the 4mC call in the modkit pileup. " if SUF.startswith("_C4") else "C positions = motif offsets 3 and 5 (NOTE: the modkit pileup places the 4mC call at offset 4; see README). ") +
       " OR = (n_both × n_neither)/(n_A_only × n_C_only); 95% CI by the "
       "log-OR normal approximation; *p* by two-sided Fisher exact test. "
       + ("MM/ML offsets are read in ALIGNED coordinates, so reverse-strand reads are handled correctly. " if SUF.endswith("_orientation_fixed") else "")
       + SRC_LINE + " "
       "")
(H / f"SuppTable9_comod_full_denominator_T1{SUF}.md").write_text(hdr + "\n\n" + "\n".join(md) + "\n")
print(open(H / f"SuppTable9_comod_full_denominator_T1{SUF}.md").read())
