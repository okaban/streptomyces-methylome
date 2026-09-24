> **SUPERSEDED 2026-09-23 — do not ship.** Built at the corrected C₄ offset but before the read-orientation fix (MM/ML offsets read in original-read rather than aligned coordinates), so its odds ratios are inflated. The shipped table is `SuppTable9_comod_full_denominator_T1_C4_orientation_fixed.md`.

**Supplementary Table 9 | Per-read probability threshold-sensitivity sweep for the AAGCCCG 4mC+6mA co-modification (full read denominator, T1).** Unit: (motif instance × read) pairs in which the read covers both an A (A0/A1) and a C (C3/C5) position of the same AAGCCCG instance (N = 200,312 pairs; 1,334 instances; T1 replicates 1-1, 1-2, 1-3). A pair is called 6mA-positive (A) / 4mC-positive (C) if the maximum per-read modification probability over its A / C positions is ≥ threshold. C position = motif offset 4 (AAGCC*CG), the cytosine carrying the 4mC call in the modkit pileup.  OR = (n_both × n_neither)/(n_A_only × n_C_only); 95% CI by the log-OR normal approximation; *p* by two-sided Fisher exact test. Source: 11_epigenome_integration/analysis/79_comod_full_denominator/comod_full_denominator_T1_C4.tsv (script 11_epigenome_integration/analysis/79_comod_full_denominator/comod_full_denominator_C4.py; table built by build_suppTable9.py).

| Threshold | n_both | n_A_only (6mA only) | n_C_only (4mC only) | n_neither | OR | 95% CI | *p* (Fisher) | P(6mA \| 4mC) | P(6mA) |
|---|---|---|---|---|---|---|---|---|---|
| ≥ 0.50 | 25,125 | 22,947 | 28,094 | 124,146 | 4.84 | [4.73, 4.95] | < 1e-300 | 0.472 | 0.240 |
| ≥ 0.55 | 21,317 | 21,213 | 29,253 | 128,529 | 4.42 | [4.32, 4.52] | < 1e-300 | 0.422 | 0.212 |
| ≥ 0.60 | 18,457 | 20,162 | 29,625 | 132,068 | 4.08 | [3.99, 4.18] | < 1e-300 | 0.384 | 0.193 |
| ≥ 0.65 | 15,871 | 19,253 | 29,425 | 135,763 | 3.80 | [3.71, 3.90] | < 1e-300 | 0.350 | 0.175 |
| ≥ 0.70 | 12,769 | 17,635 | 29,451 | 140,457 | 3.45 | [3.36, 3.54] | < 1e-300 | 0.302 | 0.152 |
| ≥ 0.75 | 10,108 | 16,052 | 28,856 | 145,296 | 3.17 | [3.08, 3.26] | < 1e-300 | 0.259 | 0.131 |
| ≥ 0.80 | 7,829 | 14,390 | 27,741 | 150,352 | 2.95 | [2.86, 3.04] | < 1e-300 | 0.220 | 0.111 |
| ≥ 0.85 | 5,239 | 11,874 | 25,986 | 157,213 | 2.67 | [2.58, 2.76] | < 1e-300 | 0.168 | 0.085 |
| ≥ 0.90 | 2,927 | 8,913 | 22,691 | 165,781 | 2.40 | [2.30, 2.51] | 1.2 × 10^-291 | 0.114 | 0.059 |
