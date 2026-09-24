> **SUPERSEDED 2026-09-23 — do not ship.** Built before the C₄ offset correction (this table places 4mC at C₃/C₅, positions that carry ~0 modification). The shipped table is `SuppTable9_comod_full_denominator_T1_C4_orientation_fixed.md`.

**Supplementary Table 9 | Per-read probability threshold-sensitivity sweep for the AAGCCCG 4mC+6mA co-modification (full read denominator, T1).** Unit: (motif instance × read) pairs in which the read covers both an A (A0/A1) and a C (C3/C5) position of the same AAGCCCG instance (N = 200,830 pairs; 1,334 instances; T1 replicates 1-1, 1-2, 1-3). A pair is called 6mA-positive (A) / 4mC-positive (C) if the maximum per-read modification probability over its A / C positions is ≥ threshold. C positions = motif offsets 3 and 5 (NOTE: the modkit pileup places the 4mC call at offset 4; see README).  OR = (n_both × n_neither)/(n_A_only × n_C_only); 95% CI by the log-OR normal approximation; *p* by two-sided Fisher exact test. Source: 11_epigenome_integration/analysis/79_comod_full_denominator/comod_full_denominator_T1.tsv (script comod_full_denominator.py; table built by build_suppTable9.py).

| Threshold | n_both | n_A_only (6mA only) | n_C_only (4mC only) | n_neither | OR | 95% CI | *p* (Fisher) | P(6mA \| 4mC) | P(6mA) |
|---|---|---|---|---|---|---|---|---|---|
| ≥ 0.50 | 1,408 | 46,866 | 1,401 | 151,155 | 3.24 | [3.01, 3.49] | 5.5 × 10^-199 | 0.501 | 0.240 |
| ≥ 0.55 | 1,060 | 41,659 | 1,277 | 156,834 | 3.12 | [2.88, 3.39] | 2.3 × 10^-150 | 0.454 | 0.213 |
| ≥ 0.60 | 803 | 37,989 | 1,170 | 160,868 | 2.91 | [2.65, 3.18] | 9.4 × 10^-107 | 0.407 | 0.193 |
| ≥ 0.65 | 585 | 34,689 | 1,056 | 164,500 | 2.63 | [2.37, 2.91] | 4.6 × 10^-69 | 0.356 | 0.176 |
| ≥ 0.70 | 410 | 30,126 | 940 | 169,354 | 2.45 | [2.18, 2.76] | 3.7 × 10^-45 | 0.304 | 0.152 |
| ≥ 0.75 | 298 | 25,977 | 787 | 173,768 | 2.53 | [2.22, 2.90] | 2.6 × 10^-36 | 0.275 | 0.131 |
| ≥ 0.80 | 193 | 22,125 | 646 | 177,866 | 2.40 | [2.04, 2.82] | 1.3 × 10^-22 | 0.230 | 0.111 |
| ≥ 0.85 | 107 | 17,076 | 528 | 183,119 | 2.17 | [1.76, 2.68] | 1.9 × 10^-11 | 0.169 | 0.086 |
| ≥ 0.90 | 40 | 11,849 | 384 | 188,557 | 1.66 | [1.20, 2.30] | 3.8 × 10^-3 | 0.094 | 0.059 |
