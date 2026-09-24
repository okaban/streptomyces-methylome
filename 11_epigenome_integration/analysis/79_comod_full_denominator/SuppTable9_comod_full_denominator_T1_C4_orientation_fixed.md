**Supplementary Table 9 | Per-read probability threshold-sensitivity sweep for the AAGCCCG 4mC+6mA co-modification (full read denominator, T1).** Unit: (motif instance × read) pairs in which the read covers both an A (A0/A1) and a C (C4) position of the same AAGCCCG instance (N = 200,312 pairs; 1,334 instances; T1 replicates 1-1, 1-2, 1-3). A pair is called 6mA-positive (A) / 4mC-positive (C) if the maximum per-read modification probability over its A / C positions is ≥ threshold. C position = motif offset 4 (AAGCC*CG), the cytosine carrying the 4mC call in the modkit pileup.  OR = (n_both × n_neither)/(n_A_only × n_C_only); 95% CI by the log-OR normal approximation; *p* by two-sided Fisher exact test. MM/ML offsets are read in ALIGNED coordinates, so reverse-strand reads are handled correctly. Source: 11_epigenome_integration/analysis/79_comod_full_denominator/comod_full_denominator_T1_C4_orientation_fixed.tsv (script 11_epigenome_integration/analysis/91_comod_stratified_and_hemi/comod_full_denominator_C4_orientation_fixed.py; table built by build_suppTable9.py). 

| Threshold | n_both | n_A_only (6mA only) | n_C_only (4mC only) | n_neither | OR | 95% CI | *p* (Fisher) | P(6mA \| 4mC) | P(6mA) |
|---|---|---|---|---|---|---|---|---|---|
| ≥ 0.50 | 51,238 | 42,182 | 57,492 | 49,400 | 1.04 | [1.03, 1.06] | 2.0 × 10^-6 | 0.471 | 0.466 |
| ≥ 0.55 | 43,452 | 39,848 | 59,830 | 57,182 | 1.04 | [1.02, 1.06] | 5.3 × 10^-6 | 0.421 | 0.416 |
| ≥ 0.60 | 37,565 | 38,236 | 60,619 | 63,892 | 1.04 | [1.02, 1.05] | 1.5 × 10^-4 | 0.383 | 0.378 |
| ≥ 0.65 | 32,298 | 36,806 | 60,131 | 71,077 | 1.04 | [1.02, 1.06] | 1.0 × 10^-4 | 0.349 | 0.345 |
| ≥ 0.70 | 26,091 | 34,017 | 60,066 | 80,138 | 1.02 | [1.00, 1.04] | 1.9 × 10^-2 | 0.303 | 0.300 |
| ≥ 0.75 | 20,596 | 31,173 | 58,883 | 89,660 | 1.01 | [0.99, 1.03] | 5.7 × 10^-1 | 0.259 | 0.258 |
| ≥ 0.80 | 15,937 | 28,131 | 56,576 | 99,668 | 1.00 | [0.98, 1.02] | 8.6 × 10^-1 | 0.220 | 0.220 |
| ≥ 0.85 | 10,672 | 23,225 | 52,945 | 113,470 | 0.98 | [0.96, 1.01] | 2.3 × 10^-1 | 0.168 | 0.169 |
| ≥ 0.90 | 5,932 | 17,479 | 46,256 | 130,645 | 0.96 | [0.93, 0.99] | 8.0 × 10^-3 | 0.114 | 0.117 |
