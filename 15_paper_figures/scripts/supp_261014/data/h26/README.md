Single-molecule concordance tables for Supplementary Figure 15, regenerated
2026-10-14 from the nine modification-tagged alignments.

Chain:
  1. extract_perread_calls.py <lib>   (../../extract_perread_calls.py)
     pysam `modified_bases`, argmax against the canonical base
     (canonical score = 255 - sum of modified scores), emitted for every call
     whose aligned reference position is in target_positions.bed:
       tgg = 2,093 strand-oriented TGGCCGGC modified cytosines
       a01 = 2,668 AAGCCCG A0/A1 adenines
     -> pr2b/{tgg,a01}_<lib>.tsv  (read, ref_pos, read_strand, prob, state)
  2. compute_h26_stats.py -> H26.csv, PC26.csv, D26.csv

The published values are reproduced exactly: bistable prediction 0.470,
within-instance excess 0.171, between-molecule excess 0.016, coherence 3.0 %
(TGGCCGGC 4mC) and 3.3 % (AAGCCCG 6mA).

USE `modified_bases`, NOT `modified_bases_forward`. In this pysam build
`modified_bases` returns query positions in the stored-SEQ orientation, which
is the coordinate system `get_aligned_pairs` uses, while
`modified_bases_forward` returns them in the original read orientation and
therefore mismatches on reverse-aligned reads. Checked on library 1-2
(first 4,000 reads, high-confidence 4mC calls, score >= 204): under
`modified_bases` 99.44 % of calls land on a reference C or G, under
`modified_bases_forward` only 85.27 %. Swapping to the _forward accessor also
moves every published value (0.452 / 0.170 / 0.015; 2.6 % / 3.2 %), so an
accidental swap is visible in the figure.
