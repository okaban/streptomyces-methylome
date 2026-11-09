"""124_occupancy_tracks_261104 — C2-1/C2-2/C2-3, I2-2/I2-3, E1-3/E1-4.
Reconstruction of the analysis run in the epi-rna kernel on 2026-11-04.
Inputs
  pileups : /Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015/*_nofilt.bed
  genome  : 11_epigenome_integration/data/NC_003888.3.fna
  TSS     : 11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv
  tracks  : raw/tracks/*, raw/tracks_hupa/*  (downloaded from GEO; see README.md)
Design
  unit of analysis = PEAK (never the site).
  null = random genomic region of IDENTICAL width with matched region GC
         (tolerance ladder 0.02/0.035/0.05/0.08/0.12/0.20), 400 replicates.
         level test additionally requires the control to contain >=1 measurable
         site of the same motif; density test does not.
  tests = (a) empirical p of mean-of-peak-means against the null distribution
          (b) per-peak percentile u of the observed value within its own 400
              controls, one-sample Wilcoxon on u-0.5  (calibrated, higher power)
          (c) BH across all tracks x motifs
  power = mde_* columns give the 2.5-97.5 percentile of the null as a ratio,
          i.e. the effect that would have been rejected at two-sided 5%.
Peak RSS of the full run: ~0.66 GB (resource.getrusage RUSAGE_SELF).

# Key numbers reproduced by this script:
#   methylatable positions: TGGCCGGC C5 = 2093, AAGCCCG C5/A1 = 1334 each
#   AL645882.2 md5 == NC_003888.3 md5 == fa44769ad49da012d251d8b1ec8672fc (8,667,507 bp)
#   22 occupancy tracks; 54 density tests, 29 level tests; 0 survive BH 5%
