# Phase variation in M145: repeat tracts in MTase/BREX genes and promoter modification clusters

Analysis ID `123_phasevar_261104` (ledger steps E2-1, E2-2, E2-3, A5-1).
Date 2026-11-04. All numbers below are from
`11_epigenome_integration/analysis/123_phasevar_261104/tables/`.

## Question

Phasevarions in *Haemophilus*, *Neisseria* and *Helicobacter* arise when a simple
sequence repeat inside an MTase gene changes length by slipped-strand mispairing,
switching the MTase ON or OFF and remodelling the methylome. This study has not
examined that possibility. *S. coelicolor* is directly relevant: Sumby & Smith
(2003, *J Bacteriol* 185:4558-4563) attributed Pgl (phage growth limitation) phase
variation at least partly to expansion/contraction of a polyguanine tract in
*pglX*, the adenine MTase of the BREX-2 locus that M145 carries at 7.35-7.37 Mb.

## E2-1 Repeat inventory (reference sequence)

Scope: the 129-gene M145 MTase catalogue plus the 15 genes of the BREX-2 locus
(union 143 genes). For each gene, the coding sequence and its promoter were
scanned (promoter = -300..TSS using the strongest Jeong 2016 TSS where one exists,
otherwise 300 bp upstream of the CDS). Homopolymer runs >= 6 bp and perfect
tandem repeats of unit 2-6 bp with >= 3 copies and >= 9 bp total were recorded.

- **316 repeats in 117 of the 143 genes**: 62 homopolymer runs and 254 SSRs.
- Homopolymer lengths are short: 56 runs of 6 bp, 5 of 7 bp, **1 of 8 bp**. No
  run reaches 9 bp in any MTase or BREX gene.
- SSRs never exceed 5 copies (longest total 18 bp, a trimer/hexamer); 47 gene-body
  repeats have a unit length not divisible by 3 and so could in principle frameshift.
- The single longest homopolymer in the whole target set is a **G8 tract at
  7,353,983-7,353,990, inside *pglX* (SCO6627 / SC_RS35335)**, 776 bp into the
  3,633 bp CDS (codon 259 of 1210). This is the expected location of the
  polyguanine tract reported for the Pgl system.
- A second PglX paralogue, SCO5331 (SC_RS28835, BREX-2 adenine MTase at 5.80 Mb),
  carries a C7 tract in its CDS and a G7 tract in its promoter.

Genome context: at 72.1% GC, long homopolymers are rare genome-wide — only 104
runs >= 7 bp, 36 >= 8 bp and 2 >= 10 bp in 8,667,507 bp, longest 12 bp. The
*pglX* G8 is therefore unusual for this genome, but it is far shorter than the
poly-G tracts (typically >= 9-17 bp) that drive characterised phasevarions.

## E2-2 Read-level tract length (nanopore)

Read level was reached. All nine nanopore libraries
(`/Users/okaban/bioinfo/methyl/260102_M145/analysis/{lib}_mapped.bam`) were used.

Design. For each tract, reads were required to be primary, MAPQ >= 20, and to
align across both flanking anchor bases with >= 25 bp of aligned flank on each
side. The read bases strictly between the two anchors were extracted by a CIGAR
walk, and the tract length was taken as the longest run of the tract base within
that segment.

**Error control.** Homopolymer length is exactly where nanopore basecalling fails,
so the spread of observed lengths cannot be read as polymorphism on its own. For
every candidate, control tracts were drawn from the same genome: same base, same
reference length, >= 10 kb away from any MTase or BREX gene, hence carrying a
single allele by assumption. 143 control tracts were genotyped alongside the 5
candidates (23-24 controls each for C7/G7, 9 for G8; 40 each at length 6 and 7 at
length 9, used to build the kernels for the +-1 bp alternative alleles).

Depth after filtering, pooled over the nine libraries: **693-1,209 reads per
candidate tract**, median 713 per control tract; 112,900 per-read measurements total (4,729 candidate, 108,171 control).

## E2-3 Multimodality

**No candidate tract is multimodal. 0 of 5.**

At every candidate the modal observed length equals the reference length, as it
does at 100% of control tracts. The statistic was the per-read G-statistic of the
observed length distribution against the pooled control kernel, and the null was
built control-against-control (each control versus the kernel of the others), so
it absorbs both basecalling error and locus-to-locus sequence context.

| tract | gene | region | ref | reads | modal f | 2nd mode f | control 2nd-mode median / max | p (vs control null) |
|---|---|---|---|---|---|---|---|---|
| G8 7,353,983 | *pglX* SCO6627 | CDS | 8 | 756 | 0.455 | 0.173 (len 7) | 0.195 / 0.301 | 0.40 |
| C7 5,804,339 | *pglX* SCO5331 | CDS | 7 | 1209 | 0.873 | 0.089 (len 6) | 0.103 / 0.285 | 0.83 |
| G7 5,805,837 | *pglX* SCO5331 | promoter | 7 | 1169 | 0.944 | 0.031 (len 6) | 0.119 / 0.320 | 0.12 |
| C7 7,347,663 | SCO6624/6625 | promoter | 7 | 693 | 0.844 | 0.117 (len 6) | 0.103 / 0.285 | 0.83 |
| C7 5,169,963 | SCO4758 | CDS | 7 | 902 | 0.860 | 0.069 (len 4) | 0.103 / 0.285 | 0.54 |

The *pglX* G8 tract has the broadest distribution of the five (modal fraction
0.455), but the breadth is a smear toward **shorter** calls (0.173 at 7 bp, 0.155
at 6 bp, 0.134 at 5 bp), which is the direction and shape of nanopore homopolymer
under-calling, not a discrete second allele. The matched G8 control tracts behave
the same way (pooled kernel 0.573 at 8 bp, 0.189 at 7, 0.116 at 6), and the
candidate sits inside the spread of those nine controls (p = 0.40).

Hartigan's dip test was computed but is **uninformative here**: it rejects
unimodality at p < 0.05 for 47.6% of the single-allele control tracts, because the
length variable is integer-valued and error-broadened. Dip p-values for the
candidates are reported in the tables but were not used for the call.

**Detection bound (what this result can and cannot exclude).** By simulating
two-allele mixtures against the control-to-control null at each tract's own depth,
with 80% power at alpha = 0.05:

| tract | minor allele at L+1 detectable from | at L-1 detectable from |
|---|---|---|
| G8 *pglX* SCO6627 | 29% of reads | 45% |
| C7/G7 tracts (4 loci) | 17-18% | 24-29% |

So a slipped-strand allele present in **more than roughly a fifth to a third of the
population** is excluded at all five tracts. A **minority switched subpopulation
below that fraction is not excluded**, and nanopore homopolymer error is what sets
the floor — this is a limit of the assay, not evidence of homogeneity at low
frequency.

## A5-1 Promoter modification clusters (pap/Dam-like architecture)

In *E. coli*, the *pap* switch works because two Dam GATC sites sit close together
in the promoter and are differentially protected. The analogous question here:
do M145 promoters carry clustered modifiable positions?

Modifiable position = an instance of a methylated recognition sequence
(TGGCCGGC, 4mC at C5; AAGCCCG, 4mC at C5 and 6mA at A1), counted on both strands
within -300..TSS. All 3,570 experimental TSS were scanned.

- 373 modifiable sites in total, 0.105 per promoter.
- **30 promoters carry >= 2 sites**; 18 of those have two sites within 50 bp;
  the maximum is 4 sites (SCO5151 and SCO5281).
- Against the null the ledger specified — random placement of the observed total
  across promoters — 18.2 promoters would be expected, so the observed 30 looks
  like a 1.65x excess.
- Against a **GC-matched null** (10,000 random 300 bp genomic windows, binned on
  promoter GC decile) the expectation is **25.5 (95% CI 16-36), p = 0.21**.
- **No promoter survives multiple testing**: minimum GC-matched p = 4.6e-4,
  minimum q = 1.0, 0 promoters at q < 0.05.

The apparent clustering is a property of GC composition, not of promoter identity.
This repeats the pattern already documented in this project for the target-gene
function bias (D119) and must be reported with the GC-matched null, not the naive one.

## What may be written

1. M145 MTase and BREX genes contain no repeat tract of the length that drives
   characterised phasevarions. The longest is a G8 in *pglX*; phase-variable MTase
   loci elsewhere typically carry >= 9-17 bp homopolymers or long tetranucleotide
   arrays.
2. At read level across nine libraries, none of the five candidate tracts, including
   the *pglX* G8, shows a second length allele above the detection floor
   (17-45% depending on tract and direction). Within that bound the cultures are
   monomorphic at these tracts.
3. No promoter shows pap-like clustering of modifiable positions beyond GC expectation.

## What must not be written

- Do **not** write that *Streptomyces* has no phase-variable MTase. Sumby & Smith
  2002/2003 establish *pglX*-based Pgl phase variation in this very strain lineage.
  The correct statement is that **we found no length polymorphism at the *pglX*
  tract in these cultures, above a detectable minor-allele fraction of ~29-45%**.
- Do **not** present this as excluding phase variation. Pgl switching is reported
  at low frequency, and a low-frequency switched subpopulation is below our floor.
  The assay limit is nanopore homopolymer error, and that limit must be stated
  whenever the negative is stated.
- Do **not** quote the 30 clustered promoters as a finding. The GC-matched
  expectation is 25.5 and p = 0.21.

## Suggested placement

Supplementary text plus one supplementary table (E2-1 inventory and E2-3 bounds).
It closes a question a referee can reasonably raise — "could a phasevarion explain
the methylome you describe?" — with a bounded negative rather than silence. It
does not belong in the main text, because nothing here is positive.

## Files

- `tables/e2_1_repeats_all_261104.tsv` - 316 repeats, with unit, copies, position, frame
- `tables/e2_1_genomewide_homopolymers_261104.tsv` - genome-wide homopolymer background
- `tables/e2_2_loci_261104.tsv` - 5 candidate + 143 control tracts
- `tables/e2_2_perread_lengths_261104.tsv.gz` - 112,900 per-read tract measurements
- `tables/e2_3_locus_stats_261104.tsv`, `e2_3_candidate_vs_control_261104.tsv`
- `tables/e2_3_candidate_length_histograms_261104.tsv`
- `tables/e2_3b_detection_bounds_261104.tsv` - minimum detectable minor-allele fraction
- `tables/a5_1_promoter_sites_261104.tsv`, `a5_1_clustered_promoters_261104.tsv`,
  `a5_1_summary_261104.tsv`
- scripts in `11_epigenome_integration/analysis/123_phasevar_261104/scripts/`

## Declared limits

- Read-level genotyping covered the 5 homopolymer tracts >= 7 bp. The 254 SSRs
  were inventoried from the reference only; none exceeds 5 copies, below the
  usual phase-variable threshold, and short SSR length cannot be read reliably
  from nanopore data either.
- Promoters without an experimental TSS used CDS-300 for the E2-1 scan; A5-1 used
  experimental TSS only (3,570 promoters).
- Repeat definitions are perfect-repeat only; imperfect/interrupted repeats were
  not scored.
- The ledger file `Writing/lit_methylation_function_approaches_261104.tsv` was not
  present on disk, so the prescribed per-step procedure could not be read. The
  designs above follow the task description; the A5-1 null was additionally run
  GC-matched because the specified random-placement null is confounded by GC in
  this genome, and the difference is decisive (18.2 vs 25.5 expected).
