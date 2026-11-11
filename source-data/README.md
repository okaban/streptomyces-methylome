# Source Data tables

*Streptomyces coelicolor* A3(2) M145 methylome and transcriptome.

These tables support statements in the manuscript but are not required to follow its argument.
They are published in the analysis repository (see Data Summary), under `source-data/`, one file per section together with this index.
The manuscript cites them as **Source Data Table SD1-SD35**; this index has 35 sections.

---

## SD1

Per-site TGGCCGGC and AAGCCCG methylation dynamics classification. This table is the record behind the dynamic-site claim withdrawn in Results §R4, retained for that purpose and not as evidence of methylation dynamics. For each dynamic site (maximum-minus-minimum modification frequency ≥ 10 percentage points; TGGCCGGC: *n* = 243 of the 855 positions called at all three timepoints; AAGCCCG: *n* = 97 of the 309 positions called at all three timepoints), modification frequencies at T1, T2, T3, range (Δ, percentage points), and trajectory type (Monotonic-Up, Monotonic-Down, Valley, Peak, and a single TGGCCGGC site classified Other). Rows: 1164.

## SD2

Sequencing depth and site yield per nanopore library. Median true depth and the number of TGGCCGGC 4mC sites called for each of the nine nanopore libraries, on the canonical pileups (`pileup_canonical_261015`, which include the 2026-02 top-up flow cells for libraries 3-2 and 3-4). Depth is the median of (Nvalid_cov + Nfail) over all cytosine rows; a site is called when pooled true depth is at least 30 and Nmod/Nvalid is at least 0.5, no separate floor being placed on Nvalid. Library 3-2 remains the lowest-yield library at 25× and 463 sites against 35–53× and 1,117–1,815 sites for the other eight; all nine were retained and every analysis is reported on per-timepoint pools of three libraries. Rows: 9.

## SD3

Sequencing data deposition. Accessions for the data generated in this study: nine Oxford Nanopore modification-tagged alignment sets and nine Illumina RNA-seq read sets in the NCBI Sequence Read Archive under BioProject PRJNA1533939, study accession SRP743765, BioSample accessions SAMN64069313–SAMN64069321 (three timepoints × three biological replicates) and run accessions SRR41105505–SRR41105522. Raw nanopore signal (POD5) is not accepted by the SRA; the per-read 4mC, 6mA and 5mC calls are carried in the MM and ML tags of the deposited BAM files. The `biological_replicate` attribute runs 1–3 at T1, 1, 3 and 4 at T2 and 2, 3 and 4 at T3: four cultures were grown at each timepoint and three taken forward, so the replicate numbers are not consecutive. *Source: BioProject PRJNA1533939; study SRP743765; BioSample SAMN64069313–321; runs SRR41105505–522.* Rows: 9.

## SD4

Occupancy on the filtered and the unfiltered pileup, per library. Per-library TGGCCGGC 4mC occupancy computed on the confidence-filtered pileup, on all reads, and on the filtered pileup with a depth floor of 10, with the resulting bias in percentage points. This is the basis for computing occupancy and level statistics on the unfiltered pileup throughout. *n* = nine libraries, three biological replicates at each of T1, T2 and T3. Rows: 9.

## SD5

Sensitivity of the three chemistries to the true-depth floor. Genome-wide counts of positions measurable at a pooled true depth (valid plus failed calls) of ≥ 10 and of ≥ 30 reads, and of those positions called at a modification level ≥ 0.5, for 6mA, 4mC and 5mC at each timepoint, with the mean level over the measurable set. The floor applies to true depth, not to the number of valid calls, so it does not act as a floor on modification state. Lowering it from 30 to 10 changes the called-site counts by at most five positions in any cell and recovers no 5mC site at any timepoint. Source: `15_paper_figures/scripts/recompute_261015/pileup_cache_261015`. Rows: 18.

## SD6

Distribution of modified positions across feature classes, by mark and recognition sequence. For each of 4mC at TGGCCGGC position 5, 4mC at AAGCCCG position 5, 6mA at AAGCCCG positions 1/2, and the two composition backgrounds (all callable cytosines, all callable adenines), the number and percentage of occurrences, measurable positions and modified positions falling in promoters, 5′ UTRs, CDSs and elsewhere, at T1–T3. Promoter and 5′ UTR follow the Jeong et al. TSS set; CDS and other follow the RefSeq annotation. Rows: 60 (5 series × 4 feature classes × 3 timepoints). Source data for the feature-class percentages quoted in the Figure 2 legend.

**Supplementary Table S3 | Same-strand dual modification of AAGCCCG.**

## SD7

Per-read probability threshold-sensitivity sweep for the AAGCCCG 4mC+6mA co-modification (full read denominator, position 5). For each per-read probability cutoff from ≥ 0.5 to ≥ 0.9, the 2 × 2 table of 6mA (A1/A2) against 4mC (position 5) over all 200,312 instance × read pairs at T1 in which the read covers at least one adenine position and position 5 of the same AAGCCCG instance. **This sweep predates the correction of the position-5 census and is superseded**; its odds ratios must not be read as current values. On the corrected census the same test gives OR 1.317 at T1 on the HAAGCCCG denominator (80,227 pairs) and OR 1.151 across all AAGCCCG instances (93,187 pairs); the superseded figures were OR 1.429 on 84,047 pairs and OR 1.174 on 199,276 pairs (Results §R3, Figure 3b). The table is retained because the sweep shows the result is not an artefact of the auto-estimated 0.64 adenine pass-threshold. | threshold | n_both | n_A_only | n_C_only | n_neither | OR | 95% CI | *p* | P(6mA \| 4mC) | P(6mA) | |---|---|---|---|---|---|---|---|---|---| | ≥ 0.50 | 51,238 | 42,182 | 57,492 | 49,400 | 1.04 | [1.03, 1.06] | 2.0 × 10⁻⁶ | 0.471 | 0.466 | | ≥ 0.55 | 43,452 | 39,848 | 59,830 | 57,182 | 1.04 | [1.02, 1.06] | 5.3 × 10⁻⁶ | 0.421 | 0.416 | | ≥ 0.60 | 37,565 | 38,236 | 60,619 | 63,892 | 1.04 | [1.02, 1.05] | 1.5 × 10⁻⁴ | 0.383 | 0.378 | | ≥ 0.65 | 32,298 | 36,806 | 60,131 | 71,077 | 1.04 | [1.02, 1.06] | 1.0 × 10⁻⁴ | 0.349 | 0.345 | | ≥ 0.70 | 26,091 | 34,017 | 60,066 | 80,138 | 1.02 | [1.00, 1.04] | 1.9 × 10⁻² | 0.303 | 0.300 | | ≥ 0.75 | 20,596 | 31,173 | 58,883 | 89,660 | 1.01 | [0.99, 1.03] | 5.7 × 10⁻¹ | 0.259 | 0.258 | | ≥ 0.80 | 15,937 | 28,131 | 56,576 | 99,668 | 1.00 | [0.98, 1.02] | 8.6 × 10⁻¹ | 0.220 | 0.220 | | ≥ 0.85 | 10,672 | 23,225 | 52,945 | 113,470 | 0.98 | [0.96, 1.01] | 2.3 × 10⁻¹ | 0.168 | 0.169 | | ≥ 0.90 | 5,932 | 17,479 | 46,256 | 130,645 | 0.96 | [0.93, 0.99] | 8.0 × 10⁻³ | 0.114 | 0.117 | Rows: 9.

## SD8

Occupancy and the 4mC-to-6mA ratio across the full range of the site-defining modification-level cut. The level cut that defines a called site was varied from 0.1 to 0.9 in steps of 0.1 and occupancy recomputed at each value, for each timepoint and for three motif classes (TGGCCGGC 4mC, AAGCCCG 6mA at position 1, and the same 6mA restricted to HAAGCCCG), giving 81 rows. Absolute occupancy is strongly cut-dependent (AAGCCCG 6mA at T1 moves from 95.41% at a cut of 0.1 to 24.29% at the census cut of 0.5 and 0.08% at 0.9; this sweep uses the measurable census of 1,330 instances, so 24.29% is the counterpart of the 24.21% reported on the canonical 1,334 in the main text), but 4mC occupancy exceeds 6mA occupancy in all 81 rows and the ratio spans 1.02 to 355.8 without inverting. At the census cut the TGGCCGGC-4mC-to-AAGCCCG-6mA ratio is 3.71, 3.44 and 3.54 at T1, T2 and T3 — a different pair of motifs from the 1.89–1.97 reported in the main text. *Source: Source Data Table SD8 (15_paper_figures/tables/supplementary/SuppTable29_level_cut_sweep_261021.tsv).* Rows: 81.

## SD9

Both-strand methylation at matched duplex units. One row per duplex unit for the 9 units at which both strands meet the depth requirement, giving the modified fraction on each strand and the difference between them. Of the 16,059 genomic GCCGGC hexamer duplexes, the number methylated on both strands, restricted to the occurrences that carry a TGGCCGGC recognition context on each strand. Hemimethylation of this hexamer is a consequence of where the recognition context lies, not a property of the enzyme: 14,021 duplexes carry no TGGCCGGC context on either strand and are never methylated, 1,983 carry it on one strand only and can therefore only ever be hemimethylated, and 55 carry it on both strands (the 10 bp palindrome TGGCCGGCCA). *The table that stood here (both-strand methylation by read set and recognition context) was withdrawn with the authors' agreement on 2026-10-25; both-strand methylation is now carried by Source Data Table SD9 alone. Later tables have been renumbered accordingly.* Rows: 9.

**Supplementary Table S4 | Determinants of the adenine mark.**

## SD10

Methylation and expression within each predicted cluster across the time course. Per-region mean 4mC level, mean 6mA level and mean transcript level at T1, T2 and T3 for all 29 antiSMASH regions. A region enters the correlation for a mark when at least five positions of that mark are measurable under the census rule, which admits 22 regions for 4mC and 21 for 6mA; the level is the mean of Nmod/Nvalid over those positions. Two abundance columns are not the same quantity: the table's TPM columns are the mean TPM of all genes whose midpoint lies inside the region, whereas the correlation below is computed on the core-biosynthetic-gene TPM of the region. Neither mark tracks cluster expression: 4mC against TPM gives Spearman *ρ* = −0.2842 (*p* = 0.0207, 66 region × timepoint observations) and 6mA against TPM *ρ* = +0.2161 (*p* = 0.0889, 63 observations); neither reaches the family threshold for the eighteen cluster-level tests, α = 0.0028. *Source: Source Data Table SD10 (15_paper_figures/tables/supplementary/SuppTable20a_cluster_methylation_expression_allread_261101.tsv).* Rows: 29.

## SD11

Motif occurrence inside clusters relative to their flanks. Observed-over-expected occurrence of TGGCCGGC and AAGCCCG inside each of the 29 regions and in the two flanks of equal length (Supplementary Figure 4g instead uses a fixed 50 kb flank, so the two interior-against-flank contrasts are not the same comparison), with the GC content of each; the expectation is the genome-wide density of the motif multiplied by the length of the interval, so it conditions on interval length and not on the local base composition that the GC column of this table records. Comparing each region with its own flanks, the observed-over-expected occurrence is indistinguishable inside and outside (TGGCCGGC median 0.994 against 1.007, Wilcoxon signed-rank *p* = 0.970; AAGCCCG 0.923 against 0.923, *p* = 0.936; 29 regions). *Source: Source Data Table SD11 (15_paper_figures/tables/supplementary/SuppTable20b_cluster_motif_occurrence_261021.tsv).* Rows: 29.

## SD12

Interior-against-flank 6mA level of biosynthetic gene clusters, with and without transcription start sites. One row per library: numbers of measurable positions and mean 6mA level in cluster interiors and flanks, the same two quantities after removing positions within the promoter window, and the interior-minus-flank differences for both sets. The promoter-excluded columns separate a cluster effect from the promoter-proximal deficit reported in R5. Positions are the strand-oriented AAGCCCG A₁/A₂ adenines (2,664–2,668 per library), not all genomic adenines, so levels here are not comparable in absolute value with the earlier genome-wide version of this table; the interior-minus-flank difference, which is the quantity used, is unaffected. Computed on the unfiltered (all-read) pileups. *n* = nine libraries, three biological replicates at each of T1, T2 and T3. Rows: 9.

## SD13

Recognition-sequence content, modified-site count and within-region positions for every predicted cluster. One row per antiSMASH region (29): region number and known-cluster name, core-gene expression class, length, GC content, and for each mark the number of recognition-sequence positions the region contains, the number of those that are modified at T1 under the canonical site rule, and the modified positions normalised to the region's own length (0 = region start, 1 = region end). 382 of the 2,909 genome-wide 4mC positions and 50 of the 416 AAGCCCG 6mA positions lie inside a region. Rows: 29.

## SD14

Gene-body methylation level by gene set. Mean gene-body level at T1 (all-read) for three sets: primary metabolism as defined by KEGG BRITE br:sco00001, A09100 Metabolism excluding B-categories 09109–09112 (1,013 genes; brite release 205, retrieved 2026-11-02), antiSMASH cluster genes (926) and all other protein-coding genes (6,057). One row per mark × set with the number of genes carrying a measurable position, the mean, median and quartiles, followed by the Mann–Whitney *U* comparison of the first two sets: 4mC *p* = 0.601 (258 against 189 genes) and 6mA *p* = 0.934 (193 against 113). *Source: Source Data Table SD14 (15_paper_figures/tables/supplementary/SuppTableS5e_geneset_level_comparison_261104.tsv).* Rows: 8.

## SD15

Cluster interior against flanks of the cluster's own length. One row per antiSMASH region × mark × timepoint with the number of measurable positions and the mean level inside the region and over two flanking intervals of the region's own length, one each side. This flank definition differs from the fixed 50 kb flank used in Supplementary Figure 4(g). The paired comparison with the cluster as the unit is null for both marks (4mC signed-rank *p* = 0.562, *n* = 21, median ratio 1.016; 6mA *p* = 0.768, *n* = 19, median ratio 0.987). *Source: Source Data Table SD15 (15_paper_figures/tables/supplementary/SuppTableS5f_cluster_interior_own_length_flank_261104.tsv).* Rows: 174.

**Supplementary Table S6 | The promoter-proximal 6mA deficit.**

## SD16

The promoter-proximal 6mA deficit under four transcription-start-site definitions and two coordinate conventions. The same ratio computed on four transcription-start-site sets (the 3,570 TSSs measured by Jeong *et al.*, the Jeong-derived subset of the comprehensive table, the 4,977 annotated gene starts used as a negative control, and the two combined) at two window widths (±100 and ±300 bp), two adenine universes (both adenines of the motif, and position 1 alone) and two coordinate conventions, giving 32 rows. The convention labelled `corrected` places both coordinates on the same 1-based footing and is the one reported in the main text; `ledger` measures the distance between a 0-based pileup start and a 1-based TSS and is therefore offset by one base. Under `corrected` at ±100 bp on the position-1 universe, the measured TSS set gives 0.9105 and the annotated gene-start control 1.0134. *Source: Source Data Table SD16 (15_paper_figures/tables/supplementary/SuppTable22b_TSS_set_controls_261021.tsv).* Rows: 32.

## SD17

Main statistical tests reported in the text. One row per test: test name, comparison, *n*, test statistic, uncorrected *p*, correction method, corrected *p*, the source of *n*, and the Results section the test belongs to. The table lists the tests whose statistics the text quotes, together with four supporting tests whose statistics it does not — the near-cognate context enrichment of the unassigned adenine calls, and the three Clark–Evans tests of whether dual-modified instances are spatially clustered. It is not an exhaustive inventory of every test run, and tests described only in the Methods or shown only as significance marks in a figure may be absent. *(The protection-zone classifier and its cross-validated AUC have been removed; see Supplementary Table S10 and the main-text non-circular definition. Three further rows — a FIMO threshold, a central-versus-flanking permutation test and a Mann–Whitney on 4mC level inside predicted transcription-factor sites — were removed on 2026-11-09 because the text no longer reports those analyses and their source tables are archived.)* *n* = 26 tests, one per row. Source: `15_paper_figures/tables/supplementary/SuppTable16_statistical_test_summary_261021.tsv`. Rows: 26.

**Supplementary Table S9 | Methyltransferase catalogue, genus-wide comparison and candidate functions.** Rows: 29.

## SD18

*Streptomyces* genus-level AAGCCCG and GCCGGC-hexamer-constituent motif O/E statistics. Per-assembly counts, densities, expected counts and observed/expected ratios for seven motifs (CCGG, AAGCCCG, GATC, GGCCGG, GCCCG, GCGG and the CGACNNNCTCC control) across 833 *Streptomyces* genome assemblies. The GCCGGC hexamer is represented by its CCGG and GGCCGG constituents. The O/E < 0.75 threshold and the REBASE comparison are applied in the text, not in the table. Rows: 833.

## SD19

Reported occupancy of bacterial 6mA methyltransferase target motifs in the published literature. Seventy-five occupancy values read from the tables or Results text of 15 independent primary methylome studies (14 contributing a 6mA value), with the organism, motif, mark, method, the study's own definition of occupancy, the system type and the exact locator for each number; 64 of the values are 6mA, five 4mC and six 5mC. Aggregates over the 63 6mA values whose methyltransferase is active (the one inactive-enzyme value, 0% at ATGCAT, is excluded): median 99.43%, IQR 95.00–99.91%, range 24.70–100%; median of the 14 per-study medians 97.53%. Six values are lower bounds from "> X%" wording and two are complements of a reported unmethylated fraction, so the true median is slightly higher. This is a survey of retrievable open-access papers rather than an exhaustive census, so the median is provisional and the differing definitions are recorded verbatim rather than harmonised. *Source: Source Data Table SD19 (15_paper_figures/tables/supplementary/SuppTable25_published_occupancy_261019.tsv).* Rows: 77.

## SD20

Intergenic 6mA by gene orientation. Non-motif (unassigned) 6mA site rate in intergenic gaps classified by the orientation of the two flanking genes: divergent (a promoter on each side), convergent (a terminator on each side), tandem (one promoter and one terminator), against within-gene positions as the reference. Mantel-Haenszel odds ratios stratified on local GC and pooled true depth, with 95% confidence intervals, at T1-T3 and under two local-GC windows (±250 bp and ±25 bp). Denominators are measurable non-motif positions. The excess is absent from divergent gaps under both windows and at all three timepoints; the stronger claim that convergent gaps are special is not robust, because convergent and tandem are not separable under the ±25 bp window. Source: `11_epigenome_integration/analysis/129_open_questions_261107/`. Rows: 48.


## SD21

AAGCCCG 6mA orientation contrast under eleven nested model specifications. Logistic regression of modification on host-gene orientation (coding versus template strand) at measurable intragenic AAGCCCG instances, at T1-T3, under a ladder running from the crude contrast to the fully adjusted model, each reported with conventional and gene-cluster-robust standard errors. Columns give the odds ratio, 95% confidence interval, p value, instance and host-gene counts, the optimiser, its convergence flag and the model formula. The ladder isolates where the published values diverge: the covariate that moves the estimate is codon phase, not replication-strand direction, which changes it by 2.1 per cent. The specification without the phase term (M5b) supplies the adjusted odds ratios quoted in the text. Source: `11_epigenome_integration/analysis/131_sixmA_round2_261108/`. Rows: 66.

## SD22

Background expectation for the sequence-space overlap of the published cytosine positions. For each of the three sequence spaces used in the reassignment — the GCCGGC/GGCCGG hexamers, contexts containing GCCCG, and their union — the table gives the number of genomic cytosines that fall inside it, the total over both strands, the percentage, the corresponding count among the 572 positions of Pisciotta et al., and the fold enrichment. Computed from the reference genome NC_003888.3 alone, with no methylation call used: cytosines are counted on both strands (C and G of the plus strand), plus-strand membership from occurrences of GCCGGC, GGCCGG and GCCCG, minus-strand membership from their reverse complements GCCGGC, CCGGCC and CGGGC. The union covers 6.712% of all 6,250,890 cytosines against 84.97% of the published positions, a 12.66-fold enrichment. Source: `15_paper_figures/tables/supplementary/SourceDataSD22_overlap_background_261108.tsv`. Rows: 3.

## SD23

Inventory of the pre-registered search space for a second recognition sequence. One row per search route, giving the class of motif architecture, whether the route was covered by an earlier scan or added by the exhaustive scan of 2026-11-07, the enumeration size (words, architectures, patterns or models as applies to that route), the significance criterion used, the highest occupancy any word of the route reached, and that word with its methylated and genomic counts. The architectures sum to 1,399 (149 from the earlier bipartite and contiguous scan, 70 for contiguous k = 12-16, 280 for bipartite with spacer N11-N20, 900 for three-block) and the highest occupancy over all routes is 32.258%, reached by a 14-mer with 31 genomic instances. The routes enumerated as words total 1.526 x 10^10; the IUPAC route is counted as 852,803,616 degenerate patterns and the position-weight route as models, so the three counts are not summed. Supplementary Table S16 carries the per-word detail of the bipartite and contiguous route only. Source: `11_epigenome_integration/analysis/126_sixmA_motif_exhaustive_261107/`. Rows: 7.

## SD24

Promoter-proximal TGGCCGGC inside predicted biosynthetic gene clusters, against the chromosome-wide rate. Supplies the counts burned into Supplementary Figure 4 panel d (127 of 2,720 and 18 of 310) and the TSS coverage stated in its footnote (3,030 of 3,570). Source: `15_paper_figures/tables/supplementary/SourceDataSD24_bgc_promoter_tggccggc_261108.tsv`. Rows: 5.

## SD25

4mC occupancy at promoter-proximal positions inside predicted biosynthetic gene clusters against the rest of the chromosome. Supplies the 19 and 261 positions burned into Supplementary Figure 4 panel e. Source: `15_paper_figures/tables/supplementary/SourceDataSD25_bgc_promoter_vs_rest_4mC_261108.tsv`. Rows: 6.

## SD26

Unassigned 6mA site rate in cluster interiors against their flanks, per timepoint. Supplies the interior-to-flank ratios burned into the Supplementary Figure 4 footnote (0.99, 0.82 and 0.91 at T1-T3). Source: `15_paper_figures/tables/supplementary/SourceDataSD26_unassigned_6mA_interior_flank_261108.tsv`. Rows: 3.

## SD27

Orientation contrast of AAGCCCG 6mA split by host-gene expression quartile. Supplies the four p values of Supplementary Figure 8 panel b and the four ratios in its footnote. Source: `15_paper_figures/tables/supplementary/SourceDataSD27_strand_orientation_expression_261108.tsv`. Rows: 4.

## SD28

Share of modified positions falling in the chromosomal core, per motif series and timepoint, against the composition-matched null. Supplies the six values of Supplementary Figure 2 panel d, the three of panel b, and the null range in the footnote (69.82-69.92%). Source: `15_paper_figures/tables/supplementary/SourceDataSD28_core_share_vs_null_261108.tsv`. Rows: 9.
## SD29

Positions called 4mC in the CCGG context, and the 5mC control, under the canonical site rule. One row per timepoint x context (three timepoints; three nested contexts). The contexts are defined by the position of the modified cytosine: CCGG is the inner C, which is the same position as C5 of TGGCCGGC; GGCCGG is the hexamer containing that C; TGGCCGGC is the octamer. GGCCGG is a subset of CCGG and TGGCCGGC a subset of GGCCGG. Columns give the measurable positions, the number called 4mC, the percentage called, and the number called 5mC. The same number is called in all three contexts at every timepoint, so every 4mC call in a CCGG context lies within TGGCCGGC; the occupancy in the last row of each timepoint reproduces the published TGGCCGGC values (90.11 / 90.40 / 89.15%), which fixes the site rule. No 5mC position was called in any library in any context. Source: canonical pileups `15_paper_figures/scripts/recompute_261015/pileup_cache_261015/*.npz`, recomputed 2026-11-09. Rows: 9.

## SD30

Exhaustive tests of chromosome three-dimensional structure against methylation. One row per axis of testing (four axes, 364 tests in total): the axis, the number of tests, what was tested, and the outcome. Hi-C contact matrices are the L-phase and M-phase data of Deng *et al.* (NCBI BioProject PRJNA891940); methylation is this study under the canonical site rule. The sequence control is the AAGCCCG occurrence count per bin, which correlates with the structural indices at least as strongly as 4mC does, so the six correlations that survive Bonferroni are attributable to sequence composition rather than to the mark. Source: `11_epigenome_integration/analysis/130_hic_6mA_261108/tables/hic_exhaustive_summary_261108.tsv`, run 2026-11-08. Rows: 4.

## SD31

Per-replicate values behind Figure 1b: dry cell weight (g/L) and actinorhodin and undecylprodigiosin concentrations (mg/L) for four independent biological replicate cultures at 12, 24 and 50 h (36 measurements). Per-replicate absorbances, blanks and the full post-hoc output are in `11_epigenome_integration/analysis/106_act_red_dcw/tables/act_red_quantification_261102.tsv` and `act_red_dcw_stats_author_261102.tsv`. Source: `15_paper_figures/tables/supplementary/SuppTable10_growth_and_production_261021.tsv`. Rows: 36.

## SD32

Enrichment of site-level methylation change in predicted biosynthetic gene clusters. One row per test (72 tests): mark, event class (loss or gain), timepoint pair, region class, the observed and expected counts under a stratified permutation null (2,000 draws, strata = local GC quintile over 500 bp x coverage-depth quintile), the permutation p and the Benjamini-Hochberg q. A change event is a position measurable at both timepoints (merged true depth >= 30, reference base matched) that crosses the 50% calling threshold; the robust subset additionally requires |delta level| >= 0.2. No test survives correction. Source: `11_epigenome_integration/analysis/135_bgc_switching_261109/tables/enrichment_permutation.tsv`, run 2026-11-09. Rows: 72.

## SD33

What explains the 3' cytosine-run response of non-motif 6mA. Two stacked tables, identified by the first column. `A1_effect_sizes`: for each of sixteen candidate features per timepoint, the logistic odds ratio adjusted for local GC over 50 bp and log coverage depth, its 95% interval, the Benjamini-Hochberg q across features, and the Mantel-Haenszel odds ratio over GC-decile x depth-tertile strata where the feature is binary. `A2_attenuation`: the residual odds ratio of the cytosine-run exposure when each candidate, and finally all seventeen together, are added to the base model, with the percentage of excess odds removed. DNA-shape terms come from the pentamer model of DNAshapeR 1.34.0 rather than a dinucleotide bendability scale; melting temperature and dG37 use the unified nearest-neighbour set of SantaLucia and Hicks (2004) as shipped in Biopython. Source: `11_epigenome_integration/analysis/134_context_depth_261109/tables/A1_feature_effects_261109.tsv` and `A2_crun_attenuation_261109.tsv`, run 2026-11-09. Rows: 105.

## SD34

How far a nested sequence context narrows the called sites. One row per rung of four ladders x three timepoints: the context, the offset of the modified base within it, the number of genomic occurrences, the number of high-confidence calls it contains, the next shorter context, and the inclusion rate - the percentage of the calls satisfying the shorter context that also satisfy the longer one. Occurrence counts are listed for reference and are never used as a denominator for comparison. Source: `11_epigenome_integration/analysis/134_context_depth_261109/tables/B1_nested_context_inclusion_261109.tsv`, run 2026-11-09. Rows: 75.

## SD35

False-positive modification calls on unmodified DNA. Two stacked tables, identified by the first column. `C_false_positive_rate`: per-replicate read-level call rates for 6mA, 4mC and 5mC at two call-probability thresholds on the ONT all-5-mer synthetic control (16,000 reads over 32 references of 155 bp, every adenine and cytosine canonical), basecalled with the same model chain as this study (`sup@v5.2.0` with `6mA@v1` and `4mC_5mC@v1`) and piled up under the same site rule. `C_level_tail`: the distribution of per-position modified fraction in the control against the unassigned 6mA positions of this study. No control position met the site rule in 597,500 depth-matched draws and no control level exceeded 0.364, against a median of 0.557 at the unassigned sites here. Limits: the control is a synthetic template of about 47% GC read essentially on one strand, so it bounds basecaller model error and not sample-side artefact; the position-count ceiling alone reaches only 2.5e-3, and the 5.0e-6 figure assumes the control positions represent this genome's; the modkit version used for the study pileups is not recorded. Source: `11_epigenome_integration/analysis/132_unmod_control_261109/tables/C_control_fp_summary_261109.tsv` and `C_level_tail_study_vs_control_261109.tsv`, run 2026-11-10. Rows: 21.
