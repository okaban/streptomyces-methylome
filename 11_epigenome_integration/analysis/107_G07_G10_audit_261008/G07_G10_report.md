# G07 / G10 — blind-spot audit, S. coelicolor M145 methylome manuscript

Date: 2026-10-08 · Analysis directory: `11_epigenome_integration/analysis/107_G07_G10_audit_261008/`
Every number below was computed this session from the pileups, the genome and the annotation.
Manuscript files were read only.

## Pipeline validation

The canonical site rule — per timepoint, pool the three libraries; pooled true depth
(ΣNvalid_cov + ΣNfail) ≥ 30; pooled ΣNmod/ΣNvalid ≥ 0.5 — applied to the 4mC (ChEBI 21839)
channel at motif offset 4 of TGGCCGGC (TGGC**C**GGC) gives **1,886 / 1,892 / 1,877** sites at
T1 / T2 / T3. The target value of 1,886 is reproduced exactly, and the T2/T3 counts match the
series quoted in R6 (Supplementary Figure 14). Note that this is not the only TGGCCGGC T1 census in
circulation on this project — `canon_values.json` carries 1,801 from a `modkit pileup --no-filtering`
route, and the manuscript supplement also carries 1,289, 1,313 and 1,594. All results below use the
1,886 pooled-rule set named in the audit brief; counts are not transferable between these censuses. The genome carries 2,093 TGGCCGGC instances, so the per-instance
methylation rate is 90.11% at T1.

Under the same rule the AAGCCCG 4mC channel (offsets C3, C5) yields **zero** called sites at any
timepoint — pooled occupancy never reaches 0.5 — and the 6mA channel (A0, A1) yields 363 / 386 / 392
methylated instances out of 1,334. R6's AAGCCCG counts (1,078 sites = 674 4mC + 404 6mA) therefore
do not come from the canonical rule; they come from the weighted high-confidence table, as the
manuscript states. All AAGCCCG results below use the canonical rule.

---

# G07 — Is the Exposed/Shielded partition confounded with recognition-sequence occurrence?

## Answer: **R4 requires retraction of its retention argument.**

### 1. The partition is a recognition-sequence partition

Within the classification window (TSS ± 293 bp), **100% of Exposed promoters carry at least one
TGGCCGGC instance and 0.9% of Shielded promoters do** (Fisher OR = ∞, p = 1.4 × 10⁻¹¹⁵,
Cliff's δ = 0.991). The confound is not partial — it is near-total collinearity. Of the 94
regulatory promoters that contain a TGGCCGGC instance in the window, 85 are Exposed (90.4%),
which is the genome-wide per-instance methylation rate (90.11%) to within 0.3 pp.

Knowing that a promoter contains the recognition sequence tells you its class. The methylation
call adds essentially nothing.

The confound is specific to the sequence used for classification. AAGCCCG shows no difference
between the classes at any window (± 150 bp: 5.9% vs 4.5%, p = 0.58; ± 293 bp: 8.2% vs 7.3%,
p = 0.67; ± 500 bp: p = 0.87; ± 2 kb: 35.3% vs 43.3%, p = 0.17). This rules out a generic
promoter-composition explanation and points at the circularity directly.

Collinearity is not an artefact of the 293 bp cut-off. At ± 2 kb — four times the window width,
where class membership is not definitional — Exposed is still 100% motif-containing against 54.9%
of Shielded (p = 2.2 × 10⁻²¹).

### 2. The conditioned null

Null model: hold each promoter's recognition-sequence count in TSS ± 293 bp fixed, and permute the
per-instance (T1, T2, T3) methylation-state vectors across all 2,093 genomic TGGCCGGC instances
(10,000 permutations). This preserves the exact marginal site counts (1,886 / 1,892 / 1,877), the
T1→T2→T3 dependence within an instance, and the positional distribution of the recognition
sequence; it destroys only the association between *which* instance is methylated and *where* it is.
The partition is rebuilt inside every permutation.

| R4 statistic | observed | null mean [95%] | p |
|---|---|---|---|
| Exposed promoters | 85 | 85.31 [80, 90] | 1.00 |
| retain a site within 293 bp at T2 | 82 | 83.06 [77, 89] | 0.82 |
| retain a site within 293 bp at T3 | 82 | 82.11 [76, 88] | 1.00 |
| Shielded gaining a site at T2 | 2 | 2.45 [0, 6] | 1.00 |
| retention rate at T2 | 0.9647 | 0.9736 [0.9390, 1.0000] | 0.59 |
| retention rate at T3 | 0.9647 | 0.9625 [0.9186, 1.0000] | 1.00 |

Every statistic sits at the centre of its null. There is no retention signal above what the
placement of the recognition sequence plus the genome-wide per-instance methylation rates already
produce.

### 3. The specific R4 numbers do not reproduce

The stored `nearest_GCCGGC_T2` and `nearest_GCCGGC_T3` columns of
`SuppTable_regulatory_gene_classification_n1051.tsv` are **stale**. They agree with a recomputation
on the canonical site set for 22.7% and 1.0% of the 1,051 genes respectively, and they give a median
nearest-site distance for Exposed promoters of 57,885 bp at T2 and 196,231 bp at T3 — against
181 bp and 178 bp when recomputed. The T1 column is also built on a smaller site set: 338 of 1,051
stored distances are larger than recomputed, only 3 smaller.

Recomputing the partition on the canonical 1,886-site set gives **85 Exposed, not 62**. Twenty-four
genes become Exposed and one leaves; 21 of the 24 are arm-localised, which is the direction that
would also soften R3's core-localisation claim (not tested here).

Holding the published n = 62 set fixed and recomputing only the distances:

| R4 as published | recomputed on canonical sites |
|---|---|
| median nearest-site distance 172 → 181 → 184 bp | 172 → 181 → 178.5 bp |
| 56 of 62 retain at T2 | 59 of 62 |
| 45 of 62 retain at T3 | 59 of 62 |
| 20 of 989 Shielded gain at T2 | 25 of 989 |

The T3 retention count (45) is the one that fails hardest; it is the number that makes the time
course look like decay.

### 4. What R4 can and cannot say

- **Cannot stand:** "Exposed promoters retain their methylation across the metabolic transition",
  read as a property of Exposed promoters. The retention rate is indistinguishable from the
  motif-conditioned null (p = 0.82 at T2, p = 1.00 at T3), and the Exposed/Shielded contrast is a
  restatement of whether a TGGCCGGC recognition sequence is present within 293 bp of the TSS.
- **Cannot stand as stated:** the counts 62 / 56 / 45 and the 20-of-989 comparator. They were
  computed from columns that do not reproduce under the canonical rule.
- **Survives:** the manuscript's own hedge — "this retention is the expected consequence of a
  largely stable core methylome rather than evidence of promoter-specific maintenance" — is now the
  measured result rather than a caution. The permissive conclusion of R4 (promoter methylation does
  not instruct the direction of expression) is untouched by G07; G07 tests the retention premise
  only.
- **Minimum repair:** state the partition for what it is — promoters that do and do not contain the
  recognition sequence near the TSS — report the conditioned null, and drop the retention sentence
  or recast it as the null result it is.

---

# G10 — Does R6 survive collapsing genes to transcription units?

## Answer: **the TU roll-up does not rescue R6; the CDS-dominance argument stays retracted, and the enrichment argument fails for a different reason.**

### 1. Transcription units

Chromosome NC_003888.3, protein_coding genes (n = 7,593), same strand, intergenic gap ≤ 50 bp →
**5,440 units** (mean 1.396 genes per unit; 4,134 single-gene units), reproducing the published
5,432 to 0.15%.

**The published CDS universe is wrong, and it is not the one used for the TU count.** Taking
protein_coding genes on *all* contigs — chromosome plus SCP1 (NC_003903.1) and SCP2 (NC_003904.1) —
and merging their coordinates reproduces the G09 audit exactly: 6,129 intervals, coding fraction
**87.85%**, TGGCCGGC motif null **93.65%**, T1 observed **93.90%**. That universe folds 403
plasmid-contig genes into the chromosome by treating their coordinates as chromosome coordinates.
M145 is plasmid-free. The effect is to call an extra 0.66 pp of the chromosome coding (87.85% vs the
correct 87.20%), and the published 93.65% and 93.90% inherit it. Meanwhile the published TU count of
5,432 reproduces only on the chromosome-only universe. The two published numbers are built on
different gene sets.

Neither error changes a conclusion — the sensitivity analysis below runs everything under both
universes — but the three percentages 87.85 / 93.65 / 93.90 should be corrected to 87.20 / 93.02 /
93.27 before they appear in print.

### 2. Pseudo-replication is real but modest for the positional claim

| set | sites | units | sites per occupied unit | effective-N loss |
|---|---|---|---|---|
| TGGCCGGC 4mC T1 | 1,886 | 1,406 | 1.341 | 25.5% |
| TGGCCGGC 4mC T2 | 1,892 | 1,409 | 1.343 | 25.5% |
| TGGCCGGC 4mC T3 | 1,877 | 1,408 | 1.333 | 25.0% |
| AAGCCCG 6mA T1 | 363 | 337 | 1.077 | 7.2% |
| AAGCCCG 6mA T2 | 386 | 356 | 1.084 | 7.7% |
| AAGCCCG 6mA T3 | 392 | 358 | 1.095 | 8.7% |

A quarter of the GCCGGC sample is clustered duplication, not independent observation.

### 3. CDS dominance: the roll-up confirms the retraction

Against the genome coding fraction (87.20%, chromosome-only) the GCCGGC CDS fraction is highly
significant at every timepoint (93.27% at T1, binomial p = 1.5 × 10⁻¹⁷). Against the motif-position
null — the CDS fraction of all 2,093 TGGCCGGC instances, 93.02% — it is not (p = 0.72 / 0.86 / 1.00
at T1 / T2 / T3). Under the published all-contig universe the same contrast gives 93.90% observed vs
93.65% null — reproducing the G09 audit percentages exactly (93.90% vs 93.65%). My two-sided
binomial p is 0.706 against that audit's 0.671; the gap is a rounding artefact of the null, not a
disagreement — the exact null is p₀ = 0.936455 (1,960 of 2,093 instances in CDS under that universe) and gives p = 0.706, while
entering the rounded 0.9364 gives p = 0.671. A two-sided binomial p is unstable to the fourth
decimal of p₀ near a discrete jump; both values support the same conclusion of no excess.

Collapsing to transcription units changes nothing. The mean per-TU CDS fraction is 92.83% at T1
against a permutation null of 92.67% [92.29, 93.11], **p = 0.46** (10,000 permutations of methylation
state across motif instances, positions and TU membership held fixed); T2 p = 0.39, T3 p = 0.45. The
CDS-majority-of-sites statistic agrees (p = 0.37 / 0.83 / 0.24). Repeating the roll-up under the
published all-contig CDS universe gives p = 0.31 / 0.64 / 0.66 — the same answer.

The TU roll-up was the remaining hope for R6's CDS argument. It does not deliver: the claim was
already dead against the correct null, and it is equally dead at the correct unit of analysis.

**One exception.** AAGCCCG 6mA does show a genuine CDS excess over its own motif-position null:
90.63% vs 86.43% at T1 (binomial p = 0.017; TU-level permutation p = 0.010), 90.67% vs 86.43% at T2
(p = 0.014; TU p = 0.006), and 89.29% at T3 (p = 0.10; TU p = 0.083). This survives the roll-up —
AAGCCCG sites are barely clustered (1.08 sites per unit). It is a small effect on a small sample
(363–392 sites) and the six tests are uncorrected; Bonferroni across them leaves T2 at p = 0.035 and
T1 at p = 0.060. Provisionally: the two systems differ, and only the dual-modified AAGCCCG system
shows CDS preference beyond sequence availability. This is the opposite assignment from the one R6
implies, where GCCGGC carries the CDS-dominance story.

### 4. The position-stratified enrichment: the unit is not the problem

Recomputing R6's CDS-internal GCCGGC gene set under the canonical rule gives 1,759 sites → 1,451
genes → 1,317 transcription units (R6 reports 1,413 sites / 1,159 genes from its own site table).

R6's three headline COG enrichments reproduce at gene level against an all-gene background —
amino acid metabolism OR 1.78 (FDR 3 × 10⁻⁴), carbohydrate metabolism OR 1.73 (FDR 0.0085), signal
transduction OR 1.45 (FDR 0.016) — directionally matching the published OR 2.03 / 1.98 / 1.76. They
**survive the TU roll-up** against an all-TU background, and get stronger (OR 2.44 / 2.24 / 1.72, all
FDR < 10⁻⁴). Operon pseudo-replication is not what breaks this analysis.

What breaks it is the background. 91.4% of genes (92.5% of units) that contain a CDS-internal
TGGCCGGC instance are methylated at T1, so "marked" and "motif-containing" are nearly the same set,
exactly as in G07. Against a motif-matched background every enrichment disappears (gene level:
FDR 1.00 / 1.00 / 1.00; TU level: FDR 0.95 / 1.00 / 0.95). That test is underpowered — the
motif-matched background is only 137 genes / 107 units — so it is **not** evidence of absence. The
defensible statement is structural rather than statistical: the marked set is the motif-containing
set, so the COG enrichment of the marked set is an enrichment of the recognition sequence in
particular functional classes, not of methylation.

### 5. What R6 can and cannot say

- **Stays retracted:** "CDS-dominant distribution supports the permissive nature of methylation."
  The TU roll-up does not restore it (p = 0.46 at T1). The CDS fraction of GCCGGC 4mC sites equals
  the CDS fraction of GCCGGC recognition sequences, at site level and at transcription-unit level,
  at all three timepoints.
- **Needs rewriting, not retracting:** the position-stratified COG/KEGG enrichments. They are robust
  to the unit of analysis; they are claims about where the recognition sequence sits, not about
  where methylation is deployed. Either re-state them with a motif-matched background and the
  honest power caveat, or re-frame them as sequence-distribution results.
- **Survives unchanged:** the descriptive claim that the positional composition of the methylome is
  constant across the time course. The site-level CDS fraction moves 93.27 → 92.92 → 93.02% for
  GCCGGC; the TU-level statistic moves 92.83 → 92.50 → 92.51%.
- **New, provisional:** AAGCCCG 6mA shows a CDS excess over its motif-position null that survives
  the TU roll-up at T1 and T2. Small n, uncorrected across six tests — needs its own treatment
  before it enters the manuscript.

---

## Files

- `tables/G07_exposed_shielded_motif_control_261008.csv` — 33 rows: validation, partition recompute, confound tests, conditioned null
- `tables/G07_per_gene_motif_counts_261008.csv` — per-gene recognition-sequence counts and recomputed distances, 1,051 regulatory genes
- `tables/G10_transcription_unit_rollup_261008.csv` — 63 rows: TU construction, pseudo-replication, CDS tests at both units, enrichment sets, COG re-test, CDS-universe diagnosis, sensitivity under the published universe
- `tables/G10_COG_enrichment_full_261008.csv` — all COG categories, unit × background
- `tables/G10_gene_to_transcription_unit_261008.csv` — 7,593 genes → 5,440 units
- `tables/canonical_motif_sites_261008.tsv` — per-instance site calls, both motifs, three timepoints
- `figures/G07_exposed_shielded_motif_control_261008.png`
- `figures/G10_transcription_unit_rollup_261008.png`

## Caveats

1. The published 87.85% coding fraction and the 93.65%/93.90% CDS fractions derive from a gene set
   that includes 403 plasmid-contig (SCP1/SCP2) genes whose coordinates are overlaid on the
   plasmid-free M145 chromosome. The corrected chromosome-only values are 87.20%, 93.02% and
   93.27%. Every test in this report was run under both universes and no conclusion changes.
   The published TU count 5,432 uses the chromosome-only universe, so the two published numbers are
   mutually inconsistent.
2. The motif-matched-background COG test has 137 genes / 107 units in the background. Its null
   results are not evidence of absence and are labelled provisional.
3. The AAGCCCG CDS-excess result rests on 363–392 sites and six uncorrected tests. Provisional.
4. R3's core-localisation claim is not tested here, but the 24 genes that become Exposed under the
   canonical rule are 21/24 arm-localised, which bears on it.
