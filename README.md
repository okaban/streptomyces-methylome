# Nanopore methylome and strand-resolved transcriptome of *Streptomyces coelicolor* A3(2) M145

Analysis code for the manuscript

> *Two recognition sequences in the Streptomyces coelicolor M145 methylome,
> one carrying 4mC and 6mA on the same DNA strand, with no transcriptome link*

(manuscript in review).

Nine biological samples across three time points of the developmental cycle
(T1 = 12 h, T2 = 24 h, T3 = 50 h; sample IDs `1-1 1-2 1-3 / 2-1 2-3 2-4 /
3-2 3-3 3-4`, the first digit being the time point). Genomic DNA and total RNA
were extracted from the **same harvested culture** for each sample, so the
methylome and the transcriptome derive from one biological specimen.

## Data availability

Sequencing data are deposited in the NCBI Sequence Read Archive under BioProject
[PRJNA1533939](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1533939)
(study SRP743765). Per-sample run accessions are listed in Supplementary Table 24.

| assay | platform | n | files |
|---|---|---|---|
| methylome | Oxford Nanopore MinION, FLO-MIN114 (R10.4.1), SQK-RBK114-96 | 9 | modification-tagged BAM (MM/ML: 6mA, 4mC, 5mC) aligned to `NC_003888.3` |
| transcriptome | Illumina NovaSeq X Plus, PE150 | 9 | raw FASTQ (Ribo-Zero Plus depletion, NEBNext Ultra II Directional) |

Raw Nanopore signal (POD5) is not accepted by the SRA and is not deposited; it is
retained by the authors and available on request.

Sequencing data are **not** in this repository. Reference genome:
`GCF_000203835.1` / `NC_003888.3`.

## Repository layout

Numbered directories follow the order of the analysis. Each holds
`scripts/` (code), `analysis/` (dated output) and sometimes `data/`.

| | |
|---|---|
| `01_qc` – `04_deseq2` | RNA-seq: fastp QC, HISAT2 alignment, featureCounts, DESeq2 |
| `05_annotation` | genome annotation |
| `11_epigenome_integration` | methylome × transcriptome integration — the core of the paper |
| `15_paper_figures` | current figure and supplementary-table generators, and the tables they emit |
| `13_TF_binding-site` | transcription-factor binding-site provenance |

## Reproducing

Conda environments are exported per stage, e.g.

```
conda env create -f 01_qc/scripts/environment_rnaseq_260127.yml
```

**Scripts hard-code absolute paths under the original author's home
directory.** Roughly 200 files contain such paths; they must be adapted before
the code will run elsewhere. Absolute paths appearing inside *tool output*
(MultiQC, featureCounts, sample tables) are left as-is deliberately: they
record what was actually executed, and rewriting them would falsify the
analysis record.

Large regenerable artifacts are excluded by `.gitignore` — sequencing data,
BLAST/aligner indexes, downloaded NCBI genomes, per-read co-modification
tables. All are reproducible from the deposited data with the scripts here.

## Scope of this archive

This archive holds the code, the per-stage conda environment exports and the
derived tables that the current manuscript rests on. Two classes of material are
deliberately absent: large intermediates that are regenerable from the deposited
BAM and FASTQ files, and the analysis lines that were withdrawn during revision.

An earlier version of this work partitioned regulatory genes into two classes by
their distance to the nearest methylated site, and then scored a classifier whose
labels were derived from that same distance. The analysis was circular; its
results were withdrawn from the manuscript and its scripts and outputs are not
part of this archive. They remain in the project history and with the authors.

## License

MIT — see [LICENSE](LICENSE). The manuscript itself is not covered by this
license.
