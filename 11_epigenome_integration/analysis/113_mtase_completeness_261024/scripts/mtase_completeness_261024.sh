#!/usr/bin/env bash
# 113_mtase_completeness_261024 — annotation-independent MTase census of S. coelicolor M145
# conda env: epi-rna.  Run from the analysis directory.
set -euo pipefail
REF="/Users/okaban/Library/CloudStorage/Dropbox-SFC-CNS/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1"
PFAM=/Users/okaban/bioinfo/data/pfam/Pfam-A.hmm
# 1. chromosome proteome (NC_003888.3 CDS carrying a protein_id) -> data/M145_chrom_proteome.faa (7,553)
# 2. full Pfam-A at the gathering threshold
hmmsearch --cpu 12 --cut_ga --domtblout data/pfam_full_chrom.domtbl "$PFAM" data/M145_chrom_proteome.faa > logs/hmmsearch_full.log
# 3. permissive pass, core DNA-MTase models, whole assembly incl. SCP1/SCP2
cp "$REF/protein.faa" data/M145_all_proteome.faa
hmmfetch --index "$PFAM"
for m in N6_Mtase N6_N4_Mtase MmeI_Mtase Eco57I DNA_methylase MethyltransfD12 Dam MT-A70 EcoRI_methylase N6-adenineMlase; do
  hmmfetch "$PFAM" "$m"
done > data/core_mtase_models.hmm
hmmsearch --cpu 12 -E 10 --domE 10 --domtblout data/core_permissive_all.domtbl \
  data/core_mtase_models.hmm data/M145_all_proteome.faa > logs/hmmsearch_permissive.log
