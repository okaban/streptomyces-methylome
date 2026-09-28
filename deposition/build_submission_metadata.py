#!/usr/bin/env python3
"""Emit the filled SRA/BioSample submission tables.

Author-supplied facts (2026-09-10 Q&A) are the CONSTANTS below; everything
else is measured from the files. Re-run after editing a constant.
"""
import csv, glob, os, subprocess

BIOPROJECT = "PRJNA1533939"      # obtained 2026-09-23 (SUB16507616)
ASSEMBLY   = "GCF_000203835.1"   # ASM20383v1, S. coelicolor A3(2); required for aligned (BAM) rows
SAMPLE_TYPE = "cell culture"     # *sample_type is MANDATORY in BioSample package Microbe.1.0

METHYL = "/Users/okaban/bioinfo/methyl/260102_M145/analysis"
RNARAW = ("/Users/okaban/Library/CloudStorage/Dropbox-SFC-CNS/Takeda Tomoki/"
          "my_projects/2026/M145_RNA-seq/data/RH250715199/01_RawData")
NANO   = "/Users/okaban/Library/CloudStorage/Dropbox-SFC-CNS/Takeda Tomoki/nanopore"
PROV   = "deposition/nanopore_run_provenance.tsv"
OUT    = "deposition"

# --- author-supplied (Q&A 2026-09-10) --------------------------------------
ORGANISM   = "Streptomyces coelicolor A3(2)"
STRAIN     = "M145"
ISOLATION  = "R5 liquid culture"
COLLECTED  = "2025-12"          # time-course sampling 2025-12-04..07
GEO        = "Japan"
LAB_HOST   = "not applicable"
# --- vendor report RH250715199 ---------------------------------------------
ILLUMINA_MODEL = "Illumina NovaSeq X Plus"
RNA_DESIGN = ("Total RNA; rRNA depletion with Illumina Ribo-Zero Plus; "
              "strand-specific library with NEBNext Ultra II Directional RNA "
              "Library Prep Kit; 150 bp paired-end")
# --- MinKNOW sample sheets / final_summary ---------------------------------
ONT_MODEL  = "MinION"
ONT_DESIGN = ("Native genomic DNA, no PCR; SQK-RBK114-96 rapid barcoding; "
              "FLO-MIN114 (R10.4.1) flow cell, 400 bps; basecalled with "
              "modification calling (6mA, 4mC, 5mC), aligned to NC_003888.3 "
              "with minimap2 -ax map-ont -y (MM/ML tags retained)")

TP = {"1": "12", "2": "24", "3": "50"}
SAMPLES = ["1-1","1-2","1-3","2-1","2-3","2-4","3-2","3-3","3-4"]

def size(p):
    try: return os.path.getsize(p)
    except OSError: return 0

def dsize(p):
    return sum(size(os.path.join(dp,f)) for dp,_,fn in os.walk(p) for f in fn)

# ---------------- BioSample (9 rows, one per biological sample) ------------
# Column names and order follow BioSample package Microbe.1.0
# (https://www.ncbi.nlm.nih.gov/biosample/docs/packages/Microbe.1.0/?format=xml):
# mandatory = collection_date, geo_loc_name, sample_type; at least one of
# strain / isolate / host / isolation_source. Asterisks match NCBI's template.
bios = []
for s in SAMPLES:
    bios.append({
        "*sample_name": f"M145_{s.replace('-','_')}",
        "sample_title": f"S. coelicolor M145, {TP[s[0]]} h, biological replicate {s.split('-')[1]}",
        "bioproject_accession": BIOPROJECT,
        "*organism": ORGANISM,
        "strain": STRAIN,
        "isolate": "not applicable",
        "host": LAB_HOST,
        "isolation_source": ISOLATION,
        "*collection_date": COLLECTED,
        "*geo_loc_name": GEO,
        "*sample_type": SAMPLE_TYPE,
        "description": (f"Time point T{s[0]} ({TP[s[0]]} h after inoculation). "
                        "Genomic DNA and total RNA were extracted from the same "
                        "harvested culture, so the nanopore methylome and the "
                        "Illumina transcriptome for this sample derive from one "
                        "biological specimen.")})

# ---------------- SRA runs (18: 9 ONT BAM + 9 Illumina fastq pairs) --------
# Columns and order are exactly the SRA_data sheet of NCBI's SRA_metadata.xlsx
# (ftp-trace.ncbi.nlm.nih.gov/sra/metadata_table/). `assembly` is required for
# aligned data, so it is filled on the BAM rows only. No extra columns: file
# sizes and checksums go to upload_manifest.tsv instead.
SRA_COLS = ["sample_name","library_ID","title","library_strategy","library_source",
            "library_selection","library_layout","platform","instrument_model",
            "design_description","filetype","filename","filename2","filename3",
            "filename4","assembly","fasta_file"]

# 3-2 and 3-4 were topped up on extra flow cells and merged before pileup
# (scripts/reanalysis_merge_pipeline.sh, 2026-02-24); the deposited BAM is the merge.
MERGED = {"3-2": ". Reads from two additional MinION flow cells (2026-02-16, 2026-02-17)"
                 " were aligned identically and merged into this BAM before pileup.",
          "3-4": ". Reads from one additional MinION flow cell (2026-02-16)"
                 " were aligned identically and merged into this BAM before pileup."}

runs, manifest = [], []
for s in SAMPLES:
    bam = f"{METHYL}/{s}_mapped.bam"
    lib = f"M145_{s.replace('-','_')}_ONT"
    runs.append(dict(
        sample_name=f"M145_{s.replace('-','_')}",
        library_ID=lib,
        title=f"Nanopore methylome of S. coelicolor M145, {TP[s[0]]} h",
        library_strategy="WGS", library_source="GENOMIC",   # random shotgun of the whole genome; WGS is the listed term, OTHER means "not listed"
        library_selection="RANDOM", library_layout="single",
        platform="OXFORD_NANOPORE", instrument_model=ONT_MODEL,
        design_description=ONT_DESIGN + MERGED.get(s, ""), filetype="bam",
        filename=os.path.basename(bam), filename2="", filename3="", filename4="",
        assembly=ASSEMBLY, fasta_file=""))
    manifest.append(dict(library_ID=lib, filename=os.path.basename(bam),
                         local_path=bam, size_bytes=size(bam),
                         size_gb=round(size(bam)/1073741824, 3), md5=""))
for s in SAMPLES:
    u = s.replace('-','_')
    lib = f"M145_{u}_RNA"
    f1, f2 = f"{RNARAW}/M145_{u}_1.fastq.gz", f"{RNARAW}/M145_{u}_2.fastq.gz"
    runs.append(dict(
        sample_name=f"M145_{u}",
        library_ID=lib,
        title=f"Strand-specific RNA-seq of S. coelicolor M145, {TP[s[0]]} h",
        library_strategy="RNA-Seq", library_source="TRANSCRIPTOMIC",
        library_selection="Inverse rRNA", library_layout="paired",
        platform="ILLUMINA", instrument_model=ILLUMINA_MODEL,
        design_description=RNA_DESIGN, filetype="fastq",
        filename=os.path.basename(f1), filename2=os.path.basename(f2),
        filename3="", filename4="", assembly="", fasta_file=""))
    for f in (f1, f2):
        manifest.append(dict(library_ID=lib, filename=os.path.basename(f),
                             local_path=f, size_bytes=size(f),
                             size_gb=round(size(f)/1073741824, 3), md5=""))
assert [k for k in runs[0]] == SRA_COLS, "SRA column order drifted from the NCBI template"

# checksums come from the two verified tables, not recomputed here
_md5 = {}
for tsv, kf, mf in ((f"{OUT}/ont_bam_md5.tsv", "filename", "md5"),
                    (f"{OUT}/raw_checksum_report.tsv", "filename", "computed_md5")):
    if os.path.exists(tsv):
        for r in csv.DictReader(open(tsv), delimiter="\t"):
            _md5[r[kf].strip()] = r[mf].strip()
for m in manifest:
    m["md5"] = _md5.get(m["filename"], "")

# ---------------- pod5 upload manifest -------------------------------------
runroot = {}
for d in os.listdir(NANO):
    p = os.path.join(NANO, d, "no_sample_id")
    if os.path.isdir(p):
        for r in os.listdir(p):
            runroot[r] = os.path.join(p, r)
pod5 = []
for r in (csv.DictReader(open(PROV), delimiter="\t") if os.path.exists(PROV) else []):
    d = os.path.join(runroot.get(r["run"], ""), "pod5_pass", r["barcode"])
    pod5.append(dict(sample_name=f"M145_{r['sample'].replace('-','_')}",
                     run=r["run"], barcode=r["barcode"],
                     pct_of_sample_reads=r["pct_of_sampled"],
                     path=d, size_gb=round(dsize(d)/1073741824, 2) if os.path.isdir(d) else 0))

for name, rows in (("biosample_attributes.tsv", bios),
                   ("sra_metadata.tsv", runs),
                   ("upload_manifest.tsv", manifest),
                   ("pod5_upload_manifest.tsv", pod5)):
    if not rows:
        print(f"{name}: SKIPPED (no input)"); continue
    with open(f"{OUT}/{name}", "w", newline="\n") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()),
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"{name}: {len(rows)} rows")

# ascp/ftp file lists — one per source directory, so no staging copy is needed
for tag, pred in (("ont", lambda m: m["filename"].endswith(".bam")),
                  ("illumina", lambda m: m["filename"].endswith(".fastq.gz"))):
    sel = [m for m in manifest if pred(m)]
    with open(f"{OUT}/upload_list_{tag}.txt", "w", newline="\n") as fh:
        fh.write("".join(m["local_path"] + "\n" for m in sel))
    print(f"upload_list_{tag}.txt: {len(sel)} files, "
          f"{sum(m['size_bytes'] for m in sel)/1073741824:.2f} GB")

# blank accession worksheet for the author to fill as the portal returns IDs
acc = [dict(sample_name=b["*sample_name"], biosample_accession="",
            ont_library_ID=f"{b['*sample_name']}_ONT", ont_run_accession="",
            rna_library_ID=f"{b['*sample_name']}_RNA", rna_run_accession="")
       for b in bios]
with open(f"{OUT}/accession_record.tsv", "w", newline="\n") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(acc[0].keys()),
                       lineterminator="\n")
    w.writeheader(); w.writerows(acc)
print(f"accession_record.tsv: {len(acc)} rows (blank, author fills)")

print(f"\nBioSamples {len(bios)} | SRA runs {len(runs)} "
      f"(ONT {sum(1 for r in runs if r['platform']=='OXFORD_NANOPORE')}, "
      f"Illumina {sum(1 for r in runs if r['platform']=='ILLUMINA')})")
print(f"run file volume: {sum(m['size_gb'] for m in manifest):.2f} GB "
      f"across {len(manifest)} files")
print(f"pod5 volume:     {sum(p['size_gb'] for p in pod5):.2f} GB "
      f"across {len(pod5)} dirs")
print("missing files:", [m['filename'] for m in manifest if m['size_bytes'] == 0] or "none")
print("missing md5:   ", [m['filename'] for m in manifest if not m['md5']] or "none")

# mandatory fields only: BioSample = *-prefixed columns; SRA = every column
# except the optional filename2-4 / assembly / fasta_file slots
OPTIONAL_SRA = {"filename2","filename3","filename4","assembly","fasta_file"}
blank  = [f"BioSample {r['*sample_name']}:{k}" for r in bios
          for k,v in r.items() if k.startswith("*") and not v]
blank += [f"SRA {r['library_ID']}:{k}" for r in runs
          for k,v in r.items() if k not in OPTIONAL_SRA and not v]
blank += [f"SRA {r['library_ID']}:assembly" for r in runs
          if r["filetype"] == "bam" and not r["assembly"]]
print("empty required fields:", blank or "none")
