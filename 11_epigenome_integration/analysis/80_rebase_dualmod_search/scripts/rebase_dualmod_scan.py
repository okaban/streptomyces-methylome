#!/usr/bin/env python3
"""D2: REBASE-wide scan for single recognition motifs carrying BOTH m4C and m6A
(the cross-genus context for the AAGCCCG same-strand dual modification).

Parses bairoch.txt (REBASE v602). For each enzyme record: organism (OS),
enzyme type (ET), recognition motif (RS, cut coords stripped), and the set of
methylation types + strand signs from MS lines. Then:
  (1) single-enzyme dual specificity (one enzyme deposits both m4C and m6A);
  (2) per-motif dual (a recognition motif has both an m4C and an m6A enzyme).
"""
import re
from pathlib import Path
from collections import defaultdict

REBASE = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/rebase/bairoch.txt")
OUT = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/80_rebase_dualmod_search/tables")

records = []
cur = {}
for line in REBASE.read_text(errors="replace").splitlines():
    tag = line[:2]
    val = line[5:].strip()
    if tag == "ID":
        cur = {"id": val, "os": "", "et": "", "rs": "", "ms": []}
    elif tag == "OS":
        cur["os"] = val
    elif tag == "ET":
        cur["et"] = val
    elif tag == "RS":
        cur["rs"] = val
    elif tag == "MS":
        cur["ms"].append(val)
    elif line.startswith("//"):
        if cur:
            records.append(cur)
        cur = {}

def motif_of(rs):
    # "GCCGGC, 2;" -> "GCCGGC" ; take first recognition token, strip coords
    if not rs:
        return None
    first = rs.split(";")[0]
    m = first.split(",")[0].strip().upper()
    return m if re.fullmatch(r"[ACGTRYSWKMBDHVN]+", m or "") else None

def mods_strands(ms_lines):
    """return (set of mod types, set of strand signs) from MS lines.
    REBASE notation: m6A, m5C, and N4-methyl-C written as 'Nm4C' — so match the
    mod label as a substring inside the parens, not anchored to '('."""
    mods, strands = set(), set()
    # 2026-10-26: 修飾ごとに位置が解決しているかを別に持つ。従来は strands 集合しか
    # 返しておらず、片方の位置が '?' の記録（M2.Pae95I: 'MS 2(Nm4C),?(m6A)'）で
    # strands == {'+'} となり same-strand と誤分類していた。'?' は未決定であって
    # プラス鎖ではない。この 1 件は REBASE 二重 51 件の 51 件目で、原稿の新規性主張
    # （短い連続モチーフで両マークを持つ系は無い）の分岐点だった。
    resolved = {}
    for ms in ms_lines:
        for pos, label in re.findall(r"(-?\??\d*)\(([^)]*)\)", ms):
            hit = None
            for canon in ("m4C", "m5C", "m6A"):
                if canon in label:
                    mods.add(canon)
                    hit = canon
            if hit is not None:
                resolved.setdefault(hit, False)
            if pos.startswith("-"):
                strands.add("-")
                if hit is not None:
                    resolved[hit] = True
            elif pos and pos != "?":
                strands.add("+")
                if hit is not None:
                    resolved[hit] = True
    return mods, strands, resolved

motif_mods = defaultdict(set)            # motif -> set of mod types (across enzymes)
motif_orgs = defaultdict(set)
single_dual = []                         # enzymes with both m4C and m6A
for r in records:
    motif = motif_of(r["rs"])
    if not motif:
        continue
    mods, strands, resolved = mods_strands(r["ms"])
    if not mods:
        continue
    motif_mods[motif] |= mods
    motif_orgs[motif].add(r["os"])
    if "m4C" in mods and "m6A" in mods:
        both_resolved = resolved.get("m4C", False) and resolved.get("m6A", False)
        if not both_resolved:
            cls = "undetermined"
        elif "-" in strands:
            cls = "opposite-strand"
        else:
            cls = "same-strand"
        single_dual.append((r["id"], r["os"], motif, sorted(mods), cls))

dual_motifs = {m: mods for m, mods in motif_mods.items()
               if "m4C" in mods and "m6A" in mods}

tot_motifs = len(motif_mods)
print(f"REBASE v602 records parsed: {len(records)}")
print(f"distinct recognition motifs with >=1 methylation annotation: {tot_motifs}")
print(f"\n(1) SINGLE-ENZYME dual 4mC+6mA specificity: {len(single_dual)} enzymes")
for eid, os, motif, mods, strand in single_dual[:40]:
    short = "SHORT<=7bp" if len(motif) <= 7 else ""
    strep = "STREPTOMYCES" if "streptomyces" in os.lower() else ""
    print(f"   {eid:18s} {motif:14s} {','.join(mods):10s} {strand:16s} {short:10s} {strep}  [{os}]")

print(f"\n(2) PER-MOTIF dual (motif has both an m4C and an m6A enzyme): "
      f"{len(dual_motifs)} of {tot_motifs} motifs "
      f"({100*len(dual_motifs)/tot_motifs:.2f}%)")
short_dual = {m: v for m, v in dual_motifs.items() if len(m) <= 7}
print(f"   of which short (<=7 bp): {len(short_dual)}")
for m in sorted(dual_motifs, key=len)[:40]:
    strep = any("streptomyces" in o.lower() for o in motif_orgs[m])
    print(f"   {m:16s} ({len(m)}bp) orgs={len(motif_orgs[m])}"
          f"{'  [has Streptomyces]' if strep else ''}")

# is AAGCCCG present in REBASE at all?
print(f"\nAAGCCCG in REBASE motif list: {'AAGCCCG' in motif_mods}")
OUT.mkdir(parents=True, exist_ok=True)
with open(OUT / "single_enzyme_dual.tsv", "w") as fh:
    fh.write("enzyme\torganism\tmotif\tmods\tstrand\n")
    for eid, os, motif, mods, strand in single_dual:
        fh.write(f"{eid}\t{os}\t{motif}\t{'+'.join(mods)}\t{strand}\n")
with open(OUT / "per_motif_dual.tsv", "w") as fh:
    fh.write("motif\tlen\tn_organisms\tmods\thas_streptomyces\n")
    for m in sorted(dual_motifs, key=len):
        strep = any("streptomyces" in o.lower() for o in motif_orgs[m])
        fh.write(f"{m}\t{len(m)}\t{len(motif_orgs[m])}\t{'+'.join(sorted(dual_motifs[m]))}\t{strep}\n")
print(f"\n[wrote] {OUT}/single_enzyme_dual.tsv, per_motif_dual.tsv")
