"""Per-genome strand-specific k-mer counting for the HAAGCCCG flank-depletion test.

Conventions (documented in the report):
  * All contigs of an assembly are concatenated with a single 'N' separator, so no
    counted word can span a contig boundary.
  * Counts are STRAND-SPECIFIC and taken on BOTH strands: each word is counted in the
    forward concatenate F and again in revcomp(F). For a non-palindromic word this
    enumerates distinct genomic loci; for a palindromic core (GCCGGC) every hexamer
    locus is counted once per strand, which is what makes a strand-specific flank
    (-1) breakdown well defined.
  * Overlapping occurrences are counted (GGC-prefixed words can self-overlap), via a
    per-length regex with a capturing lookahead.
  * Ambiguity codes other than ACGT are left in place; they simply never match.
"""
import re, sys, gzip, json
from collections import Counter

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# ---- word sets ------------------------------------------------------------
BASES = "ACGT"
# 6mA family: core 7-mer AAGCCCG, variable -1 flank -> 8-mer
MID6_A, CORE7_A = "AAGCCC", "AAGCCCG"
# 4mC family: core 7-mer GGCCGGC, variable -1 flank -> 8-mer (TGGCCGGC is the methylated form)
MID6_C, CORE7_C = "GGCCGG", "GGCCGGC"
# secondary: bare hexamer GCCGGC and its -1 breakdown
HEX6_C = "GCCGGC"

WORDS = set()
WORDS |= {MID6_A, CORE7_A, MID6_C, CORE7_C, HEX6_C}
for b in BASES:
    WORDS |= {b + CORE7_A, b + MID6_A, b + CORE7_C, b + MID6_C, b + HEX6_C}
# controls / composition reference words used by the manuscript
WORDS |= {"GATC", "CCGG"}
WORDS = sorted(WORDS)

BY_LEN = {}
for w in WORDS:
    BY_LEN.setdefault(len(w), []).append(w)
PATTERNS = {L: re.compile("(?=(" + "|".join(ws) + "))")
            for L, ws in BY_LEN.items()}


def count_words(seq, tally):
    for L, pat in PATTERNS.items():
        for m in pat.finditer(seq):
            tally[m.group(1)] += 1


def read_fasta_concat(path):
    op = gzip.open if path.endswith(".gz") else open
    chunks, n_contigs = [], 0
    with op(path, "rt") as fh:
        cur = []
        for line in fh:
            if line[0] == ">":
                if cur:
                    chunks.append("".join(cur)); cur = []
                n_contigs += 1
            else:
                cur.append(line.strip())
        if cur:
            chunks.append("".join(cur))
    return ("N".join(chunks)).upper(), n_contigs


def scan(path, accession):
    F, n_contigs = read_fasta_concat(path)
    R = F.translate(COMP)[::-1]
    tally = Counter()
    count_words(F, tally)
    count_words(R, tally)
    # base composition from the forward concatenate only (both-strand composition is
    # then symmetric by construction: p_A=p_T, p_C=p_G)
    bc = Counter(F)
    nA, nC, nG, nT = bc["A"], bc["C"], bc["G"], bc["T"]
    acgt = nA + nC + nG + nT
    out = {"accession": accession, "n_contigs": n_contigs, "seq_len": len(F),
           "acgt_len": acgt, "n_A": nA, "n_C": nC, "n_G": nG, "n_T": nT,
           "n_other": len(F) - acgt - F.count("N") + 0}
    out.update({"w_" + k: v for k, v in tally.items()})
    return out


if __name__ == "__main__":
    import time
    p = sys.argv[1]
    t0 = time.time()
    r = scan(p, "TIMING")
    print(json.dumps(r, indent=0))
    print("elapsed_s", round(time.time() - t0, 2))
