"""Supplementary Table 28 — AAGCCCG occupancy ladder, both marks on the same instances.

生成元のなかった表 28 に生成元を与え、同時に 4mC 列（同一インスタンス上）を追加する。
本文 L105/L115 の分解（最深層で 6mA − 4mC = -3.53 / 0.00 / -5.88 pp）は、
この表の 2 列を引き算すれば読者が検算できる。

規則（正準）:
  - インスタンス = 参照配列上の AAGCCCG 出現（両鎖）。位置 1 = 先頭 A、位置 5 = 3 番目の C。
  - 6mA は位置 1 のコール、4mC は位置 5 のコールで数える。
  - コールの定義は正準部位一覧 canonical_sites_genomewide_261015.tsv に従う
    （3 ライブラリをプールし、真の深度 >= 30、Nmod/Nvalid >= 0.5、参照塩基一致）。
  - 分母 n = その時点で位置 1 が測定可能だったインスタンス数（既発表の梯子と同一）。

検証: 6mA 列は parent_reverify_ladder_261016.tsv を 3 時点 x 4 段階すべてで再現する。

実行: cd 15_paper_figures/scripts && python3 recompute_261015/build_supptable28_261021.py
"""
import pathlib
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
FASTA = pathlib.Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/data/NC_003888.3.fna")
SITES = HERE / "canonical_sites_genomewide_261015.tsv"
PREV = HERE / "parent_reverify_ladder_261016.tsv"
OUT = pathlib.Path("/Users/okaban/bioinfo/rna-seq/15_paper_figures/tables/supplementary/SuppTable28_occupancy_ladder_261021.tsv")
MOTIF = "AAGCCCG"
TPS = ("T1", "T2", "T3")
COMP = str.maketrans("ACGT", "TGCA")
C1 = {"A": "T", "C": "G", "G": "C", "T": "A"}


def genome():
    return "".join(l.strip() for l in open(FASTA) if not l.startswith(">")).upper()


def instances(seq):
    """(pos1, pos5, strand_is_plus, base-1, base-3, base-5) を両鎖から集める。"""
    out, rc = [], MOTIF.translate(COMP)[::-1]
    i = seq.find(MOTIF)
    while i >= 0:
        if i >= 5:
            out.append((i, i + 4, 1, seq[i - 1], seq[i - 3], seq[i - 5]))
        i = seq.find(MOTIF, i + 1)
    j = seq.find(rc)
    while j >= 0:
        if j + 11 < len(seq):
            out.append((j + 6, j + 2, 0, C1.get(seq[j + 7], "N"),
                        C1.get(seq[j + 9], "N"), C1.get(seq[j + 11], "N")))
        j = seq.find(rc, j + 1)
    return out


def main():
    seq = genome()
    inst = instances(seq)
    H = [x for x in inst if x[3] in "ACT"]
    H3 = [x for x in H if x[4] == "G"]
    H35 = [x for x in H3 if x[5] == "C"]
    tiers = [("AAGCCCG (all instances)", inst), ("HAAGCCCG (5' base not G)", H),
             ("+ G at -3", H3), ("+ C at -5 and G at -3", H35)]

    s = pd.read_csv(SITES, sep="\t")
    call = {(tp, m): set(zip(g.pos, g.strand))
            for (tp, m), g in s.groupby(["timepoint", "mod"])}
    prev = pd.read_csv(PREV, sep="\t")
    nmap = {(r.timepoint, r.definition): int(r.n) for r in prev.itertuples()}
    keymap = {"AAGCCCG (all instances)": "AAGCCCG (all)", "HAAGCCCG (5' base not G)": "HAAGCCCG",
              "+ G at -3": "+ -3=G", "+ C at -5 and G at -3": "+ -5=C & -3=G"}

    rows = []
    for tp in TPS:
        for label, xs in tiers:
            n = nmap[(tp, keymap[label])]
            c6 = sum(1 for x in xs if (x[0], "+" if x[2] else "-") in call[(tp, "6mA")])
            c4 = sum(1 for x in xs if (x[1], "+" if x[2] else "-") in call[(tp, "4mC")])
            rows.append(dict(timepoint=tp, tier=label, n_measurable=n,
                             called_6mA_pos1=c6, occupancy_6mA_pct=round(100 * c6 / n, 2),
                             called_4mC_pos5=c4, occupancy_4mC_pct=round(100 * c4 / n, 2),
                             difference_6mA_minus_4mC_pp=round(100 * (c6 - c4) / n, 2)))
    out = pd.DataFrame(rows)

    chk = prev[prev.definition != "(ref) same-instance 4mC C5"]
    for r in chk.itertuples():
        got = out[(out.timepoint == r.timepoint) &
                  (out.tier.map(keymap) == r.definition)].called_6mA_pos1.iloc[0]
        assert got == int(r.called), f"6mA mismatch {r.timepoint} {r.definition}: {got} vs {r.called}"
    ref = prev[prev.definition == "(ref) same-instance 4mC C5"]
    for r in ref.itertuples():
        got = out[(out.timepoint == r.timepoint) &
                  (out.tier == "AAGCCCG (all instances)")].called_4mC_pos5.iloc[0]
        assert got == int(r.called), f"4mC mismatch {r.timepoint}: {got} vs {r.called}"

    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
