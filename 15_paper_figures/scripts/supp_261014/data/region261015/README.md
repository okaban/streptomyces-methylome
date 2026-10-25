# region261015 — 修飾部位の領域組成（4 区分 x 3 マーク x 3 時点）

    source  : /Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015/<lib>_pileup.bed
              （scripts/supp_261014/data/pileup_{GCCGGC,AAGCCCG_4mC,AAGCCCG_6mA}/ 経由）
    tss     : 11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv
              category が P または S かつ gene_id が空でない 3,030 件（本文記載と一致）
    genes   : 11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/data/
              GCF_000203835.1_ASM20383v1_genomic.gff の CDS 8,199 件（ゲノムの 89.1% を被覆）
    computed_on : 2026-10-15

## 区分規則（本文の記載どおり、優先順）

1. promoter   — **遺伝子体の外**かつ最近接 TSS の上流 1–500 bp（鎖を考慮）
2. 5'UTR      — TSS の 0 〜 +100 bp
3. CDS-internal
4. intergenic

旧版は「遺伝子体に関係なく promoter を優先」で、百分率は版間で比較できない。
補足図 18 のパネル a は旧版の規則で計算されているため、そちらの promoter 割合
（TGGCCGGC 8.4% など）とこの表（1.54%）は直接比較してはいけない。

## 検証

TGGCCGGC は本文記載 93.6 / 93.5 / 93.6% に対し 93.58 / 93.45 / 93.73% を再現（T3 は
正準パイルアップ統一で部位集合が 1,877 → 1,866 に変わったぶん移動）。
AAGCCCG 6mA も 91.4 / 91.1 に対し 91.35 / 91.11 と一致（T3 は 89.7 → 90.29 へ移動）。
**AAGCCCG 4mC のみ本文 87.3 / 87.8 / 87.5% に対し 86.61 / 87.10 / 86.95% と 3 時点とも
約 0.7 ポイント低い。** 部位数（1,023 / 977 / 973）は一致するので部位集合の違いではなく、
本文値が別世代の区分割当で算出されたものと判断し、再計算値を採る。
