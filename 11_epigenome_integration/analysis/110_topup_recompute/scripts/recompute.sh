#!/bin/bash
# 2026-02 追加フローセル分を含むパイルアップで正準値を数え直す。
# 対象は T3 のみ（追加ランがあるのは 3-2 と 3-4 だけ）。T1/T2 は対照として同条件で再計算する。
# 原稿・レジストリには触れない。差分表を出すところまで。
set -u
P=/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup
C=/Users/okaban/bioinfo/methyl/260102_M145/analysis/cursor_results/20260108/pileup
OUT=/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/110_topup_recompute

# 正準規則: 時点ごとに 3 ライブラリをプール / プール後 (Nvalid+Nfail) >= 30 / Nmod/Nvalid >= 0.5
census () {  # $1 = modcode, $2 = out, $3.. = pileups
  local code=$1 out=$2; shift 2
  awk -v C="$code" 'BEGIN{OFS="\t"}
    $4==C {k=$2"_"$6; nv[k]+=$10; nm[k]+=$12; nf[k]+=$16}
    END{for(k in nv) if((nv[k]+nf[k])>=30 && nv[k]>0 && nm[k]/nv[k]>=0.5){
          split(k,a,"_"); print a[1], a[2], nm[k]/nv[k]*100 }}' "$@" > "$out"
}

echo "[$(date +%H:%M:%S)] start"
# --- 4mC (code 21839)
census 21839 $OUT/tables/raw_4mC_T1.tsv  $P/1-1_pileup_restored_261009.bed $P/1-2_pileup.bed $P/1-3_pileup.bed &
census 21839 $OUT/tables/raw_4mC_T2.tsv  $P/2-1_pileup.bed $P/2-3_pileup.bed $P/2-4_pileup.bed &
census 21839 $OUT/tables/raw_4mC_T3old.tsv $P/3-2_pileup.bed $P/3-3_pileup.bed $P/3-4_pileup.bed &
wait
echo "[$(date +%H:%M:%S)] 4mC T1/T2/T3old done"
census 21839 $OUT/tables/raw_4mC_T3new.tsv $C/3-2_pileup.bed $P/3-3_pileup.bed $C/3-4_pileup.bed
echo "[$(date +%H:%M:%S)] 4mC T3new done"
# --- 6mA (code a)
census a $OUT/tables/raw_6mA_T1.tsv  $P/1-1_pileup_restored_261009.bed $P/1-2_pileup.bed $P/1-3_pileup.bed &
census a $OUT/tables/raw_6mA_T2.tsv  $P/2-1_pileup.bed $P/2-3_pileup.bed $P/2-4_pileup.bed &
census a $OUT/tables/raw_6mA_T3old.tsv $P/3-2_pileup.bed $P/3-3_pileup.bed $P/3-4_pileup.bed &
wait
census a $OUT/tables/raw_6mA_T3new.tsv $C/3-2_pileup.bed $P/3-3_pileup.bed $C/3-4_pileup.bed
echo "[$(date +%H:%M:%S)] 6mA all done"
# --- 5mC (code m) : T3 のみ、0 件が保たれるかの確認
census m $OUT/tables/raw_5mC_T3new.tsv $C/3-2_pileup.bed $P/3-3_pileup.bed $C/3-4_pileup.bed
echo "[$(date +%H:%M:%S)] 5mC T3new done"
wc -l $OUT/tables/raw_*.tsv
echo "[$(date +%H:%M:%S)] finished"
