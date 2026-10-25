#!/bin/bash
# 正準パイルアップ（pileup_canonical_261015）で全時点・全マークを数え直す。
set -u
D=/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015
O=/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/110_topup_recompute/tables
census () {
  local code=$1 out=$2; shift 2
  awk -v C="$code" 'BEGIN{OFS="\t"}
    $4==C {k=$2"_"$6; nv[k]+=$10; nm[k]+=$12; nf[k]+=$16}
    END{for(k in nv) if((nv[k]+nf[k])>=30 && nv[k]>0 && nm[k]/nv[k]>=0.5){
          split(k,a,"_"); print a[1], a[2], nm[k]/nv[k]*100 }}' "$@" > "$out"
}
echo "[$(date +%H:%M:%S)] canonical recompute start"
for code in 21839 a m; do
  tag=$([ $code = 21839 ] && echo 4mC || { [ $code = a ] && echo 6mA || echo 5mC; })
  census $code $O/canon_${tag}_T1.tsv $D/1-1_pileup.bed $D/1-2_pileup.bed $D/1-3_pileup.bed &
  census $code $O/canon_${tag}_T2.tsv $D/2-1_pileup.bed $D/2-3_pileup.bed $D/2-4_pileup.bed &
  census $code $O/canon_${tag}_T3.tsv $D/3-2_pileup.bed $D/3-3_pileup.bed $D/3-4_pileup.bed &
  wait
  echo "[$(date +%H:%M:%S)] $tag done"
done
wc -l $O/canon_*.tsv
echo "[$(date +%H:%M:%S)] finished"
