#!/bin/sh
REP="$1"; OUT="$2"; P="$3"
awk -v KF="$OUT/data/keep_positions.txt" '
BEGIN{ while((getline p < KF)>0) K[p]=1 }
($4=="21839"||$4=="a") && ($2 in K) { print $2"\t"$6"\t"$4"\t"$10"\t"$12"\t"$16 }
' "$P/${REP}_pileup.bed" > "$OUT/data/${REP}_motif.tsv"
