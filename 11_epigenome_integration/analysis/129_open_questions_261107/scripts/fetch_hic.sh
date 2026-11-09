#!/bin/bash
set -u
OUT=$1; CAP=$2; NCH=4
mkdir -p $OUT; CH=$((CAP/NCH))
grab () {  # $1 key, $2 url
  for i in $(seq 0 $((NCH-1))); do
    A=$((i*CH)); B=$((A+CH-1))
    curl -s --retry 3 -r ${A}-${B} "$2" -o $OUT/$1.part${i} &
  done
}
grab SRR23618920_1 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR236/020/SRR23618920/SRR23618920_1.fastq.gz
grab SRR23618920_2 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR236/020/SRR23618920/SRR23618920_2.fastq.gz
grab SRR23618922_1 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR236/022/SRR23618922/SRR23618922_1.fastq.gz
grab SRR23618922_2 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR236/022/SRR23618922/SRR23618922_2.fastq.gz
wait
for k in SRR23618920_1 SRR23618920_2 SRR23618922_1 SRR23618922_2; do
  cat $OUT/$k.part* > $OUT/$k.fastq.gz && rm -f $OUT/$k.part*
done
ls -l $OUT
