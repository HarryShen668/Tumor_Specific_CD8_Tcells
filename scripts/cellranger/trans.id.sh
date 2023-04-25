#!/usr/bin/bash

if [ $# -eq 0 ]; then
  echo "please give srr2sample.txt file"
  exit 1
fi

input=$1
while read key value
do
  mv ${key}_[12].fastq.gz ${value}_S1_L001_R[12]_001.fastq.gz
done < "$input"

