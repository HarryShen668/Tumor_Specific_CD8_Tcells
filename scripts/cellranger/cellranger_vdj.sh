#!/usr/bin/bash

# check input parameters
if [ $# -ne 2 ]; then
  echo "Usage: $0 sample.txt fq_dir"
  exit 1
fi

input=$1  ##sample name file
fq_dir=$2
while read id
do
/home/liuzipei/bin/cellranger vdj --id=${id} \
--fastqs=$fq_dir \
--reference=/picb/lilab5/liuzipei/cellranger_ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
--localcores=30
--localmem=100
--sample=${id} > ${id}.log 
done < "$input"
