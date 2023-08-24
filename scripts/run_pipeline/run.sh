#!/bin/bash
##SampleID,dataset,nCount_RNA,nFeature_RNA must include in the metadata


data_file="path/to/seurat_object.rds"
cancer_type="breast cancer"
#quality_control="1000,10,10000,10,5,10,90,50,200"
use_dynamic_threshold=1
output_file="path/to/outputdir"
scefile='/sce.rds'



Rscript QC.R -d "$data_file" -t "$cancer_type" -q "$quality_control" --dynamic -o "$output_file"
Rscript exhcell.R -d "$output_file$scefile" -s "path/to/signature/file" -o "$output_file"

