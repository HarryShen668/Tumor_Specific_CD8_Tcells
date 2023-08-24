#!/bin/bash
tumorTCRfile='/home/shenhaoyu/lilab5/xinda/dataset/Tcr.select.raw/NPC/tumor_type.rds'
referenceTCRTCRfile='/home/shenhaoyu/lilab5/xinda/dataset/Tcr.select.raw/NPC/blood_type.rds'
clonetype_col="cdr3"
output="/home/shenhaoyu/lilab5/xinda/dataset/Tcr.select.raw/NPC"
matchfile='/tcr.matched.rds'
Rscript ClonalExpansiontest.R  -t "$tumorTCRfile" -r "$referenceTCRTCRfile" -c "$clonetype_col" -o "$output"
Rscript TCR.bionomalTest.R -t "$output$matchfile" -o "$output"




