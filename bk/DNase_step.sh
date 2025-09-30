#!/bin/bash
# input: read_prefix, 

raw_align=""
# step 1: 
bwa_align.py $1 human -outdir $raw_align -n 24 
# step 2: 
markdup.sh $raw_align/$1.bam $mark_prefix 
# step 3: 
bam_qc.sh $mark_prefix.bam $qc_bam 
# merge file should be clarified 

