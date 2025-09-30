#!/bin/bash
# input: bam, outbam_prefix, # thread 

ml picard 
ml java 

picard="/usr/local/apps/picard/3.1.0/picard.jar"

java -jar $picard MarkDuplicates -I $1 -O $2.bam -M $2.txt --REMOVE_DUPLICATES true # && samtools index -@ $3 $2.bam 

