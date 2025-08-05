#!/bin/bash
# input: bam, outdir, threads

ml samtools 
ml bedtools 

name=$(basename $1)
name=${name%.bam}

samtools sort -@ $3 -n $1 | bedtools bamtobed -i stdin -bedpe > $2/$name.bed

