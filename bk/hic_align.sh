#!/bin/bash
# align hic paired-read to genome 
# input: read_prefix, reference, outdir, threads
ml bwa 
ml samtools 
r1=${1}_R1.fastq.gz
r2=${1}_R2.fastq.gz

# $2 is reference 
# $3 is thread 
name=$(basename $1)
name=$3/$name
bwa mem -5SP $2 $r1 $r2 -t $4 | samtools view -@ $4 -Su - | samtools sort -@ $4 - -o $name.bam && samtools index -@ $4 $name.bam
