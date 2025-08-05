#!/bin/bash
# input: bam, outdir
ml macs/3

name=$(basename $1)
name=${name%.bam}

macs3 callpeak --keep-dup all -t $1 -f BAMPE --outdir $2 -q 0.01 --nomodel --shift -75 --extsize 150 -n $name

