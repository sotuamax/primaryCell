#!/bin/bash
# transform bam to bigwig file
# input: bam, bigwig output, threads
ml deeptools

bamCoverage -b $1 -o $2 -of bigwig -bs 10 -p $3 --normalizeUsing CPM --smoothLength 20
