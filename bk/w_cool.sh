#!/bin/bash
# input: pairs, Bin size (bp), processor

source myconda 
mamba activate bio 

name=${1%.pairs.gz}
chrsize="/data/jim4/Seq/primary_cell_project/data/GRCh38_chrsize.txt"
cooler cload pairix $chrsize:$2 $1 $name.cool --assembly hg38 -p $3

