#!/bin/bash
# input: bedpe, outdir, threads

source myconda 
mamba activate bio 

name=$(basename $1)
name=${name%.bed}
echo $1
fseq2 callpeak $1 -f 0 -l 600 -pe -v -o $2 -t 4.0 -name $name -cpus $3 -standard_narrowpeak
