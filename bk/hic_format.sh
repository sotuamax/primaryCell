#!/bin/bash
# input: bam, outdir, thread

ml samtools 
ml bedtools 
ml java 
ml juicebox 

source myconda 
mamba activate bio 

name=$(basename $1)
name=${name%.bam}
out=$2/$name
chrsize="/data/jim4/Seq/primary_cell_project/data/GRCh38_chrsize.txt"
juicertools="/usr/local/apps/juicer/juicer-1.6/scripts/juicer_tools.jar"

# samtools sort -@ $3 -n $1 | bedtools bamtobed -i stdin -bedpe > $out.bedpe
# bed2pairs.py -bed $out.bedpe -chrsize $chrsize -assembly hg38 -o $out 
# juicer_tools pre [options] <infile> <outfile> <genomeID>
java -Xmx10000m -Djava.awt.headless=true -jar $juicertools pre $out.pairs.gz $out.hic --threads $3 hg38 -j $3/2 -k KR,VC,ICE
