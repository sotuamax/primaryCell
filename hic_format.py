#!/usr/bin/env python3
"""
Transform BAM input file into pairs/cool/hic.

Input: 
    required: 
    - bam: bam file (recommend in name sorted)
    - ouput: output name 
    optional:
    - assembly: genome assembly name 
    - threads: threads to use 
    - chrsize: chromosome size file 

===============
Example: 
source myoncda 
mamba activate bio 
ml pairtools juicer
hic_format.py -n 12 -o <sample> <sample>.bam

"""
# import cooler 
#import numpy as np 
import subprocess 
import argparse 
#import pandas as pd 
#import logging 
import os 
import time 
from utilities.misc import timeit

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = True, usage="hic_format.py -o <out> input.bam ")
    parser.add_argument("bam", help = "input bam file (either coordinate or name sorted, provide name sorted file to save time)")
    parser.add_argument("-o", "--output", required = True, help = "output file name (prefix)")
    parser.add_argument("-n", "--thread", default = 2, type = int, help = "thread to process bam")
    parser.add_argument("-chrsize", "--chrsize", required = False, default = "/data/jim4/Reference/human/GRCh38.p14/GRCh38_chrsize.bed", help = "file for ordered chromosome and size")
    parser.add_argument("-assembly", "--assembly", default = "hg38", help = "genome assembly name")
    args = parser.parse_args()
    return args

def main():
    start = time.time()
    args = args_parser()
    bam = args.bam
    from utilities.bam_tools import examine_sort
    try:
        b_sort = examine_sort(bam)
        print("BAM is sorted by", b_sort)
    except Exception as e: 
        print(e)
        exit(1)
    n = args.thread
    if b_sort == "coordinate":
        print("Sort BAM by name ...")
        # when bam is not sorted, sort it and store in default directory 
        sort_dir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/sortn"
        new_bam = os.path.join(sort_dir, os.path.basename(bam))
        command = f"samtools sort -@ {n} -n -o {new_bam} {bam}"
        print(command)
        subprocess.call(command, shell=True)
        bam = new_bam
    # output name 
    out = args.output
    chrsize = args.chrsize 
    assembly = args.assembly
    # generate pairs 
    if not os.path.exists(out + ".pairs.gz"):
        # use standard pairtools to generate pairs file 
        print("Generate compressed pairs ...")
        # BAM must be name sorted 
        print(f"Run pairtools w/ {n} threads ...")
        pair_command = f"pairtools parse {bam} -c {chrsize} --assembly {assembly} --nproc-in {n} --nproc-out {n} --drop-sam | pairtools sort --nproc {n} -o {out}.pairs.gz --memory 20g"
        # drop-sam: do not add sams to the output 
        print(pair_command)
        subprocess.call(pair_command, shell = True)
    # if not os.path.exists(out + ".pairs.gz"):
    #     subprocess.call(f"pairtools sort {out}.pairs.gz -o {out}_sort.pairs.gz --nproc {n} --memory 40g", shell = True)
    # generate cool
    if not os.path.exists(out + ".cool"):
        print("Load pairs to generate cool at 1k ...")
        # https://cooler.readthedocs.io/en/latest/cli.html#cooler-cload-pairs
        # pairs do NOT be sorted. Accept compressed file. 
        cool_command = f"cooler cload pairs --assembly {assembly} {chrsize}:1000 {out}.pairs.gz {out}.cool -c1 2 -p1 3 -c2 4 -p2 5"
        print(cool_command)
        subprocess.call(f"cooler cload pairs --assembly {assembly} {chrsize}:1000 {out}.pairs.gz {out}.cool -c1 2 -p1 3 -c2 4 -p2 5", shell = True)
    if not os.path.exists(out + ".mcool"):
        print("Generate mcool at 1,5,10,25,50,100 kb resolutions ...")
        mcool_command = f"cooler zoomify {out}.cool -n {n} -r 1000,5000,10000,25000,40000,50000,100000 -o {out}.mcool"
        print(mcool_command)
        print(f"Run cooler w/ {n} threads ...")
        subprocess.call(mcool_command, shell = True)
    # generate hic
    if not os.path.exists(out + ".hic"):
        # although cool no need to sort pairs, juicer must sort pairs by chromosome 
        print("Load pairs to generate HiC at 1,5,10,25,50,100 kb resolutions ... ")
        juicertools="/usr/local/apps/juicer/juicer-1.6/scripts/juicer_tools.jar"
        juicer_command = f"java -Xmx48g -jar {juicertools} pre -r 1000,5000,10000,25000,40000,50000,100000 -k KR {out}.pairs.gz {out}.hic --threads {n} {assembly}"
        print(juicer_command)
        print(f"Run juicer w/ {n} threads ...")
        subprocess.call(juicer_command, shell = True)
    print(timeit(start))

if __name__ == "__main__":
    main()
