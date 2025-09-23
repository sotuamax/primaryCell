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
source myconda 
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
    parser.add_argument("input", help = "input file (either BAM file coordinate or name sorted, provide name sorted file to save time; or pairs file sorted and compressed for UU reads)")
    parser.add_argument("-of", "--output_format", help = "output file format", choices = ["pairs", "cool", "hic"])
    parser.add_argument("-o", "--output", required = True, help = "output file name (prefix)")
    parser.add_argument("-n", "--thread", default = 2, type = int, help = "thread to process bam")
    parser.add_argument("-chrsize", "--chrsize", required = False, default = "/data/jim4/Reference/human/GRCh38.p14/GRCh38_chrsize.bed", help = "file for ordered chromosome and size")
    parser.add_argument("-assembly", "--assembly", default = "hg38", help = "genome assembly name")
    args = parser.parse_args()
    return args

def main():
    start = time.time()
    args = args_parser()
    input_file = args.input
    n = args.thread
    out = args.output
    chrsize = args.chrsize 
    assembly = args.assembly
    if args.output_format == "pairs":
        if not input_file.endswith(".bam"):
            raise ValueError("Input file is not in BAM format!")
        bam = input_file 
        from utilities.bam_tools import examine_sort
        try:
            b_sort = examine_sort(bam)
            print("BAM is sorted by", b_sort)
        except Exception as e: 
            print(e)
            exit(1)
        if b_sort == "coordinate":
            pair_command = f"samtools sort -@ {n} -n {bam} | pairtools parse -c {chrsize} --assembly {assembly} --nproc-in {n} --nproc-out {n} --drop-sam | pairtools select --nproc-in {n} --nproc-out {n} '(pair_type==\"UU\")' | pairtools sort --nproc-in {n} --nproc-out {n} -o {out}.pairs.gz"
        if b_sort == "name":
            pair_command = f"pairtools parse {bam} -c {chrsize} --assembly {assembly} --nproc-in {n} --nproc-out {n} --drop-sam | pairtools select --nproc-in {n} --nproc-out {n} '(pair_type==\"UU\")' | pairtools sort --nproc-in {n} --nproc-out {n} -o {out}.pairs.gz"
        if not os.path.exists(out + ".pairs.gz"):
            # use standard pairtools to generate pairs file 
            print("Generate pairs ...")
            print(pair_command)
            subprocess.call(pair_command, shell = True)
    if args.output_format == "cool":
        if not input_file.endswith(".pairs.gz"):
            raise ValueError("Input file is not in pairs format!")
        if not os.path.exists(out + ".cool"):
            print("Load pairs to generate cool at 1k ...")
            # https://cooler.readthedocs.io/en/latest/cli.html#cooler-cload-pairs
            # pairs do NOT to be sorted. Accept compressed file. 
            cool_command = f"cooler cload pairs --assembly {assembly} {chrsize}:1000 {out}.pairs.gz {out}.cool -c1 2 -p1 3 -c2 4 -p2 5"
            print(cool_command)
            subprocess.call(cool_command, shell = True)
        if not os.path.exists(out + ".mcool"):
            print("Generate mcool at 1,5,10,25,50,100 kb resolutions ...")
            mcool_command = f"cooler zoomify {out}.cool -n {n} -r 1000,5000,10000,25000,40000,50000 -o {out}.mcool"
            print(mcool_command)
            print(f"Run cooler w/ {n} threads ...")
            subprocess.call(mcool_command, shell = True)
    # generate hic
    if not os.path.exists(out + ".hic"):
        # although cool no need to sort pairs, juicer must sort pairs by chromosome 
        print("Load pairs to generate HiC at 1,5,10,25,50,100 kb resolutions ... ")
        juicertools="/usr/local/apps/juicer/juicer-1.6/scripts/juicer_tools.jar"
        juicer_command = f"java -Xmx48g -jar {juicertools} pre -r 1000,5000,10000,25000,40000,50000 -k KR {out}.pairs.gz {out}.hic --threads {n} {assembly}"
        print(juicer_command)
        print(f"Run juicer w/ {n} threads ...")
        subprocess.call(juicer_command, shell = True)
    print(timeit(start))

if __name__ == "__main__":
    main()
