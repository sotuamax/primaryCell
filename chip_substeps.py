#!/usr/bin/env python3
import argparse 
import subprocess
import os 
import logging

from align_tools import *
from bam_tools import bam2bed
from fastq_tools import read_count, read_len

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # master parameters 
    parser.add_argument("-ref", "--reference", required = False, help = "reference genome")
    parser.add_argument("-o", "--output", help="output directory name")
    parser.add_argument("-n", "--thread", help = "threads to use")
    # separate argumentparser for sub-parsers
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    # initiate sub-parser 
    sub_parsers = parser2.add_subparsers(dest = 'command', help = "mode to run")
    # alignment
    align = sub_parsers.add_parser("align", help = "perform alignment", parents = [parser], add_help = False)
    align.add_argument("-fastq", "--fastq", nargs = 2, required = True, help="fastq input file (for paired-end reads, order in R1 R2)")
    # transform to bigwig format
    bw = sub_parsers.add_parser("bigwig", help = "transform to bigwig format", parents = [parser], add_help = False)
    bw.add_argument("-bam", "--bam", required = True, help = "bam alignment file")
    # markdup
    markdup = sub_parsers.add_parser("markdup", help = "perform feature count", parents = [parser], add_help = False)
    markdup.add_argument("-bam", "--bam", required = True, help = 'bam alignment file')
    # parse arguments
    args=parser2.parse_args()
    return args

def main():
    args = args_parser()
    ref = args.reference
    # gtf = args.gtf 
    output = args.output
    thread = args.thread
    # perform fastq alignment
    if args.command == "align":
        logging.basicConfig(filename = os.path.basename(output)+".log", filemode = "w", format = '%(message)s', datefmt = None, level = logging.DEBUG)
        # bwa_index = args.index_dir
        r1 = args.fastq[0]; r2 = args.fastq[-1]
        read_c = read_count(r1)
        logging.info(f"read count\t{read_c}")
        read_l = read_len(r1)
        logging.info(f"read length\t{read_l}")
        if not os.path.exists(output + ".bam"):
            BWA_PE(r1, r2, ref, output, thread)
        BWA_log(output + ".bam")
    if args.command == "markdup":
        logging.basicConfig(filename = os.path.basename(output)+".log", filemode = "a", format = '%(message)s', datefmt = None, level = logging.DEBUG)
        bam = args.bam 
        if not os.path.exists(output + ".bam"):
            markdup(bam, output)
        mark_log_df = markdup_log(output + ".txt")
        read_pairs = mark_log_df.loc["READ_PAIRS_EXAMINED"].item()
        logging.info(f"Read pairs\t{read_pairs}")
        dup_rate = mark_log_df.loc["PERCENT_DUPLICATION"].item()
        logging.info(f"Duplicate rate\t{dup_rate}")
    if args.command == "bigwig":
        bam = args.bam
        if not os.path.exists(output + ".bw"):
            bw = output + ".bw"
            b2bw(bam, bw, thread)

if __name__ == "__main__":
    main()
