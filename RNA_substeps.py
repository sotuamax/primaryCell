#!/usr/bin/env python3
import argparse 
import subprocess
import os 

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # important parameters 
    # parser.add_argument("-ref", "--reference", required = False, help = "reference genome")
    parser.add_argument("-gtf", "--gtf", help = "reference gtf file used for alignment")
    parser.add_argument("-o", "--output", help="output name")
    parser.add_argument("-n", "--thread", help = "threads to use")
    parser.add_argument("-index_dir", "--index_dir", required = False, help = "alignment index file")
    # separate argumentparser for sub-parsers
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    # initiate sub-parser 
    sub_parsers = parser2.add_subparsers(dest = 'command', help = "mode to run")
    # alignment
    align = sub_parsers.add_parser("align", help = "perform alignment", parents = [parser], add_help = False)
    align.add_argument("-fastq", "--fastq", nargs = 2, required = True, help="fastq input file (for paired-end reads, order in R1 R2)")
    # add step for bigwig transformation
    bw = sub_parsers.add_parser("bigwig", help = "transform bam to bigwig file", parents = [parser], add_help = False)
    bw.add_argument("-bam", "--bam", required = True, help = 'bam alignment file')
    # count
    count = sub_parsers.add_parser("count", help = "perform feature count", parents = [parser], add_help = False)
    count.add_argument("-bam", "--bam", required = True, help = 'bam alignment file')
    # RSEM
    quantify = sub_parsers.add_parser("quantify", help = "quantify feature in a normalized manner", parents = [parser], add_help = False)
    quantify.add_argument("-fastq", "--fastq", nargs = 2, required = True, help="fastq input file (for paired-end reads, order in R1 R2)")
    # parse arguments
    args=parser2.parse_args()
    return args

def main():
    args = args_parser()
    # ref = args.reference
    gtf = args.gtf 
    output = args.output
    thread = args.thread
    if args.command == "align":
        from utilities.align_tools import STAR_PE
        star_index = args.index_dir
        r1 = args.fastq[0]; r2 = args.fastq[-1]
        STAR_PE(r1, r2, star_index, gtf, output, thread)
    if args.command == "count":
        bam = args.bam 
        from utilities.align_tools import featurecount, parse_featurecount
        featurecount(bam, gtf, output, thread)
        output_clean = parse_featurecount(output)
        output_clean.to_csv(output + ".featurecount", sep = "\t", header = True, index = False)
    if args.command == "quantify":
        from utilities.align_tools import rsem_quantify, parse_rsem
        r1 = args.fastq[0]; r2 = args.fastq[-1]
        rsem_index = args.index_dir
        if not os.path.exists(output + ".genes.results") or not os.path.exists(output + ".isoforms.results"):
            rsem_quantify(r1, r2, rsem_index, output, thread)
        output_clean = parse_rsem(output)
        output_clean.to_csv(output + ".rsem", sep = "\t", header = True, index = False)
    if args.command == "bigwig":
        bam = args.bam 
        bw = output + ".bw"
        from utilities.bam_tools import bam2bw 
        if not os.path.exists(bw):
            bam2bw(bam, bw, thread)

if __name__ == "__main__":
    main()
    