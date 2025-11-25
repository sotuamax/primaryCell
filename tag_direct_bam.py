#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import argparse 
import pysam 
import os 
import subprocess

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # important parameters 
    parser.add_argument("bam", help="bam input")
    parser.add_argument("-flag", "--flag", help = "library barcode and its flag gene in a table, used to split BAM based on its flag gene (no header).")
    parser.add_argument("-n", "--threads", default = 1, type = int, help="threads to use to parse BAM")
    parser.add_argument("-split", "--split", action = "store_true", help = "Split BAM based on @RG")
    parser.add_argument("-bw", "--bw", action = "store_true", help = "generate matched Bigwig file on the splitted BAM file")
    parser.add_argument("-o", "--output", help = "output bam file name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    bam_file = args.bam
    flag = args.flag 
    n = args.threads 
    if not os.path.exists(bam_file + ".bai"):
        index_command = f"samtools index -@ {args.threads} {bam_file}"
        subprocess.call(index_command, shell = True)
    flag_df = pd.read_table(flag, sep = "\t", header = None, names = ["barcode", "label"])
    # flag_df.columns = [c.strip("_") for c in flag_df.columns]
    bam_handle = pysam.AlignmentFile(bam_file, "rb", threads = n)
    barcode_gene_dict = {row.barcode:row.label for row in flag_df.itertuples()}
    if args.output is None:
        new_bam = bam_file.replace(".bam", ".RG.bam")
    else:
        new_bam = args.output
    if not new_bam.endswith(".bam"):
        new_bam = new_bam + ".bam"
    header = str(bam_handle.header)
    # remove old RG if necessary 
    header = "\n".join([line for line in header.split("\n") if not line.startswith("@RG")])
    print("Add @RG to header ...")
    for gene in sorted(set(flag_df["label"])):
        header += f"@RG\tID:{gene}\n"
    print("Write new BAM with @RG ...")
    with pysam.AlignmentFile(new_bam, "wb", text = header, threads = n) as newbam:
        for read in bam_handle.fetch():
            if read.reference_name.startswith("chr") and read.reference_name not in ["chrY", "chrM"]:
                try:
                    read.set_tag("RG", barcode_gene_dict[read.query_name.split(":")[-1]], value_type = "Z")
                    newbam.write(read)
                except:
                    pass
    subprocess.call(f"samtools index -@ {n} {new_bam}", shell = True)
    # split bam file based on its RG
    if args.split:
        print("Generate split BAM ...")
        subprocess.call(f"samtools split {new_bam} -f '%*_%!.%.' -@ {n}", shell = True)
        if args.bw:
            from utilities.bam_tools import bam2bw
            print("Generate bigwig ...")
            all_flag = sorted(set(flag_df["label"]))
            for flag in all_flag:
                if os.path.exists(f"{args.output}_{flag}.bam"):
                    n_bam = f"{args.output}_{flag}.bam"
                    subprocess.call(f"samtools index -@ {n} {n_bam}", shell = True)
                    bam2bw(n_bam, n_bam.replace(".bam", ".bw"), n, smooth=100) # by default, no normalization 

if __name__ == "__main__":
    main()
