#!/usr/bin/env python3
import os 
import subprocess 
import argparse
import pandas as pd 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-dir", "--directory", help = "directory for sequencing data")
    parser.add_argument("-tsv", "--tsv", help = "GC and name for softlink in a tab-seprated file (no header).")
    parser.add_argument("-outdir", "--outdir", help = "the directory for softlink file using new name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    tsv = args.tsv 
    df = pd.read_table(tsv, sep = "\t", header = None, names = ["GC", "name"])
    GC_name = {row.GC:row.name for row in df.itertuples()}
    for GC in GC_name:
        r1 = [os.path.join(args.directory, f) for f in os.listdir(args.directory) if f.startswith(GC + "_") and "R1" in f][0]
        r2 = [os.path.join(args.directory, f) for f in os.listdir(args.directory) if f.startswith(GC + "_") and "R2" in f][0]
        # r1 = os.path.abspath(os.path.join(args.directory, GC+"_R1.fastq.gz"))
        # r2 = os.path.abspath(os.path.join(args.directory, GC+"_R2.fastq.gz"))
        r1_ln = os.path.abspath(os.path.join(args.outdir, GC_name[GC] + "_R1.fastq.gz"))
        r2_ln = os.path.abspath(os.path.join(args.outdir, GC_name[GC] + "_R2.fastq.gz"))
        if os.path.exists(r1) and os.path.exists(r2):
            rename_command = f"ln -s {r1} {r1_ln} & ln -s {r2} {r2_ln}"
            print(rename_command)
            subprocess.call(rename_command, shell = True)
        else:
            print(r1 + " not exist in the folder")
            print(r2 + " not exist in the folder")

if __name__ == "__main__":
    main()
