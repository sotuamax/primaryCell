#!/usr/bin/env python3
import pandas as pd 
import argparse 
import gzip 
from utilities.misc import timeit
import time 
import numpy as np 
from utilities.misc import ignore_warning
ignore_warning()
import os 

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = True, 
                                     usage="hic_format.py -o <out> inputfile -n 10")
    parser.add_argument("-pairs", "--pairs", nargs = "+", help = "input pairs file (support multiple pairs)")
    parser.add_argument("-o", "--output", required = True, help = "output bed name")
    parser.add_argument("-cis", "--cis", action = "store_true", help = "only intra-chromosomal PETs retained")
    #parser.add_argument("-r1", "--r1", action = "store_true", help = "only R1 (that is in forward orientation) retained")
    parser.add_argument("-rlen", "--rlen", default = 50, type = int, help = "read length used to extend position")
    args = parser.parse_args()
    return args

def main():
    start = time.time()
    args = args_parser()
    pairs_files = args.pairs 
    bed_output = args.output if args.output.endswith(".bed") else args.output + ".bed"
    readlen = args.rlen
    print("Write bed file ...")
    if os.path.exists(bed_output):
        print(bed_output, " exists!")
        exit(0)
    try:
        os.remove(bed_output)
    except:
        pass 
    #with gzip.open(bed_output + ".gz", "wt") as fw:
    for p in pairs_files:
        print(f"Read {p} ...")
        for chunk in pd.read_table(p, sep = "\t", comment="#", chunksize=5_000_000, header = None, names = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2", "pair_type"]):
            chunk = chunk.query("pair_type == 'UU'")
            if args.cis:
                chunk = chunk.query("chrom1 == chrom2")
            chunk["pos1"] = np.where(chunk["strand1"] == "-", chunk["pos1"]-readlen, chunk["pos1"])
            chunk["pos1_1"] = chunk["pos1"] + readlen
            chunk["pos2"] = np.where(chunk["strand2"] == "-", chunk["pos2"]-readlen, chunk["pos2"])
            chunk["pos2_1"] = chunk["pos2"] + readlen
            chunk[["chrom1", "pos1", "pos1_1","strand_1"]].to_csv(bed_output, sep = "\t", mode = "a", header = False, index = False)
            chunk[["chrom2", "pos2", "pos2_1", "strand_2"]].to_csv(bed_output, sep = "\t", mode = "a", header = False, index = False)
            # fw.write("\n".join([f"{row.chrom1}\t{row.pos1}\t{row.pos1_1}\n{row.chrom2}\t{row.pos2}\t{row.pos2_1}" for row in chunk.itertuples()]))
            #fw.write("\n".join([f"{row.chrom1}\t{row.pos1-readlen}\t{row.pos1}\n" for row in chunk.query("strand1 == '-'")[['chrom1', "pos1"]].itertuples()]))
        print(f"Finished Writing {p}!")
    print(timeit(start))

if __name__ == "__main__":
    main()
