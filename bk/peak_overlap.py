#!/usr/bin/env python3
import sys 
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

import pandas as pd 
import bioframe as bf 
import os
import numpy as np 
import argparse 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("peak1", help="peak 1 in bed format")
    parser.add_argument("peak2", help = "peak 2 in bed format")
    parser.add_argument("-o", "--output", help = "output name")
    args=parser.parse_args()
    return args

def overlap(peak1_df, peak2_df):
    peak_overlap = bf.overlap(peak1_df, peak2_df, how = "left")
    peak_overlap["overlap"] = np.where(pd.isna(peak_overlap["chrom_"]), 0, 1)
    peak_overlap.drop(["chrom_", "start_", "end_", "name_"], axis = 1, inplace = True)
    return peak_overlap

def main():
    args = args_parser()
    peak1 = args.peak1
    peak2 = args.peak2
    output = args.output
    peak1_df = bf.read_table(peak1, schema = "bed4")
    peak2_df = bf.read_table(peak2, schema = "bed4")
    # 
    peak1_overlap = overlap(peak1_df, peak2_df)
    peak2_overlap = overlap(peak2_df, peak2_df)
    # 
    peak1_overlap.drop_duplicates(keep = "first", inplace = True)
    peak2_overlap.drop_duplicates(keep = "first", inplace = True)
    # 
    peak1_overlap.to_csv(output + "1.bed", header = False, index = False, sep = "\t")
    peak2_overlap.to_csv(output + "2.bed", header = False, index = False, sep = "\t")

if __name__ == "__main__":
    main()

