#!/usr/bin/env python3

import bioframe as bf 
import argparse 
import sys
import numpy as np 

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("bed1", help = "input bed1")
    parser.add_argument("bed2", help = "input bed2")
    parser.add_argument("-o", "--output", required = False, help = "output overlap regions")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    bed1 = args.bed1; bed2 = args.bed2 
    bed_df1 = bf.read_table(bed1, schema = "bed5", dtype = {"name":str})
    bed_df2 = bf.read_table(bed2, schema = "bed5", dtype = {"name":str})
    bed12_over = bf.overlap(bed_df1, bed_df2, how = "inner", return_overlap = True, suffixes = ("_1", "_2"))
    over1_unique = round(len(np.unique(bed12_over["name_1"]))/len(bed_df1)*100, 2)
    over2_unique = round(len(np.unique(bed12_over["name_2"]))/len(bed_df2)*100, 2)
    print("overlap1, ", over1_unique, "; overlap2, ", over2_unique)
    if args.output != None:
        # bed12_over.to_csv(args.output + "_overlap.txt", sep = "\t", header = True, index = False)
        # bed1_unique = bed_df1[~bed_df1["name"].isin(np.unique(bed12_over["name_1"]))]
        bed_df1["label"] = np.where(bed_df1["name"].isin(np.unique(bed12_over["name_1"])), 1, 0)
        # bed2_unique = bed_df2[~bed_df2["name"].isin(np.unique(bed12_over["name_2"]))]
        bed_df2["label"] = np.where(bed_df2["name"].isin(np.unique(bed12_over["name_2"])), 1, 0)
        bed_df1.to_csv(args.output + "_1.bed", sep = "\t", header = True, index = False)
        bed_df2.to_csv(args.output + "_2.bed", sep = "\t", header = True, index = False)

if __name__ == "__main__":
    main()