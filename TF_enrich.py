#!/usr/bin/env python
"""
This script is to identify TF motifs that overlaps regions in a given bed file. 
Input: 
    required: 
    - bed: bed file for regions to look at; 
    - TF: one or multiple TF bed files (where TF motifs identified);
    optional:
    - query: a subset of bed file, for query regions used to test if TF enriched compared to all regions; 

Output: 
    a matrix-like file for which each row is bed region, and each column match to a TF where 0 for no-overlap, and 1 for overlap. 
"""
import pandas as pd 
import bioframe as bf
import numpy as np
# from scipy.stats import fisher_exact
import os 
import argparse 
import time 
from utilities.misc import timeit 
from utilities.misc import ignore_warning
ignore_warning()
import os 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="TF_enrich.py -query query.bed -TF TF1.bed TF2.bed TF3.bed -bg bg.bed -o bg_TF")
    parser.add_argument('bed', help = "bed input to search for ovelapping with TFs")
    parser.add_argument("-TF", "--TF", nargs = "+", help = "Occupancy of transcription factor in bed format")
    parser.add_argument("-query", "--query", required = False, nargs = "*", default = [], help = "query region in bed format")
    parser.add_argument("-o", "--output", help = "output name (prefix)")
    args=parser.parse_args()
    return args

def main():
    start = time.time()
    # get parameters 
    args = args_parser()
    query = args.query 
    tf_list = args.TF # TF is a list with 1 or > 1 TF bed files 
    bg = args.bed
    output = args.output
    bg_df = bf.read_table(bg, schema = "bed4")
    # add query matrix
    if len(query) > 0:
        for q in query:
            query_df = bf.read_table(q, schema = "bed5") # minimum bed5 
            # add background df with columns that w/wo overlapping query 
            query_bg = pd.merge(bg_df, query_df, on = ["chrom", "start", "end"], how = "inner")
            if len(query_bg) != len(query_df):
                print("Update background file ....")
                # when background not include all query, update background to include query 
                # bg_df.set_index(["chrom", "start", "end"]).index
                qindex = query_df.set_index(["chrom", "start", "end"]).index
                bg_add_df = query_df[~qindex.isin(bg_df.set_index(["chrom", "start", "end"]).index)]
                bg_df = pd.concat([bg_df, bg_add_df[bg_df.columns.tolist()]], axis = 0)
            query_df[os.path.basename(q)] = 1
            bg_df = pd.merge(bg_df, query_df[["chrom", "start", "end", os.path.basename(q)]], on = ["chrom", "start", "end"], how = "left")
            bg_df[os.path.basename(q)] = np.where(pd.isna(bg_df[os.path.basename(q)]), 0, 1)
    # add TF matrix
    for tf in tf_list:
        factor = os.path.basename(tf).split(".bed")[0]
        tf_df = bf.read_table(tf, schema = "bed5")
        # add background df with column that w/wo overlapping TF region 
        bg_df = bf.overlap(bg_df, tf_df, how = "left")
        bg_df[factor] = np.where(pd.isna(bg_df["name_"]), 0, 1)
        bg_df.drop([c for c in bg_df.columns if c.endswith("_")], axis = 1, inplace = True)
        bg_df.drop_duplicates(keep = "first", ignore_index = True, inplace = True)
        if len(bg_df) == len(bg_df):
            print(f"{factor} annotation passed!")
        #tab_freq = pd.crosstab(bg_df[f"query"], bg_df[factor])
        #tab_fisher = fisher_exact(tab_freq)
        #print(factor, ":", tab_fisher)
    bg_df.to_csv(output + ".txt", sep = "\t", header = True, index = False)
    print(timeit(start))

if __name__ == "__main__":
    main()
