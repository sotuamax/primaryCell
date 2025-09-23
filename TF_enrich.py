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

Note: 
    to save time, for bg overlap TF factor, only run one time and add other queries to save time. 

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
from tqdm.auto import tqdm

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="TF_enrich.py -query query.bed -TF TF1.bed TF2.bed TF3.bed -bg bg.bed -o bg_TF")
    parser.add_argument('-bg', "--bg", help = "bed input to search for ovelapping with TFs (as background)")
    parser.add_argument("-TF", "--TF", nargs = "+", help = "Occupancy of transcription factor in bed format")
    parser.add_argument("-query", "--query", nargs = "+", help = "bed input for query region to test for enrichment")
    parser.add_argument("-o", "--output", help = "output name (prefix)")
    args=parser.parse_args()
    return args

def main():
    start = time.time()
    # get parameters 
    args = args_parser()
    query = args.query 
    tf_list = args.TF # TF is a list with 1 or > 1 TF bed files 
    bg = args.bg
    output = args.output
    bg_df = bf.read_table(bg, schema = "bed4")
    bg_df = bg_df.query("chrom != 'chrY'").copy()
    # add query matrix
    print("Add query overlap ...")
    mat_dict = dict()
    for q in query:
        # query not necessarity use the same coordinates as bg 
        query_df = bf.read_table(q, schema = "bed5") # minimum bed5 
        # add background df with columns that w/wo overlapping query 
        bg_query = bf.overlap(bg_df, query_df, how = "left") 
        # query_bg = pd.merge(bg_df, query_df, on = ["chrom", "start", "end"], how = "inner")
        mat_dict[os.path.basename(q)] = np.where(pd.isna(bg_query["chrom_"]), 0, 1)
        bg_query.drop([c for c in bg_query.columns if c.endswith("_")], axis = 1, inplace = True)
    # add TF matrix
    print("Add TF overlap ...")
    for tf in tqdm(tf_list):
        try:
            factor = os.path.basename(tf).split(".bed")[0]
            tf_df = bf.read_table(tf, schema = "bed3")
            # add background df with column that w/wo overlapping TF region 
            bg_query = bf.overlap(bg_query, tf_df, how = "left")
            mat_dict[factor] = np.where(pd.isna(bg_query["chrom_"]), 0, 1)
            bg_query.drop([c for c in bg_query.columns if c.endswith("_")], axis = 1, inplace = True)
            # bg_query.drop_duplicates(keep = "first", ignore_index = True, inplace = True)
        except Exception as e:
            print(e)
    mat_df = pd.DataFrame.from_dict(mat_dict)
    bg_query_mat = pd.merge(bg_query, mat_df, left_index = True, right_index = True)
    bg_query_mat.to_csv(output + ".txt", sep = "\t", header = True, index = False)
    print(timeit(start))

if __name__ == "__main__":
    main()
