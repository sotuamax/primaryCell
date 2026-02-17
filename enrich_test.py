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

Updated 11/10/2025: enable multi-threading process for searching of overlaps for TF peaks; 

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
import os 
from tqdm.auto import tqdm
from joblib import Parallel, delayed
ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="TF_enrich.py -bed region.bed -query specific1.bed specific2.bed -target GATA3.bed RUNX3.bed -o <output_name>")
    parser.add_argument('-bed', "--bed", help = "bed input as background")
    parser.add_argument("-query", "--query", required = False, nargs = "+", help = "bed input for query region to test for enrichment")
    parser.add_argument("-target", "--target", nargs = "+", help = "bed file for target region of interest (e.g., TF motif)")
    parser.add_argument("-n", "--threads", type = int, default = 1, help = "number of threads to use.")
    parser.add_argument("-o", "--output", help = "output name (prefix)")
    args=parser.parse_args()
    return args

def overlap_fun(bed1:str, bed2:str, coln:str):
    # no need to get strand information for TF (only interested in its coordinates)
    bed1_df = bf.read_table(bed1, schema = "bed3").drop_duplicates(keep = "first", ignore_index = True)
    bed2_df = bf.read_table(bed2, schema = "bed3").drop_duplicates(keep = "first", ignore_index = True)
    bed_overlap = bf.overlap(bed1_df,bed2_df,how="left")[["chrom", "start","end", "chrom_"]].drop_duplicates(keep = "first", ignore_index = True)
    bed_overlap[coln] = np.where(pd.isna(bed_overlap["chrom_"]), 0, 1)
    bed_overlap.drop("chrom_", axis = 1, inplace = True)
    return bed_overlap.set_index(["chrom", "start", "end"])

def main():
    start = time.time()
    # get parameters 
    args = args_parser()
    bed = args.bed 
    target_list = args.target; 
    query_list = args.query
    n = args.threads
    output = args.output
    print("Read bed file")
    bed_df = bf.read_table(bed, schema = "bed3").drop_duplicates(keep = "first", ignore_index = True)
    print("Report bed overlap with target bed files ...")
    overlap_target = Parallel(n_jobs=n, prefer = "threads")(
                delayed(overlap_fun)(bed, target, os.path.basename(target).split(".bed")[0]) for target in tqdm(target_list)
            )
    overlap_target_df = overlap_target[0].join(overlap_target[1:]).reset_index()
    assert len(overlap_target_df) == len(bed_df), print("Overlap target is not match to bed input file")
    # to add query 
    if query_list is not None:
        pd.DataFrame([os.path.basename(q).split(".bed")[0] for q in query_list]).to_csv(f"{output}_query.txt", sep = "\t", header = False, index = False)
        pd.DataFrame([os.path.basename(t).split('.bed')[0] for t in target_list]).to_csv(f"{output}_target.txt", sep = "\t", header = False, index = False)
        print("Report bed overlap with query bed file ...")
        overlap_query = Parallel(n_jobs=n)( # query included in background 
            delayed(overlap_fun)(bed, query, os.path.basename(query).split(".bed")[0]) for query in query_list
        )
        overlap_query_df = overlap_query[0].join(overlap_query[1:]).reset_index()
        assert len(overlap_query_df) == len(bed_df), print("Overlap query is not match to bed input file")
        overlap_df = pd.merge(overlap_query_df, overlap_target_df, on=["chrom", "start", "end"])
    else:
        overlap_df = overlap_target_df
    overlap_df.to_csv(f"{output}_bed.txt", sep = "\t", header = True, index = False)
    if query_list is not None:
        from scipy.stats import fisher_exact
        # make contingency table 
        from itertools import product
        query_col = pd.read_table(f"{output}_query.txt", sep = "\t", header = None, names = ["query"])["query"].tolist()
        target_col = pd.read_table(f"{output}_target.txt", sep = "\t", header = None, names = ["target"])["target"].tolist()
        query_target_combination = list(product(query_col, target_col))
        fisher_res = Parallel(n_jobs = n)(
            delayed(fisher_exact)(pd.crosstab(overlap_df[qc[0]],overlap_df[qc[-1]]), alternative = "greater") for qc in query_target_combination
        )
        qc_list = list()
        for qc in query_target_combination:
            qc_freq = pd.crosstab(overlap_df[qc[0]],overlap_df[qc[-1]]).reset_index().melt(id_vars = qc[0], var_name = qc[1], value_name = "freq")
            qc0_total = qc_freq[qc_freq[qc[0]] == 1]["freq"].sum()
            qc1_total = qc_freq[qc_freq[qc[1]] == 1]["freq"].sum()
            overlap_total = qc_freq[(qc_freq[qc[0]] == 1) & (qc_freq[qc[1]] == 1)]["freq"].sum()
            qc_list.append((qc[0], qc[1],qc0_total, qc1_total, overlap_total))
        # add query/target enrichment results 
        fisher_res_df = pd.merge(pd.DataFrame(qc_list, columns = ["query", "target", "query_freq", "target_freq", "overlap"]), pd.DataFrame(fisher_res), left_index = True, right_index = True)
        fisher_res_df.to_csv(output + ".txt", sep = "\t", header = True, index = False)
        # fisher_res_df.sort_values(["statistic", "pvalue"], ascending = [False, True], inplace = True, ignore_index = True)
        print("Write significant enrichment ...")
        # fisher_res_df.to_csv(output + ".stat", sep = "\t", header = True, index = False)
        fisher_res_df.query("pvalue < 0.001 and statistic > 1").to_csv(output + "_sig.txt", sep = "\t", header = True, index = False)
    print(timeit(start))

if __name__ == "__main__":
    main()
