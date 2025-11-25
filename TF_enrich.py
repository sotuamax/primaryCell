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


def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="TF_enrich.py -bg query.bed -TF TF1.bed TF2.bed TF3.bed -bg bg.bed -o bg_TF")
    parser.add_argument('-bed', "--bed", help = "bed input as background")
    parser.add_argument("-query", "--query", required = False, nargs = "+", help = "bed input for query region to test for enrichment")
    parser.add_argument("-target", "--target", nargs = "+", help = "bed file for target region of interest (e.g., TF motif)")
    parser.add_argument("-n", "--threads", type = int, default = 1, help = "number of threads to use.")
    parser.add_argument("-o", "--output", help = "output name (prefix)")
    args=parser.parse_args()
    return args

def overlap_fun(bed_df, target_file, reverse = False):
    # no need to get strand information for TF (only interested in its coordinates)
    if reverse: 
        tf_df = bf.read_table(bed_df, schema = "bed3").drop_duplicates(keep = "first", ignore_index = True)
        bed_target_overlap = bf.overlap(tf_df,target_file,how="left")[["chrom", "start","end", "chrom_"]].drop_duplicates(keep = "first", ignore_index = True)
        assert len(bed_target_overlap) == len(tf_df), print("Report overlap is not match to bed shape.")
        return np.array(np.where(pd.isna(bed_target_overlap["chrom_"]), 0, 1))
    else:
        tf_df = bf.read_table(target_file, schema = "bed3").drop_duplicates(keep = "first", ignore_index = True)
        bed_target_overlap = bf.overlap(bed_df,tf_df,how="left")[["chrom", "start","end", "chrom_"]].drop_duplicates(keep = "first", ignore_index = True)
        assert len(bed_target_overlap) == len(bed_df), print("Report overlap is not match to bed shape.")
        return np.array(np.where(pd.isna(bed_target_overlap["chrom_"]), 0, 1))

def main():
    ignore_warning()
    start = time.time()
    # get parameters 
    args = args_parser()
    bed = args.bed 
    target_list = args.target
    query_list = args.query
    n = args.threads
    output = args.output
    print("Read bed file")
    bed_df = bf.read_table(bed, schema = "bed3").drop_duplicates(keep = "first", ignore_index = True)
    print("Report bed overlap with target bed files ...")
    overlap_target = Parallel(n_jobs=n, prefer = "threads")(
                delayed(overlap_fun)(bed_df, target) for target in tqdm(target_list)
            )
    target_col = [os.path.basename(tf).split(".bed")[0] for tf in target_list]
    overlap_mat = np.vstack(overlap_target).transpose()
    overlap_df = pd.DataFrame(overlap_mat, columns = target_col)
    assert len(overlap_df) == len(bed_df), print("Overlap target is not match to bed input file")
    bed_target = pd.merge(bed_df, overlap_df, left_index = True, right_index = True)
    # to add query 
    if query_list is not None:
        query_col = [os.path.basename(q).split(".bed")[0] for q in query_list]
        
        qtotal = Parallel(n_jobs=n)( # total input query 
            delayed(len)(bf.read_table(q, schema = "bed3").drop_duplicates(keep = "first", ignore_index = True)) for q in query_list
        )
        ttotal = Parallel(n_jobs=n)( # total input target 
            delayed(len)(bf.read_table(t, schema = "bed3").drop_duplicates(keep = "first", ignore_index = True)) for t in target_list
        )
        qsample = Parallel(n_jobs=n)( # query included in background 
            delayed(overlap_fun)(q, bed_df, reverse = True) for q in query_list
        )
        tsample = Parallel(n_jobs=n)( # target included in background
            delayed(overlap_fun)(t, bed_df, reverse = True) for t in target_list
        )
        query_ = pd.DataFrame(list(zip(query_col, np.vstack(qsample).sum(axis =1))), columns = ["query", "qinput"]); query_["qtotal"] = qtotal
        target_ = pd.DataFrame(list(zip(target_col, np.vstack(tsample).sum(axis = 1))), columns = ["target", "tinput"]); target_["ttotal"] = ttotal
        # 
        print("Report bed overlap with query bed file ...")
        overlap_query = Parallel(n_jobs=n)(
            delayed(overlap_fun)(bed_df, q) for q in query_list
        )
        # 
        overlap_query_mat = np.vstack(overlap_query).transpose()
        overlap_query_df = pd.DataFrame(overlap_query_mat, columns = query_col)
        assert len(overlap_query_df) == len(bed_df), print("Overlap query is not match to bed input file")
        bed_query = pd.merge(bed_df, overlap_query_df, left_index = True, right_index = True)
        bed_target = pd.merge(bed_query, bed_target, on = ["chrom", 'start', 'end'])
        print("Perform fisher enrichment test ...")
        from scipy.stats import fisher_exact
        # make contingency table 
        from itertools import product
        query_target_combination = list(product(query_col, target_col))
        fisher_res = Parallel(n_jobs = n)(
            delayed(fisher_exact)(pd.crosstab(bed_target[qc[0]],bed_target[qc[-1]]), alternative = "greater") for qc in query_target_combination
        )
        # add query/target enrichment results 
        fisher_res_df = pd.merge(pd.DataFrame(query_target_combination, columns = ["query", "target"]), pd.DataFrame(fisher_res), left_index = True, right_index = True)
        # add query/target input count 
        fisher_res_df = pd.merge(pd.merge(fisher_res_df, query_, on = "query", how = "left"), target_, on = "target", how = "left")
        # query_total = [bed_target[qc[0]].sum() for qc in query_target_combination]; 
        # target_total = [bed_target[qc[-1]].sum() for qc in query_target_combination]; 
        overlap_total = [len(bed_target[bed_target[list(qc)].sum(axis = 1) == 2]) for qc in query_target_combination]
        fisher_res_df["overlap"] = overlap_total
        fisher_res_df["query_count"] = fisher_res_df["qinput"].astype(str) +"/" + fisher_res_df["qtotal"].astype(str)
        fisher_res_df["target_count"] = fisher_res_df["tinput"].astype(str) + "/" + fisher_res_df["ttotal"].astype(str)
        fisher_res_df.drop(["qtotal", "ttotal", "qinput", "tinput"], axis = 1, inplace = True)
        fisher_res_df["total"] = len(bed_target)
        fisher_res_df.sort_values(["statistic", "pvalue"], ascending = [False, True], inplace = True, ignore_index = True)
        print("Write query and target enrichment statistics ...")
        fisher_res_df.to_csv(output + ".stat", sep = "\t", header = True, index = False)
        fisher_res_df.query("pvalue < 0.001 and statistic > 1").to_csv(output + "_sig.stat", sep = "\t", header = True, index = False)
    print("Write overlap dataframe ...")
    bed_target.to_csv(output + ".txt", sep = "\t", header = True, index = False)
    print(timeit(start))

if __name__ == "__main__":
    main()
