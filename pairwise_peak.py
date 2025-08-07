#!/usr/bin/env python3
"""
given peaks files, to report 
pairwise files comparison by 
shared, unique peaks for each file. 
"""
import pandas as pd 
import bioframe as bf 
import argparse 
import os 
from itertools import combinations 
import numpy as np 
from utilities.misc import ignore_warning
ignore_warning()

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", 
                                     formatter_class = argparse.RawDescriptionHelpFormatter, 
                                     description="pairwise_peak.py -peak p1.bed p2.bed p3.bed -tol 20 -o p_compare")
    parser.add_argument("-peak", nargs="+", help = "peaks files to compare")
    parser.add_argument("-tol", "--tolerance", default = 1, type = int, help = "tolerance for minimum distance when forming peak cluster")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def main():
    args = args_parser()
    peak_files = args.peak
    tol = args.tolerance
    df_list = list()
    for f in peak_files:
        df = bf.read_table(f, schema = "bed3")
        df["sample"] = os.path.basename(f)
        df_list.append(df)
    # combine all peak data, and report cluster when peaks within min_dist 
    peak_all = pd.concat(df_list, axis = 0)
    peak_cluster = bf.cluster(peak_all, min_dist = tol)
    peak_cluster.sort_values("cluster", inplace = True)
    # unique clusters 
    cluster_uni = peak_cluster[["chrom", "cluster_start", "cluster_end", "cluster"]].copy().drop_duplicates(keep = "first")
    # get all samples and paired samples 
    samples = sorted(set(peak_cluster["sample"]))
    paired_samples = combinations(samples, 2)
    cluster_mat = [[cluster] + [sum(cluster_df["sample"].str.contains(s)) for s in samples] for cluster, cluster_df in peak_cluster.groupby("cluster")]
    cluster_mat_df = pd.DataFrame(cluster_mat, columns = ["cluster"]+samples).set_index("cluster")
    # add cluster coordinates to sample mat 
    cluster_mat_sample = pd.merge(cluster_uni, cluster_mat_df.reset_index(), on="cluster")
    cluster_mat_sample.to_csv(f"{args.output}.mat", sep = "\t", header = True, index = False)
    # report in pairwise format 
    p_comparison = list()
    for p1,p2 in paired_samples:
        p1p2_shared = sum((cluster_mat_df[[p1,p2]] != 0).sum(axis = 1) == 2)
        p1_unique = np.where((cluster_mat_df[p1] != 0) & (cluster_mat_df[p2] == 0))[0].shape[0]
        p2_unique = np.where((cluster_mat_df[p1] == 0) & (cluster_mat_df[p2] != 0))[0].shape[0]
        p_comparison.append((p1,p2,p1p2_shared,p1_unique,p2_unique))
    pairwise_summary = pd.DataFrame(p_comparison, columns = ["sample1", "sample2", "shared", "unique1", "unique2"])
    pairwise_summary["ratio1"] = pairwise_summary["shared"]/(pairwise_summary["shared"] + pairwise_summary["unique1"])
    pairwise_summary["ratio2"] = pairwise_summary["shared"]/(pairwise_summary["shared"] + pairwise_summary["unique2"])
    pairwise_summary.to_csv(f"{args.output}.txt", sep = "\t", header = True, index = False)

if __name__ == "__main__":
    main()

