#!/usr/bin/env python 
"""
This script is to identify super interactive loci (SIL) based on significant interactions across chromatin. 

major steps: 
- identify highly interactive anchor (report anchor frequencies); 
- tether interactive anchor within a distance cutoff; 

SIL: 
span > 50 k region; or 
loop > 100; or 
top 5% interval ranked by loops; 

Input: 
    required: 
    - loop: loop interaction file 
    - output: output prefix 
    optional:
    - min_interaction: 
    - min_dist: 
    - report_loop: report loops for stitched interval

"""

import pandas as pd 
import os 
import numpy as np 
import bioframe as bf 
import argparse 
import matplotlib.pyplot as plt

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help=True, usage="\nSIL.py <input.bedpe> -o <output>")
    parser.add_argument("loop", help = "input bedpe/longrange format file (loop interactions)")
    # parser.add_argument("-bam", "--bam", help = "BAM file used to report for signal value at candidate sites.")
    parser.add_argument("-min_interaction", "--min_interaction", type = int, default = 3, help = "minimum interactions for anchor")
    parser.add_argument("-min_dist", "--min_distance", type = int, default = 10000, help = "minimum distance to stitch anchors")
    parser.add_argument("-report_loop", "--report_loop", action = "store_true", help = "report loops that overlap with Super interactive loci")
    parser.add_argument("-o", "--output", help = "output name")
    args=parser.parse_args()
    return args

def anchor_freq(loop_df):
    """
    Identify candidates of super interactive anchors, that is loop anchors that involved in considerable more interactions. 
    Inputs:
    bedpe: file to read in the format of bedpe 
    relative: filtering on relative scale (e.g., involved in loops number more than other 90% anchors); 
    absolute: filtering on absolute sclae (e.g., involved in loops number of at least 3);
    Return: 
    anchor_df: pandas df that pass the filtering (w/ columns: chrom, start, end, freq)
    """
    loop_anchor = pd.concat([loop_df[["chrom1", 'start1', "end1"]].rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"}), loop_df[["chrom2", "start2", "end2"]].rename(columns = {"chrom2":"chrom", "start2":"start", "end2":"end"})], axis = 0) 
    anchor_freq = loop_anchor.groupby(["chrom", "start", "end"]).size().reset_index(name = "freq")
    return anchor_freq 

def tether_bed(anchor_df, distance):
    """
    Tether anchors when the distance between is less than a distance cutoff
    Inputs: 
    anchor_df: pandas df in bed format (minimum: chrom, start, end);
    distance: distance cutoff use to tether anchors;
    Return: 
    df: tethered anchor in bed format (labeled by: its length, and number of anchors contributed to the tethering)
    """
    anchor_cluster = bf.cluster(anchor_df, min_dist = distance)
    cluster_anchor_freq = anchor_cluster.groupby(["chrom", "cluster_start", "cluster_end"]).size().reset_index(name = "anchor_count") # ["start"].transform("nunique")
    cluster_loop_freq = anchor_cluster.groupby(["chrom", "cluster_start", "cluster_end"])["freq"].sum().reset_index(name = "loop_count")
    anchor_cluster = pd.merge(cluster_loop_freq, cluster_anchor_freq, on = ["chrom", "cluster_start", "cluster_end"])
    anchor_cluster.rename(columns = {"cluster_start":"start", "cluster_end":"end"}, inplace = True)
    return anchor_cluster

def add_value(anchor_tethered):
    anchor_tethered["per_anchor"] = anchor_tethered["loop_count"]/anchor_tethered["anchor_count"]
    # anchor_tethered["value"] = anchor_tethered["score"]*(anchor_tethered["loop_count"]/anchor_tethered["anchor_count"])
    anchor_tethered["name"] = "loop:" + anchor_tethered["loop_count"].astype(str) + ";anchor:" + anchor_tethered["anchor_count"].astype(str) + ";per_anchor:" + anchor_tethered["per_anchor"].astype(str)

def main():
    args = args_parser()
    file = args.loop
    from utilities.bed_tools import open_loop 
    print(f"Read loop data {file} ...")
    loop_df = open_loop(file)
    output = args.output
    #loop_df = loop_df[["chrom1", 'start1', "end1", "chrom2", "start2", "end2"]].copy()
    print("Generate loop anchor frequencies ...")
    anchor_freq_df = anchor_freq(loop_df)
    q80 = np.quantile(anchor_freq_df["freq"], 0.8)
    min_interaction = q80 if q80 > args.min_interaction else args.min_interaction
    print("Plot anchors ranked by loops per anchor involved in ...")
    plt.figure(figsize = (6,5))
    plt.title(os.path.basename(output))
    plt.plot(range(len(anchor_freq_df)), sorted(anchor_freq_df["freq"]), marker = "o", markersize = 0.3)
    plt.text(x = 1, y = min_interaction+1, s = f"Threshold={min_interaction}", color = "red")
    plt.axhline(y = min_interaction, color = "r", linestyle = "--", linewidth = 0.8)
    plt.xlabel("Rank")
    plt.ylabel("Loops per anchor")
    plt.savefig(f"{output}_anchor.pdf")
    plt.close()
    print(f"Filter top anchors w/ > {min_interaction} loops ...")
    top_anchor = anchor_freq_df.query("freq > @min_interaction")
    print(f"Stitch anchors based on {args.min_distance} bp ...")
    anchor_tethered_top = tether_bed(top_anchor, args.min_distance)
    add_value(anchor_tethered_top)
    anchor_min = np.quantile(anchor_tethered_top["anchor_count"], 0.8)
    anchor_min = anchor_min if anchor_min > 3 else 3 
    print(f"Min anchor per cluster > {anchor_min} ...")
    anchor_tethered_top_min_anchor = anchor_tethered_top.query("anchor_count > @anchor_min")
    print("Write clustered-anchor as candidate file ...")
    anchor_tethered_top_min_anchor[["chrom", "start", "end", "name"]].to_csv(f"{output}_candidate.bed", sep = "\t", header = False, index = False)
    per_anchor_min = np.quantile(anchor_tethered_top_min_anchor["per_anchor"], 0.8)
    per_anchor_min = per_anchor_min if per_anchor_min > 3 else 3
    print("Plot clustered-anchors ranked by loops per-anchor involved in ...")
    plt.figure(figsize = (6,5))
    plt.title(os.path.basename(output))
    plt.plot(range(len(anchor_tethered_top_min_anchor)), sorted(anchor_tethered_top_min_anchor["per_anchor"]), marker = "o", markersize = 0.3)
    plt.axhline(y = per_anchor_min, color = "r", linestyle = "--", linewidth = 0.8)
    plt.text(x = 1, y = per_anchor_min+1, s = f"Threshold={per_anchor_min}", color = "red")
    plt.xlabel("Rank")
    plt.ylabel("Loops per anchor within clustered-anchor")
    plt.savefig(f"{output}_clusteranchor.pdf")
    plt.close()
    print(f"Min loops per anchor > {per_anchor_min} ...")
    print("Write super interaction into file ...")
    anchor_tethered_filtered = anchor_tethered_top_min_anchor.query("per_anchor > @per_anchor_min")
    anchor_tethered_filtered[["chrom", "start", "end", "name"]].to_csv(f"{output}_clean.bed", sep = "\t", header = False, index = False)
    if args.report_loop:
        print("Report loops ... ")
        from utilities.bed_tools import bedpe_retrieve
        loop_count = list()
        overlap_loops = list()
        for c, c_df in anchor_tethered.groupby(["chrom", "start", "end"]):
            c_loop = bedpe_retrieve(c, loop_df)
            c_loop["sil"] = f"{c[0]}:{c[1]}-{c[2]}"
            overlap_loops.append(c_loop)
            c_df["loops"] = len(c_loop)
            loop_count.append(c_df)
        # add SIL and its overlaying loop number 
        anchor_tethered = pd.concat(loop_count, axis = 0)
        #anchor_tethered.drop("count", axis = 1, inplace = True)
        anchor_tethered["name"] = "A:" + anchor_tethered["count"].astype(str) + ";" + "L:" + anchor_tethered['loops'].astype(str) 
        # A: for anchor; L for loops 
        # add SIL and its loop interactions 
        S_loops = pd.concat(overlap_loops, axis = 0)
        S_loops.to_csv(args.output + ".bedpe", sep = "\t", header = False, index = False)

if __name__ == "__main__":
    main()
