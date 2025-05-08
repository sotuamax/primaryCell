#!/usr/bin/env python 
"""
This script is to identify super interactive loci (SIL)

major steps: 
- identify highly interactive anchor (report anchor frequencies); 
- tether interactive anchor within a distance cutoff; 
- count loops that relates to the tethered anchor loci (report loops that overlaps SILs);

Input: 
    required: 
    - bedpe: loop interaction file 
    - 
"""

import pandas as pd 
import os 
import numpy as np 
import bioframe as bf 
import argparse 
from bed_tools import bedpe_retrieve

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("bedpe", help = "input bedpe file (loop interactions)")
    parser.add_argument("-min_interaction", "--min_interaction", type = int, default = 3, help = "minimum interactions used to filter highly interactive anchor")
    parser.add_argument("-tether_dist", "--tether_distance", type = int, default = 10000, help = "minimum distance cutoff to tether genomic intervals together")
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
    anchor_freq = pd.DataFrame(loop_anchor.groupby(["chrom", "start", "end"]).size(), columns = ["freq"]).reset_index()
    return anchor_freq 

def tether_bed(anchor_df, distance):
    """
    Tether anchor intervals when the distance between is less than a cutoff
    Inputs: 
    anchor_df: pandas df in bed format (minimum: chrom, start, end);
    distance: distance cutoff use to tether anchors;
    Return: 
    df: tethered anchor in bed format (labeled by: its length, and number of anchors contributed to the tethering)
    """
    anchor_stitch_df = bf.cluster(anchor_df, min_dist = distance)
    anchor_stitch_df["length"] = anchor_stitch_df["cluster_end"] - anchor_stitch_df["cluster_start"]
    anchor_stitch_df["count"] = anchor_stitch_df.groupby(["chrom", "cluster_start", "cluster_end"])["start"].transform("nunique")
    # select columns of interest
    anchor_stitch_df = anchor_stitch_df[["chrom", "cluster_start", "cluster_end", "length", "count"]].rename(columns = {"cluster_start":"start", "cluster_end":"end"})
    # remove duplicates 
    anchor_return = anchor_stitch_df.drop_duplicates(keep = "first")
    return anchor_return

def main():
    args = args_parser()
    bedpe = args.bedpe 
    loop_df = pd.read_table(bedpe, sep = "\t", header = None, names = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2"])
    anchor_freq_df = anchor_freq(loop_df)
    dist_list = list()
    for c, chr_df in anchor_freq_df.groupby("chrom"):
        i = 0
        for row in chr_df.itertuples():
            if i == 0:
                start = row.start 
                end = row.end 
                dist_list.append(0)
            else: 
                dist_list.append(row.start - end)
                start = row.start 
                end = row.end 
            i += 1 
    anchor_freq_df["distance"] = dist_list 
    # 
    anchor_freq_df.to_csv(args.output + ".anchor", sep = "\t", header = True, index = False)
    high_anchor = anchor_freq_df.query("freq >= @args.min_interaction")
    anchor_tethered = tether_bed(high_anchor, args.tether_distance)
    anchor_tethered.to_csv(args.output + ".tethered", sep = "\t", header = True, index = False)
    # find loops with either end overlap the tethered super anchor 
    if args.report_loop:
        loop_count = list()
        overlap_loops = list()
        for c, c_df in anchor_tethered.groupby(["chrom", "start", "end"]):
            c_loop = bedpe_retrieve(c, loop_df)
            c_loop["chrom"] = c[0]; c_loop["start"] = c[1]; c_loop["end"] = c[2]
            overlap_loops.append(c_loop)
            c_df["loops"] = len(c_loop)
            loop_count.append(c_df)
        # 
        anchor_tethered = pd.concat(loop_count, axis = 0)
        anchor_tethered.to_csv(args.output + ".tethered", sep = "\t", header = True, index = False)
        S_loops = pd.concat(overlap_loops, axis = 0)
        S_loops.to_csv(args.output + ".loop", sep = "\t", header = False, index = False)

if __name__ == "__main__":
    main()
