#!/usr/bin/env python 
"""
This script is to identify super interactive loci (SIL) based on 

major steps: 
- identify highly interactive anchor (report anchor frequencies); 
- tether interactive anchor within a distance cutoff; 

SIL: 
span > 50 k region; or 
loop > 100; or 
top 5% interval ranked by loops; 

Input: 
    required: 
    - bedpe: loop interaction file 
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
from utilities.bed_tools import bedpe_retrieve, tether_bed

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help=True, usage="\nSIL.py <input.bedpe> -o <output>")
    parser.add_argument("bedpe", help = "input bedpe file (loop interactions)")
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
    anchor_freq = pd.DataFrame(loop_anchor.groupby(["chrom", "start", "end"]).size(), columns = ["freq"]).reset_index()
    return anchor_freq 

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
    print(f"Filter loop anchors w/ {args.min_interaction} ...")
    high_anchor = anchor_freq_df.query("freq >= @args.min_interaction")
    print(f"Stitch anchor based on {args.min_distance}")
    anchor_tethered = tether_bed(high_anchor, args.min_distance)
    # bedGraph with its value for the number of stitched anchors 
    anchor_tethered["name"] = "A:" + anchor_tethered["count"].astype(str)
    # anchor_tethered.to_csv(args.output + ".bed", sep = "\t", header = False, index = False) 
    # find loops with either end overlap the tethered super anchor 
    if args.report_loop:
        print("Report loops ... ")
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
    anchor_tethered[["chrom", "start", "end", "name"]].to_csv(args.output + ".bed", sep = "\t", header = False, index = False) 

if __name__ == "__main__":
    main()
