#!/usr/bin/env python 
"""
This script is for retrieving the bigwig value using a step size in a chrom region

"""
import pandas as pd 
import argparse 
from utilities.chrom_func import chrom_arms
from bw_tools import bw_value

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="bw2bedGraph sample.bw -step 5000 -o sample")
    parser.add_argument("bw", help = "input bw file")
    parser.add_argument("-step", "--step_size", type = int, help = "step size used to profile")
    parser.add_argument("-o", "--output", required = False, help = "output name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    bw = args.bw 
    step = args.step_size 
    # get chrom arms 
    arms = chrom_arms("hg38")
    # bw_value return a list of lists (per list container matches to per row in bed input)
    bw_score = bw_value(bw, arms, step)
    arm_score_list = list()
    # 
    for row in arms.itertuples():
        arm_score = pd.DataFrame(bw_score[row.Index], columns = ["score"])
        arm_score["chrom"] = row.chrom 
        arm_score["start"] = range(row.start, row.end, step)
        arm_score["end"] = arm_score["start"] + step
        arm_score_list.append(arm_score[["chrom", "start", "end", "score"]])
    # combine all score 
    arm_score_all = pd.concat(arm_score_list, axis = 0)
    arm_score_all.to_csv(f"{args.output}.bedGraph", sep = "\t", header = False, index = False)

if __name__ == "__main__":
    main()
