#!/usr/bin/env python3
"""
Docstring for loop_union
Given multiple loops files, and join all loops by chromosome positions. 

Report loop contact (count) when requested. 

"""
from utilities.bed_tools import open_loop
import pandas as pd 
import os 
import numpy as np 
import bioframe as bf 
import argparse 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help=True, usage="")
    parser.add_argument("-loops", "--loops", nargs = "+", help = "input bedpe file (loop interactions)")
    parser.add_argument("-count", "--count", action = "store_true", help = "add count data into loop file")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args=parser.parse_args()
    return args

def group_fun(df):
    group_ids = []
    current_group = 0
    for i, row in df.iterrows():
        if i == 0:
            # first row starts group 0
            group_ids.append(current_group)
        else:
            prev = df.iloc[i-1]
            # "touching" if current start - previous start <= 1; 
            if row['bin1']-prev['bin1'] <= 1 and row['bin2']-prev["bin2"] <= 1:
                group_ids.append(current_group)
            else:
                current_group += 1
                group_ids.append(current_group)
    return group_ids 

def main():
    args = args_parser()
    loop_files = args.loops
    all_loop = list()
    for file in loop_files:
        print(f"Read {file} ...")
        loop_df = open_loop(file)
        if args.count:
            loop_df.rename(columns = {"contact":os.path.basename(file).split("_")[0]}, inplace = True)
        else:
            try:
                loop_df.drop("contact", axis = 1, inplace = True)
            except:
                pass 
        all_loop.append(loop_df.set_index(["chrom1", "start1", "end1", "chrom2", "start2", "end2"]))
    print("Join all loops ...")
    all_loop_join = all_loop[0].join(all_loop[1:], how = "outer")
    all_loop_join = all_loop_join.fillna(0)
    all_loop_join.reset_index(inplace = True)
    all_loop_join.index = all_loop_join["chrom1"].astype(str) + ":" + all_loop_join["start1"].astype(str) + "-" + all_loop_join["end1"].astype(str) + "|" + all_loop_join["chrom2"].astype(str) + ":" + all_loop_join["start2"].astype(str) + "-" + all_loop_join["end2"].astype(str)
    all_loop_join.drop(["chrom1", "start1", "end1", "chrom2", "start2", "end2"], axis = 1, inplace = True)
    all_loop_join.to_csv(f"{args.output}.txt", sep = "\t", header = True, index = True)

if __name__ == "__main__":
    main()
