from utilities.bed_tools import open_loop
import pandas as pd 
import os 
import numpy as np 
import bioframe as bf 
import argparse 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help=True, usage="\nSIL.py <input.bedpe> -o <output>")
    parser.add_argument("loop", help = "input bedpe file (loop interactions)")
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
    file = args.loop
    loop_df = open_loop(file)
    loop_df["bin1"] = loop_df["start1"]//1000
    loop_df["bin2"] = loop_df["start2"]//1000 
    loop_connect = pd.DataFrame(list(zip(loop_df["chrom1"], loop_df["bin1"], loop_df["bin2"])), columns = ["chrom", "bin1", "bin2"])
    for chr, chr_sub in loop_connect.groupby("chrom"):
        chr_group = group_fun(chr_sub.sort_values(["bin1", "bin2"]))
        chr_sub["group"] = chr_group
