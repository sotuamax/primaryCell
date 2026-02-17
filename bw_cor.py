#!/usr/bin/env python3
import pandas as pd 
import bioframe as bf 
import os 
import numpy as np 
import argparse
import pyBigWig 

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, usage="")
    parser.add_argument("-bw1", "--bw1", nargs = "+", help = "bigwig files for correlation analysis") 
    parser.add_argument("-bw2", "--bw2", nargs = "+", help = "bigwig files used for correlation analysis")
    parser.add_argument("-bed", "--bed", help = "bed file")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def main():
    args = args_parser()
    bw_files1 = args.bw1 
    bw_files2 = args.bw2 
    bed = args.bed 
    output = args.output 
    bed_df = bf.read_table(bed, schema = "bed3")
    # scan bed regions with bigwig files to get its score 
    # bigwig files in 1st input 
    score1_list = list(); score1_name = list()
    for bw1 in bw_files1:
        print(f"Process {bw1} ...")
        bw1_name = os.path.basename(bw1).split(".bw")[0] + "_1"
        score1_name.append(bw1_name)
        bw_handle1 = pyBigWig.open(bw1, "r")
        bw1_score = np.array([bw_handle1.stats(row.chrom, row.start, row.end) for row in bed_df.itertuples()])
        score1_list.append(bw1_score)
    score1_all= np.hstack(score1_list)
    score1_df = pd.DataFrame(score1_all, columns=score1_name)
    score1_df.to_csv(f"{output}_bw1.txt", sep = "\t", header = True, index = False)
    # bigwig files in 2nd input 
    score2_list = list(); score2_name = list()
    for bw2 in bw_files2:
        print(f"Process {bw2} ...")
        bw2_name = os.path.basename(bw2).split(".bw")[0] + "_2"
        score2_name.append(bw2_name)
        bw_handle2 = pyBigWig.open(bw2, "r")
        bw2_score = np.array([bw_handle2.stats(row.chrom, row.start, row.end) for row in bed_df.itertuples()])
        score2_list.append(bw2_score)
    score2_all = np.hstack(score2_list)
    score2_df = pd.DataFrame(score2_all, columns = score2_name)
    score2_df.to_csv(f"{output}_bw2.txt", sep = "\t", header = True, index = False)
    # combine all scores 
    all_score_df = pd.concat([score1_df, score2_df], axis = 1)
    # score correlation in spearman 
    score_spearman = all_score_df.corr(method = "spearman")
    # filter the correlation matrix so that bigwig 1st input as row and 2nd input as column 
    score_i = [i for i in score_spearman.index if i in score1_name]
    score_c = [c for c in score_spearman.columns if c in score2_name]
    # write input output file
    score_spearman.loc[score_i][score_c].to_csv(f"{output}.txt", sep = "\t", header = True, index = True, float_format="%.2f")

if __name__ == "__main__":
    main()
