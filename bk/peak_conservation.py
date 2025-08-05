#!/usr/bin/env python3
"""
Input bigwig (conservation score) and bed file
Conservation score for regions around peak summits
"""

import pyBigWig
import numpy as np 
import pandas as pd 
import argparse 
import os 
from bed_tools import bed_center,orphan_bed
import matplotlib.pyplot as plt 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="this program is for retriving the conservative score for assigned bed region")
    parser.add_argument("bed", help = "bed input")
    parser.add_argument("bigwig", help = "bigwig file with nucleo-level conserve scores")
    parser.add_argument("-flank", "--flank", default = 500, type = int, help = "flanking side size centered at bed region pospoint (default: 500)")
    parser.add_argument("-o", "--output", required = False, help = "output directory")
    args=parser.parse_args()
    return args

def main():
    # all ranks run for argument parsing 
    args = args_parser()
    bw = args.bigwig
    bed = args.bed
    flank = args.flank
    bw_handle = pyBigWig.open(bw)
    bed_df = bed_center(bed) # 2 columns: chrom, pos
    bed_df = orphan_bed(bed_df, flank)
    # to get values within a range 
    if args.output == None:
        output = os.path.basename(bed).split(".bed")[0]
    else:
        output = args.output
    # bed_pos_df_chr = bed_pos_df[bed_pos_df["contig"].isin(chr_pick)]
    bed_score = np.array([bw_handle.values(row.chrom, row.pos-flank, row.pos+flank+1) for row in bed_df.itertuples()])
    # pd.DataFrame(bed_pos_score, columns = list(range(-flank, flank+1, 1))).to_csv(os.path.join(outdir, bed_name+".txt"), sep = "\t", header = True, index = False)
    bed_score = pd.DataFrame(bed_score)
    bed_score = bed_score[~bed_score.isna().any(axis = 1)].copy()
    score_mean = bed_score.mean(axis = 0)
    score_df = pd.DataFrame(score_mean, columns = ["score"])
    score_df["pos"] = range(-flank, flank+1, 1)
    score_df[["pos", "score"]].to_csv(f"{output}.txt", sep = "\t", index = False, header = True)
    # plot 
    fig, axis = plt.subplots(1,1)
    axis.plot(score_df["pos"], score_df["score"])
    # axis.axvline(x = 0, ymin = score_df["scaled_coverage_mean"].min(), ymax = score_df["scaled_coverage_mean"].max(), linestype = "--", size = 0.5, )
    # axis.set_xticks(np.array(score_df["pos"][::tick_break]))
    axis.set_xlabel("distance from centroid")
    bwname = os.path.basename(bw).split(".bw")[0]
    axis.set_ylabel(f"score for {bwname}")
    plt.savefig(f"{output}.pdf", format = "pdf")

if __name__ == "__main__":
    main()
