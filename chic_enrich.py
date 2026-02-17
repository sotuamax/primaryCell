#!/usr/bin/env python3


"""
Enrichment is foldchange between the center vs. surrounding region.
"""

import pandas as pd 
import numpy as np 
import argparse 
import os
import pysam
import bioframe as bf 
from utilities.misc import ignore_warning
ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # parser.add_argument("-sample", "--sample", help="sample prefix")
    parser.add_argument("-bed", "--bed", help = "bed file for feature enrichment")
    parser.add_argument("-midbed", "--midbed", action = "store_true", help = "use the midpoint of bed peak file as center")
    parser.add_argument("-extend", "--extend", default = 0, type = int, help = "extend regions for bed.")
    parser.add_argument("-win", "--window", default = 50, type = int, help = "window size in screen")
    parser.add_argument("-bam", "--bam", help = "bam file with barcode data available.")
    # parser.add_argument("-barcode", "--barcode", action = "store_true", help = "cell barcode used to generate matrix")
    # parser.add_argument("-gene", "--gene", help = "gene set used to filter barcode (symbol, set)")
    parser.add_argument("-tag", "--tag", help = "GFP tag gene with matched barcode for ChiC reads (barcode, symbol)")
    parser.add_argument("-output", "--output", help="output file")
    parser.add_argument("-n", "--threads", type = int, default = 1, help = "threads used to process bam file")
    args = parser.parse_args()
    return args

def main():
    args = args_parser()
    bed = args.bed
    bam = args.bam
    # gene = args.gene 
    tag = args.tag
    output = args.output 
    n = args.threads
    bed_df = bf.read_table(bed, schema = "bed3")
    # gene_df = pd.read_table(gene, sep = "\t", header = 0); gene_df.columns = ["symbol", "set"]
    tag_df = pd.read_table(tag, sep = "\t", header = None, names = ["barcode", "symbol"])
    # when extend, span surrounding region on the input bed interval region
    if args.midbed: 
        print("Find midpoint of bed file ...")
        bed_df["pos"] = (bed_df["start"] + bed_df["end"])//2 
        bed_df["start"] = bed_df["pos"]; bed_df["end"] = bed_df["pos"]+1
    if args.extend != 0:
        print(f"extend bed region by {args.extend} bp")
        bed_df = bf.expand(bed_df, pad = args.extend)
        #bed_df = bf.cluster(bed_df)[["chrom", "cluster_start", "cluster_end"]].drop_duplicates(keep = "first")
        #bed_df.columns = ["chrom", "start", "end"]
        #bed_df["start"] = np.where(bed_df["start"] < 0, 0, bed_df["start"])
    bed_df = bed_df[(bed_df["chrom"].str.startswith("chr")) & (bed_df["chrom"] != 'chrM') & (bed_df["chrom"] != 'chrY')].copy()
    # bed_df.to_csv(bed.replace(".bed", ".tmp.bed"), sep = "\t", header = False, index = False)
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    # tag_df = pd.merge(tag_df, gene_df, on = "symbol", how = "left")
    # tag_df["set"] = tag_df["set"].fillna("other")
    print("Identify ChiC barcode ... ")
    tag_barcode = tag_df["barcode"].tolist()
    # total number of reads on chrom region 
    # bc_total = {bc:0 for bc in tag_barcode}
    # for c in sorted(set(bed_df["chrom"])):
    #     for read in bam_handle.fetch(c):
    #         try:
    #             bc_total[read.query_name.split(":")[-1]] += 1
    #         except:
    #             pass 
    # bc_total_df = pd.DataFrame.from_dict(bc_total, orient="index", columns = ["count"]).reset_index().rename(columns = {"index":"ChiC_barcode"})
    # each row by row
    win = args.window
    for row in bed_df.itertuples():
        size = row.end-row.start
        print(f"Screen region size {size//1000} kb")
        break 
    # count_list = list()
    bc_mat = {bc:np.zeros((len(bed_df), len(range(0, size, win))), dtype = int) for bc in tag_barcode}
    r = 0 # row number 
    for row in bed_df.itertuples():
        for e,i in enumerate(range(0, size, win)):
            for read in bam_handle.fetch(row.chrom, row.start+i, row.start+i+win):
                try:
                    bc_mat[read.query_name.split(":")[-1]][r, e] += 1
                except:
                    pass
        r += 1
        # count_list.append(row_count)
    print("Save array matrix in binary file ...")
    np.savez(f"{output}.npz", **bc_mat)
    #bc_score = np.array(bc_array).sum(axis = 0)+1
    #all_score.append(bc_score)
    #count_df = pd.DataFrame(all_score, index = [bc for bc in tag_barcode])
    # count_df.to_csv(output + ".mat", sep = "\t", header = False, index = True)
    # save all read count at cell barcode level to a matrix file 
    # calculate enrichment score by its center (+/- 20 windows around the center window)
    # and its side (+/- 10 windows around its left/right end)
    bc_list = list()
    for bc in bc_mat: 
        count_df = bc_mat[bc]
        print("Calculate fc between center / edge ...")
        count_sum = count_df.sum(axis = 0)
        center_sum = count_sum[size//win//2 - 20:size//win//2 + 20].sum()
        side_sum = count_sum[:20].sum()+count_sum[-20:].sum() + 1
         # +1 to avoid 0 in division
        bc_enrich = center_sum/side_sum
        bc_list.append((bc, bc_enrich))
    bc_enrich_df = pd.DataFrame(bc_list, columns=["barcode", "score"])
    bc_enrich_df = pd.merge(tag_df, bc_enrich_df, on = "barcode")
    bc_enrich_df.to_csv(output + ".score", sep = "\t", header = True, index = False)
    bc_score_list = list()
    for bc in tag_barcode:
        bc_c = pd.DataFrame(count_df[bc].tolist(), columns = ['up', "mid", "down"])
        # remove any peak and its surround with no reads (as 0)
        bc_c = bc_c[bc_c.sum(axis = 1) > 0].copy()
        # get the mid mean as signal value; two ends average as background (add 1 to avoid 0 as devidend)
        if len(bc_c) > 0:
            sig = bc_c["mid"].mean(); bg = (bc_c["up"].sum() + bc_c["down"].sum() + 1)/(len(bc_c)*2)
            bc_score = sig/bg
            bc_score_list.append((bc, bc_c.to_numpy().sum(), bc_c["mid"].sum(), bc_score, len(bc_c)))
    #
    score_df = pd.DataFrame(bc_score_list, columns = ["barcode", "score_read", "peak_read", "score", "score_peak"])
    score_df = pd.merge(tag_df, score_df, on = "barcode", how = "inner")
    score_df = pd.merge(bc_total_df, score_df, on = "barcode", how = "inner")
    score_df.to_csv(output + ".txt", sep = "\t", header = True, index = False)

if __name__ == "__main__":
    main()


