#!/usr/bin/env python3

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
    parser.add_argument("-extend", "--extend", default = 0, type = int, help = "extend regions for bed.")
    # parser.add_argument("-step", "--step", default = 50, type = int, help = "step size in screen")
    parser.add_argument("-bam", "--bam", help = "bam file with barcode data available.")
    # parser.add_argument("-barcode", "--barcode", action = "store_true", help = "cell barcode used to generate matrix")
    # parser.add_argument("-gene", "--gene", help = "gene set used to filter barcode (symbol, set)")
    parser.add_argument("-tag", "--tag", help = "GFP tag gene with matched barcode (barcode, symbol, label)")
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
    tag_df = pd.read_table(tag, sep = "\t", header = 0); tag_df = tag_df.query("label == 1")
    try:
        tag_df.drop(["total", "overlap_peak", "enrichment", "label"], axis = 1, inplace = True)
    except:
        pass 
    # when extend, span surrounding region on the input bed interval region
    if args.extend != 0:
        bed_df = bf.expand(bed_df, pad = args.extend)
        bed_df = bf.cluster(bed_df)[["chrom", "cluster_start", "cluster_end"]].drop_duplicates(keep = "first")
        bed_df.columns = ["chrom", "start", "end"]
        bed_df["start"] = np.where(bed_df["start"] < 0, 0, bed_df["start"])
    bed_df = bed_df[(bed_df["chrom"].str.startswith("chr")) & (bed_df["chrom"] != 'chrM') & (bed_df["chrom"] != 'chrY')].copy()
    bed_df.to_csv(bed.replace(".bed", ".tmp.bed"), sep = "\t", header = False, index = False)
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    # tag_df = pd.merge(tag_df, gene_df, on = "symbol", how = "left")
    # tag_df["set"] = tag_df["set"].fillna("other")
    tag_barcode = tag_df["ChiC_barcode"].tolist()
    # each row by row
    count_list = list()
    # total number of reads on chrom region 
    bc_total = {bc:0 for bc in tag_barcode}
    for c in sorted(set(bed_df["chrom"])):
        for read in bam_handle.fetch(c):
            try:
                bc_total[read.query_name.split(":")[-1]] += 1
            except:
                pass 
    bc_total_df = pd.DataFrame.from_dict(bc_total, orient="index", columns = ["count"]).reset_index().rename(columns = {"index":"ChiC_barcode"})
    
    for row in bed_df.itertuples():
        size = row.end-row.start
        row_count = {bc:[0,0,0] for bc in tag_barcode}
        for read in bam_handle.fetch(row.chrom, row.start, row.end):
            read_code = read.query_name.split(":")[-1]
            try:
                row_count[read_code][1] += 1
            except:
                pass
        for read in bam_handle.fetch(row.chrom, row.start-size, row.start-20):
            read_code = read.query_name.split(":")[-1]
            try:
                row_count[read_code][0] += 1
            except:
                pass
        for read in bam_handle.fetch(row.chrom, row.end + 20, row.end+size):
            read_code = read.query_name.split(":")[-1]
            try:
                row_count[read_code][2] += 1
            except:
                pass
        count_list.append(row_count)
    count_df = pd.DataFrame(count_list)
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
    score_df = pd.DataFrame(bc_score_list, columns = ["ChiC_barcode", "score_read", "peak_read", "score", "score_peak"])
    score_df = pd.merge(tag_df, score_df, on = "ChiC_barcode", how = "inner")
    score_df = pd.merge(bc_total_df, score_df, on = "ChiC_barcode", how = "inner")
    score_df.to_csv(output + ".txt", sep = "\t", header = True, index = False)

if __name__ == "__main__":
    main()


