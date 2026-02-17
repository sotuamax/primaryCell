#!/usr/bin/env python3
"""
Given bed file and BAM file, to evaluate reads overlapping bed region at single cell level or bulk cell level. 
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
    parser.add_argument("-extend", "--extend", default = 0, type = int, help = "extend regions for bed.")
    parser.add_argument("-bam", "--bam", help = "bam file used to sort for read overlap")
    parser.add_argument("-barcode", "--barcode", action = "store_true", help = "cell barcode used to generate matrix")
    parser.add_argument("-output", "--output", help="output file")
    parser.add_argument("-n", "--threads", type = int, default = 1, help = "threads used to process bam file")
    args = parser.parse_args()
    return args

def main():
    # read sample GFP/ChiC results 
    args = args_parser()
    #sample = args.sample
    #output_dir = args.outdir
    bed = args.bed
    bam = args.bam
    output = args.output 
    n = args.threads 
    # bam = "output/GC760805/GC760805.dedup.qc.RG.bam"
    # bed = "ChiC_peaks.qc.bed"
    #res_file = os.path.join(output_dir, sample + "_clean_barcode.stat")
    # if not os.path.exists(res_file):
    #     raise IOError(f"Cannot find result table for {sample}") 
    # df = pd.read_table(res_file, sep = "\t", header = 0)
    # df = df.query("label == 'good'").copy()
    # chic_dict = {chic_barcode:list() for chic_barcode in df["ChiC_barcode"].tolist()}
    # locate BAM reads for cells identified as quality good
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    bed_df = bf.read_table(bed, schema = "bed3")
    if args.extend != 0:
        bed_df = bf.expand(bed_df, pad = args.extend)
        bed_df = bf.cluster(bed_df)[["chrom", "cluster_start", "cluster_end"]].drop_duplicates(keep = "first")
        bed_df.to_csv(bed.replace(".bed", ".tmp.bed"), sep = "\t", header = False, index = False)
    bed_df.columns = ["chrom", "start", "end"]; bed_df["start"] = np.where(bed_df["start"] < 0, 0, bed_df["start"])
    barcode_enrich_list = list()
    from utilities.bam_tools import bam_header
    bam_seq = bam_header(bam, seq = True)
    for c in bam_seq:
        read_list = [(read.reference_start, read.reference_end, read.query_name.split(":")[-1]) for read in bam_handle.fetch(c)]
        read_df = pd.DataFrame(read_list, columns = ["start", "end", "barcode"])
        read_df["chrom"] = c
        c_df = bed_df.query("chrom == @c")
        if len(c_df) > 0:
            read_overlap_bed = bf.overlap(read_df, c_df, how = "left")
            barcode_enrich = pd.DataFrame([(bc, len(bc_df), len(bc_df[~pd.isna(bc_df["start_"])])) for bc, bc_df in read_overlap_bed.groupby("barcode")], columns = ["barcode", "total", "overlap_peak"])
        else:
            barcode_enrich = pd.DataFrame([(bc, len(bc_df), 0) for bc, bc_df in read_df.groupby("barcode")], columns = ["barcode", "total", "overlap_peak"])
        barcode_enrich["chrom"] = c
        barcode_enrich_list.append(barcode_enrich)
    all_overlap = pd.concat(barcode_enrich_list, axis = 0)
    all_overlap.to_csv(output + ".txt", sep = "\t", header = True, index= False)

    if args.barcode:
        # barcode_df = pd.read_table(args.barcode, header = None, sep = "\t", names = ["barcode"])
        barcode_dict = {b:0 for b in sorted(set(all_overlap["barcode"]))}
        row_container = list()
        for row in bed_df.itertuples():
            row_value = barcode_dict.copy()
            for read in bam_handle.fetch(row.chrom, row.start, row.end):
                read_code = read.query_name.split(":")[-1]
                row_value[read_code] += 1
            row_container.append(np.array(list(row_value.values())))
        bed_barcode_mat = np.vstack(row_container)
        from scipy.io import mmwrite 
        from scipy.sparse import coo_matrix 
        # coo_matrix(bed_barcode_mat)
        mmwrite(output + ".mtx", coo_matrix(bed_barcode_mat))
    # # find the overlaps between ChiC reads and feature 
    # overlap_list = list()
    # for chic_barcode in chic_dict:
    #     chic_barcode_bed = pd.DataFrame(chic_dict[chic_barcode], columns = ["chrom", "start", "end"])
    #     bed_tmp = list()
    #     for bed in bed_list:
    #         feature_bed = bf.read_table(bed, schema = "bed3")
    #         chic_overlap_feature = bf.overlap(chic_barcode_bed, feature_bed, how = "inner")
    #         overlap_ratio = round(len(chic_overlap_feature)/(len(chic_barcode_bed))*100, 1)
    #         bed_tmp.append(overlap_ratio)
    #     overlap_list.append([chic_barcode] + bed_tmp)
    # overlap_df = pd.DataFrame(overlap_list, columns = ["ChiC_barcode"] + [os.path.basename(b) for b in bed_list])
    # overlap_df = pd.merge(overlap_df, df[["ChiC_barcode", "gene", "symbol"]], on = ["ChiC_barcode"], how = "left")

if __name__ == "__main__":
    main()
