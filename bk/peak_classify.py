#!/usr/bin/env python3
"""
Input GTF and peaks
To classify peaks into different features based on gtf annotation. 
Briefly, 
gene body reside (promoter, non-coding, cds, intron)
outside of gene body region (intergenic region)
Note: 
peak is filtered for chromosome region (default)


use homer instead, for faster processing 

"""

from mpi4py import MPI
import bioframe as bf 
import pandas as pd 
import pysam 
import sys 
import argparse
from gtf_tools import * 
from peak_tools import * 
import logging 

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("peak", help="peak call in narrowPeak format")
    parser.add_argument("feature", help="feature file in gtf or bed format (assign using -format)")
    parser.add_argument("-o", "--output", help="output name")
    parser.add_argument("-format", "--format", choices = ["bed", "gtf"], help = "feature file format", default = "gtf")
    parser.add_argument("-TSS_range", "--TSS_range", help = "TSS range for upstream and downstream of TSS (default 500)", default = 500, type = int)
    args=parser.parse_args()
    return args

def peakassign(peak_onside):
    peak_list = list()
    for label, peak_group in peak_onside.groupby(["chrom", "start", "end", "name"]):
        if len(set(peak_group["label_"])) > 1:
            label_dict = {"cds":0, "promoter":1, "ncoding":2, "intron":3}
            peak_group = peak_group.sort_values(by = "label_", ignore_index = True, key = lambda x:x.map(label_dict))
            peak_range = np.arange(label[1],label[2], 1)
            overlap = [len(np.intersect1d(peak_range, np.arange(row.start_, row.end_, 1)))/len(np.union1d(peak_range, np.arange(row.start_, row.end_, 1))) for row in peak_group.itertuples()]
            # use iloc to find the max overlap value index (in its local index range)
            peak_pick = peak_group.iloc[[np.argmax(overlap)]]
        else:
            peak_pick = peak_group.iloc[[0]]
        peak_list.append(peak_pick)
    peak_filtered = pd.concat(peak_list, axis = 0)
    peak_filtered.drop(["chrom_", "start_", "end_"], axis = 1, inplace = True)
    return peak_filtered

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # 
    args = args_parser()
    peak = args.peak
    feature = args.feature
    output = args.output
    logging.basicConfig(filename = output+".log", filemode = "a", format = '%(asctime)s  %(message)s', datefmt = "%H:%M:%S", level = logging.DEBUG)
    if rank == 0:
        peak_df = bf.read_table(peak, schema = "bed4")
        # if q != None:
        #    peak_df = peak_q(peak_df, args.q)
        logging.info(f"Run on {size} cores")
        chrom_set = sorted(set(peak_df["chrom"]))
        worker_tasks = {w:[] for w in range(size)}
        i = 0
        for chrom in chrom_set:
            worker_tasks[i].append(chrom)
            i = (i+1)%size
    else:
        worker_tasks = None 
        peak_df = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    peak_df = comm.bcast(peak_df, root = 0)
    # retrieve each feature
    if args.format == "gtf":
        gtf_list = list()
    rank_list = list()
    if rank == 0:
        if args.format == "gtf":
            logging.info("Input feature format as GTF")
            logging.info("Start parsing GTF file ....")
        elif args.format == "bed":
            logging.info("Input feature format as bed")
        else:
            logging.info("Feature format is not supported (gtf or bed)")
    for chr in worker_tasks[rank]:
        peak_chr = peak_df[peak_df["chrom"] == chr]
        if args.format == "gtf":
            feature_chr = combine_features(feature, chr, args.TSS_range)
            gtf_list.append(feature_chr)
        elif args.format == "bed":
            feature_df = pd.read_table(feature, sep = "\t", header = None, names = ["chrom", "start", "end", "transcript", "label"])
            feature_chr = feature_df[feature_df["chrom"] == chr]
        else:
            exit(1)
        peak_overlap = bf.overlap(peak_chr, feature_chr, how = "left")
        peak_outside = peak_overlap[peak_overlap.isna().any(axis = 1)].copy()
        peak_onside = peak_overlap[~peak_overlap.isna().any(axis = 1)].copy()
        peak_filtered_onside = peakassign(peak_onside)
        peak_outside["label_"] = "intergenic"
        peak_outside.drop(["chrom_", "start_", "end_"], axis = 1, inplace = True)
        peak_examined = pd.concat([peak_filtered_onside, peak_outside], axis = 0)
        rank_list.append(peak_examined)
    if args.format == "gtf":
        gtf_rank = pd.concat(gtf_list, axis = 0)
        gtf_gather = comm.gather(gtf_rank, root = 0)
    rank_peak = pd.concat(rank_list, axis = 0)
    rank_gather = comm.gather(rank_peak, root = 0)
    if rank == 0:
        peak_assigned = pd.concat(rank_gather, axis = 0)
        peak_assigned.to_csv(f"{output}.txt", sep = "\t", index = False, header = True)
        if args.format == "gtf":
            logging.info("write GTF feature to file")
            feature_assigned = pd.concat(gtf_gather, axis = 0)
            feature_assigned.to_csv(f"gtf_feature.bed", sep = "\t", header = False, index = False)
    exit(0)

if __name__ == "__main__":
    main()
