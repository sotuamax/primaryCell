#!/usr/bin/env python3
"""
This program is to generate unified peaks from multiple peaks data (summited based clustering). 

Input:  
    required: 
    - peak (in narrowPeak format); 
    - output prefix;
    optional: 
    - min_dist: distance used for summit clustering; 
    - min_score: used to filter peak on score; 
    - min_fc: used  to fileter peak on fc; 
    - saf: to generate saf format output

"""

import sys 
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

import pandas as pd 
import bioframe as bf 
import os 
import argparse 
import numpy as np 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-peaks", "--peaks", nargs = "+", help="peaks file in narrowPeak format (contain summits)")
    parser.add_argument("-blacklist", "--blacklist", help = "blacklist region", default = "/data/jim4/data/blacklist/hg38_blacklist_jointed.bed")
    parser.add_argument("-min_dist", "--min_dist", default = 50, type = int, help = "minimum distance for peak summits, used for clustering")
    parser.add_argument("-min_score", "-min_score", default = 0, type = int, help = "minimum score used to filter peaks (10xqValue). ")
    parser.add_argument("-min_fc", "--min_fc", default = 0, type = float, help = "minimum foldchange for peak signal")
    parser.add_argument("-saf", "--saf", action = "store_true", help = "output SAF format for the consensus peaks")
    parser.add_argument("-o", "--output", help = "output name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    peaks = args.peaks
    output= args.output
    # peak_df.to_csv("test.narrowPeak", header = False, index = False, sep = "\t")
    # exit(1)
    peak_list = list()
    for peak in sorted(peaks):
        print(f"Read {peak} ...")
        peak_df = bf.read_table(peak, schema = "narrowPeak")
        peak_count_raw = len(peak_df)
        print("Total peaks: ", len(peak_df))
        if args.min_score > 0 or args.min_fc > 0:
            from utilities.chrom_func import chrom_size
            print("Filtering peaks on chromosome only ...")
            chroms = sorted(set(chrom_size()["chrom"]))
            print(f"Filtering peaks on minimum score {args.min_score}")
            peak_df = peak_df.query("(chrom in @chroms) and (score >= @args.min_score or fc >= @args.min_fc)").copy()
            if len(peak_df) != peak_count_raw:
                print("Write filtered peak into file ...")
                peak_df.to_csv(peak.replace(".narrowPeak", "_qc.narrowPeak"), sep = "\t", header = False, index = False)
            print("Filtered peaks: ", len(peak_df), round(len(peak_df)/peak_count_raw*100, 2))
        peak_list.append(peak_df)
    peak_df = pd.concat(peak_list, axis = 0)
    # make sure peaks name are unique
    peak_df["name"] = peak_df["name"].str.split("_peak", expand = True)[0]
    new_peak_df = list()
    for na, na_df in peak_df.groupby("name"):
        na_df = na_df.sort_values(["chrom", "start"])
        na_df["name"] = range(1, 1+len(na_df))
        na_df["name"] = na + "_"+ na_df["name"].astype(str)
        new_peak_df.append(na_df)
    peak_df = pd.concat(new_peak_df, axis = 0)
    print("Peak names are unique !")
    narrowPeak_col = peak_df.columns.tolist()
    # extend peaks given min dist
    # peak_df["name"] = peak_df["name"].astype(str)
    # peak_df["summit"] = peak_df["start"] + peak_df["relSummit"]
    # peak_df["start"] = np.where(peak_df["summit"] - peak_df["start"] < args.min_dist, peak_df["summit"]-args.min_dist, peak_df["start"])
    # peak_df["end"] = np.where(peak_df["end"] - peak_df["summit"] < args.min_dist, peak_df["summit"]+args.min_dist, peak_df["end"])
    # peak_df["relSummit"] = peak_df["summit"] - peak_df["start"]
    # peak_df.drop("summit", axis=1, inplace = True)
    # get summit bed 
    summit_df = peak_df[["chrom", "start", "relSummit", "name"]].copy()
    summit_df["summit"] = summit_df["start"] + summit_df["relSummit"]
    summit_df["start"] = summit_df["summit"]; summit_df["end"] = summit_df["summit"]+1
    summit_df = summit_df[["chrom", "start", "end", "name", "summit"]]
    # cluster peaks/summits based on distance cutoff (min_dist)
    from utilities.peak_tools import cluster_
    print(f"Use min distance at {args.min_dist} bp to cluster summits ...")
    peak_cluster = cluster_(peak_df, min_dist = 3) # cluster of peak 
    summit_cluster = cluster_(summit_df, min_dist = args.min_dist) # cluster of summits
    # cluster is for peak_cluster, cluster_ is for summit_cluster
    clustered_df = pd.merge(peak_cluster, summit_cluster[["name", "summit", "cluster"]].rename(columns = {"cluster":"cluster_"}), on = "name", how = "outer")
    del peak_df, summit_df, peak_cluster, summit_cluster
    # update new summit as the midpoint of summit clusters 
    print("For summit cluster, update new summit to the median position ...")
    clustered_df["summit"] = clustered_df.groupby(["cluster", "cluster_"])["summit"].transform("median").astype(int)
    # update peak interval to contain the range of all peaks for clustered summits (name update to refer to all the peaks source)
    clustered_df["start"] = clustered_df.groupby("cluster")["start"].transform("min").astype(int)
    clustered_df["end"] = clustered_df.groupby("cluster")["end"].transform("max").astype(int)
    clustered_df["name"] = clustered_df.groupby(["cluster", "cluster_"])["name"].transform(lambda x:",".join(x))
    # update score and fc as the mean when clustering applied 
    print("Update score value as the mean score value of clustered peaks ...")
    clustered_df["score"] = clustered_df.groupby(["cluster", "cluster_"])["score"].transform("mean").astype(float).round(1)
    print("Update fc value to the mean fc value of clustered peaks ...")
    clustered_df["fc"] = clustered_df.groupby(["cluster", "cluster_"])["fc"].transform("mean").astype(float).round(1)
    clustered_df["-log10p"] = clustered_df.groupby(["cluster", "cluster_"])["-log10p"].transform("mean").astype(float).round(2)
    clustered_df["-log10q"] = clustered_df.groupby(["cluster", "cluster_"])["-log10q"].transform("mean").astype(float).round(2)
    clustered_df.drop(["relSummit"], axis = 1, inplace = True)
    #clustered_df["relSummit"] = clustered_df["summit"]-clustered_df["start"]
    clustered_df.drop_duplicates(keep = "first", inplace = True, ignore_index = True)
    # clustered_df.to_csv("tmp.txt", sep = "\t", header = True, index = False)
    clustered_df["n"] = clustered_df.groupby(["chrom", "start", "end"])["summit"].transform('nunique')
    # clustered_df.query("n == 1")
    # for row in .iterrows():
    # for peaks with > 1 summits, renew its peak intervals
    print("Update peaks with more than 1 summits (split peaks based on summit ...")
    print("Separate by 50 bp ...")
    multi_summit = clustered_df.query("n != 1").copy().sort_values(by = ["chrom", "start", "end", "summit"])
    new_c_df = list()
    for c,c_df in multi_summit.groupby(["chrom", "start", "end"], as_index = False):
        summit_list = c_df["summit"].tolist()
        for s in range(len(summit_list)-1):
            summit_mid = int(np.mean(summit_list[s:s+2]))
            c_df.iat[s, 2] = summit_mid - 25 
            c_df.iat[s+1, 1] = summit_mid + 25 # separated by 40 bp distance 
        new_c_df.append(c_df)
    multi_summit_new = pd.concat(new_c_df, axis = 0)
    clustered_df_new = pd.concat([clustered_df.query("n == 1"), multi_summit_new], axis = 0)
    clustered_df_new["relSummit"] = clustered_df_new["summit"]-clustered_df_new["start"]
    clustered_df_new.sort_values(by = ["chrom", "start"], inplace = True)
    if args.blacklist is not None:
        print("Remove peaks overlappying blacklist site ...")
        from utilities.peak_tools import rmblacklist
        blacklist_df = bf.read_table(args.blacklist, schema = "bed3")
        clustered_df = rmblacklist(clustered_df_new, blacklist_df)
    # write into output file
    print(f"Total consensus peaks: {len(clustered_df_new[narrowPeak_col])}")
    consensus_all = clustered_df_new[narrowPeak_col].copy()
    consensus_all.to_csv(output + ".narrowPeak", sep = "\t", header = False, index = False)
    if args.saf:
        print("report in SAF & bed format ...")
        saf = clustered_df_new[["name", "chrom", "start", "end"]].copy()
        saf["name_new"] = range(len(saf))
        saf["name_new"] = "peak_" + saf["name_new"].astype(str)
        saf["strand"] = "."
        saf[["name", "name_new"]].to_csv(output+".saf.peak", sep = "\t", header = False, index = False)
        saf[["name_new", "chrom", "start", "end", "strand"]].to_csv(output + ".saf", sep = "\t", header = False, index = False)
        saf_new_bed = saf[["chrom", "start", "end", "name_new"]]
        saf_new_bed.to_csv(output + ".bed", sep = "\t", header = False, index = False) # standard bed4 format
        saf_new_bed = bf.read_table(output + ".bed", schema = "bed4")
        for p in args.peaks:
            print(p)
            pdf = bf.read_table(p, schema = "narrowPeak")
            saf_new_bed = bf.overlap(saf_new_bed, pdf, how = "left")
            saf_new_bed[os.path.basename(p).split("_peaks.narrowPeak")[0]] = np.where(saf_new_bed["chrom_"].isna(), 0, 1)
            saf_new_bed.drop([c for c in saf_new_bed.columns if c.endswith("_")], axis = 1, inplace = True)
            saf_new_bed.drop_duplicates(keep = "first", ignore_index = True, inplace = True)
        saf_new_bed.drop(["chrom", "start", "end"], axis = 1).to_csv(output + "_01.txt", sep = "\t", header = True, index = False)
        for c in saf_new_bed.columns:
            saf_new_bed[saf_new_bed[c] == 1][[c]].to_csv(f"{c}_consensus.bed", sep = "\t", header = False, index = False)
            

if __name__ == "__main__":
    main()
