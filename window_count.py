#!/usr/bin/env python
"""
Scan the genome using 1 kb window and count reads per window. KMeans with 2 clusters used to model the window count data and identify windows with more counts. Then the identified counts cluster was fit to a negative binomial distribution, so that each count can report a p-value. To take consideration of local signal, fold change per window was calculated by comparing its count value with surrounding two windows each site. 

Final enrichment regions reported as: 

1) adjusted p-value for the window count < 0.05 (control at genome wide); 
2) fold change > 1.5 (control at local site). 

"""
"""
        count values: 
            - used for kmeans to remove very small count bins; 
            - used to fit negative binomial distribution; 
            - used to calculate foldchange of bin against neighboring bins; 
        max values (to make sure reads are relatively concentrated within 1 k bin): 
        """

import pandas as pd 
import pysam 
import bioframe as bf 
import argparse 
from utilities.chrom_func import chrom_arms
import numpy as np
from scipy import stats
from sklearn.cluster import KMeans
from tqdm import tqdm 
from statsmodels.stats.multitest import multipletests
from scipy.optimize import minimize
from utilities.peak_tools import rmblacklist
import pyBigWig 
from utilities.misc import timeit
import time 
from joblib import Parallel, delayed
from utilities.misc import ignore_warning
ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-bw", "--bw", help = "BigWig input to examine (only support alignment hg38 for now)")
    parser.add_argument("-bam", "--bam", help = "Bam input to examine (only support alignment hg38 for now)")
    # parser.add_argument("-mode", "--mode", choices=["bw", "bam"], help = "input file mode to run for score (for input as bigwig, the score must be integer)")
    parser.add_argument("-w", "--window", default = 1000, type = int, help = "window size used to screen genome")
    parser.add_argument("-n", "--threads", default = 1, type = int, help = "threads to use to process bam file (applicable for BAM input)")
    parser.add_argument("-blacklist", "--blacklist", required = False, help = "blacklist region in bed file")
    # parser.add_argument("-padj", "--padj", type = float, help = "cutoff at adjusted p-value")
    parser.add_argument("-o", "--output", help = "output name")
    args = parser.parse_args()
    return args

def bw_max(chr, start, end, bw_handle):
    return bw_handle.stats(chr, start, end, type = "max")

def main():
    starttime = time.time()
    args = args_parser()
    bw = args.bw; bam = args.bam
    n = args.threads
    w = args.window  
    i = 0
    #padj_cutoff = args.padj 
    #print(f"Cutoff for padj {padj_cutoff}!")
    blacklist = args.blacklist 
    output = args.output
    arm_df = chrom_arms()
    enrich_list = list()
    for c, chr_df in arm_df.groupby("chrom"):
        print(f"Process {c} ...")
        chr_region_list = list()
        chr_max_list = list()
        chr_count_list = list()
        for row in chr_df.itertuples():
            chr = row.chrom; start = row.start; end = row.end
            # print(chr, start, end)
            # all region as array 
            region_ar = np.array(range(start, end-w, w)).reshape(-1)
            chr_region_list.append(region_ar); del region_ar
            # from bam to get read count across a window region
            bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
            region_count = np.array([bam_handle.count(chr, region, region + w) for region in range(row.start, row.end-w, w)]).reshape(-1)
            # 
            bw_handle = pyBigWig.open(bw)
            # get window score into 1-D array 
            # region_count = np.array([bw_handle.stats(row.chrom, region, region+w) for region in range(row.start, row.end-w, w)]).reshape(-1) # by default, get average value
            # from bigwig to get max read depth within a window region
            # 1. Read in parallel with threads
            region_max_count = Parallel(n_jobs=n, backend="threading")(
                delayed(bw_max)(chr, region, region + w, bw_handle) for region in range(row.start, row.end-w, w)
            )
            region_max_count = np.array(region_max_count).reshape(-1)
            # region_max_count = np.array([bw_handle.stats(row.chrom, region, region+w, type = "max") for region in range(row.start, row.end-w, w)]).reshape(-1)
            # replace any None as 0
            region_count[region_count == None] = 0 
            region_max_count[region_max_count == None] = 0
            # maximum array should be integer only 
            region_max_count = np.round(region_max_count).astype(int)
            chr_count_list.append(region_count); del region_count
            chr_max_list.append(region_max_count); del region_max_count
        #print(f"Concatenate values for {chr} ...")
        chr_region = np.hstack(chr_region_list).reshape(-1,1)
        chr_count = np.hstack(chr_count_list).reshape(-1,1)
        chr_max = np.hstack(chr_max_list).reshape(-1,1)
        # for the scanned value, assuming two types of read counts (more reads, less reads)
        assert len(chr_count) == len(chr_max) and len(chr_region) == len(chr_count), print("Retrieve value for mean and max differ in length.")
        #print("Perform KMeans clustering for two components ...")
        kmeans = KMeans(n_clusters=2, n_init = 10)
        kmeans.fit(np.log10(chr_count+1))
        # find the window have more reads 
        if kmeans.cluster_centers_[0] > kmeans.cluster_centers_[-1]:
            condition = (kmeans.labels_ == 0).reshape(-1,1)
        else:
            condition = (kmeans.labels_ == 1).reshape(-1,1)
        #print("Find cluster of counts that is larger ...")
        chr_region_select = chr_region[condition]
        chr_max_select = chr_max[condition]
        chr_count_select = chr_count[condition]
        # fit a negative binomial distribution (use maximum array, all are integer)
        #print("Fit negative binomial distribution on maximum window count ...")
        fit_nb = fit_negative_binomial(chr_count_select)
        size = fit_nb["n"]; p = fit_nb["p"]; mu = fit_nb["mu"] # nbinomial parameters 
        uniq_count = np.unique(chr_count_select)
        # report the p-value for the read count within window 
        #print("Collect p-value for maximum count per window ...")
        p_val = stats.nbinom.sf(uniq_count-1, size, p)
        cnt_p_df = pd.DataFrame({"count":uniq_count, "p":p_val})
        #print("Perform multiple test correction ...")
        cnt_p_df["padj"] = multipletests(p_val, alpha=0.05, method='fdr_bh')[1]
        # get the count cutoff at padj < 0.05 
        #print("Set cutoff at padj < 0.1 ...")
        # cnt_min = cnt_p_df.query("padj < @padj_cutoff")["count"].min()
        # make array container same length as itself fill by nan 
        #print("Get neighbor value for each window ...")
        neighbor_cnt = np.full_like(chr_count, np.nan); neighbor_max = np.full_like(chr_count, np.nan) #; epsilon = 1e-8
        neighbor_cnt[1:-1] = ((chr_count[:-2] + chr_count[2:])/2); neighbor_max[1:-1] = ((chr_max[:-2] + chr_max[2:])/2) # import tiny amount value to avoid 0 signal in neighborning region
        #region_fc[1:-1] = chr_count[1:-1]/neighbor_cnt; region_fc_max[1:-1] = chr_max[1:-1]/neighbor_max
        #print("Write into new file ...")
        neighbor_cnt = neighbor_cnt.reshape(-1,1); neighbor_max = neighbor_max.reshape(-1,1)
        if i == 0:
            pd.merge(pd.DataFrame({"chrom":chr, "start":chr_region_select, "end":chr_region_select + w, "count":chr_count_select, "max":chr_max_select, "neighbor_count":neighbor_cnt[condition], "neighbor_max":neighbor_max[condition]}), cnt_p_df, on = "count", how = "left").to_csv(output+ ".txt", sep = "\t", header = True, index = False, float_format="%.2e", mode = "w")
        else:
            pd.merge(pd.DataFrame({"chrom":chr, "start":chr_region_select, "end":chr_region_select + w, "count":chr_count_select, "max":chr_max_select, "neighbor_count":neighbor_cnt[condition], "neighbor_max":neighbor_max[condition]}), cnt_p_df, on = "count", how = "left").to_csv(output+ ".txt", sep = "\t", header = False, index = False, float_format="%.2e", mode = "a")
        i += 1 
        # filter condition (count using maximum value, fc using mean value 
        # filter_condition = (chr_max >= cnt_min) & (region_fc > 1.5)
        #print(region_fc.shape, filter_condition.shape, region_ar.shape)
        # enrich_site = chr_region[filter_condition]
        #enrich_mean_cnt = chr_count[filter_condition]
        #enrich_max_cnt = chr_max[filter_condition]
        #enrich_p_dict = {row.count:(row.p, row.padj) for row in cnt_p_df.itertuples()}
        # score as the -log10 transformed padj
        # enrich_p = np.array([enrich_p_dict[cnt] for cnt in chr_max[condition]])
        # add raw data into temp file 
        # , "p":[enrich_p_dict[cnt][0] for cnt in chr_max[condition]], "padj":[enrich_p_dict[cnt][1] for cnt in chr_max[condition]]
        # add enrichment df using artifical cutoffs 
    #     enrich_fc = region_fc[filter_condition]
    #     # enrich_score = enrich_p * enrich_fc 
    #     enrich_df = pd.DataFrame({"chrom":chr, "start":enrich_site, "end":enrich_site + w, "max":enrich_max_cnt, "mean":enrich_mean_cnt, "fc":enrich_fc})
    #     enrich_df = pd.merge(enrich_df, cnt_p_df, left_on = "max", right_on = "count").drop("count", axis = 1)
    #     enrich_list.append(enrich_df)
    # enrich_all_df = pd.concat(enrich_list, axis = 0)
    # if blacklist is not None: 
    #     # remove blacklist 
    #     blacklist_df = bf.read_table(blacklist, schema = "bed3")
    #     enrich_all_df = rmblacklist(enrich_all_df, blacklist_df)
    # enrich_all_df.to_csv(output+ ".txt", sep = "\t", header = True, index = False, float_format="%.2e")
    # print("Add score to enrichment regions ...")
    # enrich_all_df["score"] = (-np.log10(enrich_all_df["padj"])) * enrich_all_df["fc"]
    # max_score = enrich_all_df[enrich_all_df["score"] < np.inf]["score"].max()
    # enrich_all_df = enrich_all_df.replace(np.inf, max_score)
    # enrich_all_df[["chrom", "start", "end", 'score']].to_csv(output+ ".bedGraph", sep = "\t", header = False, index = False)
    print(timeit(starttime))

def fit_negative_binomial(data):
    """
    Fit Negative Binomial distribution to count data using MLE.
    Returns estimated (n, p, mu).
    """
    data = np.asarray(data)
    # Negative log-likelihood
    def nbinom_nll(params):
        n, p = params
        if n <= 0 or not (0 < p < 1):
            return np.inf
        return -np.sum(stats.nbinom.logpmf(data, n, p))
    # Initial guesses (method of moments)
    mean_guess = np.mean(data)
    var_guess = np.var(data)
    n0 = (mean_guess**2) / (var_guess - mean_guess) if var_guess > mean_guess else 10
    p0 = n0 / (n0 + mean_guess)
    init_params = [max(n0, 1e-2), min(max(p0, 1e-2), 0.99)]
    # Optimize likelihood
    res = minimize(nbinom_nll, init_params, bounds=[(1e-6, None), (1e-6, 1-1e-6)])
    n_est, p_est = res.x
    mu_est = n_est * (1 - p_est) / p_est
    return {"n": n_est, "p": p_est, "mu": mu_est}


if __name__ == "__main__":
    main()

