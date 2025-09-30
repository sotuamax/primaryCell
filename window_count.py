#!/usr/bin/env python
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

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", help = "file input to examine (only support alignment hg38 for now)")
    parser.add_argument("-mode", "--mode", choices=["bw", "bam"], help = "input file mode to run for score (for input as bigwig, the score must be integer)")
    parser.add_argument("-w", "--window", default = 1000, type = int, help = "window size used to screen genome")
    parser.add_argument("-n", "--threads", default = 1, type = int, help = "threads to use to process bam file")
    parser.add_argument("-blacklist", "--blacklist", required = False, help = "blacklist region in bed file")
    parser.add_argument("-o", "--output", help = "output name")
    args = parser.parse_args()
    return args

def main():
    start = time.time()
    args = args_parser()
    input_file = args.input
    m = args.mode 
    n = args.threads
    w = args.window  
    blacklist = args.blacklist 
    output = args.output
    arm_df = chrom_arms()
    enrich_list = list()
    for row in tqdm(arm_df.itertuples(index = False), total = len(arm_df)):
        chr = row.chrom; start = row.start; end = row.end
        print(chr, start, end)
        if m == "bam":
            bam_handle = pysam.AlignmentFile(input_file, "rb", threads = n)
            region_count = np.array([bam_handle.count(row.chrom, region, region + w) for region in range(row.start, row.end, w)]).reshape(-1,1)
        if m == "bw":
            bw_handle = pyBigWig.open(input_file)
            region_count = np.array([bw_handle.stats(row.chrom, region, region+w, type = "max") for region in range(row.start, row.end-w, w)])
            region_count[region_count == None] = 0 
            region_count = region_count.astype(int).reshape(-1, 1)
        kmeans = KMeans(n_clusters=2, n_init = 10)
        kmeans.fit(np.log10(region_count+1))
        if kmeans.cluster_centers_[0] > kmeans.cluster_centers_[-1]:
            count_cluster = region_count[kmeans.labels_ == 0]
        else:
            count_cluster = region_count[kmeans.labels_ == 1]
        fit_nb = fit_negative_binomial(count_cluster)
        size = fit_nb["n"]; p = fit_nb["p"]; mu = fit_nb["mu"]
        uniq_count = np.unique(count_cluster)
        p_val = stats.nbinom.sf(uniq_count-1, size, p)
        cnt_p_df = pd.DataFrame({"count":uniq_count, "p":p_val})
        try:
            cnt_p_df["padj"] = multipletests(p_val, alpha=0.05, method='fdr_bh')[1]
            cnt_min = cnt_p_df.query("padj < 0.05")["count"].min()
            # all region as array 
            region_ar = np.array(range(start, end, w))
            # region_count reshape to 1 dimension (-1)
            region_count = region_count.reshape(-1)
            region_fc = np.full_like(region_count, np.nan)
            neighbor_cnt = ((region_count[:-2] + region_count[2:] + 1)/2)
            region_fc[1:-1] = region_count[1:-1]/neighbor_cnt
            # filter condition 
            filter_condition = np.where((region_count > cnt_min) & (region_fc > 1.5))
            enrich_site = region_ar[filter_condition]
            enrich_cnt = region_count[filter_condition]
            enrich_p_dict = {row.count:row.padj for row in cnt_p_df.itertuples()}
            # score as the -log10 transformed padj
            enrich_p = (-np.log10(np.array([enrich_p_dict[cnt] for cnt in enrich_cnt])))
            enrich_fc = region_fc[filter_condition]
            enrich_score = enrich_p * enrich_fc 
            enrich_df = pd.DataFrame({"chrom":chr, "start":enrich_site, "end":enrich_site + w, "name":enrich_cnt, "score":enrich_score, "strand":"."})
            enrich_list.append(enrich_df)
        except Exception as e:
            print(cnt_p_df)
            print(e)
    enrich_all_df = pd.concat(enrich_list, axis = 0)
    if blacklist is not None: 
        # remove blacklist 
        blacklist_df = bf.read_table(blacklist, schema = "bed3")
        enrich_all_df = rmblacklist(enrich_all_df, blacklist_df)
    enrich_all_df.to_csv(output+ ".bed", sep = "\t", header = False, index = False)
    enrich_all_df[["chrom", "start", "end", 'score']].to_csv(output+ ".bedGraph", sep = "\t", header = False, index = False)
    timeit(start)

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

