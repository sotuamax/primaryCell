#!/usr/bin/env python
import bioframe as bf 
import numpy as np 
import pandas as pd
import argparse 
import logging 
import time 
from utilities.misc import timeit 
import sys
from tqdm.auto import tqdm

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="\
    peak_interaction.py \
    -ref hg38\
    -cool sample.mcool \
    -peak sample.narrowPeak\
    -chrom chr8,chr19\
    -o sample")
    parser.add_argument("-ref", "--reference", choices = ["hg38", "mm39"], help = "reference genome")
    parser.add_argument("-cool", "--cool", help = "cool input file")
    parser.add_argument("-peak", "--peak", help = "peak file in narrowPeak/bed format. Pre-filter peaks based on fc and score is recommended.")
    parser.add_argument("-chrom", "--chromosome", required = False, nargs = "+", help = "chromosome to examine")
    parser.add_argument("-blacklist", "--blacklist", required = False, help = "blacklist region of the reference genome.")
    parser.add_argument("-region", "--region", required = False, help = "chromosome region to examine in a bed file (no header line)")
    parser.add_argument("-r", "--resolution", default = 1000, type = int, help = "resolution in bp for cool file (default: 1000)")
    parser.add_argument("-w", "--window", type = int, nargs = "+", help = "number of windows (resolution is per window size) used for expected value calculation on local scale")
    parser.add_argument("-screen_0s", "--screen_0s", default = 50, type = int, help = "number of 0s encountered before give up (default: 50)")
    parser.add_argument("-global", "--global", action = "store_true", help = "Caculate global expected count")
    # parser.add_argument("-min_contact", "--min_contact", default = 3, type = int, help = "minimum contacts for each pair combination")
    parser.add_argument("-fdr", "--fdr", default = 0.1, type = float, help = "false discovery rate (FDR) used for BH multiple testing correction")
    parser.add_argument("-pairs", "--pairs", required = False, help = "pairs file used to refine the interaction site")
    parser.add_argument("-o", "--output", help = "output prefix to use")
    parser.add_argument("-t", "--threads", default = 8, type = int, help = "number of threads to use")
    args=parser.parse_args()
    return args

def obs_filter(bp_obs, min_contact):
    """
    filter observed contact datafrmae by minimum contacts
    Input: 
    bp_obs: dataframe of observed contact 
    Return: 
    dataframe: filtered contacts with minimum contact given"""
    bp_obs_filter = bp_obs.query(f"obs >= {min_contact}").copy()
    return bp_obs_filter

def poisson_p(obs_exp_df, window:list):
    """
    Perform poisson test given observed counts and the expected counts
    input: 
    """
    BASE_pos = 2**(1/3)
    BASE_neg = 2**(-1/3)
    from scipy.stats import poisson 
    for w in window:
        exp_max = obs_exp_df[f"exp_{w}"].max();bin_end = np.ceil(3*(np.log2(exp_max))).astype(int)
        ledges = np.concatenate(([0], np.logspace(start = 5, stop = 1, num = 5, base = BASE_neg, dtype = np.float64), np.logspace(start=0, stop= bin_end, num = bin_end, base=BASE_pos, dtype = np.float64)))
        lbins = pd.cut(x = obs_exp_df[f"exp_{w}"], bins = ledges)
        obs_exp_df[f"lambda_{w}"] = [b.right if not pd.isna(b) else 0 for b in lbins]
        obs_exp_pair = set(zip(obs_exp_df["obs"], obs_exp_df[f"lambda_{w}"]))
        o, e = zip(*obs_exp_pair)
        oe_p = poisson.sf(np.array(o), np.array(e))
        df_eop = pd.DataFrame(); df_eop["obs"] = o; df_eop[f"lambda_{w}"] = e; df_eop[f"p_{w}"] = oe_p
        obs_exp_df = pd.merge(obs_exp_df, df_eop, on = ["obs", f"lambda_{w}"])
    return obs_exp_df

def bh_correction(obs_exp_p, window:list, fdr:float):
    # from scipy.stats import false_discovery_control 
    for w in window:
        obs_exp = list()
        for lbin, lbin_df in obs_exp_p.groupby(f"lambda_{w}"):
            lbin_df[f"q_{w}"] = lbin_df[f"p_{w}"].rank()*fdr/len(lbin_df)
            obs_exp.append(lbin_df)
        obs_exp_p = pd.concat(obs_exp, axis = 0)
    for w in window:
        obs_exp_p[f"fdr_{w}"] = np.where(obs_exp_p[f"p_{w}"] <= obs_exp_p[f"q_{w}"], 0, 1)
    return obs_exp_p

def bin2peak(obs_exp_enrich, peak_df):
    """
    Combine bins with original peak info
    Inputs: 
    obs_exp_enrich: at bin scale, the observed/expected contacts;
    peak_df: original input narrowPeak file 
    Returns: 
    dataframe of contact peaks that matched to the bin
    """
    obs_exp_enriched = obs_exp_enrich.drop([col for col in obs_exp_enrich.columns if "lambda" in col or "q_" in col or "p_" in col], axis = 1).copy()
    if "fc" in peak_df.columns:
        left_merge = pd.merge(obs_exp_enriched, peak_df[['chrom', "start", "end", "bin", "fc"]], left_on = ["chrom", "bin1"], right_on = ["chrom", "bin"]).drop("bin", axis = 1)
        right_merge = pd.merge(left_merge, peak_df[['chrom', "start", "end", "bin", "fc"]], left_on = ["chrom", "bin2"], right_on = ["chrom", "bin"]).drop("bin", axis = 1)
        right_merge.rename(columns = {"start_x":"start1", "end_x":"end1", "start_y":"start2", "end_y":"end2", "fc_x":"fc1", "fc_y":"fc2"}, inplace = True)
    else:
        left_merge = pd.merge(obs_exp_enriched, peak_df[["chrom", "start", "end", "bin"]], left_on = ["chrom", "bin1"], right_on = ["chrom", "bin"]).drop("bin", axis = 1)
        right_merge = pd.merge(left_merge, peak_df[['chrom', "start", "end", "bin"]], left_on = ["chrom", "bin2"], right_on = ["chrom", "bin"]).drop("bin", axis = 1)
        right_merge.rename(columns = {"start_x":"start1", "end_x":"end1", "start_y":"start2", "end_y":"end2"}, inplace = True)
    return right_merge

def longrange(peak_contact):
    """
    transfer significant contacts into long-range format (chr1 1000 2000 chr1:5000-6000,30)
    Input: 
    sig_contact: contacts pass statistical test
    Return:
    standard long-range format of significant contacts (a pair of interacting loc has two records)
    """
    peak_contact["B"] = peak_contact["chrom"].astype(str) + ":" + peak_contact["start2"].astype(str) + "-" + peak_contact["end2"].astype(str) + "," + peak_contact["obs"].astype(str)
    peak_contact["A"] = peak_contact["chrom"].astype(str) + ":" + peak_contact["start1"].astype(str) + "-" + peak_contact["end1"].astype(str) + "," + peak_contact["obs"].astype(str)
    peak_B = peak_contact[["chrom", "start1", "end1", "B"]].copy()
    peak_A = peak_contact[['chrom', "start2", "end2", "A"]].copy()
    peak_B.rename(columns = {"start1":"start", "end1":"end", "B":"link"}, inplace = True)
    peak_A.rename(columns = {"start2":"start", "end2":"end", "A":"link"}, inplace = True)
    peak_longrange = pd.concat([peak_A, peak_B], axis = 0)
    return peak_longrange

def main():
    start = time.time()
    args = args_parser()
    logging.basicConfig(level = logging.INFO, filename = f"{args.output}.log", filemode = "w", format = '%(message)s', datefmt = "")
    logging.info(f"Reference: {args.reference}")
    logging.info(f"Resolution: {args.resolution}")
    logging.info(f"Stop sign on #0s: {args.screen_0s}")
    ref = args.reference 
    cool = args.cool
    peak = args.peak
    resolution = args.resolution
    stop_0s = args.screen_0s
    window = args.window
    #
    if args.region is None:
        from utilities.chrom_func import chrom_arms 
        ref_arm = chrom_arms(ref)
    else:
        ref_arm = bf.read_table(args.region, schema = "bed3")
    if args.chromosome is not None:
        ref_arm = ref_arm.query("chrom in @args.chromosome")
        all_chrom = args.chromosome
    else:
        all_chrom = sorted(set(ref_arm["chrom"]))
    ### 
    #logging.info(ref_arm)
    # parse cool file 
    from utilities.cool_tools import parser_cool
    clr = parser_cool(cool, resolution = args.resolution)
    # parse peak file
    from utilities.peak_tools import peak_parse
    peak_df = peak_parse(peak, chrom = all_chrom)
    # update chrom list and ref_arm based on chroms in peak_df
    all_chrom = sorted(set(peak_df["chrom"]))
    ref_arm = ref_arm.query("chrom in @all_chrom").copy()
    ref_arm.index = range(len(ref_arm))
    # remove blacklist region
    if args.blacklist is not None:
        print("Remove peaks overlapping blacklist region")
        from utilities.peak_tools import rmblacklist
        blacklist_df = bf.read_table(args.blacklist, schema = "bed3")
        peak_df = rmblacklist(peak_df, blacklist_df)
    logging.info(f"Peaks #: {len(peak_df)}")
    if len(peak_df) == 0:
        print('no peak detected.')
        exit(1)
    # split the job by chromosome arms
    chrom_paired_peak_dict = {}
    for row in ref_arm.itertuples():
        region = (row.chrom, row.start, row.end) 
        from utilities.peak_tools import generate_paired_peak 
        # generate dict of (chr, start, end):{bin1:[bin2_list]}
        chrom_paired_peak_dict[region] = generate_paired_peak(peak_df, resolution, region)
    for arm_region in chrom_paired_peak_dict:
        clr_matrix = clr.matrix(balance = False).fetch(arm_region, arm_region)
        clr_offset = row.start//resolution
        for b1 in chrom_paired_peak_dict[arm_region]:
            obs_interaction(clr_matrix, clr_offset, arm_region[0], b1, chrom_paired_peak_dict[arm_region][b1], stop_0s, resolution)



    worker_tasks = {w:{(row.chrom, row.start, row.end):[] for row in ref_arm.itertuples()} for w in range(args.threads)}
    w_idx = 0
    for chrom in chrom_paired_peak_dict:
        for bin1 in chrom_paired_peak_dict[chrom]:
            worker_tasks[w_idx][chrom].append(bin1)
            w_idx = (w_idx + 1) % args.threads

    Process(target = obs_interaction, args = ())
    for row in tqdm(ref_arm.itertuples(), desc = "chromosomes", position = 0):
        arm_region = (row.chrom, row.start, row.end)
        # print(arm_region)
        # logging.info(f"Processing {arm_region} ")
        clr_matrix = clr.matrix(balance = False).fetch(arm_region, arm_region)
        # global expected counts (chrom-arm count mean)
        # ref_exp_global[arm_region] = global_exp(clr_matrix)
        # print(clr_matrix.shape)
        clr_offset = row.start//resolution
        ### get observed contacts 
        logging.info("Retrieve obs. counts ...")
        from utilities.cool_tools import obs_interaction
        # rank_chr_obs is a list of dataframes
        rank_chr_obs = [obs_interaction(clr_matrix, clr_offset, row.chrom, bin1_pos, chrom_paired_peak_dict[arm_region][bin1_pos], stop_0s, resolution) for bin1_pos in tqdm(worker_tasks[rank][arm_region], desc = f"Obs:core{rank}/{size}", position = 1)]
        # print(f"RRR {rank}", rank_chr_obs)
        # print(f"Retrieve obs/exp counts on {arm_region} with {len(rank_chr_obs)} rows")
        if len(rank_chr_obs) > 0:
            # print("++++++++++++")
            obs_chr_rank = pd.concat(rank_chr_obs, axis = 0, ignore_index = True)
            # logging.info("Estimate exp. counts ...")
            for w in window:
                from utilities.cool_tools import local_exp_interaction
                exp_chr_rank = local_exp_interaction(clr_matrix, clr_offset, obs_chr_rank, w)
                obs_chr_rank = pd.merge(obs_chr_rank, exp_chr_rank, on = ["bin1", "bin2"])
            # print(obs_chr_rank)
            rank_list.append(obs_chr_rank)
    del clr_matrix
    obs_exp_rank_df = pd.concat(rank_list, axis = 0, ignore_index = True)
    obs_exp_rank_df.to_csv(f"{args.output}_{rank}.tmp", sep = "\t", header = True, index = False)
    if rank == 0:
        test = 1
    else:
        test = None
    test = comm.bcast(test, root = 0)
    # all_rank_df = comm.gather(obs_exp_rank_df, root = 0)
    # if rank == 0:
    #     pd.concat(all_rank_df, axis = 0).to_csv(f"combined.tmp", sep = "\t", header = True, index = False)
    # exit(1)
    # all_rank_df.to_csv(f"{rank}.tmp", sep = "\t", header = True, index = False)
    if rank == 0:
        df_list = [pd.read_table(f"{args.output}_{i}.tmp", sep = "\t", header = 0) for i in range(size)]
        obs_exp_df = pd.concat(df_list, axis = 0, ignore_index = True)
        try:
            for i in range(size):
                os.remove(f"{args.output}_{i}.tmp")
        except:
            pass
        obs_exp_p = poisson_p(obs_exp_df, window)
        # all_rank_df = [a for a in all_rank_df if a is not None]
        # obs_exp_p = pd.concat(all_rank_df, axis = 0, ignore_index = True)
        # order obs_exp_p dataframe by 
        obs_exp_p.to_csv(args.output + "_p.txt", sep = "\t", header = True, index = False)
        obs_exp_padj = bh_correction(obs_exp_p, window, args.fdr)
        obs_exp_padj.to_csv(args.output + "_padj.txt", sep = "\t", header = True, index = False)
        fdr_w = [f"fdr_{w}" for w in window]
        obs_exp_enrich = obs_exp_padj[obs_exp_padj[fdr_w].sum(axis = 1) == 0].query("obs >= 3")
        obs_exp_enrich.to_csv(args.output + "_padj_sig.txt", sep = "\t", header = True, index = False)
        # obs_exp_fail = obs_exp_padj[obs_exp_padj[["fdr_10", "fdr_20", "fdr_50"]].sum(axis = 1) != 0]
        peak_df["bin"] = peak_df["summit"]//resolution
        obs_exp_peak_dup = bin2peak(obs_exp_enrich, peak_df)
        if "fc" in peak_df.columns:
            obs_exp_peak_dup["fc_total"] = obs_exp_peak_dup[[col for col in obs_exp_peak_dup.columns if "fc" in col]].sum(axis = 1)
            obs_exp_peak_df = obs_exp_peak_dup.loc[obs_exp_peak_dup.groupby(["bin1", "bin2"])["fc_total"].transform("max") == obs_exp_peak_dup["fc_total"]]
        else:
            # replace order of columns 
            obs_exp_peak_df = obs_exp_peak_dup
        obs_exp_peak_df[["chrom", "start1", "end1", "start2", "end2"] + [col for col in obs_exp_peak_df.columns if col not in ["chrom", "start1", "end1", "start2", "end2"]]].to_csv(args.output + "_peak_interaction.txt", sep = "\t", header = True, index = False)
        longrange(obs_exp_peak_df).to_csv(args.output + ".longrange", sep = "\t", index = False, header = False)
    print(timeit(start), f"on core {rank}")
    exit(0)
    

if __name__ == "__main__":
    main()
