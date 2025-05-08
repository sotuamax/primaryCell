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
from datetime import datetime 
from utilities.misc import ignore_warning
ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="\
    peak_interaction.py \
    -ref hg38\
    -cool sample.mcool \
    -peak sample.narrowPeak\
    -o sample")
    parser.add_argument("-ref", "--reference", choices = ["hg38", "mm39"], default = "hg38", help = "reference genome")
    parser.add_argument("-cool", "--cool", required = True, help = "cool input file")
    parser.add_argument("-chrom", "--chromosome", required = False, nargs = "+", help = "chromosome to examine")
    #parser.add_argument("-region", "--region", required = False, help = "chromosome region to examine in a bed file (no header line)")
    parser.add_argument("-r", "--resolution", default = 1000, type = int, help = "resolution in bp for cool file (default: 1000)")
    parser.add_argument("-o", "--output", help = "output prefix to use")
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    # in the mode of dectection, input peak file to generated paired peaks and report significant interactions
    sub_parsers = parser2.add_subparsers(dest = "command", help = "subcommand to run")
    discovery = sub_parsers.add_parser("detect", help = "identify significant interactions", add_help = False, parents = [parser]) 
    discovery.add_argument("-peak", "--peak", help = "peak file in narrowPeak or bed format.")
    discovery.add_argument("-screen_0s", "--screen_0s", default = 50, type = int, help = "number of 0s encountered before give up (default: 50)")
    discovery.add_argument("-w", "--window", type = int, nargs = "+", help = "number of windows (resolution is per window size) used for expected value calculation on local scale")
    discovery.add_argument("-min_contact", "--min_contact", default = 8, type = int, help = "minimum contacts for each pair combination")
    discovery.add_argument("-fdr", "--fdr", default = 0.05, type = float, help = "false discovery rate (FDR) used for BH multiple testing correction (default: 0.05)")
    # in the mode of count, perform read count at the queried interaction site
    count = sub_parsers.add_parser("count", help = "read count at the interaction site", add_help = False, parents = [parser])
    count.add_argument("-interaction", "--interaction", required = False, help = "paired anchors in bedpe format (minimum columns: chrom1, start1, end1, chrom2, start2, end2), no header")
    args=parser2.parse_args()
    return args

def lambda_group(obs_exp_df, window:list):
    """
    Based on expected count, to group them into lambda group using a base-value based logspace. """
    BASE_pos = 2**(1/3)
    for w in window: 
        exp_max = obs_exp_df[f"exp_{w}"].max();bin_end = np.ceil(3*(np.log2(exp_max))).astype(int)
        ledges = np.logspace(start=0, stop= bin_end, num = bin_end, base=BASE_pos, dtype = np.float64)
        lbins = pd.cut(x = obs_exp_df[f"exp_{w}"], bins = ledges)
        obs_exp_df[f"lambda_{w}"] = [b.right if not pd.isna(b) else 1 for b in lbins]
    return obs_exp_df 

def poisson_p(obs_exp_df, window:list):
    """
    Perform poisson test given observed counts and the expected counts
    input: 
    """
    from scipy.stats import poisson
    for w in window:
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

def bin2peak(obs_exp_enrich, peak_df, resolution):
    """
    Combine bins with original peak info
    Inputs: 
    obs_exp_enrich: at bin scale, the observed/expected contacts;
    peak_df: original input narrowPeak file 
    Returns: 
    dataframe of contact peaks that matched to the bin
    """
    # obs_exp_enriched = obs_exp_enrich.drop([col for col in obs_exp_enrich.columns if "lambda" in col or "q_" in col or "p_" in col], axis = 1).copy()
    peak_df["start1"] = peak_df["start"]//resolution
    peak_df["end1"] = peak_df["end"]//resolution+1
    peak_df["start2"] = peak_df["start"]//resolution
    peak_df["end2"] = peak_df["end"]//resolution+1
    peak_df["chrom1"] = peak_df["chrom"]
    peak_df["chrom2"] = peak_df["chrom"]
    # 
    left_merge = pd.merge(obs_exp_enrich, peak_df[["chrom1", "start1", "end1", "start", "end"]], on = ["chrom1", "start1", "end1"], how = "left")
    left_merge["start1"] = left_merge["start"]
    left_merge["end1"] = left_merge["end"]
    left_merge.drop(["start", "end"], axis = 1, inplace = True)
    right_merge = pd.merge(left_merge, peak_df[['chrom2', "start2", "end2", "start", "end"]], on = ["chrom2", "start2", "end2"], how = "left")
    right_merge["start2"] = right_merge["start"]
    right_merge["end2"] = right_merge["end"]
    right_merge.drop(["start", "end"], axis = 1, inplace = True)
    return right_merge

def longrange(peak_contact):
    """
    transfer significant contacts into long-range format (chr1 1000 2000 chr1:5000-6000,30)
    Input: 
    sig_contact: contacts pass statistical test
    Return:
    standard long-range format of significant contacts (a pair of interacting loc has two records)
    """
    peak_contact = peak_contact.copy()
    peak_contact["B"] = peak_contact["chrom2"].astype(str) + ":" + peak_contact["start2"].astype(str) + "-" + peak_contact["end2"].astype(str) + "," + peak_contact["obs"].astype(str)
    peak_contact["A"] = peak_contact["chrom1"].astype(str) + ":" + peak_contact["start1"].astype(str) + "-" + peak_contact["end1"].astype(str) + "," + peak_contact["obs"].astype(str)
    peak_longrange = peak_contact[["chrom1", "start1", "end1", "B"]].copy()
    # peak_A = peak_contact[['chrom', "start2", "end2", "A"]].copy()
    peak_longrange.columns = ["chrom", "start", "end", "link"] 
    # .rename(columns = {"start1":"start", "end1":"end", "B":"link"}, inplace = True)
    # peak_A.rename(columns = {"start2":"start", "end2":"end", "A":"link"}, inplace = True)
    # peak_longrange = pd.concat([peak_A, peak_B], axis = 0)
    return peak_longrange

def main():
    start = time.time()
    from mpi4py import MPI
    from mpi4py.util import pkl5
    comm = pkl5.Intracomm(MPI.COMM_WORLD) # to overcome the overflow error when comm data > 2 GB
    rank = comm.Get_rank()
    size = comm.Get_size()
    # arguments 
    args = args_parser()
    ref = args.reference 
    cool = args.cool
    resolution = args.resolution
    output = args.output
    if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Parse cool at assigned resolution ...")
    from utilities.cool_tools import parser_cool
    clr = parser_cool(cool, resolution = resolution)
    if args.command == "detect":
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Run in detection mode ...")
        peak = args.peak
        stop_0s = args.screen_0s
        window = args.window
        # parse chromosome arms for reference 
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Process ref chrom ...")
            #if args.region is None:
            from utilities.chrom_func import chrom_arms 
            ref_arm = chrom_arms(ref)
            #else:
            #    ref_arm = bf.read_table(args.region, schema = "bed3")
            if args.chromosome is not None:
                ref_arm = ref_arm[ref_arm["chrom"].isin(args.chromosome)]
            from utilities.peak_tools import peak_parse
            all_chrom = sorted(set(ref_arm["chrom"]))
            print(datetime.now().strftime("%H:%M:%S") + " - Collect peaks ...")
            peak_df = peak_parse(peak, chrom = all_chrom)
            all_chrom = sorted(set(peak_df["chrom"]))
            ref_arm = ref_arm.query("chrom in @all_chrom").copy()
            ref_arm.index = range(len(ref_arm))
            # print(ref_arm)
            ref_arm = [(row.chrom, row.start, row.end) for row in ref_arm.itertuples()] # to a list
        else:
            ref_arm = None 
            peak_df = None 
        ref_arm = comm.bcast(ref_arm, root = 0)
        peak_df = comm.bcast(peak_df, root = 0)
        assert len(peak_df) > 0, 'No peak detected.'
        # split the job by chromosome arms
        if rank == 0:
            w_idx = 0
            worker_tasks1 = {w:[] for w in range(size)}
            for region in ref_arm:
                worker_tasks1[w_idx].append(region)
                w_idx = (w_idx + 1) % size
        else:
            worker_tasks1 = None
        worker_tasks1 = comm.bcast(worker_tasks1, root = 0)
        from utilities.peak_tools import generate_paired_peak 
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Generate paired anchors per chrom ...")
        chrom_region_peak_pair = {region:generate_paired_peak(peak_df, region = region) for region in worker_tasks1[rank]}
        chrom_region_peak_pair_gather = comm.gather(chrom_region_peak_pair, root = 0)
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Gather paired anchors ...")
            peak_pair = dict()
            for c in chrom_region_peak_pair_gather:
                peak_pair.update(c)
            # peak_pair = pd.concat(chrom_region_peak_pair_gather, axis = 0)
        else:
            peak_pair = None 
        peak_pair = comm.bcast(peak_pair, root = 0)
        del worker_tasks1; del chrom_region_peak_pair
        
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Distribute workload ...")
            worker_tasks = {w:{region:[] for region in ref_arm} for w in range(size)}
            w_idx = 0
            for region in peak_pair:
                for anchor1 in peak_pair[region]:
                    worker_tasks[w_idx][region].append((anchor1, peak_pair[region][anchor1]))
                    w_idx = (w_idx + 1) % size
            del peak_pair
        else:
            worker_tasks = None
        try:
            worker_tasks = comm.bcast(worker_tasks, root = 0) 
            # be aware Overflow error. when comm object is over 2 Gb size limit.
        except Exception as e:
            print(e)
            exit(1)
        #
        rank_list = list()
        for arm_region in tqdm(worker_tasks[rank], desc = f"Core:{rank+1}/{size}", position = 0):
            # arm_region stores (chrom, region-start, region-end) # it is also a key for a dictionary
            clr_matrix = clr.matrix(balance = False).fetch(arm_region, arm_region)
            clr_offset = arm_region[1]//resolution
            from utilities.cool_tools import obs_interaction
            # get observation counts 
            # b2_anchor_array range object
            rank_chr_obs = [obs_interaction(clr_matrix, peak_df, clr_offset, arm_region[0], bin1_anchor, b2_anchor_array, stop_0s, resolution) for bin1_anchor,b2_anchor_array in tqdm(worker_tasks[rank][arm_region], desc = f"Obs:{arm_region}", position = 1,leave=False)]
            if len(rank_chr_obs) > 0:
                obs_chr_rank = pd.concat(rank_chr_obs, axis = 0, ignore_index = True)
                if len(obs_chr_rank) > 0:
                    from utilities.cool_tools import local_exp_interaction
                    # get expected counts for window size 'w'
                    exp_chr_rank = pd.concat([local_exp_interaction(clr_matrix, clr_offset, group_region = (a1, b1, a2, b2), wl = window) for (a1, b1, a2, b2), group_df in tqdm(obs_chr_rank.groupby(["start1", "end1", "start2", "end2"]), desc = f"Exp:{arm_region}", position = 2, leave = False)], axis = 0, ignore_index = True)
                    obs_chr_rank = pd.merge(obs_chr_rank, exp_chr_rank, on = ["start1", "end1", "start2", "end2"])
                rank_list.append(obs_chr_rank)
        del clr_matrix
        obs_exp_rank_df = pd.concat(rank_list, axis = 0, ignore_index = True)
        all_rank_df = comm.gather(obs_exp_rank_df, root = 0)
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Collect interactions ...")
            obs_exp_df = pd.concat(all_rank_df, axis = 0, ignore_index = True) # collect all counts (obs/exp)
            obs_exp_df = lambda_group(obs_exp_df, window) # group exp counts into lambda group
            obs_exp_p = poisson_p(obs_exp_df, window) # perform poisson test on obs/lambda pair 
            obs_exp_padj = bh_correction(obs_exp_p, window, args.fdr)
            obs_exp_padj.sort_values(by = ["chrom1", "start1", "start2"], ascending = [True, True, True], inplace = True, ignore_index = True)
            peak_padj = bin2peak(obs_exp_padj, peak_df, resolution)
            peak_padj.to_csv(output + "_padj.txt", sep = "\t", header = True, index = False, float_format='%.2f')
            fdr_w = [f"fdr_{w}" for w in window]
            obs_exp_enrich = obs_exp_padj[obs_exp_padj[fdr_w].sum(axis = 1) == 0].query("obs >= @args.min_contact")
            # obs_exp_enrich.to_csv(output + "_padj_sig.txt", sep = "\t", header = True, index = False, float_format='%.2f')
            print(datetime.now().strftime("%H:%M:%S") + " - Finished!")
            obs_exp_peak = bin2peak(obs_exp_enrich, peak_df, resolution)
            obs_exp_peak.to_csv(output + "_padj_sig.txt", sep = "\t", header = True, index = False, float_format='%.2f')
            longrange(obs_exp_peak).to_csv(output + ".longrange", sep = "\t", index = False, header = False)
        print(timeit(start), f"on core {rank}")
        exit(0)
    if args.command == "count":
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Read interaction sites ...")
        interaction = args.interaction 
        ### parsing in 10M matrix 
        interaction_df = pd.read_table(interaction, sep = "\t", header = None)
        interaction_df = interaction_df.iloc[:, :6] 
        interaction_df.columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
        segment = 10_000_000
        interaction_df["bin1"] = interaction_df["start1"]//segment 
        interaction_df["bin2"] = interaction_df["start2"]//segment+1
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Distribute workload ...")
            worker_tasks = {w:[] for w in range(size)}
            w_idx = 0
            for (chr,b1,b2), b_df in interaction_df.groupby(["chrom1", "bin1", "bin2"]):
                worker_tasks[w_idx].append((chr, b1, b2))
                w_idx = (w_idx + 1) % size
        else:
            worker_tasks = None 
        worker_tasks = comm.bcast(worker_tasks, root = 0)
        # 
        from utilities.cool_tools import interaction_count
        rank_list = list()
        for (chr, b1, b2) in tqdm(worker_tasks[rank], desc = f"Core:{rank+1}/{size}"):
            interaction_query = interaction_df.query("chrom1 == @chr & bin1 == @b1 & bin2 == @b2")[["chrom1", "start1", "end1", "chrom2", "start2", "end2"]].copy()
            rstart = interaction_query["start1"].min()//resolution*resolution
            rend = (interaction_query["end2"].max()//resolution+1)*resolution
            arm_region = (chr, rstart, rend)
            # print(arm_region)
            clr_matrix = clr.matrix(balance = False).fetch(arm_region, arm_region)
            # print(clr_matrix)
            clr_offset = arm_region[1]//resolution
            # print(clr_offset)
            # print(interaction_query)
            # exit(1)
            interaction_cnt = interaction_count(clr_matrix, clr_offset, resolution, interaction_query)
            interaction_query["count"] = interaction_cnt
            rank_list.append(interaction_query)
        rank_df = pd.concat(rank_list, axis = 0)
        rank_all = comm.gather(rank_df, root = 0)
        if rank == 0:
            print(datetime.now().strftime("%H:%M:%S") + " - Collect obs. count ...")
            rank_df_all = pd.concat(rank_all, axis = 0)
            rank_df_all.to_csv(output + ".txt", sep = "\t", header = True, index = False)
            print(datetime.now().strftime("%H:%M:%S") + " - Finished. ")
        print(timeit(start), f"on core {rank}")
        exit(0)

if __name__ == "__main__":
    main()
