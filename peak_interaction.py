#!/usr/bin/env python
"""
This program performs two modes. One for detection of significant interactions; One for count reads supporting interactions.

mode - detect: 
This mode takes input contact data in cool format, and examine chromosome arm by arm (supports hg38/mm39 for now). 
For interaction candidates, given a peak file and paired peaks are examined for significance as compared to its background contact signals. 
Background region is controlled by window size (for example, examine 5/10/20/50 windows at resolution 1000 surrounding the candidate region), and the expected count is estimated by background average; 
for the estimated expected contact frequencies, they can be grouped into different lambda groups (see lambda_group function); update: 09/30/2025 - further grouped for the range between 0-1 into 0, 0.01, 0.02, etc.
The p-value for significance between observed vs. expected is estimated using poisson distribution (that is at certain distance, an 'lambda' contacts can happen where lambda is the average number of contacts at this distance, we can get the p-value for the observed contacts); 

Paired peaks for examination are screened one by one util certain number of 0s encountered before giving up.  

mode - count: 
This mode takes input contact data in cool format and interaction site in bedpe format for examination, and it reports the total number of reads supporting the interaction site. 

Example (trigger multiple cores using MPI): 
mpiexec -n 10 peak_interaction.py detect -w 10 20 50 100 -cool sample.mcool -r 1000 -peak sample.bed -o sample -min_dist 3000 -min_contact 3 
mpiexec -n 10 peak_interaction.py count -cool sample.mcool -r 1000 -interaction sample.bedpe -o sample 

# when multiple criteria setup, each criteria can be more relaxed. 

# update 10/20/2025: local_exp_interaction (center contacts removed from background calculation)

"""
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
    peak_interaction.py detect -w 5 10 20 50 100 -cool sample.mcool -peak sample.narrowPeak -o sample \n \
    peak_interaction.py count -cool sample.mcool -r 1000 -interaction sample.bedpe -o sample \n")
    parser.add_argument("-ref", "--reference", choices = ["hg38", "mm39"], default = "hg38", help = "reference genome")
    parser.add_argument("-cool", "--cool", required = True, help = "cool input file")
    parser.add_argument("-chrom", "--chromosome", required = False, nargs = "+", help = "chromosome to examine")
    parser.add_argument("-r", "--resolution", default = 1000, type = int, help = "resolution in bp for cool file (default: 1000)")
    parser.add_argument("-o", "--output", help = "output prefix to use")
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    # in the mode of dectection, input peak file to generated paired peaks and report significant interactions
    sub_parsers = parser2.add_subparsers(dest = "command", help = "subcommand to run")
    discovery = sub_parsers.add_parser("detect", help = "identify significant interactions", add_help = False, parents = [parser]) 
    discovery.add_argument("-peak", "--peak", help = "peak file in narrowPeak or bed format." \
    "Note: peaks used as loop anchors have to be pre-processed to make sure no multiple peaks at the same bin region.")
    discovery.add_argument("-screen_0s", "--screen_0s", default = 100, type = int, help = "number of 0s encountered before give up (default: 100)")
    discovery.add_argument("-w", "--window", type = int, nargs = "+", help = "number of windows (resolution is per window size) used for expected value calculation on local scale")
    discovery.add_argument("-obs", "--observed", required = False, help = "table of observed+expected values if available")
    discovery.add_argument("-min_contact", "--min_contact", default = 3, type = int, help = "minimum contacts for each pair combination")
    discovery.add_argument("-score", "--score", default = 20, type = float, help = "minimum score used to filter significant interactions.")
    discovery.add_argument("-min_dist", "--min_distance", default = 0, type = int, help = "minimum distance for the identified loop interactions")
    discovery.add_argument("-max_dist", "--max_distance", default = 10_000_000, type = int, help = "maximum distance for the identified loop interactions")
    # in the mode of count, perform read count at the queried interaction site
    count = sub_parsers.add_parser("count", help = "read count at the interaction site", add_help = False, parents = [parser])
    count.add_argument("-interaction", "--interaction", required = False, help = "paired anchors in bedpe format (minimum columns: chrom1, start1, end1, chrom2, start2, end2), no header")
    args=parser2.parse_args()
    return args

def lambda_group(obs_exp_df, window:list):
    """
    Based on expected count, to group them into lambda group using a base-value based logspace. 
    see: https://github.com/open2c/cooltools/blob/master/cooltools/api/dotfinder.py (function dots)
    """
    # logspace start from negative value can help: for loop enriched with weak background, smaller lambda group can help on its background statistical test (e.g., background as 0.01 can be tested in a group range 0.2). 
    BASE_pos = 2**(1/3)
    # although not all samples share the same exp max value, to keep the logspace the same for all test, use 1000 as default.
    bin_end = 30 # np.ceil(3*(np.log2(1000))).astype(int)
    logspace = np.logspace(start=-5, stop= bin_end, num = bin_end+5+1, base=BASE_pos, dtype = np.float64)
    # add start as -6 (tried -10, may be too small, and value < 1 should be refined).
    logspace = np.insert(logspace, 0, 0) # add 0 at position 0
    for w in window: 
        lbins = pd.cut(x = obs_exp_df[f"exp_{w}"], bins = logspace, include_lowest=True)
        obs_exp_df[f"lambda_{w}"] = [b.right if not pd.isna(b) else logspace[-1] for b in lbins]
    # exp_lambda = pd.cut(x = obs_exp_df["obs"], bins = logspace, include_lowest=True)
    # obs_exp_df["obs_lambda"] = [e.right if not pd.isna(e) else logspace[-1] for e in exp_lambda]
    return obs_exp_df 

def poisson_p(obs_exp_df, window:list):
    """
    Perform poisson test given observed counts and the expected counts
    """
    from scipy.stats import poisson
    for w in window:
        obs_exp_pair = set(zip(obs_exp_df["obs"], obs_exp_df[f"lambda_{w}"]))
        o, e = zip(*obs_exp_pair)
        oe_p = poisson.sf(np.round(np.array(o)), np.array(e))
        df_eop = pd.DataFrame({"obs":o, f"lambda_{w}":e, f"p_{w}":oe_p})
        obs_exp_df = pd.merge(obs_exp_df, df_eop, on = ["obs", f"lambda_{w}"])
    return obs_exp_df

def bh_correction(obs_exp_p, window:list, fdr:float):
    # from scipy.stats import false_discovery_control 
    """
    for each of the expected lambda group, control false discovery rate, and 
    caculate its q value.
    
    Note: 
    when the size of each lambda group is different, for the same p-value, its p-adj can be different. 
    The differences arised from lambda group size can import biases. 
    No correction may be more reasonable. (10/22/2025)
    Use score: enrichment * log10(p-value);
    enrichment is the fold of observed/expected contacts. 
    """
    from statsmodels.stats.multitest import multipletests
    obs_exp = list()
    for _, lbin_df in obs_exp_p.groupby("obs"):
        for w in window:
            _, lbin_df[f"bh_{w}"], _, _ = multipletests(lbin_df[f"p_{w}"], alpha = fdr, method = "fdr_bh")
        obs_exp.append(lbin_df)
    obs_exp_p = pd.concat(obs_exp, axis = 0)
    # for w in window:
    #     obs_exp_p[f"fdr_{w}"] = np.where(obs_exp_p[f"p_{w}"] <= obs_exp_p[f"q_{w}"], 0, 1)
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
    left_merge = pd.merge(obs_exp_enrich, peak_df[["chrom1", "start1", "end1", "start", "end"]], on = ["chrom1", "start1", "end1"], how = "left")
    left_merge["start1"] = left_merge["start"]
    left_merge["end1"] = left_merge["end"]
    left_merge.drop(["start", "end"], axis = 1, inplace = True)
    right_merge = pd.merge(left_merge, peak_df[['chrom2', "start2", "end2", "start", "end"]], on = ["chrom2", "start2", "end2"], how = "left")
    right_merge["start2"] = right_merge["start"]
    right_merge["end2"] = right_merge["end"]
    right_merge.drop(["start", "end"], axis = 1, inplace = True)
    return right_merge

def bedpe(peak_contact):
    """
    Docstring for bedpe
    
    transfer peak contact txt into bedpe format
    """
    peak_contact = peak_contact.copy()
    peak_contact["name"] = "."
    return peak_contact[["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "obs"]]

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
    if rank == 0:
        print(datetime.now().strftime("%H:%M:%S") + " - Process arguments ...")
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
        peak = args.peak
        stop_0s = args.screen_0s
        window = args.window
        max_dist = args.max_distance
        if args.observed is None:
            if rank == 0:
                print(datetime.now().strftime("%H:%M:%S") + " - Run in detection mode ...")
                print(datetime.now().strftime("%H:%M:%S") + " - Process ref chrom ...")
            from utilities.chrom_func import chrom_arms 
            ref_arm = chrom_arms(ref)
            from utilities.peak_tools import peak_parse
            if args.chromosome is not None:
                ref_arm = ref_arm[ref_arm["chrom"].isin(args.chromosome)]
            all_chrom = sorted(set(ref_arm["chrom"]))
            peak_df = peak_parse(peak, chrom = all_chrom)
            assert len(peak_df[peak_df["start"]//resolution != peak_df["end"]//resolution]) == 0, f"peak anchor occupy more than one bin at {resolution} bp"
            # get peak bin occupancy 
            peak_bin_df = pd.DataFrame(zip(peak_df["chrom"].tolist(), (peak_df["start"]//resolution).tolist()), columns = ["chrom", "bin"]) 
            assert len(peak_bin_df.drop_duplicates(keep = "first")) == len(peak_df), f"peak overlapped on bins at {resolution} bp."
            # parse chromosome arms for reference 
            if rank == 0:
                print(datetime.now().strftime("%H:%M:%S") + " - Collect peaks ...")
                all_chrom = sorted(set(peak_df["chrom"]))
                ref_arm = ref_arm.query("chrom in @all_chrom").copy()
                # split the genome region into segments of max_dist*2 (so that each time, loops within max_dist are screened)
                # for the end of the genome, only present once
                ref_arm_list = list()
                for row in ref_arm.itertuples():
                    start = row.start
                    while start <= row.end-max_dist*2:
                        end = start + max_dist*2
                        ref_arm_list.append((row.chrom, start, end))
                        start += max_dist 
                    ref_arm_list.append((row.chrom, start, row.end))
                ref_arm_new = pd.DataFrame(ref_arm_list, columns=["chrom", "start","end"]).sort_values(["chrom", "start"], ignore_index=True)
                #ref_arm.index = range(len(ref_arm))
                ref_arm_list = [(row.chrom, row.start, row.end) for row in ref_arm_new.itertuples()] # to a list
            else:
                ref_arm_list = None 
            ref_arm_list = comm.bcast(ref_arm_list, root = 0)
            assert len(peak_df) > 0, 'No peaks detected.'
            assert ref_arm_list is not None, "chromosome region is not valid."
            # split the job by chromosome arms
            if rank == 0:
                w_idx = 0
                worker_tasks = {w:[] for w in range(size)}
                for region in ref_arm_list:
                    worker_tasks[w_idx].append(region)
                    w_idx = (w_idx + 1) % size
                print(datetime.now().strftime("%H:%M:%S") + " - Distribute chromosome to each core ...")
            else:
                worker_tasks = None
            worker_tasks = comm.bcast(worker_tasks, root = 0)
            from utilities.peak_tools import generate_paired_peak 
            if rank == 0:
                print(datetime.now().strftime("%H:%M:%S") + " - Generate paired anchors per chrom ...")
            try:
                # paired anchors as index format (index of peak_df)
                # generate_paired_peak: return {peak_index_1: range(peak_index_2:max_distance)}
                # for region in worker_tasks[rank]:
                #     t = generate_paired_peak(peak_df, max_dist, region = region)
                chrom_region_peak_pair = {region:generate_paired_peak(peak_df, max_dist, region = region) for region in worker_tasks[rank]}
            except Exception as e:
                print("Error generating paired anchors!")
                print(region)
                exit(1)
            rank_list = list()
            for arm_region in tqdm(worker_tasks[rank], desc = f"Core:{rank+1}/{size}", position = 0):
                # arm_region stores (chrom, region-start, region-end) # it is also a key for a dictionary
                clr_matrix = clr.matrix(balance = False).fetch(arm_region, arm_region)
                clr_offset = arm_region[1]//resolution
                from utilities.cool_tools import obs_interaction
                # get observation counts 
                # b2_anchor_array range object
                # return data frame start/end is resolution normalized
                rank_chr_obs = [obs_interaction(clr_matrix, peak_df, clr_offset, arm_region[0], bin1_anchor, chrom_region_peak_pair[arm_region][bin1_anchor], stop_0s, resolution) for bin1_anchor in tqdm(chrom_region_peak_pair[arm_region], desc = f"Obs:{arm_region}", position = 1,leave=False)]
                if len(rank_chr_obs) > 0:
                    obs_chr_rank = pd.concat(rank_chr_obs, axis = 0, ignore_index = True)
                    if len(obs_chr_rank) > 0:
                        from utilities.cool_tools import local_exp_interaction
                        # get expected counts for window size 'w'
                        # taking inside window average interaction frequencies as an expected interaction
                        # what if the bin range is > 1 bin (cross multiple bins)? 
                        exp_chr_rank = pd.concat([local_exp_interaction(clr_matrix, clr_offset, group_region = (a1, b1, a2, b2), wl = window) for (a1, b1, a2, b2), _ in tqdm(obs_chr_rank.groupby(["start1", "end1", "start2", "end2"]), desc = f"Exp:{arm_region}", position = 2, leave = False)], axis = 0, ignore_index = True)
                        obs_chr_rank = pd.merge(obs_chr_rank, exp_chr_rank, on = ["start1", "end1", "start2", "end2"])
                    rank_list.append(obs_chr_rank)
            del clr_matrix
            obs_exp_rank_df = pd.concat(rank_list, axis = 0, ignore_index = True)
            all_rank_df = comm.gather(obs_exp_rank_df, root = 0)
        if rank == 0:
            try:
                if args.observed is None:
                    print(datetime.now().strftime("%H:%M:%S") + " - Collect all interactions ...")
                    obs_exp_df = pd.concat(all_rank_df, axis = 0, ignore_index = True) # collect all counts (obs/exp)
                    print(datetime.now().strftime("%H:%M:%S") + " - From bin match to peak coordinates ...")
                    obs_exp_df = bin2peak(obs_exp_df, peak_df, resolution)
                    print(datetime.now().strftime("%H:%M:%S") + " - Write obs/exp contact counts ...")
                    obs_exp_df.to_csv(output + ".txt", sep = "\t", header = True, index = False, float_format='%.3f')
                else:
                    print(datetime.now().strftime("%H:%M:%S") + " - Observed/expected contacts provided ...")
                    obs_exp_df = pd.read_table(args.observed, sep = "\t", header = 0)
                # add lambda group for expected value
                print(datetime.now().strftime("%H:%M:%S") + " - Lambda group expected values ...")
                obs_exp_df = lambda_group(obs_exp_df, window) # group expected counts into lambda group
                print(datetime.now().strftime("%H:%M:%S") + " - Get p-value of observed contacts compare to background ...")
                obs_exp_p = poisson_p(obs_exp_df, window) # perform poisson test on obs/lambda pair 
                # print(datetime.now().strftime("%H:%M:%S") + " - Adjust p-value ...")
                # obs_exp_p = bh_correction(obs_exp_p, window, args.fdr)
                print(datetime.now().strftime("%H:%M:%S") + " - Add score for each test w/ obs & exp ...")
                lambda_max = obs_exp_p[[f"lambda_{w}" for w in window]].max(axis = 1); enrichment = obs_exp_p["obs"]/lambda_max; del lambda_max
                p_max = obs_exp_p[[f"p_{w}" for w in window]].max(axis = 1);logp = -np.log10(p_max); del p_max 
                obs_exp_p["score"] = logp*enrichment; del logp, enrichment
                obs_exp_p.sort_values(by = ["chrom1", "start1", "start2"], ascending = [True, True, True], inplace = True, ignore_index = True)
                obs_exp_p.drop_duplicates(keep = "first", ignore_index=True, inplace=True)
                print(datetime.now().strftime("%H:%M:%S") + " - Write scored file ...")
                obs_exp_p.drop([c for c in obs_exp_p.columns if "lambda" in c or "p_" in c], axis = 1).to_csv(output + "_score.txt", sep = "\t", header = True, index = False, float_format='%.2f')
                print(datetime.now().strftime("%H:%M:%S") + f" - Filter by minimum contact >= {args.min_contact} & score >= {args.score} ...")
                # filter by contact frequency
                contact_filter = obs_exp_p["obs"] >= args.min_contact
                # filter by p-value (adjusted)
                # bh_filter = (obs_exp_p[[c for c in obs_exp_p.columns if c.startswith("bh")]] > args.fdr).sum(axis = 1) == 0
                score_filter = obs_exp_p["score"] >= args.score
                # get the filtered peak interactions 
                obs_exp_filter = obs_exp_p[score_filter & contact_filter].copy()
                print(datetime.now().strftime("%H:%M:%S") + f" - Write filtered contacts ...")
                obs_exp_filter.to_csv(output + "_filtered.txt", sep = "\t", header = True, index = False, float_format='%.3e')
                longrange(obs_exp_filter).to_csv(output + ".longrange", sep = "\t", index = False, header = False)
                bedpe(obs_exp_filter).to_csv(output + ".bedpe", sep = "\t", index = False, header = False)
                if args.min_distance > 0:
                    print(datetime.now().strftime("%H:%M:%S") + f" - Filter by minimum distance >= {args.min_distance} ...")
                    longrange(obs_exp_filter.query("end2 - start1 > @args.min_distance")).to_csv(output + f"_{args.min_distance}.longrange", sep = "\t", index = False, header = False)
                    bedpe(obs_exp_filter.query("end2 - start1 > @args.min_distance")).to_csv(output + f"_{args.min_distance}.bedpe", sep = "\t", index = False, header = False)
                    obs_exp_filter.query("end2 - start1 > @args.min_distance").to_csv(output + f"_filtered_{args.min_distance}.txt", sep = "\t", header = True, index = False, float_format='%.3e')
                print(datetime.now().strftime("%H:%M:%S") + " - Finished!")
            except Exception as e:
                print(e)
                exit(1)
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
            for (chr,b1,b2), _ in interaction_df.groupby(["chrom1", "bin1", "bin2"]):
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
