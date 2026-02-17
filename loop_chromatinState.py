#!/usr/bin/env python 
import pandas as pd 
import bioframe as bf
from utilities.bed_tools import open_loop, pe2bed 
import argparse

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, usage="")
    parser.add_argument("-loop", "--loop", help = "loop data in bedpe format (preferred), also support longrange format")
    parser.add_argument("-peak", "--peak", help = "peak data in narrowPeak format (this is before reporting consensus)")
    parser.add_argument("-state", "--chromstate", help = "chromotin states for the same cell type.")
    parser.add_argument("-o", "--output", help = "output prefix")
    args=parser.parse_args()
    return args

def category_dict():
    TSS_dict = {"1_TssA":"TSS", "2_TssAFlnk":"TSS"}; 
    BIV_dict = {"10_TssBiv":'BIV', "11_BivFlnk":'BIV', "12_EnhBiv":'BIV'}; 
    TX_dict = {"3_TxFlnk":"TX", "4_Tx":"TX", "5_TxWk":"TX"}; 
    ENH_dict = {"6_EnhG":"ENH", "7_Enh":"ENH"}; 
    REPRESS_dict = {"13_ReprPC":"REPRESS", "14_ReprPCWk":"REPRESS"}; 
    REPEAT_dict = {"8_ZNF/Rpts":"REPEAT"}; 
    HET_dict = {"9_Het":"HET"}; 
    QUIES_dict = {"15_Quies":"QUIES"}
    combined_dict = TSS_dict | BIV_dict | TX_dict | ENH_dict | REPRESS_dict | REPEAT_dict | HET_dict | QUIES_dict 
    return combined_dict 

def main():
    args = args_parser()
    loop = args.loop; peak = args.peak; state = args.chromstate; output = args.output
    state_df = bf.read_table(state, schema = "bed9", skiprows=1)
    loop_df = open_loop(loop)
    loop_anchor = pe2bed(loop_df)
    loop_anchor_initial = bf.read_table(peak, schema = "narrowPeak").query("score > 100 or fc > 2.5")
    loop_anchor_initial["summit"] = loop_anchor_initial["start"] + loop_anchor_initial["relSummit"]
    # for each loop consensus anchor, identify its cell type peak summit 
    consensus_overlap_initial = bf.overlap(loop_anchor, loop_anchor_initial, how = "inner")
    consensus_summit = consensus_overlap_initial[["chrom", "start", "end", "summit_"]].drop_duplicates(keep = "first")
    loop_df = pd.merge(loop_df, consensus_summit, left_on = ["chrom1", "start1", "end1"], right_on = ["chrom", "start", "end"], how = "inner").drop(["chrom", "start", "end"], axis = 1).rename(columns = {"summit_":"summit1"}) # add summit1 for upstream anchor 
    loop_df = pd.merge(loop_df, consensus_summit, left_on = ["chrom2", "start2", "end2"], right_on = ["chrom", "start", "end"],how = "inner").drop(["chrom", "start", "end"], axis = 1).rename(columns = {"summit_":"summit2"}) # add summit2 for downstream anchor 
    loop_df["summit1_1"] = loop_df["summit1"] + 1
    loop_df["summit2_1"] = loop_df["summit2"] + 1
    # overlap of loop_df and chromotin states 
    left_join = bf.overlap(loop_df.rename(columns = {"chrom1":"chrom", "summit1":"start", "summit1_1":"end"}), state_df[["chrom", "start", "end", "name"]], how = "inner").rename(columns = {"chrom":"chrom1", "name_":"state1"}) 
    right_join = bf.overlap(loop_df.rename(columns = {"chrom2":"chrom", "summit2":"start", "summit2_1":"end"}), state_df[["chrom", "start", "end", "name"]], how = "inner").rename(columns = {"chrom":"chrom2", "name_":"state2"}) 
    # join left and right anchor chromotin states
    loop_state = pd.merge(left_join[["chrom1", "start1", "end1", "chrom2", "start2", "end2", "state1"]], right_join[["chrom1", "start1", "end1", "chrom2", "start2", "end2", "state2"]], on = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]).drop_duplicates(keep = "first")
    # write loop state into a file 
    loop_state.sort_values(["chrom1", "start1", "chrom2", "start2"], inplace = True, ignore_index=True)
    loop_state["name"] = loop_state["state1"].astype(str) + "|" + loop_state["state2"].astype(str)
    loop_state[["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name"]].to_csv(f"{output}_state.bedpe", sep ="\t", header = True, index = False)
    print(f"Used loop number:{len(loop_df)}; Evaluated loop number:{len(loop_state)}")
    # summarize state pairs statistics 
    print("report chromatin state for loop connection ...")
    state_pair = pd.crosstab(loop_state["state1"], loop_state["state2"]).reset_index().melt(id_vars = "state1", var_name = "state2", value_name = "freq")
    pair_list = list()
    for _, row in state_pair.iterrows():
        pair_list.append("-".join([i for i in sorted([row.state1, row.state2])]))
    state_pair["pair"] = pair_list
    state_pair_freq = state_pair.groupby("pair").sum("freq").reset_index().sort_values("freq", ascending=False, ignore_index=True)
    state_pair_freq["ratio"] = state_pair_freq["freq"]/sum(state_pair_freq["freq"])
    state_pair_freq.to_csv(f"{output}_freq.txt", sep = "\t", header = True, index = False)
    print("add pair category ...")
    state_df = pd.DataFrame.from_dict(category_dict(), orient="index", columns = ["category"]).reset_index().rename(columns = {"index":"state"})
    state_pair_freq[["state1", "state2"]] = state_pair_freq["pair"].str.split("-", expand = True)
    state_temp_ = pd.merge(state_pair_freq, state_df, left_on = "state1", right_on = "state").rename(columns = {"category":"category1"}).drop("state", axis = 1)
    state_pair = pd.merge(state_temp_, state_df, left_on = "state2", right_on = "state").rename(columns = {"category":"category2"}).drop("state", axis = 1)
    state_pair["category"] = state_pair["category1"].astype(str) + "-" + state_pair["category2"].astype(str)
    category_freq = state_pair.groupby("category").sum("freq").reset_index().sort_values("freq", ascending=False, ignore_index=True)
    category_freq.to_csv(f"{output}_category_freq.txt", sep = "\t", header = True, index = False)


if __name__ == "__main__":
    main()
