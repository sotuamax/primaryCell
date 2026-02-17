#!/usr/bin/env python
import pandas as pd 
import bioframe as bf
import numpy as np
from statsmodels.stats.multitest import multipletests
import os 
import argparse 
from utilities.misc import ignore_warning
import os 
from utilities.cal_tools import t_statistic, z_statistic
ignore_warning()


def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="group_specific_feature.py -count count.mat -group cell_group -o specific_count ")
    parser.add_argument("-count", "--count", help = "count data in a matrix (e.g., expression matrix, loop contact matrix). \n " \
    "Note that among samples (columns), data values must be normalized to be comparable. ")
    parser.add_argument("-group", "--group", help = "group of cell-types (no header)")
    parser.add_argument("-o", "--output", help = "output name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    count = args.count
    group = args.group
    output = args.output
    count_mat = pd.read_table(count, sep = "\t", header = 0, index_col = 0)
    group_df = pd.read_table(group, sep = "\t", header = None, names = ["celltype", "group"], comment="#")
    group_df.sort_values("group", inplace = True, ignore_index=True)
    try:
        print("Retrieve score for each cell-type included in the group ...")
        count_mat = count_mat[group_df["celltype"].tolist()].copy()
        count_mat.to_csv(f"{output}_mat.tsv", sep = "\t", header = True, index = True)
    except Exception as e:
        print(e)
        exit(0)
    g_list = list()
    for g, g_df in group_df.groupby("group"):
        print(f"Specific analysis for {g} ...")
        if len(g_df) > 1:
            g_statistic = t_statistic(count_mat, group_df["group"], g)
            
        else:
            g_statistic = z_statistic(count_mat, group_df["group"], g)
        _, g_statistic["padj"], _, _ = multipletests(g_statistic["p"], alpha = 0.05, method = "fdr_bh")
        g_statistic.sort_values("statistic", ascending=False).to_csv(f"{output}_{g}.tsv", sep = "\t", header = True, index = True)
        g_list.append(g_statistic[["statistic", "A", "p"]].copy().rename(columns = {"statistic":f"statistic_{g}", "A":g, "p":f"p_{g}"}))
    g_statistic = g_list[0].join(g_list[1:])
    print("Collect all specific score ...")
    g_statistic.to_csv(f"{output}.tsv", sep = "\t", header = True, index = True)

if __name__ == "__main__":
    main()

