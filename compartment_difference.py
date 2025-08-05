#!/usr/bin/env python3
"""

To explore shift in compartment status between different cell types

"""

import pandas as pd 
import bioframe as bf 
import numpy as np 
import os 
import sys 
from sklearn.preprocessing import StandardScaler
import argparse 

### read the pair of celltype bedGraph 
def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help=True, 
                                     usage="")
    parser.add_argument("sample1", help = "sample1 bedgraph")
    parser.add_argument("sample2", help = "sample2 bedgraph")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

args = args_parser()

celltype1 = args.sample1
celltype2 = args.sample2
celltype1_df = bf.read_table(celltype1, schema = "bedGraph")
celltype2_df = bf.read_table(celltype2, schema = "bedGraph")

### normalize PC values 
def value_scale(value):
    value = np.array(value)
    scaled_value = StandardScaler().fit_transform(value.reshape(-1,1))
    return scaled_value.reshape(-1)

def bin_scale(value):
    b = pd.cut(value, bins = np.quantile(value, np.arange(0,1.1,0.1)), include_lowest=True, retbins = True, labels = list(range(0,10)))
    return b[0]

def df_group(df):
    df["sign"] = np.sign(df["value"])
    df["z"] = df.groupby("chrom")["value"].transform(value_scale)
    df["bin"] = df.groupby("chrom")["z"].transform(bin_scale)
    return df

pairwise_celltype = bf.overlap(df_group(celltype1_df), df_group(celltype2_df), how = "inner")
pairwise_celltype['diff'] = abs(pairwise_celltype["bin"].astype(int) - pairwise_celltype["bin_"].astype(int))

diff_pairwise = pairwise_celltype.query("diff > 1").copy()
diff_pairwise.sort_values(["chrom", "start"], inplace = True)
diff_pairwise.to_csv(args.output + ".txt", sep = "\t", header = True, index = False)
