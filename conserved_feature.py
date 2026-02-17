#!/usr/bin/env python
"""
There are different ways to get most conserved features:
1. absolution deviation from mean (mean distance)
2. absolution deviation from median (median distance)
When values no negative values: 
3. entropy + static score 
4. coeffient Variance

"""
import pandas as pd 
import bioframe as bf
import numpy as np
import os 
import argparse 
from utilities.misc import ignore_warning
import os 
ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="group_specific_feature.py -count count.mat -group cell_group -o specific_count ")
    parser.add_argument("-count", "--count", help = "count data in a matrix (e.g., expression matrix, loop contact matrix). \n " \
    "Note that among samples (columns), data values must be normalized to be comparable. ")
    parser.add_argument("-negative", "--negative", action = "store_true", help = "if count includes negative values")
    parser.add_argument("-o", "--output", help = "output name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    count = args.count 
    output = args.output 
    count_mat = pd.read_table(count, sep = "\t", header = 0, index_col = 0)
    from utilities.cal_tools import distance_from_mean, distance_from_median
    mean_distance = distance_from_mean(count_mat)
    median_distance = distance_from_median(count_mat)
    mean_distance_df = pd.DataFrame(mean_distance.mean(axis = 1), columns=["distance_mean"])
    median_distance_df = pd.DataFrame(np.median(median_distance, axis = 1), columns = ["distance_median"], index=count_mat.index)
    count_min = pd.DataFrame(count_mat.min(axis = 1), columns = ["min"])
    count_max = pd.DataFrame(count_mat.max(axis = 1), columns = ["max"])
    count_average = pd.DataFrame(count_mat.mean(axis = 1), columns = ["average"])
    count_median = pd.DataFrame(np.median(count_mat, axis = 1), columns = ["median"], index = count_mat.index)
    count_score = count_average.join([count_median, count_min, count_max, mean_distance_df, median_distance_df])
    # distance_sum_df.sort_values("distance", inplace = True)
    if not args.negative:
        from utilities.cal_tools import entropy, CoeVar
        count_entropy = entropy(count_mat, static=True)
        count_score = count_score.join([count_entropy])
        count_coevar = CoeVar(count_mat)
        count_score = count_score.join([count_coevar])
    print("Write conservation score into file ...")
    count_score.to_csv(output + ".tsv", sep = "\t", header = True, index = True)

if __name__ == "__main__": 
    main()
