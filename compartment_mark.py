#!/usr/bin/env python3
"""
This script is to report overlaps between compartment and histone modification marks; 

Input: 
    required: 
    - compartment (from dchic, in bedGraph format);
    - histone modification marks (in bed/narrowPeak format, can be multiple files); 

Output compartment overlaps with histone mark per file (for A/B count, see comment line)
"""

import pandas as pd
import bioframe as bf 
import argparse 
import numpy as np 

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = False)
    parser.add_argument("compartment", help = "compartment score in bedGraph format")
    parser.add_argument("-mark", nargs = "+", help = "histone modification mark in bed format")
    parser.add_argument("-mark_labels", nargs = "+", help = "mark labels used to distinguish between marks input (1 by 1 match to mark)")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def main():
    args = args_parser()
    compartment = args.compartment 
    mark = args.mark 
    output = args.output
    # read compartment file 
    compartment_df = bf.read_table(compartment, schema = "bedGraph")
    mark_dict = dict(zip(mark, args.mark_labels))
    for m in mark:
        if m.endswith(".narrowPeak"):
            mark_df = bf.read_table(m, schema = "narrowPeak")
        if m.endswith(".bed"):
            mark_df = bf.read_table(m, schema = "bed4")
        # find overlap between compartment and marks 
        compartment_mark = bf.overlap(compartment_df, mark_df, how = 'inner')
        AB_compt = np.where(compartment_mark["score"] > 0, "A", "B")
        AB_count_dict = dict(np.unique(AB_compt, return_counts = True))
        with open(f"{output}_{mark_dict[m]}.txt", "w") as fw:
            fw.write(f"#{mark_dict[m]}_A:{AB_count_dict["A"]}\n")
            fw.write(f"#{mark_dict[m]}_B:{AB_count_dict["B"]}\n")
            compartment_mark.to_csv(f"{output}_{mark_dict[m]}.txt", header = True, index = False)


if __name__ == "__main__":
    main()
