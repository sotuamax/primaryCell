import bioframe as bf 
import numpy as np 
import pandas as pd

def args_parser():
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-cool", "--cool", help = "cool input file")
    parser.add_argument("-r", "--resolution", help = "cool resolution to work with")
    parser.add_argument("-interaction", "--interaction", help = "interaction site to examine count for")
    args=parser.parse_args()
    return args

def cool_retrieve(cool, r):
    b1_relative = b1-offset; b2_relative = b2-offset
    left = b1_relative-w if b1_relative-w > 0 else 0
    right = b2_relative-w if b2_relative-w > 0 else 0 
    left_limit = b1_relative+w+1 if b1_relative+w+1 < clr_matrix.shape[0] else clr_matrix.shape[0]
    right_limit = b2_relative+w+1 if b2_relative+w+1 < clr_matrix.shape[0] else clr_matrix.shape[0]
    wb_mat = clr_matrix[left:left_limit, right:right_limit]


def main():
    args = args_parser()
    cool = args.cool
    interaction = args.interaction
    r = args.resolution
    from utilities.cool_tools import parser_cool
    clr = parser_cool(cool, resolution = r)
    # processing interaction file 
    interaction_df = pd.read_table(interaction, sep = "\t", header = None)
    # 
    max_dist = (interaction_df["start2"] - interaction_df["start1"]).max()
    split_Mb = max_dist//1_000_000
    interaction_df["bin1"] = interaction_df["start1"]//split_Mb
    interaction_df["bin2"] = interaction_df["start2"]//split_Mb
    for c,b1,b2, bdf in interaction_df.groupby(["chrom", "bin1", "bin2"]):
        arm_region = (c, b1*split_Mb, b2*split_Mb)
        clr_matrix = clr.matrix(balance = False).fetch(arm_region, arm_region)
        clr_offset = arm_region[1]//r
        for a,b in list(zip(bdf["start1"], bdf["end2"])):
            clr_value = clr_matrix[a-clr_offset, b-clr_offset]
            print(c,a,b, clr_value)
    



