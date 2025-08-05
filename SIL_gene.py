#!/usr/bin/env python3
"""
To identify SIL target genes; 

SIL is a genomic interval where siginificant chromotin interactions observed.

Input: 
required: 
- annotated loops; 


"""
import pandas as pd 
import os 
import numpy as np 
import bioframe as bf 
import argparse 
from utilities.bed_tools import bedpe_retrieve, tether_bed

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help=True, usage="\nSIL.py <input.bedpe> -o <output>")
    parser.add_argument("ann", help = "annotated loops")
    parser.add_argument("-o", "--output", help = "output name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    ann = args.ann 
    ann_df = pd.read_table()
    
if __name__ == "__main__":
    main()



SIL_loop_gene = gene_loop_count(SIL_loop_ann, ["SG", "IG", "GG"], disease_gene = disease_g["gene"].tolist())
SIL_loop_gene.to_csv("/data/jim4/Seq/primary_cell_project/analysis/Loops/output/disease/SIL_gene_loop_disease.txt", sep = "\t", header = True, index = False)

SIL_loop_gene.query("loops > 20").to_csv("/data/jim4/Seq/primary_cell_project/analysis/Loops/output/SIL/SIL_highly_targeted_gene.txt", sep = "\t", header = True, index = False)

