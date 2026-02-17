#!/usr/bin/env python 
"""
Generate the correlation between loop-contact and gene-expression.

gene expression filtered: 
    1. remove lowly expressed expression; 
    2. remove genes without variation (CoVar < 0.1); 

example: 

"""
import pandas as pd 
import numpy as np 
# from scipy.stats import spearmanr 
# from joblib import Parallel, delayed
from scipy.stats import rankdata
import argparse 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, usage="")
    parser.add_argument("-loop", "--loop_count", help = "loop count data in the format of a matrix (must be normalized). " \
    "Loop data should be filtered beforehand.")
    parser.add_argument("-exp", "--expression", help = "gene expression data in the format of a matrix (when multiple replicates for one celltype, calculate its average, must be normalized).")
    parser.add_argument("-exp_coevar", "--expression_coevar", type = float, help = "cutoff of coefficient variance for gene expression data.")
    parser.add_argument("-loop_ann", "--loop_annotation", help = "loop annotation to connect from loop to its target genes")
    parser.add_argument("-gg", "--gg", action = "store_true", help = "Use Gene-Gene interaction loop. ")
    parser.add_argument("-shuffle", "--shuffle_annotation", action = "store_true", help = "shuffle the paired relation between loop and gene expression data. ")
    parser.add_argument("-o", "--output", help = "output prefix")
    args=parser.parse_args()
    return args

def exp_filter(exp_df, exp_coevar):
    """
    must filter on gene expression data 
    filter out lowly expressed genes 
    genes should show variations among cell types to be able to do rank for correlation (spearman's)
    """
    # exp_df["gene"] = exp_df.index.to_series().str.split("|", expand = True)[1]
    # exp_df.set_index("gene", inplace = True)
    # drop genes with low expression across all cell-types of interested 
    print("Remove genes with expression < 1 across all samples ...")
    exp_df = exp_df[exp_df.max(axis = 1) > 1].copy()
    print(f"Remove genes with low expression variation (coefficient variance > {exp_covar}) ... ")
    from utilities.cal_tools import CoeVar
    exp_covar = CoeVar(exp_df)
    exp_df = exp_df[exp_covar > exp_coevar].copy()
    return exp_df

def pair_contact_exp(loop2gene_df, loop_count_df, exp_df):
    """
    Docstring for pair_contact_exp
    
    :param loop2gene_df: loop2gene relation in a dataframe (index in loop, gene in column)
    :param loop_count_df: loop_count_df (index in loop, and column for each cell type)
    :param exp_df: exp_df (index in gene, and collumn for each cell type)
    
    :return paired data between loop_count and exp_count 

    """
    loop2gene_count = pd.merge(loop2gene_df, loop_count_df, left_index = True, right_index = True, how = "left").reset_index()
    loop2gene_count = pd.merge(loop2gene_count, exp_df, left_on = "gene", right_index = True, how = "left", suffixes = ["_contact", "_exp"])
    loop2gene_count.dropna(inplace=True, ignore_index=True)
    return loop2gene_count

def corr_contact_exp(loop2gene_count):
    loop_ar = np.array(loop2gene_count[[c for c in loop2gene_count if c.endswith("_contact")]])
    exp_ar = np.array(loop2gene_count[[c for c in loop2gene_count if c.endswith("_exp")]])
    return rowwise_spearman(loop_ar, exp_ar)
    
def rowwise_spearman(X, Y):
    """Row-wise Spearman correlation using average ranks (tie-aware)"""
    # Rank each row with ties handled as average
    Xr = np.apply_along_axis(rankdata, 1, X, method='average')
    Yr = np.apply_along_axis(rankdata, 1, Y, method='average')
    # Compute row-wise Pearson on ranks
    Xc = Xr - Xr.mean(axis=1, keepdims=True)
    Yc = Yr - Yr.mean(axis=1, keepdims=True)
    num = np.sum(Xc * Yc, axis=1)
    den = np.sqrt(np.sum(Xc**2, axis=1) * np.sum(Yc**2, axis=1))
    return num / den

def main():
    args = args_parser()
    loop_count = args.loop_count; loop_ann = args.loop_annotation; exp_count = args.expression
    output = args.output
    print("read loop count matrix ...")
    loop_count_df = pd.read_table(loop_count, sep = "\t", header = 0, index_col = 0) # rows as loop lable, columns as cell types 
    assert len(set(loop_count_df.index)) == len(loop_count_df), print("loops are not unique!")
    print("read expression count ...")
    exp_df = pd.read_table(exp_count, sep = "\t", header = 0, index_col = 0) # each row is a gene, and each column show cell type 
    print("Filter on gene expression data ...")
    exp_df = exp_filter(exp_df, args.expression_coevar)
    print(f"There are {len(set(exp_df.index))} genes retained ...")
    assert len(set(exp_df.index)) == len(exp_df), print("gene expression are not unique!")
    print("Pair cell-type in loop and expression dataset ...")
    loop_celltypes = set(loop_count_df.columns); exp_celltypes = set(exp_df.columns)
    # only cell-types present both in expression and loop data are considered 
    pick_celltypes = sorted(set(loop_celltypes).intersection(exp_celltypes))
    # update expression and loop count data on matched cell-types 
    exp_df = exp_df[pick_celltypes].copy(); loop_count_df = loop_count_df[pick_celltypes].copy()
    assert list(loop_count_df.columns) == list(exp_df.columns), print("Celltypes in loop data and expression does not match.")
    print("read loop annotation ...")
    loop_ann_df = pd.read_table(loop_ann, sep = "\t", header = 0, index_col = 0) 
    # each row is a loop, each column shows the annotated types 
    print(f"Retrive loop & target gene information (Consider GG:{args.gg})...")
    from utilities.bed_tools import parse_loop_ann, bed2index
    loop2gene = parse_loop_ann(loop_ann_df, gg = args.gg); 
    loop2gene_df = pd.merge(loop_ann_df[["chrom1", "start1", "end1", "chrom2", "start2", "end2"]], loop2gene, left_index = True, right_index = True)
    print("Pair loop-count & gene-expression count ...")
    bed2index(loop2gene_df, pe_format=True)
    contact2exp = pair_contact_exp(loop2gene_df, loop_count_df, exp_df)
    contact2exp["cor"] = corr_contact_exp(contact2exp)
    contact2exp.to_csv(f"{output}.txt", sep = "\t", header = True, index = False)
    if args.shuffle_annotation:
        print("Shuffle loop and gene annotations ...")
        ### generate a randomized loop-gene pairs
        all_loop = loop2gene_df.index.tolist() #; all_gene = loop2gene_df["gene1"].tolist()
        np.random.shuffle(all_loop) # ; random.shuffle(all_gene)
        loop2gene_df.index = all_loop
        shuffled_contact2exp = pair_contact_exp(loop2gene_df, loop_count_df, exp_df)
        shuffled_contact2exp["cor"] = corr_contact_exp(shuffled_contact2exp)
        shuffled_contact2exp.to_csv(f"{output}_shuffled.txt", sep = "\t", header = True, index = False)

if __name__ == "__main__":
    main()
