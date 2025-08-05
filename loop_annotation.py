#!/usr/bin/env python 
"""
Annotate loops based on gene annotation. 

annotate loop anchors first, and mirror anchor annotation 
onto loop and assign loop annotation group. 

anchor annotation sequence: 
- promoter anchor annotated; 
- anchor not overlapping promoter but with gene body; 
- anchor not overlapping with either promoter or gene body; 

Note: for loops overlapping with multiple genes, it can be annotated multiple times, one annotation per gene. 
for example: 
chr7	26857036	26858582	chr7	27146890	27147195	GG=[SKAP2]-[HOXA3]
chr7	26857036	26858582	chr7	27146890	27147195	GG=[SKAP2]-[HOXA6]


required: 
- bedpe: bedpe input;
- output: output prefix 
optional: 
- gene: gene 

"""
import pandas as pd 
import os 
import numpy as np 
import bioframe as bf 
import pysam 
import csv 
import time 
import argparse 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help=True, 
                                   usage="\nloop_annotation.py sample.bedpe -o <out>")
    parser.add_argument("bedpe", help = "input bedpe file (loop interactions)")
    parser.add_argument("-gene", "--gene", default = "/data/jim4/Reference/human/GRCh38.p14/GTF/ncbiRefSeqSelect.bed", help = "Input gene bed")
    parser.add_argument("-o", "--output", required = True, help = "output name")
    args=parser.parse_args()
    return args

def find_TSS(transcript_df):
    """
    Given transcript dataframe in bed format, to find transcript TSS in a bed format
    """
    transcript_df_new = transcript_df.copy()
    transcript_df_new["pos"] = np.where(transcript_df_new["strand"] == "-", transcript_df_new["end"], transcript_df_new["start"])
    transcript_df_new['start'] = transcript_df_new["pos"] - 500; transcript_df_new["end"] = transcript_df_new["pos"] + 500
    return transcript_df_new.drop("pos", axis = 1)

def annotate_anchor(loop_anchor, transcript_df):
    """
    Annotate loop based on input gene transcript. 
    Anchor annotation: 
    - Use transcript info to identify promoter region (500 bp around TSS), and loop anchors overlapping promoter indicated as P anchor labeled with its target gene name; 
    - Remaining loop anchors are then annotated when they are overlapping gene body region, labeled as its target gene name; 
    - For anchors not overlapping either P or gene body region, labeled as intergenic; 
    Loop annotation: 
    - fill anchor annotation onto its loop side; 
    - when both loop anchors are on intergenic region, genes the loop across are indicated by its strand; 
    Input: 
    loop_df: loop df with minimum columns "chrom1", 'start1', "end1", "chrom2", "start2", "end2"
    transcript_df: transcript df with minimum columns "chrom", 'start', "end", "gene", "strand" 
    Output: 
    loop_ann_df: loop df with annotation in columns "gene1" for anchor1, "gene2" for anchor2 
    
    Annotation label: 
    genes are wrapped by bracket (e.g., [MYC])
    P:[MYC] for anchor overlapping MYC Promoter; 
    N for intergenic region; 
    N:[ABC1][ABC2] for anchor across multiple genes ("+" in "gene1", "-" in "gene2"); 
    """
    # to avoid overwrite on the original df 
    t_df = transcript_df.copy()
    # update gene name, contained within "[]"
    t_df["gene"] = "[" + t_df["gene"].astype(str) + "]"
    # get loop anchors 
    # get TSS
    transcript_tss = find_TSS(t_df)
    # find anchor overlapping TSS 
    TSS_anchor = bf.overlap(loop_anchor, transcript_tss[["chrom", "start", "end", "gene"]], how = "inner")
    TSS_anchor = TSS_anchor[["chrom", 'start', "end", "gene_"]].drop_duplicates(keep = "first").rename(columns = {"gene_":"gene"})
    # loop anchor overlapping gene body but not TSS 
    query_anchor_1 = pd.merge(loop_anchor, TSS_anchor, on = ["chrom", "start", "end"], how = "left")
    nTSS_anchor = query_anchor_1[pd.isna(query_anchor_1["gene"])].drop("gene", axis = 1)
    # find non-TSS anchor overlapping gene body 
    body_anchor = bf.overlap(nTSS_anchor, t_df[["chrom", "start", "end", "gene"]], how = "inner")# .sort_values(["chrom", "start", 'end'])
    body_anchor = body_anchor[["chrom", "start", "end", "gene_"]].drop_duplicates(keep = "first").rename(columns = {"gene_":"gene"})
    body_anchor["gene"] = body_anchor["gene"]
    # when anchor is not overlapping either gene body or gene promoter, it is assigned as intergenic anchor
    query_anchor_2 = pd.merge(loop_anchor, pd.concat([TSS_anchor, body_anchor], axis = 0), on =["chrom", "start", "end"], how = "left")
    intergenic_anchor = query_anchor_2[pd.isna(query_anchor_2["gene"])].drop("gene", axis = 1)
    # 
    TSS_anchor["gene"] = "P:" + TSS_anchor["gene"].astype(str)
    intergenic_anchor["gene"] = "N"
    # anchor annotation based on gene
    anchor_ann = pd.concat([TSS_anchor, body_anchor, intergenic_anchor], axis = 0) 
    return anchor_ann 

def annotate_loop(loop_df, anchor_ann, transcript_df):
    t_df = transcript_df.copy()
    t_df["gene"] = "[" + t_df["gene"].astype(str) + "]"
    l_df = loop_df.copy()
    # add annotation to loop data 
    loop_ann1 = pd.merge(l_df, anchor_ann.rename(columns = {"chrom":"chrom1", "start":"start1", "end":"end1", "gene":"gene1"}), on = ["chrom1", "start1", "end1"])
    loop_ann2 = pd.merge(l_df, anchor_ann.rename(columns = {"chrom":"chrom2", "start":"start2", "end":"end2", "gene":"gene2"}), on = ["chrom2", "start2", "end2"])
    loop_ann = pd.merge(loop_ann1, loop_ann2, on = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    loop_ann_final = pd.merge(l_df, loop_ann, on = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"], how = "left").sort_values(["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    # add gene info when loop not overlapping genes, but may contain gene(s)
    loop_ann_final_N = loop_ann_final.query("gene1 == 'N' and gene2 == 'N'").copy()
    z_gene = list()
    for z in zip(loop_ann_final_N["chrom1"], loop_ann_final_N["start1"], loop_ann_final_N["end2"]):
        z_select = bf.select(t_df, z)
        if len(z_select) > 0: 
            plus_gene = z_select.query("strand == '+'")["gene"]
            minus_gene = z_select.query("strand == '-'")["gene"]
        else:
            plus_gene = ""
            minus_gene = ""
        z_gene.append((z[0], z[1], z[2], "N:" + "".join(plus_gene), "N:" + "".join(minus_gene)))
    intergenic_loop = pd.DataFrame(z_gene, columns = ["chrom1", "start1", "end2", "gene1", "gene2"])
    # add re-annotated intergenic loop into the annotation df 
    loop_ann_final_N_update = pd.merge(loop_ann_final_N[["chrom1", "start1", "end1", "chrom2", "start2", "end2"]], intergenic_loop, on = ["chrom1", "start1", "end2"])
    loop_ann_update = pd.concat([loop_ann_final.query("gene1 != 'N' or gene2 != 'N'"), loop_ann_final_N_update], axis = 0).sort_values(["chrom1", "start1", "start2"], ignore_index=True)
    return loop_ann_update

def group_loop(loop_ann_input):
    """
    Classify loops based on its pattern related to genes; 
    Input: 
    loop_ann_df: loop_df with annotation in "gene1", "gene2" 
    Output: 
    loop_ann_group: loop_df with additional columns of "class", "sublass"
    
    Loop class: 
    SG - single gene; 
    GG - gene gene; 
    IG - inside gene; 
    AG - across gene; 
    NG - none gene; 
    """
    # prevent overwrite on original input df 
    loop_ann_df = loop_ann_input.copy()
    target = loop_ann_df["gene1"].astype(str) + "_" + loop_ann_df["gene2"].astype(str)
    ### put into classes 
    single_target_P = loop_ann_df[target.str.contains("P:") & (target.str.contains("^N_") | target.str.contains("_N$"))] # single target gene on its promoter, 44916
    single_target = loop_ann_df[~target.str.contains("P:") & (target.str.contains("^N_") | target.str.contains("_N$"))] # single target gene on its gene body, 80421
    gene_gene_PP = loop_ann_df[target.str.contains("^P:") & target.str.contains("_P:")] # gene promoter - gene promote, 15531
    gene_gene_P = loop_ann_df[(target.str.contains("^P:") & ~target.str.contains("_P:") & (target.str.replace("P:", "").str.split("_", expand = True)[0] != target.str.replace("P:", "").str.split("_", expand = True)[1])) | (~target.str.contains("^P:") & target.str.contains("_P:") & (target.str.replace("P:", "").str.split("_", expand = True)[0] != target.str.replace("P:", "").str.split("_", expand = True)[1]))].query("gene1.str.contains('\[') and gene2.str.contains('\[')") # gene body - gene promoter, 40938
    gene_gene = loop_ann_df[(target.str.contains("^\[")) & (target.str.contains("\]$")) & (~target.str.contains(":"))].query("gene1 != gene2") # gene body - gene body, 31178
    inside_gene_P = loop_ann_df[(target.str.contains("P:")) & (target.str.replace("P:", "").str.split("_", expand = True)[0] == target.str.replace("P:", "").str.split("_", expand = True)[1])] # gene promoter within it self, 26838
    inside_gene = loop_ann_df.query("gene1 != 'N:' and gene1 == gene2") # within it self, 136844
    N_single_gene = loop_ann_df[target.str.contains("^N:") & target.str.contains("_N:")].query("(gene1 == 'N:' or gene2 == 'N:') and gene1 != gene2") # constrain single gene, 13396 
    N_genes = loop_ann_df[target.str.contains("^N:") & target.str.contains("_N:")].query("gene1 != 'N:' and gene2 != 'N:' and gene1 != gene2") # constrain multiple genes, 4260 
    N_N = loop_ann_df[target.str.contains("^N:") & target.str.contains("_N:")].query("gene1 == 'N:' and gene1 == gene2") # constrain no gene, 100342
    # revalidate loop annotation group is match to the original input (groups are exclusive)
    loop_ann_class_validate = pd.concat([single_target_P, single_target, gene_gene, gene_gene_P, gene_gene_PP, inside_gene_P, inside_gene, N_single_gene, N_genes, N_N], axis = 0)
    if len(loop_ann_class_validate) == len(loop_ann_df):
        print("Grouped annotation is passed !")
    # add annotation to loops 
    loop_ann_df["subcategory"] = np.where(loop_ann_df.index.isin(single_target_P.index), "SG_P", 
                                    np.where(loop_ann_df.index.isin(single_target.index), "SG", 
                                                np.where(loop_ann_df.index.isin(gene_gene_PP.index), "GG_PP", 
                                                        np.where(loop_ann_df.index.isin(gene_gene_P.index), "GG_P", 
                                                                np.where(loop_ann_df.index.isin(gene_gene.index), "GG", 
                                                                        np.where(loop_ann_df.index.isin(inside_gene_P.index), "IG_P", 
                                                                                    np.where(loop_ann_df.index.isin(inside_gene.index), "IG", 
                                                                                            np.where(loop_ann_df.index.isin(N_single_gene.index), "AG_S", 
                                                                                                    np.where(loop_ann_df.index.isin(N_genes.index), "AG_M", 
                                                                                                            np.where(loop_ann_df.index.isin(N_N.index), "NG", "other"))))))))))
    # add class 
    loop_ann_df["category"] = loop_ann_df["subcategory"].str.split("_", expand = True)[0]
    return loop_ann_df 

def main():
    args = args_parser()
    print("Read bedpe file ...")
    loop_df = pd.read_table(args.bedpe, sep = "\t", header = None, names = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2"])
    # "/data/jim4/Reference/human/GRCh38.p14/GTF/select.transcript.GCF_000001405.40_GRCh38.p14_genomic.bed"
    print("Read gene transcript file ...")
    transcript_df = pd.read_table(args.gene, sep = "\t", header = None, names = ["chrom", "start", "end", "gene", "strand"])
    # get anchor from bedpe 
    print("Get loop anchors ...")
    from utilities.bed_tools import pe2bed
    anchor_df = pe2bed(loop_df)
    print("Annotate anchors ...")
    anchor_ann = annotate_anchor(anchor_df, transcript_df)
    # 
    print("Add annotation to loops ...")
    loop_ann_df = annotate_loop(loop_df, anchor_ann, transcript_df)
    loop_ann_grouped = group_loop(loop_ann_df)
    loop_ann_grouped.to_csv(args.output + ".txt", sep = "\t", header = True, index = False)
    # get loop category count 
    loop_category = loop_ann_grouped[["chrom1", 'start1', "end1", "chrom2", "start2", "end2", "subcategory", "category"]].drop_duplicates(keep = "first")
    class_freq = pd.DataFrame(loop_category.groupby(["subcategory", "category"]).size(), columns=["freq"])
    class_freq["prop"] = round(class_freq["freq"]/sum(class_freq["freq"]) * 100, 1)
    class_freq.to_csv(args.output + "_category.txt", sep = "\t", header = True, index = True)

if __name__ == "__main__":
    main()

