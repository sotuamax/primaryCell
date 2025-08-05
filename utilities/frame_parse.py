import numpy as np 
import  pandas as pd 
import logging 

def read_frame(read_qualified, side_end = False):
    """
    Input:
    read_qualified:pandas.DataFrame\treads alignment bed file 
    Return:
    read supporting 0/1/2 frames in three separate dataframes 
    """
    if side_end:
        read_qualified = read_qualified.query("infuse_site == start_ or infuse_site == end_")
    # strand is "+", feature strand is "-"; and vice versa
    read_qualified["cds_1aa"] = np.where(read_qualified["strand"] == "+", read_qualified["end_"]-read_qualified["frame_"], read_qualified["start_"]+read_qualified["frame_"])
    read_qualified["inframe"] = abs(read_qualified["cds_1aa"] - read_qualified["infuse_site"])%3
    read_frame0 = read_qualified.query("inframe == 0")
    read_frame1 = read_qualified.query("inframe == 2")
    read_frame2 = read_qualified.query("inframe == 1")
    return (read_frame0, read_frame1, read_frame2)

def gene_read_count(frame_df):
    """
    Input: 
    frame_df (dataframe): genes match to a frame 
    Output: 
    tuple: (num of genes, num of reads) 
    """
    frame_gene_num = len(set(frame_df["gene_"]))
    frame_read_num = len(set(frame_df["name"]))
    return (frame_gene_num, frame_read_num)


def frame_log(frame_gene_num, frame_read_num, total, frame):
    """"
    Write frame of 
    genes#
    read support genes#"""
    logging.info(f"Frame{frame} gene#\t{frame_gene_num}")
    logging.info(f"Frame{frame} Read#\t{frame_read_num}({round(frame_read_num/total*100,1)}%)")

def gene_frame(gene_frame_list, frame, barcode = None):
    """
    For read overlaying gene features, report as gene/symbol/read_count and the first aligned read position to tracing back
    Input: 
    gene_frame_list: list of dataframes for read overlaying gene features 
    frame: designated frame 
    barcode: if enable bacord sorting 
    Return: 
    gene and its supporting read_count with reference position.
    """
    frame_df = gene_frame_list[frame]
    c_list = list()
    if barcode == None:
        for g, gdf in frame_df.groupby(["gene_", "symbol_"]):
            c_list.append((g[0], g[1], len(set(gdf["name"])), list(set(gdf["chrom"]))[0], gdf["start"].min()))
        inframe_out = pd.DataFrame(c_list, columns = ["gene", "symbol", "read_count", "chrom", "pos"])
    else:
        for g, gdf in frame_df.groupby("barcode"):
            barcode_target_symbol_num = len(set(gdf["symbol_"]))
            for gg, ggdf in gdf.groupby(["gene_", "symbol_"]):
                c_list.append((g, gg[0], gg[1], barcode_target_symbol_num, len(set(ggdf["name"])), list(ggdf["chrom"])[0], ggdf["start"].min()))
        inframe_out = pd.DataFrame(c_list, columns = ["barcode", "gene", "symbol", "barcode_target_gene_num", "read_count", "chrom", "pos"])
    inframe_out["frame"] = frame
    return inframe_out

def target_log(expected, observed, frame):
    """"
    When target genes indicated, report actual target vs. intended target
    """
    on_target = expected.intersection(observed)
    out_of_target = observed.difference(expected)
    logging.info(f"On/Out-target Frame{frame}#\t{len(on_target)}/{len(out_of_target)}")
