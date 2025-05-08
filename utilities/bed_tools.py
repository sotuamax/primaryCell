import pandas as pd 
import numpy as np
import bioframe as bf
import os 


def bed2fasta(bed, ref, out = None):
    from Bio import SeqIO 
    ref = SeqIO.parse(ref, "fasta")
    ref_dict = SeqIO.to_dict(ref)
    bed_df = bf.read_table(bed, schema = "bed4")
    if out is None:
        out = bed.split(".bed")[0]
    with open(out + ".fa", "w") as fw:
        for row in bed_df.itertuples(): 
            bed_seq = str(ref_dict[row.chrom].seq)[row.start:row.end]
            fw.write(f">{row.chrom}:{row.start}-{row.end} {row.name}\n{bed_seq}\n")

def pe2bed(loop_df):
    """
    Transform bedpe to bed file (remove duplicated intervals)
    """
    l_df = loop_df.copy()
    loop_anchor = pd.concat([l_df[["chrom1", 'start1', "end1"]].rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"}), l_df[["chrom2", "start2", "end2"]].rename(columns = {"chrom2":"chrom", "start2":"start", "end2":"end"})], axis = 0).drop_duplicates(keep = "first", ignore_index=True)
    return loop_anchor 

def bedpe_retrieve(region, loop_df):
    """
    Given region, to search loops with either anchor side overlap it. 
    Input: 
    region: tuple with chrom, start, end 
    loop_df: pandas df in bedpe format 
    Return: 
    selected loop_df that overlaps the given region. 
    """
    # check if loop_df index duplicated or not 
    if len(set(loop_df.index)) == len(loop_df):
        loop_df.index = range(len(loop_df))
    # get loop (left, chrom1, start1, end1) overlap region
    left_loop = bf.select(loop_df.rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"}), region)
    # get loop (right, chrom2, start2, end2)
    right_loop = bf.select(loop_df.rename(columns = {"chrom2":"chrom", "start2":"start", "end2":"end"}), region)
    index_all = np.unique(np.concatenate((left_loop.index, right_loop.index)))
    return loop_df.loc[index_all, :].copy()

def bedpe_colname(bedpe_df):
    """
    Based on the shape of bedpe, to name column names for bedpe file.
    """
    if len(bedpe_df.columns) == 6:
        bedpe_df.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2"]
    if len(bedpe_df.columns) == 7:
        bedpe_df.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2", "name"]
    if len(bedpe_df.columns) == 8:
        bedpe_df.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2", "name", "score"]
    if len(bedpe_df.columns) == 9:
        bedpe_df.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2", 'name', "strand1", "strand2"]
    if len(bedpe_df.columns) == 10:
        bedpe_df.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2", 'name', "score", "strand1", "strand2"]
    return bedpe_df

def is_bedpe(bedpe_df):
    """
    Verify bedpe format as correct. 
    """
    if len(bedpe_df.query("start1 > end1")) != 0:
        return False
    if len(bedpe_df.query("start2 > end2")) != 0:
        return False 
    if len(bedpe_df.query("start1 > start2")) != 0:
        return False 
    return True

def bed_overlap(bed1, bed2):
    bed_df1 = bf.read_table(bed1, schema = "bed3")
    bed_df2 = bf.read_table(bed2, schema = "bed3")
    print(f"Input bed {len(bed_df1)} vs. {len(bed_df2)}")
    bed_overlap = bf.overlap(bed_df1, bed_df2, how = "inner")
    bed_overlap1 = bed_overlap[["chrom", "start", "end"]].drop_duplicates(keep = "first")
    bed_overlap2 = bed_overlap[["chrom_", "start_", "end_"]].drop_duplicates(keep = "first")
    bed_overlap1.columns = bed_df1.columns 
    bed_overlap2.columns = bed_df2.columns
    print(f"Overlap bed {len(bed_overlap1)} vs. {len(bed_overlap2)}")
    return (bed_overlap1, bed_overlap2)

def loop_overlap(loop1, loop2, extend = None):
    """
    compare two loop-interaction data, and return loops that overlapped between them in overlap_loop1 and overlap_loop2
    To return the overlapped loops, set loop_return True.
    input loop file format:
    chrom1, start1, end1, chrom2, start2, end2
    """
    loop_df1 = pd.read_table(loop1, sep = "\t", header = None)
    loop_df2 = pd.read_table(loop2, sep = "\t", header = None)
    loop_df1 = loop_df1.iloc[:, :6]
    loop_df2 = loop_df2.iloc[:, :6]
    loop_df1.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2"]
    loop_df2.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2"]
    if not is_bedpe(loop_df1):
        print("filter performed on loop1")
        loop_df1 = loop_df1.query("start1 < start2")
    if not is_bedpe(loop_df2):
        print("filter performed on loop2")
        loop_df2 = loop_df2.query("start1 < start2")
    if extend is not None:
        loop_df1["start1"] = loop_df1["start1"]-extend
        loop_df1["end1"] = loop_df1["end1"]+extend
        loop_df1["start2"] = loop_df1["start2"]-extend
        loop_df1["end2"] = loop_df1["end2"]+extend
        loop_df2["start1"] = loop_df2["start1"]-extend
        loop_df2["end1"] = loop_df2["end1"]+extend
        loop_df2["start2"] = loop_df2["start2"]-extend
        loop_df2["end2"] = loop_df2["end2"]+extend
    # rename left anchor for loop1/2, and find overlaps between them for left anchor
    print(f"Total {len(loop_df1)} vs. {len(loop_df2)}")
    left_overlap = bf.overlap(loop_df1.rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"}), loop_df2.rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"}), how = "inner")
    right_interval1 = [pd.Interval(row.start2, row.end2) for row in left_overlap.itertuples()]
    right_interval2 = [pd.Interval(row.start2_, row.end2_) for row in left_overlap.itertuples()]
    # find overlaps for right anchor (using Interval)
    rows = [True if p[0].overlaps(p[-1]) else False for p in list(zip(right_interval1, right_interval2))]
    overlap_df = left_overlap[rows]
    # get loop1 end for overlapped loops 
    loop_df1_overlap = overlap_df[[c for c in overlap_df.columns if not c.endswith("_")]].copy(); loop_df1_overlap.columns = loop_df1.columns; loop_df1_overlap.drop_duplicates(keep = "first", inplace = True)
    # get loop2 end for overlapped loops 
    loop_df2_overlap = overlap_df[[c for c in overlap_df.columns if c.endswith("_")]].copy(); loop_df2_overlap.columns = loop_df2.columns; loop_df2_overlap.drop_duplicates(keep = "first", inplace = True)
    print(f"Overlap {len(loop_df1_overlap)} vs. {len(loop_df2_overlap)}")
    return loop_df1_overlap, loop_df2_overlap
 
def bedpe_strand_ann(bedpe1, bed2):
    """
    Annotate strand for loop-interaction data given bed file strand information. 
    To identify overlaps between bedpe1 and bed2, and 
    return it in the original bedpe1 format with "strand1/2" labeled by bed.
    """
    bedpe_df = bf.read_table(bedpe1, schema = "bedpe")
    bed_df = bf.read_table(bed2, schema = "bed6") # contain strand info
    if len(set(bed_df["name"])) != len(bed_df):
        bed_df["name"] = range(len(bed_df))
    left_bedpe = bedpe_df[[col for col in bedpe_df.columns if col.endswith("1")] + ["name"]].copy()
    left_bedpe.columns = [c.replace("1", "") for c in [col for col in bedpe_df.columns if col.endswith("1")] + ["name"]]
    left_bedpe_overlap = bf.overlap(left_bedpe, bed_df).sort_values(by = ["chrom", "start", "end"])
    left_single_overlap = left_bedpe_overlap[left_bedpe_overlap.groupby(["chrom", "start", "end"])["name_"].transform('count') == 1].copy()
    left_single_overlap["strand"] = left_single_overlap["strand_"]
    # 
    right_bedpe = bedpe_df[[col for col in bedpe_df.columns if col.endswith("2")] + ["name"]].copy()
    right_bedpe.columns = [c.replace("2", "") for c in [col for col in bedpe_df.columns if col.endswith("2")] + ["name"]]
    right_bedpe_overlap = bf.overlap(right_bedpe, bed_df).sort_values(by = ["chrom", "start", "end"])
    right_single_overlap = right_bedpe_overlap[right_bedpe_overlap.groupby(["chrom", "start", "end"])["name_"].transform('count') == 1].copy()
    right_single_overlap["strand"] = right_single_overlap["strand_"]
    # 
    left_strand = left_single_overlap[['name', "strand"]].rename(columns = {"strand":"strand1"})
    right_strand = right_single_overlap[['name', "strand"]].rename(columns = {"strand":"strand2"})
    strand_df = pd.merge(left_strand, right_strand, on = "name", how = 'inner')
    bedpe_df = pd.merge(bedpe_df.drop(["strand1", "strand2"], axis = 1), strand_df, on = "name", how = "left")
    return bedpe_df 

def bed_center(bed, schema = "bed3"):
    """
    Return the midpoint of the start and end region in the bed file
    Input: 
    bed: bed file (chrom, start, end)
    Output: 
    pd.DataFrame: midpoint of start and end of given bed file (chrom, pos)
    """
    bed_df = bf.read_table(bed, schema = schema)
    bed_df["pos"] = (bed_df["start"]+bed_df["end"])//2
    bed_pos_df = bed_df[["chrom", "pos"]]
    return bed_pos_df

def bed_filter(bed, schema = "bed5"):
    """
    Filter bed file for chromosome regions (not mitochondrial)
    Input: 
    bed: bed file (default schema bed5: chrom, start, end, name, strand)
    output: 
    bed: filtered bed file (chromosomes retained)
    """
    bed_pos_df = bf.read_table(bed, schema = schema)
    bed_pos_df = bed_pos_df[(bed_pos_df["chrom"].str.startswith("chr")) & (bed_pos_df["chrom"] != "chrM")].copy()
    bed_pos_df.index = range(len(bed_pos_df))
    return bed_pos_df

def rm_blacklist(bed, blacklist, schema = "bed5"):
    """
    Remove 
    Input: 
    bed: bed file 
    blacklist: blacklist file 
    Output: 
    bed: bed file after removing region overlaying blacklist region
    See: https://bioframe.readthedocs.io/en/latest/guide-intervalops.html#subtract-set-difference
    """
    bed_df = bf.read_table(bed, schema = schema)
    black_df = bf.read_table(blacklist, schema = "bed3")
    bed_no_black = bf.setdiff(bed_df, black_df)
    # setdiff is to remove any intervals in bed_df that overlaps black_df 
    # bed_overlap_black = bf.overlap(bed_df, black_df, how = "left")
    # bed_no_black = bed_overlap_black[bed_overlap_black.isna().any(axis = 1)].drop(["chrom_", "start_", "end_"], axis = 1)
    # bed_no_black.index = range(len(bed_no_black))
    return bed_no_black

def orphan_bed(bed_centered, window):
    """
    Input: 

    """
    bed_centered["start"] = bed_centered["pos"]-window
    bed_centered["end"] = bed_centered["pos"]+window
    bed_cluster = bf.cluster(bed_centered[["chrom", "start", "end", "pos"]], min_dist = 0)
    cluster, cluster_cnt = np.unique(bed_cluster["cluster"], return_counts = True)
    orphan_cluster = np.argwhere(cluster_cnt == 1).reshape(-1, ).tolist()
    orphan_bed_df = bed_cluster[bed_cluster["cluster"].isin(orphan_cluster)]
    return orphan_bed_df[["chrom", "pos"]]

def bigwig_parse(bw, chrom, start, end):
    import pyBigWig 
    bw_handle = pyBigWig.open(bw)
    bw_handle.values(chrom, start, end)
    
