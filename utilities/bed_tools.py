import pandas as pd 
import numpy as np
import bioframe as bf


def bed2index(df:pd.DataFrame, reverse = False, pe_format = False):
    """
    input df, to change bed coordinates into dataframe index or retrive information in index back to bed format. 

    """
    if not pe_format:
        if not reverse: # change from bed column to index
            df.index = df["chrom"].astype(str) + ":"+ df["start"].astype(str) + "-" + df["end"].astype(str)
            df.drop(["chrom", "start", "end"], axis = 1, inplace = True)
        if reverse: # change from index to bed columns 
            index_df = pd.DataFrame(df.index)
            df.index = range(len(df))
            chrom_ = index_df[0].str.split(":", expand = True)[0]
            pos_ = index_df[0].str.split(":", expand = True)[1]
            start_ = pos_.str.split("-",expand = True)[0].astype(int)
            end_ = pos_.str.split("-",expand = True)[1].astype(int)
            df["chrom"] = chrom_; df["start"] = start_; df["end"] = end_
    if pe_format:
        if not reverse:
            df.index = df["chrom1"].astype(str) + ":"+ df["start1"].astype(str) + "-" + df["end1"].astype(str) + "|" + df["chrom2"].astype(str) + ":"+ df["start2"].astype(str) + "-" + df["end2"].astype(str)
            df.drop(["chrom1", "start1", "end1", "chrom2", "start2", "end2"], axis = 1, inplace = True)
        if reverse:
            index_df = pd.DataFrame(df.index)
            df.index = range(len(df))
            if "_" in index_df.iloc[0].values[0]:
                pos1_ = index_df[0].str.split("_", expand = True)[0]
                pos2_ = index_df[0].str.split("_", expand = True)[1]
            if "|" in index_df.iloc[0].values[0]:
                pos1_ = index_df[0].str.split("|", expand = True)[0]
                pos2_ = index_df[0].str.split("|", expand = True)[1]
            chrom1_ = pos1_.str.split(":", expand = True)[0]
            chrom2_ = pos2_.str.split(":", expand = True)[0]
            anchor1_ = pos1_.str.split(":", expand = True)[1]
            anchor2_ = pos2_.str.split(":", expand = True)[1]
            start1_ = anchor1_.str.split("-", expand = True)[0].astype(int)
            end1_ = anchor1_.str.split("-", expand = True)[1].astype(int)
            start2_ = anchor2_.str.split("-", expand = True)[0].astype(int)
            end2_ = anchor2_.str.split("-", expand = True)[1].astype(int)
            df["chrom1"] = chrom1_; df["start1"] = start1_; df["end1"] = end1_
            df["chrom2"] = chrom2_; df["start2"] = start2_; df["end2"] = end2_

def parse_loop_ann(loop_ann_df, gg = False):
    """
    Input loop_ann_df must use loop as the index. 
    return: 
    loop - gene pair in a dataframe (only consider single-gene loop, internal-gene loop, gene-gene loop)
    returned loop2gene dataframe with its index match to input annotation dataframe
    """
    if gg:
        loop_ann_df = loop_ann_df[loop_ann_df["category"].isin(["SG", "IG", "GG"])].copy()
    else:
        loop_ann_df = loop_ann_df[loop_ann_df["category"].isin(["SG", "IG"])].copy()
    import re 
    g_ = dict()
    for (i,row) in enumerate(loop_ann_df.itertuples()):
        g_[i] = []
        g1 = re.search(r'\[(.*?)\]', row.gene1)
        g2 = re.search(r'\[(.*?)\]', row.gene2)
        if g1:
            g1 = g1.group(1)
        if g2:
            g2 = g2.group(1)
        if g1 == g2:
            g_[i].append(g1)
        else:
            if g1:
                g_[i].append(g1)
            if g2:
                g_[i].append(g2)
    loop2gene = pd.DataFrame.from_dict(g_, orient="index").rename(columns = {0:"gene1", 1:"gene2"})
    loop2gene.index = loop_ann_df.index 
    return loop2gene

def open_loop(file, format = "bedpe"):
    """
    open loop file by automatic detecting its format
    supporting format: 
    - longrange:
    - bedpe: (not standard bedpe file, no strand information)
    """
    if format == "longrange":
        print("read longrange format")
        df = pd.read_table(file, sep = "\t", header = None)
        df.columns = ["chrom1", "start1", "end1", "connect"]
        df[["connect", "contact"]] = df["connect"].str.split(",", expand = True)
        df[["chrom2", "pos2"]] = df["connect"].str.split(":", expand = True)
        df[["start2", "end2"]] = df["pos2"].str.split("-", expand = True)
        df.drop(["connect", "pos2"], axis = 1, inplace = True)
        df["start2"] = df["start2"].astype(int)
        df["end2"] = df["end2"].astype(int)
    if format == "bedpe":
        print("read bedpe format")
        df = pd.read_table(file, sep = "\t", header = None)
        df.columns = bedpe_colname(df)
    if "contact" in df.columns: 
        df["contact"] = df["contact"].astype(int)
        return df[["chrom1", "start1", "end1", "chrom2", "start2", "end2", "contact"]]
    else:
        return df[["chrom1", "start1", "end1", "chrom2", "start2", "end2"]]

def concate_chrom(df):
    """
    operate on the original input df. 
    """
    df["name"] = df["chrom1"].astype(str) + ":" + df["start1"].astype(str) +  "-" + df["end1"].astype(str) + "_" + df["chrom2"].astype(str) +  ":" + df["start2"].astype(str) + "-" + df["end2"].astype(str)

def loop2anchor(loop_df):
    """
    Docstring for loop2anchor
    
    :param loop_df: loop information in a df (loop stored in column named "index")
    :return: df with chrom1/start1/end1/chrom2/start2/end2
    :rtype: Any
    """
    loop_df[["p1", 'p2']] = loop_df["index"].str.split("_", expand = True)
    loop_df[["chrom1", "p1"]] = loop_df["p1"].str.split(":", expand = True)
    loop_df[["start1", "end1"]] = loop_df["p1"].str.split("-", expand = True)
    loop_df[["chrom2", "p2"]] = loop_df["p2"].str.split(":", expand = True)
    loop_df[["start2", "end2"]] = loop_df["p2"].str.split("-", expand = True)
    loop_df["start1"] = loop_df["start1"].astype(int); loop_df["end1"] = loop_df["end1"].astype(int)
    loop_df["start2"] = loop_df["start2"].astype(int); loop_df["end2"] = loop_df["end2"].astype(int)
    loop_anchor = pd.concat([loop_df.rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"})[["chrom", "start", "end"]], loop_df.rename(columns = {"chrom2":"chrom", "start2":"start", "end2":"end"})[["chrom", "start", "end"]]], axis = 0)
    return loop_anchor 

def union_loop(file_list):
    df_list = list()
    for file in file_list:
        df = open_loop(file).drop_duplicates(keep = "first")
        new_index = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
        df_list.append(df[new_index].set_index(new_index))
    join_df = df_list[0].join(df_list[1:], how = "outer").reset_index()
    return join_df

def intersect_loops(file_list):
    df_list = list()
    for file in file_list:
        df = open_loop(file)
        new_index = [c for c in df.columns if c != "contact"]
        df_list.append(df.set_index(new_index).rename(columns = {"contact":file})) 
    return df_list[0].join(df_list[1:], how = "inner")

def bedpe_dist(bedpe_f):
    dist_list = list()
    for chunk in pd.read_table(bedpe_f, sep = "\t", header = None, names = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"], chunksize = 10000):
        for row in chunk.itertuples():
            if row.chrom1 == row.chrom2:
                pe_dist = row.end2 - row.start1 
                dist_list.append(pe_dist)
    return dist_list 

def cis_trans_bed(bed):
    """
    parse bedpe file for cis/trans reads stat
    """
    for chunk in pd.read_table(bed, header = None, chunksize=100):
        df_col = bedpe_colname(chunk)
        break 
    cis_pair = 0
    trans_pair = 0
    close_pair = 0
    bed_df_chunk = pd.read_table(bed, header = None, names = df_col, chunksize=500000)
    for chunk_df in bed_df_chunk:
        cis_pair += len(chunk_df.query("chrom1 == chrom2"))
        trans_pair += len(chunk_df.query("chrom1 != chrom2"))
        close_pair += len(chunk_df.query("chrom1 == chrom2 and start2 - start1 < 1000"))
    valid_pair_ratio = (cis_pair - close_pair)/(cis_pair + trans_pair)
    log_dict = {"cis":int(cis_pair), "trans":int(trans_pair), "cis_yield (cis/(cis+trans))":int(cis_pair)/(cis_pair + trans_pair), "short_distance (<1 kb)":int(close_pair), "long_distance (>1 kb)":(cis_pair-close_pair)/int(cis_pair), "cis_valid_yield (>1k cis/(cis+trans))":float(valid_pair_ratio)}
    return log_dict 

def bed2fasta(bed, ref, out = None):
    """
    Reformat bed into fasta file
    """
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

def pe2bed(l_df:pd.DataFrame, count_anchor_freq = False):
    """
    Transform bedpe to bed file (remove duplicated bed regions)
    
    """
    if not count_anchor_freq:
        loop_anchor = pd.concat([l_df[["chrom1", 'start1', "end1"]].rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"}), l_df[["chrom2", "start2", "end2"]].rename(columns = {"chrom2":"chrom", "start2":"start", "end2":"end"})], axis = 0).drop_duplicates(keep = "first", ignore_index=True)
    else:
        loop_anchor = pd.concat([l_df[["chrom1", 'start1', "end1"]].rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"}), l_df[["chrom2", "start2", "end2"]].rename(columns = {"chrom2":"chrom", "start2":"start", "end2":"end"})], axis = 0)
        loop_anchor["freq"] = loop_anchor.groupby(["chrom", "start", "end"]).transform("size")
        loop_anchor.drop_duplicates(keep = "first", ignore_index=True, inplace = True)
    return loop_anchor 

def bedpe_retrieve(bedpe_df:pd.DataFrame,region_bed_df:pd.DataFrame):
    """
    Given region, to search loops with either anchor side overlap it. 
    Input: 
    region: tuple with chrom, start, end 
    loop_df: pandas df in bedpe format 
    Return: 
    selected loop_df that overlaps the given region. 
    """
    # check if loop_df index duplicated or not 
    # assert len(set(loop_df.index)) == len(loop_df), print("loop not uniquely indexed")
    # get loop (left, chrom1, start1, end1) overlap region
    # left_loop = bf.select(loop_df.rename(columns = {"chrom1":"chrom", "start1":"start", "end1":"end"}), region)
    # # get loop (right, chrom2, start2, end2)
    # right_loop = bf.select(loop_df.rename(columns = {"chrom2":"chrom", "start2":"start", "end2":"end"}), region)
    # index_all = np.unique(np.concatenate((left_loop.index, right_loop.index)))
    # return loop_df.loc[index_all, :].copy()
    # update on retrieve bedpe file in a faster way by using pybedtools (11/25/25)
    import pybedtools
    if not "chrom1" in bedpe_df.columns:
        bedpe_df.columns = bedpe_colname(bedpe_df)
    anchor1 = bedpe_df.iloc[:, :3].drop_duplicates(keep = "first", ignore_index = True)
    anchor2 = bedpe_df.iloc[:, 3:6].drop_duplicates(keep = "first", ignore_index = True)
    overlap1 = pd.DataFrame([tuple(or1) for or1 in pybedtools.BedTool.from_dataframe(anchor1).intersect(pybedtools.BedTool.from_dataframe(region_bed_df), wa = True)], columns = ["chrom1", 'start1', "end1"]); overlap1["start1"] = overlap1["start1"].astype(int); overlap1["end1"] = overlap1["end1"].astype(int) # region r1 that overlap region of interest
    overlap2 = pd.DataFrame([tuple(or2) for or2 in pybedtools.BedTool.from_dataframe(anchor2).intersect(pybedtools.BedTool.from_dataframe(region_bed_df), wa = True)], columns = ["chrom2", "start2", "end2"]); overlap2["start2"] = overlap2["start2"].astype(int); overlap2["end2"] = overlap2["end2"].astype(int) # region r2 that overlap region of interest 
    pe1 = pd.merge(overlap1, bedpe_df, how = "left", on = ["chrom1", 'start1', "end1"]) # any bedpe that on anchor1 overlap regions of interest 
    pe2 = pd.merge(overlap2, bedpe_df, how = "left", on = ["chrom2", "start2", "end2"]) # any bedpe that on anchor2 overlap regions of interest 
    overlap_pe = pd.concat([pe1[["chrom1", 'start1', "end1", "chrom2", "start2", "end2"]], pe2[["chrom1", 'start1', "end1", "chrom2", "start2", "end2"]]], axis = 0).drop_duplicates(keep = "first", ignore_index = True).sort_values(["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    return overlap_pe

def bedpe_colname(bedpe_df):
    """
    Based on the shape of bedpe, return its column name 
    """
    bedpe_col = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2", 'name', "score", "strand1", "strand2"]
    return bedpe_col[:len(bedpe_df.columns)]

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
    """
    Find overlaps between two bed files 
    """
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

def loop_overlap(loop1:str, loop2:str, extend1 = None, extend2 = None):
    """
    compare two loop-interaction data, and return loops that overlapped between them in overlap_loop1 and overlap_loop2
    To return the overlapped loops
    input loop file format (bedpe):
    chrom1, start1, end1, chrom2, start2, end2

    """
    loop_df1 = pd.read_table(loop1, sep = "\t", header = None)
    loop_df2 = pd.read_table(loop2, sep = "\t", header = None)
    loop_df1 = loop_df1.iloc[:, :6]
    loop_df2 = loop_df2.iloc[:, :6]
    loop_df1.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2"]
    loop_df2.columns = ["chrom1", 'start1', "end1", "chrom2", "start2", "end2"]
    if len(loop_df1.query("start1 > start2")) != 0: 
        print("input loop1 pos1 and pos2 is not ordered")
        loop_df1 = loop_df1.query("start1 < start2")
    if len(loop_df2.query("start1 > start2")) != 0:
        print("input loop2 pos1 and pos2 is not ordered")
        loop_df2 = loop_df2.query("start1 < start2")
    if extend1 is not None:
        loop_df1["start1"] = loop_df1["start1"]-extend1
        loop_df1["end1"] = loop_df1["end1"]+extend1
        loop_df1["start2"] = loop_df1["start2"]-extend1
        loop_df1["end2"] = loop_df1["end2"]+extend1
    if extend2 is not None:
        loop_df2["start1"] = loop_df2["start1"]-extend2
        loop_df2["end1"] = loop_df2["end1"]+extend2
        loop_df2["start2"] = loop_df2["start2"]-extend2
        loop_df2["end2"] = loop_df2["end2"]+extend2
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

def rm_blacklist(bed:str, blacklist:str, schema = "bed5"):
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
    
def rmblacklist(peak_df:pd.DataFrame, black_df:pd.DataFrame):
    """
    Remove peaks that overlap blacklist regions 
    Inputs: 
    peak_df: peak df contains minimum columns of chrom, start, end
    blacklist: blacklist df contains minimum columns of chrom, start, end 
    Returns: 
    df: peaks that not overlapping with blacklist region
    """
    peak_overlap_black = bf.overlap(peak_df, black_df, how = "left")
    col2del = [c for c in peak_overlap_black.columns if c.endswith("_")]
    peak_no_black = peak_overlap_black[peak_overlap_black.isna().any(axis = 1)].drop(col2del, axis = 1)
    return peak_no_black