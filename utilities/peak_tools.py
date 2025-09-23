import bioframe as bf 
import pandas as pd 
import numpy as np 
import subprocess

def peakcall(file):
    """
    default genome: hs (hg38)
    """
    if file.endswith(".bam"):
        output_name = file.split(".bam")[0]
        command = f"macs3 callpeak -q 0.01 --keep-dup all -g hs -t {file} -f BAM --nomodel --extsize 200 --shift -100 -n {output_name}"
    if file.enswith(".bed"):
        output_name = file.split(".bed")[0]
        command = f"macs3 callpeak -q 0.01 --keep-dup all -g hs -t {file} -f BED --nomodel --extsize 200 --shift -100 -n {output_name}"
    subprocess.call(command, shell = True)


def peak_overlap_label(peak, anno):
    """
    Find overlap between peak and compartment. 
    Note that the overlap was evaluated on peak summit, as peak may overlay more than one compartment. 
    Return:
    peak dataframe annotated with its overlapping compartment and compartment score. 
    """
    peak_df = bf.read_table(peak, schema = "narrowPeak")
    if len(set(peak_df["name"])) != len(peak_df):
        peak_df["name"] = range(len(peak_df))
    summit_df = peak_df[["chrom", "start", "relSummit", "name"]].copy()
    summit_df["start"] = summit_df["start"] + summit_df["relSummit"]
    summit_df["end"] = summit_df["start"] + 1
    summit_df = summit_df[["chrom", "start", "end", "name"]]
    # summit_all = len(summit_df["name"])
    # read annotation file
    anno_df = bf.read_table(anno, schema = "bed4", comment = "#")
    anno_df.rename(columns = {"name":"label"}, inplace = True)
    # compartment_df["compartment"] = np.where(compartment_df["score"] > 0, "A", "B")
    # find overlap 
    summit_overlap = bf.overlap(summit_df, anno_df, how = "left")
    peak_anno = pd.merge(peak_df, summit_overlap[["name", "label_"]], on = ["name"])
    return peak_anno.rename(columns = {"label_":"label"})

def peak_parse(peak:str, chrom = None):
    """
    Reformat peak file by adding coordinates of summits, and drop columns (strand, p, q, name, fc)
    Input: 
    peak (dataframe): peak file name (in narrowPeak format)
    chrom (list): list of chromosomes to filter
    Return: 
    dataframe: re-format peak input dataframe, containing columns chrom/start/end/summit
    """
    import bioframe as bf
    if peak.endswith(".narrowPeak"):
        peak_df = bf.read_table(peak, schema = "narrowPeak")
        #peak_df["summit"] = peak_df["start"] + peak_df["relSummit"]; peak_df.drop(["relSummit", "strand", "name", "score"], axis = 1, inplace = True)
    if peak.endswith(".bed"):
        peak_df = pd.read_table(peak, names=["chrom", "start", "end"])
    if chrom != None:
        peak_df = peak_df.query("chrom in @chrom").copy()
    peak_df.sort_values(["chrom", "start"], inplace = True)
    peak_df.index = range(len(peak_df))
    return peak_df

def generate_paired_peak(peak_df:pd.DataFrame, region = None):
    """
    Based on the bin-interactions, to assign peak interactions that overlaying bins.
    Note: the distance between the the paired anchors should be within range (min_dist, max_dist). 
    default: min_dist - 5000; max_dist - 10 Mb 
    """
    if region is not None: 
        view_peak = bf.select(peak_df, region) # index maintained 
    else:
        view_peak = peak_df
    if len(set(view_peak['chrom'])) != 1:
        return {}
    # view_peak.sort_values(by = ["chrom", "start", "end"], ignore_index = True, inplace = True)
    # bins = sorted(set(zip(view_peak["start"]//resolution, view_peak["end"]//resolution+1)))
    # bins = sorted(set(view_peak["summit"]//resolution))
    # if do not use summit, use peak range, and within the range, the maximum observed count used. 
    last_peak = view_peak.index.max(); first_peak = view_peak.index.min()
    pairs = {i:range(i+1, last_peak+1) for i in range(first_peak, last_peak)}
    # paired_anchor = {bins[i]:np.array([b2 for b2 in bins[i+1:] if (b2[0]-bins[i][0] < max_dist//resolution) and (b2[0]-bins[i][0] > min_dist//resolution)]) for i in range(len(bins)-1)}
    # from itertools import combinations 
    # paired_anchor = pd.DataFrame([(a1,a2) for a1,a2 in combinations(bins, 2) if (abs(a1[0]-a2[0]) > min_dist//resolution) and (abs(a1[0]-a2[0]) < max_dist//resolution)], columns = ["anchor1", "anchor2"])
    # paired_anchor["chrom"] = region[0]
    return pairs

def rmblacklist(peak_df, black_df):
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

def peak2saf(narrowPeak, rename = True):
    """
    Prepare narrowPeak format into SAF format for featureCount input
    https://subread.sourceforge.net/featureCounts.html
    """
    peak_df = bf.read_table(narrowPeak, schema = "narrowPeak")
    peak_df.sort_values(by = ["chrom", "start", "end"], inplace = True)
    if rename:
        peak_saf = peak_df[["chrom", "start", "end", "strand"]].drop_duplicates(keep = "first").copy()
        peak_saf["ID"] = range(len(peak_saf))
        peak_saf["ID"] = "peak_" + peak_saf["ID"].astype(str)
        return peak_saf[["ID", "chrom", "start", "end", "strand"]]
    else:
        return peak_df[["name", "chrom", "start", "end", "strand"]]
 
def cluster_(peak_df, min_dist=None):
    """
    Return peak interval cluster
    Input: 
    peak: dataframe
    min_dist: (int) minimum distance used for clustering
    Return:
    peak intervals with cluster ID labeled when its distance less than min_dist
    """
    peak_cluster = bf.cluster(peak_df, min_dist = min_dist)
    peak_cluster.drop(["cluster_start", "cluster_end"], inplace = True, axis = 1)
    peak_cluster.sort_values(by = "cluster", inplace = True, ignore_index = True)
    return peak_cluster 

def consensus_range(cluster_df, pprop):
    # return consensus start, end, centroid, core start, core end
    min_sample = len(cluster_df)*pprop
    cpeak_array = np.concatenate([np.arange(row.start, row.end) for row in cluster_df.itertuples()], axis = 0)
    cpeak_array, cpeak_coverage = np.unique(cpeak_array, return_counts = True)
    cpeak_region = cpeak_array[cpeak_coverage >= min_sample].reshape(-1,)
    # 
    if len(cpeak_region) > 0:
        if sum(np.diff(cpeak_region) > 10) == 0:
            return (cpeak_region.min(), cpeak_region.max(), int(np.median(cluster_df["summit"])), cluster_df["summit"].min().astype(int), cluster_df["summit"].max().astype(int))
        else:
            while pprop > 0:
                pprop -= 0.1
                recursion_return = consensus_range(cluster_df, pprop)
                if len(recursion_return) > 0:
                    return recursion_return
                    break 
    if len(cpeak_region) == 0:
        return (cluster_df["start"].min().astype(int), cluster_df["end"].max().astype(int), int(np.median(cluster_df["summit"])), cluster_df["summit"].min().astype(int), cluster_df["summit"].max().astype(int))

def consensusCall(cluster_chr, pprop):
    peak_nest = cluster_chr.copy()
    # late on, I need to add extra step to examine summit cluster nested with different peaks (not cluster)
    cluster_dict = dict()
    for cluster1, cluster_1 in peak_nest.groupby(["cluster"]):
        for cluster2, cluster_2 in cluster_1.groupby("cluster_"):
            if len(cluster_2) > 1:
                # this is when summits form cluster
                cluster_dict[cluster2] = consensus_range(cluster_2, pprop)
            else:
                # summit is orphan
                cluster_dict[cluster2] = (cluster_2["start"].item(), cluster_2["end"].item(), cluster_2["summit"].item(), cluster_2["summit"].item(), cluster_2["summit"].item()+1)
    cluster_consensus = pd.DataFrame.from_dict(cluster_dict, orient = "index", columns = ["consensus_start", "consensus_end", "centroid", "core_start", "core_end"], dtype = int)
    cluster_consensus["cluster_"] = cluster_consensus.index
    return cluster_consensus

def refine_consensus(cs_df):
    # refine is performed only on overlapping peak regions (peaks not overlap with other peaks are unchanged)
    cs_sub = cs_df[["chrom", "consensus_start", "consensus_end", "core_start", "core_end", "cluster_"]].copy().rename(columns = {"consensus_start":"start", "consensus_end":"end"}).drop_duplicates(keep = "first", ignore_index = True)
    cs_cluster = bf.cluster(cs_sub, min_dist = 0) 
    cluster_arr, cluster_cnt = np.unique(cs_cluster["cluster"], return_counts = True)
    # find overlapping consensus peaks (at least 2)
    css_clustered = cs_cluster[cs_cluster["cluster"].isin(cluster_arr[np.where(cluster_cnt != 1)].tolist())].copy()
    css_clustered["refined_start"] = np.where(css_clustered["start"]<(css_clustered["core_start"]-75), css_clustered["core_start"]-75, css_clustered["start"])
    css_clustered["refined_end"] = np.where(css_clustered["end"]>(css_clustered["core_end"]+75), css_clustered["core_end"]+75, css_clustered["end"])
    # to add single peak region (unchanged) to df
    cs_single = cs_cluster[cs_cluster["cluster"].isin(cluster_arr[np.where(cluster_cnt == 1)].tolist())].copy()
    cs_single["refined_start"] = cs_single["start"]
    cs_single["refined_end"] = cs_single["end"]
    refined_bed = pd.concat([css_clustered[["refined_start", "refined_end", "cluster_"]], cs_single[["refined_start", "refined_end", "cluster_"]]], axis = 0)
    refined_bed.sort_values(["cluster_"], inplace = True, ignore_index= True)
    return refined_bed

def peak_overlay(consensus_chrom):
    cc_df = list()
    for c, consensus_df in consensus_chrom.groupby(["chrom", "cluster"]):
        if len(consensus_df) > 1:
            chrom = c[0]
            # clustt = c[-1]
            c_range, c_coverage = np.unique(np.concatenate([np.arange(row.start,row.end) for row in consensus_df.itertuples()], axis = 0), return_counts = True)
            c_multi_coverage = c_range[c_coverage > 1]
            # get the very start and end of the overlay peaks region 
            c_start = c_multi_coverage.min().reshape(-1,)
            c_end = c_multi_coverage.max().reshape(-1,)
            multi_index = np.argwhere(np.diff(c_multi_coverage) != 1)
            # when the peak region is not continuous 
            if len(multi_index) > 0:
                c_upper = c_multi_coverage[:-1][multi_index].reshape(-1,)
                c_lower = c_multi_coverage[1:][multi_index].reshape(-1,)
                c_mid = np.concatenate([c_upper, c_lower], axis = 0)
                c_over_order = np.sort(np.concatenate([c_start, c_mid, c_end], axis = 0))
            else:
                c_over_order = np.concatenate([c_start, c_end])
            over_start = c_over_order[::2]
            over_end = c_over_order[1::2]
            c_df = pd.DataFrame({"chrom":chrom, "start":over_start, "end":over_end})
            cc_df.append(c_df)
        else:
            pass 
    cc_df = pd.concat(cc_df, axis = 0)
    return cc_df 

def overlay2peaks(peak_over_chrom, consensus_chrom, min_dist):
    # overlay region and its contributing peaks 
    overlay_peaks = bf.overlap(peak_over_chrom, consensus_chrom[["chrom", "start", "end", "name"]], how = "inner")
    # resummit based on center and range 
    summit_tmp = consensus_chrom[["chrom", "centroid", "name"]].copy()
    summit_tmp["start"] = summit_tmp["centroid"]-min_dist
    summit_tmp["end"] = summit_tmp["centroid"]+min_dist+1
    overlay_summits = bf.overlap(peak_over_chrom, summit_tmp[["chrom", "start", "end", "name"]], how = "inner")
    # overlay-summit as reference for peaks merging 
    merged_ps = pd.merge(overlay_peaks, overlay_summits, on = ["chrom", "start", "end"], suffixes = ("peak", "summit"), how = "inner")
    # peak_dict = dict()
    peak_combine = list()
    peak_drop = list()
    for c, groupdf in merged_ps.groupby(["chrom", "start", "end"]):
        peak_n = np.unique(groupdf["name_peak"]); summit_n = list(np.unique(groupdf["name_summit"]))
        if len(peak_n) == len(summit_n):
            # combine peaks 
            peak_combine.append(tuple(peak_n))
        if len(peak_n) > 1 and len(summit_n) == 1:
            peak_drop = peak_drop + summit_n
        if len(summit_n) > 1:
            peak_min = min(peak_n); peak_max = max(peak_n)
            if peak_min in summit_n:
                summit_n.remove(peak_min)
            if peak_max in summit_n:
                summit_n.remove(peak_max)
            peak_drop = peak_drop + summit_n
    return peak_combine, peak_drop



