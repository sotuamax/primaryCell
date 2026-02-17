import cooler 
import numpy as np 
import pandas as pd 


def total_counts(cool, resolution):
    """
    Count total number of reads in cool file. 
    """
    if cool.endswith(".mcool") and resolution != None:
        cool_handle = cooler.Cooler(f"{cool}::resolutions/{resolution}")
    else:
        cool_handle = cooler.Cooler(cool)
    total_reads = 0
    for b in cool_handle.matrix(balance = False, sparse = True):
        total_reads += b.data.sum()
    return total_reads 

def cool_matrix(cool, resolution, region, balance = False):
    """
    region: tuple (three-element w/ chrom, start, end); 
    resolution: int; 
    cool: str for cool input file; 
    """
    if cool.endswith(".mcool") and resolution != None:
        cool_handle = cooler.Cooler(f"{cool}::resolutions/{resolution}")
    else:
        cool_handle = cooler.Cooler(cool)
    clr_matrix = cool_handle.matrix(balance = balance).fetch(region, region)
    return clr_matrix 

def global_exp(clr_matrix:np.ndarray):
    return get_diag(clr_matrix, i = 0).mean()

def local_exp_interaction(clr_matrix:np.ndarray, offset:int, group_region:tuple, wl:list):
    """
    retrieve expected counts for given distance of peak pairs. (average for screened region as expected value)
    Inputs:
    bp_obs_filtered: dataframe of observed count after filtering 
    clr: cooler handle 
    resolution: cool resolution 
    window: number of bins for expected value calculation (at local scale, default: 0.5 M)
    Returns:
    dataframe: observed and expected contacts of a window region at a given distance

    Note: mean is for any diagonal region that overlap the requested site (it may cross multiple bins).
    """
    bg_mean = list()
    a1, b1, a2, b2 = group_region # 
    ar1, br1, ar2, br2 = np.array([a1,b1,a2,b2]) - offset
    for w in wl:
        ar1_w, ar2_w = np.where(np.array([ar1, ar2]) - w < 0, 0, np.array([ar1, ar2]) - w)
        br1_w, br2_w = np.where(np.array([br1, br2]) + w > clr_matrix.shape[0], clr_matrix.shape[0], np.array([br1, br2]) + w)
        obs_mat = clr_matrix[ar1_w:br1_w, ar2_w:br2_w] # matrix[row_start:row_end, col_start:col_end]
        obs_diag = get_diag(obs_mat, i = 0) # the obs_mat is in a oddxodd matrix
        bg_mean_w = np.concatenate([obs_diag[:w],obs_diag[w+1:]]).mean()
        bg_mean.append(bg_mean_w)
    exp_df = pd.DataFrame([[a1,b1,a2,b2] + bg_mean], columns = ["start1", "end1", "start2", "end2"] + [f"exp_{w}" for w in wl])
    return exp_df

def interaction_count(clr_matrix, offset:int, resolution:int, interaction_site:pd.DataFrame):
    """
    Count reads at query interaction sites
    """
    interaction_pos = interaction_site[["start1", "end1", "start2", "end2"]]//resolution - offset 
    interaction_cnt = [clr_matrix[row.start1:row.end1+1, row.start2:row.end2+1].sum() for row in interaction_pos.itertuples()]
    # print(clr_matrix.shape)
    # interaction_site["count"] = interaction_cnt
    return interaction_cnt

def obs_interaction(clr_matrix, peak_df, offset:int, chrom:str, bin1_pos:int, bin2:range, stop_sign_of_0:int, resolution: int):
    """
    Take fragments of chromosomes and peak data, search for the observed contacts at specified bins in the cool_matrix.
    Inputs:
    clr: cooler handle storing observation matrix of binned windows across genome
    region: tuple of three values (chrom, start, end)
    peak_df: dataframe of peaks
    stop_sign_of_0: the maximum number of 0s reached before stopping the paired peak search 
    resolution: resolution of input cool
    bin1_pos: upstream anchor position
    bin2: list of all possible downstream anchors to examine (stop based on the stop_sign_of_0)
    Returns:
    dataframe: observed counts for paired peak anchors, with columns chrom/bin1/bin2/obs (bin1 < bin2, and the position = bin*resolution)
    """
    bp_list = list()
    bp_0_cutoff = 0
    a = (peak_df.loc[bin1_pos, "start"]//resolution, peak_df.loc[bin1_pos, "end"]//resolution+1)
    a2 = peak_df.iloc[bin2, :]
    #
    for bpv in zip(a2["start"]//resolution, a2["end"]//resolution + 1):
    # a = np.array(bin1_pos)-offset
    #for bpv in bin2:
        #b = bpv-offset
        #if abs(a-b) < 10:
        try:
            bp_obs = clr_matrix[a[0]-offset:a[-1]-offset, bpv[0]-offset:bpv[-1]-offset][0][0] # in python, it's always left close, right open (exclusive)
            #else: # matrix slicing: matrix[row_start:row_end, col_start:col_end]
            #bp_obs = clr_matrix[a-1:a+2, b-1:b+2].max()
            if bp_obs > 0: 
                bp_list.append((bpv[0], bpv[-1], bp_obs))
            if bp_obs == 0:
                bp_0_cutoff += 1
            if bp_obs != 0:
                bp_0_cutoff = 0
            if bp_0_cutoff >= stop_sign_of_0:
                break
        except Exception as e:
            print(e)
            pass
    bp_df = pd.DataFrame(bp_list, columns = ["start2", "end2", "obs"]); bp_df["start1"] = a[0];bp_df["end1"] = a[-1]; bp_df["chrom1"] = chrom;bp_df["chrom2"] = chrom
    return bp_df[["chrom1", "start1", "end1", "chrom2", "start2", "end2", "obs"]]

def get_diag(arr, i=0):
    """
    adapted from cooltools; 
    Get the i-th diagonal of a matrix.
    This solution was borrowed from
    http://stackoverflow.com/questions/9958577/changing-the-values-of-the-diagonal-of-a-matrix-in-numpy
    """
    arr = np.array(arr)
    return arr.ravel()[
        max(i, -arr.shape[1] * i) : max(0, (arr.shape[1] - i))
        * arr.shape[1] : arr.shape[1]
        + 1
    ]

def parser_cool(cool, resolution = None, balance = False):
    """
    Parse cool file to get the handle for random access 
    """
    if cool.endswith(".mcool") and resolution != None:
        cool_handle = cooler.Cooler(f"{cool}::resolutions/{resolution}")
    else:
        cool_handle = cooler.Cooler(cool)
    binsize = cool_handle.binsize 
    # print(f"resultion: {binsize}")
    if balance == True:
        cooler.balance_cooler(cool_handle, cis_only = True, chunksize = 100000000, store = True)
    return cool_handle

def raw_pixel(cool, resolution = None, filter = None):
    """
    Parse cool to get chromosome bin and bin-bin pixel. 
    When filter = 'cis', intra-chromosomal contacts are retained. 
    """
    if resolution != None:
        cool_handle = cooler.Cooler(f"{cool}::resolutions/{resolution}")
    else:
        cool_handle = cooler.Cooler(cool)
    bin = cool_handle.bins()[:][["chrom", "start", "end"]] # the returned bins are always sorted 
    #bin["i"] = bed_df.index
    # bed_out = bed_df[bed_df["chrom"].isin(["chrY", "chrM"])]["i"].tolist()
    pixel = cool_handle.pixels()[:] # upper triangle, both cis and trans included (intra- and inter- chromosome)
    if filter == "cis":
        chrom_range = pd.concat([bin.groupby("chrom")["end"].idxmin().rename("min"), bin.groupby("chrom")["end"].idxmax().rename('max')], axis = 1)
        # bin1_id <= bin2_id (always valid)
        pixel = pd.concat([pixel.query("bin1_id >= @row.min and bin2_id <= @row.max") for row in chrom_range.itertuples()], axis = 0)
    return bin, pixel

def pixel_chrom(bin, pixel):
    """
    split pixel file based on chromosome in bin into a dictionary 
    """
    chrom_range = pd.concat([bin.groupby("chrom")["end"].idxmin().rename("min"), bin.groupby("chrom")["end"].idxmax().rename("max")], axis = 1)
    pixel_dict = {row.Index:pixel.query("@row.max >= bin1_id >= @row.min") for row in chrom_range.itertuples()}
    return pixel_dict

