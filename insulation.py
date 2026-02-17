#!/usr/bin/env python3
"""
Slide a triangle-shape window along the genome (one corner on the diagonal of the matrix), and score for the sum of contacts within the window for each position. At certain location, the score is significantly lower than its surrounding region (reflecting lowered contact frequencies between upstream and downstream loci), this position is referred to as boundary position.

Example: 
mpiexec -n 4 insulation.py sample.cool -r 25000 -w 10 -o sample

output: 
size: size of the square submatrix used to calculate the contact mean (score); 
score: mean of contact in submatrix; 
slope: the derivative at the position for the score (by taking surrounding 5 scores and get the differences between upstream 5 and downstream 5); 
valley: the point where derivative as zero and where the local minima reached; 
strength: difference between the valley point derivative and its neighboring point derivative; 

update: 
derivative: 3 bins on each side -> as an argument (-dv, default : 3)
minima: the score is minima in surrounding 5 bin on each side -> as an argument (-vs, default : 10); minma screen should examine regions beyond the shifting site (-1, +1)
surrounding region should be at least >fd times of the local minima; 

the final output in bedGraph format: 
score represents the foldchange between the valley point and its surrounding region's local maxima contact frequencies.

"""

# import cooler 
import bioframe as bf 
import numpy as np
import pandas as pd
import argparse
# from utilities.cool_tools import raw_pixel
from utilities.cool_tools import cool_matrix
from utilities.chrom_func import chrom_size
from utilities.misc import ignore_warning, timeit
from tqdm.auto import tqdm
import os 
from scipy.stats import zscore
import time

ignore_warning()

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", 
                                     formatter_class = argparse.RawDescriptionHelpFormatter, 
                                     description="insulation.py input.cool -r 25000 -w 10 -o output")
    parser.add_argument("cool", help = "cool input file")
    parser.add_argument("-balance", "--balance", action = "store_true", help = "use balanced matrix data stored in cool.")
    parser.add_argument("-ref", "--reference", help = "reference genome", default = "hg38")
    parser.add_argument("-r", "--resolution", help = "resolution", type = int)
    parser.add_argument("-w", "--window", help = "window size used to calculate submatrix score (mean of contacts)", type = int)
    parser.add_argument("-ignore_diags", "--ignore_diags", default = 1, type = int, help = "ignore diagon region by this number.")
    parser.add_argument("-boundary", "--boundary", action = "store_true", help = "report boundary based on insulation profile")
    parser.add_argument("-smooth", "--smooth", action = "store_true", help = "smooth contact signal")
    parser.add_argument("-q", "--quantile", type = float, default = 0.5, help = "quantile range used to filter boundaries")
    parser.add_argument("-dv", "--derivative", type = int, default = 2, help = "bin number used to calculate derivative for center position")
    #parser.add_argument("-vs", "--valley_screen", type = int, default = 10, help = "To find minima of a screening size (per side)")
    #parser.add_argument("-fc", "--foldchange", default = 1.5, type = float, help = "within the screen range, the maximum value is *fd* times of the minima score. ")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def valley_point(arr_slope, arr_value, dev):
    """
    pick for two consecutive derivatives crossing zero (which can be peak/valley); 
    trend from negative derivative to positive derivative as the region surrounding local minima (valley);

    flank: bin number, used to calculate the contact mean upstream/downstream the valley point
    """
    # find shifting signs for neighbor bins, return index
    zero_crossing = np.where((np.sign(arr_slope)[:-1] * np.sign(arr_slope)[1:] < 0))[0] 
    # find shifting derivatives (from negative to positive) where match to a valley, and return index
    valley_zero = zero_crossing[arr_slope[zero_crossing] < arr_slope[zero_crossing+1]] 
    # valley_diff = np.abs(arr_slope[valley_zero]-arr_slope[valley_zero+1])
    # valley_Q1 = np.quantile(valley_diff, 0.25)
    # valley_zero = valley_zero[valley_diff > valley_Q1]
    # find neighboring points surrounding derivative at 0, and then returning the one point closest to 0 
    valley_index = np.unique(np.array([i - dev + np.argmin(arr_value[i-dev:i+dev+1]) for i in valley_zero if i > 2])).astype(int) # key improve to expand region for searching of minima around the derivative shifting site 
    valley_index = valley_index[(valley_index >= 0) & (valley_index < len(arr_value))]
    # to avoid same index occurs twice, add unique for the index value 
    valley = np.zeros(len(arr_slope)) #valley_min[valley_index[1:-1]] = 1 # replace 1 at index for minima 
    valley[valley_index] = 1
    return valley 

def main():
    from mpi4py import MPI
    from mpi4py.util import pkl5
    # to overcome the overflow error when comm data > 2 GB
    comm = pkl5.Intracomm(MPI.COMM_WORLD) 
    rank = comm.Get_rank()
    size = comm.Get_size()
    args = args_parser()
    cool = args.cool 
    resolution = args.resolution;  win = args.window; ignore_diag = args.ignore_diags
    ref = args.reference 
    out = args.output
    ref_size = chrom_size(ref); ref_size.sort_values("chrom", inplace = True); ref_size.index = range(len(ref_size))
    if rank == 0:
        #start = time.time()
        win1 = np.ones((win,win))
        drange = args.ignore_diags+1
        for d in range(drange):
            for dd in range(d+1):
                win1[win-1-d+dd, dd] = 0
        print("Distribute work ...")
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        for i in ref_size.index.tolist():
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
        print("Contact resolution: ", resolution)
    else:
        worker_tasks = None
        win1 = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    win1 = comm.bcast(win1, root = 0)
    pixel_list = list()
    for i in tqdm(worker_tasks[rank], desc = f"Parallel: {rank}"): # each core process one chrom-arm 
        chr,start,end = ref_size.loc[i,"chrom"], ref_size.loc[i,"start"], ref_size.loc[i,"end"]
        i_offset = start//resolution
        region_mtx = cool_matrix(cool, resolution, region = (chr, start, end), balance=args.balance)
        for center in tqdm(range(start//resolution+win-i_offset, end//resolution-win-i_offset, 1), desc = f"Process:{chr}:{start}-{end}", position = 1, leave = False): # center is the diagonal position on the matrix 
            region_submtx = region_mtx[(center - win):(center),(center+1):(1+center + win)]*win1
            #region_submtx = region_mtx[(center-win-ignore_diag+1):center-ignore_diag+1, (center+ignore_diag):(center+win+ignore_diag)] 
            # diamond region with one end at the diagonal pos (left include, right exclude, the diagonal site is not included)
            # q1_contact, q2_contact = np.quantile(region_submtx, q1), np.quantile(region_submtx, q2)
            #region_submtx = region_submtx[(region_submtx > q1_contact) & (region_submtx < q2_contact)]
            count = region_submtx.sum()
            #size = (region_submtx > 0).sum()
            pixel_list.append((chr, (center+i_offset)*resolution, (center+i_offset+1)*resolution, count))
    pixel_df = pd.DataFrame(pixel_list, columns = ["chrom", "start", "end", "value"])
    pixel_all = comm.gather(pixel_df)
    if rank == 0:
        print("Process window: ", f"{win}x{win}")
        print("Collecting all insulation score ...")
        pixel_all_df = pd.concat(pixel_all, axis = 0)
        pixel_all_df.sort_values(["chrom", "start"], inplace = True, ignore_index = True)
        #print("Write raw insulation data (contacts) ...")
        #pixel_all_df[["chrom", "start", "end", "value"]].to_csv(f"{out}.contact.bedGraph", sep = "\t", header = False, index = False)
        print("Write normalized insulation data (z-scale) ...")
        pixel_all_df["z"] = pixel_all_df.sort_values(["chrom", "start"], ignore_index = True).groupby("chrom")["value"].transform(zscore)
        pixel_all_df[["chrom", "start", "end", "z"]].to_csv(f"{out}.contact.z.bedGraph", sep = "\t", header = False, index = False)
        from scipy.signal import savgol_filter, find_peaks
        ww = args.derivative*2+1
        pixel_all_df["smooth_z"] = pixel_all_df.sort_values(["chrom", "start"], ignore_index = True).groupby("chrom")["z"].transform(lambda x: savgol_filter(x.values, window_length=ww, polyorder = 3))
        pixel_all_df["dx"] = pixel_all_df.sort_values(["chrom", "start"], ignore_index = True).groupby("chrom")["z"].transform(lambda x: savgol_filter(x.values, window_length=ww, polyorder = 3, deriv = 1))
        pixel_all_df["ddx"] = pixel_all_df.sort_values(["chrom", "start"], ignore_index = True).groupby("chrom")["z"].transform(lambda x: savgol_filter(x.values, window_length=ww, polyorder = 3, deriv = 2))
        print("Write 1st derivative of contact value ...")
        pixel_all_df[["chrom", "start", "end", "dx"]].to_csv(f"{out}.dx.bedGraph", sep = "\t", header = False, index = False)
        print("Write 2nd derivative of contact value ...")
        pixel_all_df[["chrom", "start", "end", "ddx"]].to_csv(f"{out}.ddx.bedGraph", sep = "\t", header = False, index = False)
        pixel_all_df["valley1"] = (pixel_all_df.sort_values(["chrom", "start"], ignore_index = True)
                                    .groupby("chrom", group_keys = False)["dx"]
                                    .transform(lambda x: valley_point(x.values, pixel_all_df.loc[x.index, "smooth_z"], args.derivative)))
        
        def peak_features(x):
            peaks, props = find_peaks(x, prominence=0, distance=args.derivative*2+1)
            mask = np.zeros(len(x), dtype=int)
            prom = np.zeros(len(x), dtype=float)
            mask[peaks] = 1
            prom[peaks] = props["prominences"]
            return mask, prom
        valley_res = pixel_all_df.sort_values(["chrom", "start"], ignore_index = True).groupby("chrom")["smooth_z"].apply(lambda x: peak_features(-x))
        peak_mask = np.concatenate([r[0] for r in valley_res.values])
        peak_prom = np.concatenate([r[1] for r in valley_res.values])
        pixel_all_df["valley2"] = peak_mask
        pixel_all_df["prom"] = peak_prom
        print("Write contact information ...")
        pixel_all_df.to_csv(f"{out}.txt", sep = "\t", header = True, index = False)
        if args.boundary:
            print("Report boundary ...")
            ddx_cut1 = np.quantile(pixel_all_df.query("valley1 == 1")["ddx"], args.quantile)
            ddx_cut2 = np.quantile(pixel_all_df.query("valley2 == 1")["ddx"], args.quantile)
            ddx_cut = np.quantile(pixel_all_df.query("valley1 == 1 and valley2 == 1")["ddx"], args.quantile)
            pixel_all_df.query("valley1 == 1 and ddx > @ddx_cut1")[['chrom', "start", "end", "ddx"]].to_csv(f"{out}.boundary1.bedGraph", sep = "\t", header = False, index = False)
            pixel_all_df.query("valley2 == 1 and ddx > @ddx_cut2")[['chrom', "start", "end", "ddx"]].to_csv(f"{out}.boundary2.bedGraph", sep = "\t", header = False, index = False)
            pixel_all_df.query("valley1 == 1 and valley2 == 1 and ddx > @ddx_cut")[['chrom', "start", "end", "ddx"]].to_csv(f"{out}.boundary.bedGraph", sep = "\t", header = False, index = False)
        #print(timeit(start))

if __name__ == "__main__":
    main()


