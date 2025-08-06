#!/usr/bin/env python3
"""
Slide a triangle-shape window along the genome (one corner on the diagonal of the matrix), and score for the sum of contacts within the window for each position. At certain location, the score is significantly lower than its surrounding region (reflecting lowered contact frequencies between upstream and downstream loci), this position is referred to as boundary position.

Example: 
insulation.py sample.cool -r 25000 -w 10 -o sample

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

"""

# import cooler 
import bioframe as bf 
import numpy as np
import pandas as pd
import argparse
# from utilities.cool_tools import raw_pixel
from utilities.cool_tools import cool_matrix
from utilities.chrom_func import chrom_arms
from utilities.misc import ignore_warning
from tqdm.auto import tqdm
import os 

ignore_warning()

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", 
                                     formatter_class = argparse.RawDescriptionHelpFormatter, 
                                     description="insulation.py input.cool -r 25000 -w 10 -o output")
    parser.add_argument("cool", help = "cool input file")
    parser.add_argument("-ref", "--reference", help = "reference genome", default = "hg38")
    parser.add_argument("-r", "--resolution", help = "resolution", type = int)
    parser.add_argument("-w", "--window", nargs = "+", help = "window size used to calculate submatrix score (mean of contacts)", type = int)
    parser.add_argument("-dv", "--derivative", type = int, default = 3, help = "bin size used to calculate derivative for center position")
    #parser.add_argument("-vs", "--valley_screen", type = int, default = 10, help = "To find minima of a screening size (per side)")
    parser.add_argument("-fc", "--foldchange", default = 1.5, type = float, help = "within the screen range, the maximum value is *fd* times of the minima score. ")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def derivative(arr, side):
    """
    Caculate derivative of array by taking the side values for their difference (centered on one value). 
    """
    new_arr = arr.copy()
    for i in range(side, len(arr)-side):
        before = arr[i-side:i]
        after = arr[i+1:i+1+side]
        new_arr[i] = np.mean(after) - np.mean(before)
    return new_arr 

def valley_point(arr_slope, arr_value):
    """
    pick for two consecutive derivatives crossing zero (which can be peak/valley); 
    trend from negative derivative to positive derivative as the region surrounding local minima (valley);
    """
    # valley_max_downstream = np.zeros(len(arr_slope))
    #strength_arr = np.zeros(len(arr_slope))
    #minima_arr = np.zeros(len(arr_slope))
    signs = np.sign(arr_slope) # get the sign of slope
    # find shifting signs for neighbor bins, return index
    zero_crossing = np.where(signs[:-1] * signs[1:] < 0)[0] 
    # find shifting derivatives (from negative to positive) where match to a valley, and return index
    valley_zero = zero_crossing[arr_slope[zero_crossing] < arr_slope[zero_crossing+1]] 
    # find neighboring points surrounding derivative at 0, and then returning the one point closest to 0 
    valley_index = np.unique(np.array([0] + [i - 2 + np.argmin(arr_value[i-2:i+4]) for i in valley_zero if i > 2] + [len(arr_slope)])).astype(int) # key improve to expand region for searching of minima around the derivative shifting site 
    # to avoid same index occurs twice, add unique for the index value 
    # for a,b in list(zip(valley_index[:-1], valley_index[1:])):
    #     print(a,b, np.max(arr_value[a:b]))
    upstream_max = np.array([np.max(arr_value[a:b]) for a,b in list(zip(valley_index[:-2], valley_index[1:-1]))])
    downstream_max = np.array([np.max(arr_value[a:b]) for a,b in list(zip(valley_index[1:-1], valley_index[2:]))])
    valley_min = np.zeros(len(arr_slope)); valley_min[valley_index[1:-1]] = 1 # replace 1 at index for minima 
    valley_upstream = np.zeros(len(arr_slope)); valley_upstream[valley_index[1:-1]] = upstream_max
    valley_downstream = np.zeros(len(arr_slope)); valley_downstream[valley_index[1:-1]] = downstream_max
    # strength_arr[valley_index2] = arr_slope[valley_zero+1] - arr_slope[valley_zero] # add strength value at index for valley (slope change)
    return valley_min, valley_upstream, valley_downstream #, strength_arr

def update_max(valley_value, arr_value):
    valley_index = np.where(valley_value == 1)[0]
    valley_upmax = np.zeros(len(arr_value)); valley_downmax = np.zeros(len(arr_value))
    if len(valley_index) > 0:
        if valley_index[0] != 0:
            valley_index = np.insert(valley_index, 0,0)
            valley_index = np.insert(valley_index, -1, len(valley_value))
            valley_index = np.unique(valley_index)
        valley_upmax[valley_index[1:-1]] = np.array([np.max(arr_value[a:b]) for a,b in list(zip(valley_index[:-2], valley_index[1:-1]))])
        valley_downmax[valley_index[1:-1]] = np.array([np.max(arr_value[a:b]) for a,b in list(zip(valley_index[1:-1], valley_index[2:]))])
    return valley_upmax, valley_downmax 

def main():
    from mpi4py import MPI
    from mpi4py.util import pkl5
    # to overcome the overflow error when comm data > 2 GB
    comm = pkl5.Intracomm(MPI.COMM_WORLD) 
    rank = comm.Get_rank()
    size = comm.Get_size()
    args = args_parser()
    cool = args.cool 
    fc = args.foldchange
    #from utilities.cool_tools import parser_cool 
    #cool_handle = parser_cool(cool, args.resolution)
    resolution = args.resolution
    ref = args.reference 
    out = args.output
    ref_arm = chrom_arms(ref)
    # clr = parser_cool(cool, resolution)
    # only intra-chromosome (cis) PETs under consideration
    win_list = args.window
    # step = args.step # per bin (resolution) is one step, how many bin to walk each step
    if rank == 0:
        print("Distribute work ...")
        # bin, pixel = raw_pixel(cool, resolution, filter = "cis")
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        for i in ref_arm.index.tolist():
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
        print("Contact resolution: ", resolution)
    else:
        worker_tasks = None
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    for win in win_list:
        pixel_list = list()
        for i in tqdm(worker_tasks[rank], desc = f"Parallel: {rank}"): # each core process one chrom-arm 
            chr,start,end = ref_arm.loc[i,"chrom"], ref_arm.loc[i,"start"], ref_arm.loc[i,"end"]
            i_offset = start//resolution
            region_mtx = cool_matrix(cool, resolution, region = (chr, start, end))
            # = ref_row.chrom, ref_row.start, ref_row.end 
            for center in tqdm(range(start//resolution+win-i_offset, end//resolution-win-i_offset, 1), desc = f"Process:{chr}:{start}-{end}", position = 1, leave = False):
                region_submtx = region_mtx[center-win:center, center+1:center+win+1]
                count = region_submtx.mean()
                size = (region_submtx > 0).sum()
                pixel_list.append((chr, (center+i_offset)*resolution, (center+i_offset+1)*resolution, size, count))
        pixel_df = pd.DataFrame(pixel_list, columns = ["chrom", "start", "end", f"size_{win}", f"score_{win}"])
        pixel_all = comm.gather(pixel_df)
        if rank == 0:
            print("Process window: ", f"{win}x{win}")
            print("Collecting all insulation score ...")
            pixel_all_df = pd.concat(pixel_all, axis = 0)
            pixel_all_df = pixel_all_df[pixel_all_df[f"size_{win}"] > win].copy()
            c_new_list = list()
            for row in ref_arm.itertuples():
                c_df = bf.select(pixel_all_df, (row.chrom, row.start, row.end))
                c_df[f"Df_{win}"] = derivative(np.array(c_df[f"score_{win}"]), args.derivative)
                c_new_list.append(c_df)
            c_new = pd.concat(c_new_list, axis = 0)
            c_new2 = list()
            for row in ref_arm.itertuples():
                c_df = bf.select(c_new, (row.chrom, row.start, row.end))
                c_df[f"valley_{win}"],c_df[f"upmax_{win}"], c_df[f"downmax_{win}"] = valley_point(np.array(c_df[f"Df_{win}"]), np.array(c_df[f"score_{win}"]))
                c_new2.append(c_df)
            c_new2_df = pd.concat(c_new2, axis = 0)
            c_new2_df[f"upfc_{win}"] = c_new2_df[f"upmax_{win}"]/c_new2_df[f"score_{win}"]; c_new2_df[f"downfc_{win}"] = c_new2_df[f"downmax_{win}"]/c_new2_df[f"score_{win}"]
            #c_new2_df.to_csv(f"{out}_{win}.tmp", sep = "\t", header = True, index = False)
            #exit(0)
            # iterative update on maxima value surrounding minima valley 
            c_validate = c_new2_df[c_new2_df[f"valley_{win}"] == 1].copy()
            while len(c_validate[(c_validate[f"upfc_{win}"] < fc) & (c_validate[f"downfc_{win}"] < fc)]) > 0:
                print("Iteration run ...")
                c_new2_df[f"valley_{win}"] = np.where((c_new2_df[f"upfc_{win}"] < fc) & (c_new2_df[f"downfc_{win}"] < fc), 0 , c_new2_df[f"valley_{win}"])
                c_list = list()
                for row in ref_arm.itertuples():
                    c_df = bf.select(c_new2_df, (row.chrom, row.start, row.end))
                    c_df[f"upmax_{win}"], c_df[f"downmax_{win}"] = update_max(np.array(c_df[f"valley_{win}"]), np.array(c_df[f"score_{win}"]))
                    c_list.append(c_df)
                c_new2_df = pd.concat(c_list, axis = 0)
                c_new2_df[f"upfc_{win}"] = c_new2_df[f"upmax_{win}"]/c_new2_df[f"score_{win}"]; c_new2_df[f"downfc_{win}"] = c_new2_df[f"downmax_{win}"]/c_new2_df[f"score_{win}"]
                c_validate = c_new2_df[c_new2_df[f"valley_{win}"] == 1]
            c_new2_df.to_csv(f"{out}_{win}.txt", sep = "\t", header = True, index = False)
    if rank == 0:
        win_file = list()
        for win in win_list:
            c_win = pd.read_table(f"{out}_{win}.txt", sep = "\t", header = 0, index_col = ["chrom", "start", "end"])
            win_file.append(c_win)
        all_win = pd.concat(win_file, axis = 1, join = "inner"); all_win.reset_index(inplace = True)
        all_win.to_csv(f"{out}.txt", sep = "\t", header = True, index = False)
        if os.path.exists(f"{out}.txt"):
            try:
                for win in win_list:
                    os.remove(f"{out}_{win}.txt")
            except:
                pass
        col = [c for c in all_win.columns if c.startswith("valley")]
        select_win = all_win[(all_win[col] == 1).mean(axis = 1) == 1].copy()
        print("Generating boundary bedGraph file (foldchange as )...")
        for win in win_list:
            select_win[f"fc_{win}"] = np.where(select_win[f"upfc_{win}"] > select_win[f"downfc_{win}"], select_win[f"upfc_{win}"], select_win[f"downfc_{win}"])
            select_win[f"strand_{win}"] = np.where((select_win[f"upfc_{win}"] > select_win[f"downfc_{win}"]) & (select_win[f"downfc_{win}"] < fc), "+", np.where((select_win[f"upfc_{win}"] < select_win[f"downfc_{win}"]) & (select_win[f"upfc_{win}"] < fc), "-", "."))
        select_win["fc"] = select_win[[c for c in select_win.columns if c.startswith("fc")]].mean(axis = 1)
        # print(select_win)
        select_win["name"] = range(len(select_win)); select_win["name"] = "b" + select_win["name"].astype(str)
        select_win[["chrom", "start", "end", "name", "fc", f"strand_{win}"]].to_csv(f"{out}.bed", sep = "\t", header = False, index = False)
    exit(0)
        
if __name__ == "__main__":
    main()


