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
derivative: 3 bins on each side -> as an argument (default : 3)
minima: derivative 0 and also the score is minima in surrounding 5 bin on each side -> as an argument (default : 10)
minma screen should examine regions beyond the shifting site (-1, +1)

"""

# import cooler 
# import bioframe as bf 
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
    parser.add_argument("-v", "--valley_screen", type = int, default = 10, help = "To find minima of a screening size (per side)")
    parser.add_argument("-fd", "--foldchange", default = 1.8, type = float, help = "within the screen range, the maximum value is *fd* times of the minima score. ")
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

def valley_point(arr_slope, arr_value, v, fd):
    """
    pick for two consecutive derivatives crossing zero (which can be peak/valley); 
    trend from negative derivative to positive derivative as the region surrounding local minima (valley);
    """
    # 0s as initiate array 
    valley_arr = np.zeros(len(arr_slope))
    #strength_arr = np.zeros(len(arr_slope))
    #minima_arr = np.zeros(len(arr_slope))
    signs = np.sign(arr_slope) # get the sign of slope
    # find shifting signs for neighbor bins, return index
    zero_crossing = np.where(signs[:-1] * signs[1:] < 0)[0] 
    # find shifting derivatives (from negative to positive) where match to a valley, and return index
    valley_zero = zero_crossing[arr_slope[zero_crossing] < arr_slope[zero_crossing+1]] 
    # find neighboring points surrounding derivative at 0, and then returning the one point closest to 0 
    valley_update = np.array([i - 2 + np.argmin(arr_value[i-2:i+4]) for i in valley_zero]) # key improve to expand region for searching of minima around the derivative shifting site 
    valley_update2 = np.array([i for i in valley_update if i > v and np.max(arr_value[i-v:i+v+1]) >= fd*(arr_value[i])])
    valley_arr[valley_update2] = 1 # replace 1 at index for valley
    # strength_arr[valley_update2] = arr_slope[valley_zero+1] - arr_slope[valley_zero] # add strength value at index for valley (slope change)
    return valley_arr #, strength_arr

def main():
    from mpi4py import MPI
    from mpi4py.util import pkl5
    # to overcome the overflow error when comm data > 2 GB
    comm = pkl5.Intracomm(MPI.COMM_WORLD) 
    rank = comm.Get_rank()
    size = comm.Get_size()
    args = args_parser()
    cool = args.cool 
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
        for i in worker_tasks[rank]: # each core process one chrom-arm 
            chr,start,end = ref_arm.loc[i,"chrom"], ref_arm.loc[i,"start"], ref_arm.loc[i,"end"]
            i_offset = start//resolution
            region_mtx = cool_matrix(cool, resolution, region = (chr, start, end))
            # = ref_row.chrom, ref_row.start, ref_row.end 
            for center in tqdm(range(start//resolution+win-i_offset, end//resolution-win-i_offset, 1), desc = f"Process:{chr}:{start}-{end}", position = 1, leave = True):
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
            # for w in win_list:
            #     pixel_all_df = pixel_all_df[pixel_all_df[f"size_{w}"] > w*w*0.3]
            # pixel_all_df.sort_values(by = ["chrom", "start"]).to_csv(f"{out}.txt", sep = "\t", header = True, index = False)
            c_new_list = list()
            for _, c_df in pixel_all_df.groupby("chrom"):
                c_df[f"slope_{win}"] = derivative(np.array(c_df[f"score_{win}"]), args.derivative)
                c_new_list.append(c_df)
            c_new = pd.concat(c_new_list, axis = 0)
            c_new2 = list()
            for _, c_df in c_new.groupby("chrom"):#,c_df[f"strength_{win}"]
                c_df[f"valley_{win}"] = valley_point(np.array(c_df[f"slope_{win}"]), np.array(c_df[f"score_{win}"]), args.valley_screen, args.foldchange)
                c_new2.append(c_df)
            c_new2_df = pd.concat(c_new2, axis = 0)
            c_new2_df.to_csv(f"{out}_{win}.txt", sep = "\t", header = True, index = False)
    if rank == 0:
        win_file = list()
        for win in win_list:
            c_win = pd.read_table(f"{out}_{win}.txt", sep = "\t", header = 0, index_col = ["chrom", "start", "end"])
            #c_valley = c_win[c_win[f"valley_{win}"] == 1]
            #strength_q5 = np.quantile(np.log10(c_valley[f"strength_{win}"]), 0.05)
            #c_valley[np.log10(c_valley[f"strength_{win}"]) > strength_q5][[f"strength_{win}"]].to_csv(f"{win}.bedGraph", header = False, index = True, sep = "\t")
            win_file.append(c_win)
        all_win = pd.concat(win_file, axis = 1, join = "inner")
        all_win.to_csv(f"{out}.txt", sep = "\t", header = True, index = True)
        if os.path.exists(f"{out}.txt"):
            try:
                for win in win_list:
                    os.remove(f"{out}_{win}.txt")
            except:
                pass
    exit(0)
        # valley_df = c_new2_df.query("valley == 1").copy()
        # valley_df["pos1"] = valley_df["start"] + (resolution*2)/valley_df["strength"] * abs(valley_df["slope"])
        # valley_df["pos2"] = valley_df["pos1"] + resolution
        # valley_df["pos1"] = valley_df["pos1"].astype(int); valley_df["pos2"] = valley_df["pos2"].astype(int)
        # valley_df.drop("valley", axis = 1, inplace = True)
        # valley_df.to_csv(f"{out}.txt", sep = "\t", header = True, index = False)
        # valley_df[["chrom", "pos1", "pos2", "strength"]].to_csv(f"{out}.bedGraph", sep = "\t", header = False, index = False)
        
if __name__ == "__main__":
    main()


