#!/usr/bin/env python3
"""
Slide a triangle-shape window along the genome (one corner on the diagonal of the matrix), and score for the sum of contacts within the window for each position. At certain location, the score is significantly lower than its surrounding region (reflecting lowered contact frequencies between upstream and downstream loci), this position is referred to as boundary position.


"""

# import cooler 
# import bioframe as bf 
# import numpy as np
import pandas as pd
import argparse
from utilities.cool_tools import raw_pixel
from utilities.chrom_func import chrom_arms
from utilities.misc import ignore_warning
from tqdm.auto import tqdm

ignore_warning()

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = True)
    parser.add_argument("cool", help = "cool input file")
    parser.add_argument("-ref", "--reference", help = "reference genome", default = "hg38")
    parser.add_argument("-r", "--resolution", help = "resolution", type = int)
    parser.add_argument("-step", "--step", help = "step size", type = int)
    parser.add_argument("-w", "--window", help = "window size", type = int)
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def main():
    from mpi4py import MPI
    from mpi4py.util import pkl5
    # to overcome the overflow error when comm data > 2 GB
    comm = pkl5.Intracomm(MPI.COMM_WORLD) 
    rank = comm.Get_rank()
    size = comm.Get_size()
    # 
    args = args_parser()
    cool = args.cool 
    resolution = args.resolution 
    ref = args.reference 
    out = args.output 
    ref_arm = chrom_arms(ref)
    # clr = parser_cool(cool, resolution)
    # only intra-chromosome (cis) PETs under consideration
    win = args.window
    step = args.step # per bin (resolution) is one step, how many bin to walk each step
    pixel_list = list()
    if rank == 0:
        print("Distribute work ...")
        bin, pixel = raw_pixel(cool, resolution, filter = "cis")
        worker_tasks = {w:[] for w in range(size)}
        w_idx = 0
        for i in ref_arm.index.tolist():
            worker_tasks[w_idx].append(i)
            w_idx = (w_idx + 1) % size
        # print("Splited work!")
    else:
        worker_tasks = None
        bin = None
        pixel = None 
    worker_tasks = comm.bcast(worker_tasks, root = 0)
    bin = comm.bcast(bin, root = 0)
    pixel = comm.bcast(pixel, root = 0)
    # 
    # tqdm(worker_tasks[rank][arm_region], desc = f"Obs:{arm_region}", position = 1,leave=False)]
    for i in worker_tasks[rank]:
        ref_row = ref_arm.loc[i,:]
        chr,start,end= ref_row.chrom, ref_row.start, ref_row.end 
        for site in tqdm(range(start//resolution+win, end//resolution-win, step), desc = f"Process:{chr}:{start}-{end}", position = 1, leave = False):
            i = bin.query("chrom == @chr and start == @site*@resolution").index
            site_pixel = pixel.query("@i >= bin1_id >= @i-@win and @i <= bin2_id <= @i+@win and bin2_id - bin1_id <= @win and bin1_id != bin2_id")["count"]
            pixel_list.append((chr, site*resolution, site_pixel.size, site_pixel.sum()))
    pixel_df = pd.DataFrame(pixel_list, columns = ["chrom", "site", "size", "count"])
    #pixel_df["end"] = pixel_df["start"] + resolution
    pixel_df["value"] = round(pixel_df["count"]/pixel_df["size"], 1)
    # pixel_df.to_csv("../../test2.bedGraph", sep = "\t", header = False, index = False)
    pixel_all = comm.gather(pixel_df)
    if rank == 0: 
        print("Collection reso")
        pixel_all_df = pd.concat(pixel_all, axis = 0)
        pixel_all_df.sort_values(by = ["chrom", "site"]).to_csv(f"{out}.txt", sep = "\t", header = True, index = False)

if __name__ == "__main__":
    main()



# cool_insulation_value = insulation(cool_handle, view_df = hg38_arms, clr_weight_name = None, window_bp = diamond_window, nproc = 4, ignore_diags = 3)
