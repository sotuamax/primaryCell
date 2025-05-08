#!/usr/bin/env python3
"""
This is a pipeline developed to utilize dchic (https://github.com/ay-lab/dcHiC) to identify (sub)compartments. 
See dchic publication: https://www.nature.com/articles/s41467-022-34626-6 
for dchic github: https://github.com/ay-lab/dcHiC 

Input: 
    required: 
    - cool file (for chromatin interaction matrix);
    - output name
    optional: 
    - resolution: assign its value when cool contains multiple resolutions
    - backlist: blacklist site
    -ref_path: golden reference file path

Tested on 05/06/2025: 
    `dchic_pipe.py -cool /data/jim4/Seq/primary_cell_project/alignment/HiTrAC/pairs/nCAEC.mcool -r 100000 -ref_path data/hg38_10000_goldenpathData/ -o test`
    
"""

import os 
import pandas as pd
import subprocess 
import argparse 
from utilities.cool_tools import parser_cool 
import glob 
import time 
from utilities.misc import timeit

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description= "dchic_pipe.py -cool <input.mcool> -blacklist <blacklist.bed> -ref hg38 -r 5000 -o <output> \n " \
    "")
    parser.add_argument("-cool", "--cool", help = "cool input file")
    parser.add_argument("-r", "--resolution", type = int, help = "resolution to use", required = False)
    parser.add_argument("-blacklist", "--blacklist", help = "blacklist in bed format", required = False)
    parser.add_argument("-ref", "--reference", required = False, default = "hg38", help = "reference genome")
    parser.add_argument("-ref_path", "--ref_path", required = False, help = "the folder path for reference genome")
    parser.add_argument("-subcomp", "--subcompartment", action = "store_true", help = "when assigned, perform subcompartment run")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def dchic_command(file, ref, ref_path, subcomp = False):
    """"
    parameters of dchicf.r
    dchicf.r 
    --file info.file 
    --pcatype cis (perform PCA on cis interaction matrix)
    --pcatype trans (perform PCA on trans interaction matrix) not suggest 
    --pcatype select (select best PC for downstream analysis, after cis/trans step)
    --pcatype analyze (perform differential PCA analysis on the selected PC's)
    --pcatype subcomp (perform subcompartment test)
    --dirovwt force overwrite (Overwrite the existing directories and files)
    --sthread (parallel sample processing)
    --cthread (parallel chromosome processing)
    --pthread (parallel pca per chromosome per sample)
    """
    command = list()
    # step1: perform PCA analysis on input data 
    # by default PCA run on cis interactions
    step1 = f"dchicf.r --file {file} --pcatype cis --pc 2 --dirovwt T --cthread 2 --pthread 4"
    command.append(step1)
    # step2: select PC for compartment assignment
    if ref_path is None:
        step2 = f"dchicf.r --file {file} --pcatype select --genome {ref}"
    else: 
        # the gfolder should contain 3 files: genome.fa genome.tss.bed genome.chrom.sizes 
        step2 = f"dchicf.r --file {file} --pcatype select --genome {ref} --gfolder {ref_path}"
    command.append(step2)
    # step3 (optional): when subcompartment call required (at least two samples are required)
    if subcomp: 
        # not suggest
        step3 = f"dchicf.r --file {file} --pcatype subcomp "
        command.append(step3)
    # combine all steps (proceed to next step only if previous step succeed)
    return command

def main():
    start = time.time()
    args = args_parser()
    cool = args.cool
    blacklist = args.blacklist
    resolution = args.resolution
    out = args.output
    # get cool handle, and from it to get the bins and its pixels per bin interaction
    print("Parse cool input ...")
    cool_input = parser_cool(cool, resolution)
    # get bin position across the genome 
    bin_bed = cool_input.bins()[:][["chrom", "start", "end"]]
    bin_bed["index"] = bin_bed.index
    bin_bed = bin_bed[~bin_bed["chrom"].isin(["chrY", "chrM"])]
    # get pixel value (bin-bin contacts) 
    pixel_mat = cool_input.pixels()[:] # pixel_mat format: <indexA> <indexB> <count>
    # index is the index in bin_bed file 
    # 
    # remove blacklist region when blacklist is not 
    if blacklist != None:
        import bioframe as bf 
        blacklist_df = bf.read_table(blacklist, schema = "bed3")
        bed_df = bf.subtract(bed_df, blacklist_df) # double check subtract function 
    # write out bed (bin) and mat (pixel) for dchic run
    print("Write bin and pixel data into file ...")
    bin_bed.to_csv(f"{out}_{resolution}.bed", sep = "\t", header = False, index = False)
    pixel_mat.to_csv(f"{out}_{resolution}.mat", sep = "\t", header = False, index = False)
    # write a file annotating bin_bed and pixel_mat 
    with open(out + ".file", "w") as fw:
        # mat bed replicate_prefix experiment_prefix
        # <mat> <bin> <replicate> <sample>
        print("Prepare file for dchic ...")
        fw.write(f"{out}_{resolution}.mat\t{out}_{resolution}.bed\t{out}_{resolution}\t{out}\n")
    # run dchic command 
    command2run = dchic_command(out + ".file", args.reference, args.ref_path, args.subcompartment)
    try: 
        print("Run dchic ...")
        total_steps = len(command2run)
        for i in range(total_steps):
            print(f"Perform step {i+1} ...")
            print(command2run[i])
            subprocess.call(f"{command2run[i]} >> {out}.log 2>&1", shell = True)
    except Exception as e: 
        print(e)
        exit(1)
    # reorganize dchic output once succeed 
    sample_file = pd.read_table(out + ".file", sep = "\t", header = None, names = ["mat", "bed", "replicate", "sample"])
    print("Gather dchic output ... ")
    for row in sample_file.itertuples():
        # note that dchic output name is in the name of the replicate 
        file = row.replicate 
        # the selected pca results 
        # pc_selected = pd.read_table(file + "_chr_pc_selected.txt", sep = "\t", header = 0)
        # pc value stored in directory 
        dir = os.path.join(file + "_pca", "intra_pca", file + "_mat")
        # collect all final pc values per chromosome 
        pc_f = sorted([f for f in glob.glob(os.path.join(dir, "*.pc.bedGraph"))])
        # combine pc values of all chromosomes 
        if len(pc_f) > 0:
            pc_all = pd.concat([pd.read_table(f, sep = "\t", header = None) for f in pc_f], axis = 0)
            pc_all.to_csv(f"{out}_{resolution}.bedGraph", sep = "\t", header = False, index = False)
        else:
            print("pc select bedGraph file cannot be identified.")
    print(timeit(start))

if __name__ == "__main__":
    main()
