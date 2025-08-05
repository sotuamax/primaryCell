#!/usr/bin/env python3
"""
This is a pipeline developed to utilize dchic (https://github.com/ay-lab/dcHiC) to identify (sub)compartments. 
See dchic publication: https://www.nature.com/articles/s41467-022-34626-6 
for dchic github: https://github.com/ay-lab/dcHiC 

Input: 
    required: 
    - cool file (for chromatin interaction matrix);
    - output prefix
    optional: 
    - resolution: assign its value when cool contains multiple resolutions
    - blacklist: blacklist site
    - ref: reference genome version 
    - ref_path: golden reference file path

Tested on 05/06/2025: 
source myconda 
mamba activate dchic 
dchic_pipe.py -cool /data/jim4/Seq/primary_cell_project/alignment/HiTrAC/pairs/nCAEC.mcool -r 10000 -ref_path/data/jim4/Seq/primary_cell_project/analysis/Compartment/data/hg38_goldenpathData/ -o test

    # the goldenpathData for ref is at 10 kb resolution 

"""
import numpy as np 
import os 
from sklearn.preprocessing import StandardScaler
import pandas as pd
import subprocess 
import argparse 
from utilities.cool_tools import parser_cool 
import glob 
import time 
from utilities.misc import timeit

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help=True, 
                                     usage="\nInput as cool: \n  dchic_pipe.py -cool <input.mcool> -r 5000 -blacklist <blacklist.bed> -ref_path <goldref> -o <output>\n\nInput as tsv: \n  dchic_pipe.py -bin bin.txt pixel.txt -blacklist <blacklist.bed> -ref_path <goldref> -o <output> ")
    group = parser.add_mutually_exclusive_group(required = True)
    group.add_argument("-cool", help = "cool input file") # either cool provided or bin provided (as exclusive parameters)
    group.add_argument("-bin_pixel", nargs = 2, help = "chromosome bin and pixel file")
    #parser.add_argument("-pixel", help = "bin pixel file")
    parser.add_argument("-r", "--resolution", type = int, help = "resolution to use", required = False)
    parser.add_argument("-blacklist", "--blacklist", help = "blacklist in bed format", required = False)
    parser.add_argument("-ref", "--reference", required = False, default = "hg38", help = "reference genome")
    parser.add_argument("-ref_path", "--ref_path", required = False, help = "the folder path for reference genome (not related to resolution)")
    # parser.add_argument("-subcomp", "--subcompartment", action = "store_true", help = "when assigned, perform subcompartment run")
    parser.add_argument("-n", "--threads", type = int, default = 1, help = "threads to use ")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def value_scale(value):
    value = np.array(value)
    scaled_value = StandardScaler().fit_transform(value.reshape(-1,1))
    return scaled_value.reshape(-1)

def bin_scale(value):
    b = pd.cut(value, bins = np.quantile(value, np.arange(0,1.1,0.1)), include_lowest=True, retbins = True, labels = list(range(0,10)))
    return b[0]

def df_group(df):
    df["sign"] = np.sign(df["value"])
    df["z"] = df.groupby("chrom")["value"].transform(value_scale)
    df["bin"] = df.groupby("chrom")["z"].transform(bin_scale)
    return df

def dchic_command(file, ref, ref_path, n = 1):
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
    step1 = f"dchicf.r --file {file} --pcatype cis --pc 2 --dirovwt T --cthread {n}"
    command.append(step1)
    # step2: select PC for compartment assignment
    if ref_path is None:
        step2 = f"dchicf.r --file {file} --pcatype select --genome {ref}"
    else: 
        # the gfolder should contain 3 files: genome.fa genome.tss.bed genome.chrom.sizes 
        step2 = f"dchicf.r --file {file} --pcatype select --genome {ref} --gfolder {ref_path}"
    command.append(step2)
    # combine all steps (proceed to next step only if previous step succeed)
    return command

def main():
    start = time.time()
    args = args_parser()
    cool = args.cool
    blacklist = args.blacklist
    resolution = args.resolution
    out = args.output
    bin = args.bin_pixel
    # get cool handle, and from it to get the bins and its pixels per bin interaction
    if cool != None:
        print("Parse cool input ...")
        cool_input = parser_cool(cool, resolution)
        # get bin position across the genome 
        bin_bed = cool_input.bins()[:][["chrom", "start", "end"]]
        bin_bed["index"] = bin_bed.index
        bin_bed = bin_bed[~bin_bed["chrom"].isin(["chrY", "chrM"])]
        # get pixel value (bin-bin contacts) 
        pixel_mat = cool_input.pixels()[:] # pixel_mat format: <indexA> <indexB> <count>
        print(f"Write bin and pixel data with {resolution} into file ...")
        bin_bed.to_csv(f"{out}_{resolution}.bed", sep = "\t", header = False, index = False)
        pixel_mat.to_csv(f"{out}_{resolution}.pixel", sep = "\t", header = False, index = False)
    # index is the index in bin_bed file 
    # remove blacklist region when blacklist is not 
    if blacklist != None:
        import bioframe as bf 
        blacklist_df = bf.read_table(blacklist, schema = "bed3")
        bin_bed = bf.subtract(bin_bed, blacklist_df) # double check subtract function 
    # write out bed (bin) and mat (pixel) for dchic run
    # write a file annotating bin_bed and pixel_mat 
    with open(out + ".file", "w") as fw:
        # mat bed replicate_prefix experiment_prefix
        # <mat> <bin> <replicate> <sample>
        print("Prepare file for dchic ...")
        if cool != None:
            fw.write(f"{out}_{resolution}.pixel\t{out}_{resolution}.bed\t{out}_{resolution}\t{out}\n")
        else:
            fw.write(f"{bin[-1]}\t{bin[0]}\t{out}\t{out}\n")
    # run dchic command 
    command2run = dchic_command(out + ".file", args.reference, args.ref_path, n = args.threads)
    try: 
        print("Run dchic ...")
        total_steps = len(command2run)
        for i in range(total_steps):
            print(f"Perform step {i+1} ...")
            print(command2run[i])
            subprocess.call(f"{command2run[i]}", shell = True)
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
            pc_all.columns = ["chrom", "start", "end", "value"]
            pc_all["z"] = pc_all.groupby("chrom")["value"].transform(value_scale)
            # pc_all["sign"] = np.sign(pc_all["value"])
            # df["bin"] = df.groupby("chrom")["z"].transform(bin_scale)
            if np.equal(np.sign(pc_all["z"]), np.sign(pc_all["value"])).sum() == len(pc_all):
                print("After normalization (z-score), sign was not changed.")
            else:
                print("After normalization (z-score), sign changed. ")
                print("Correcting sign ...")
                pc_all["z"] = np.where(np.equal(np.sign(pc_all["z"]), np.sign(pc_all["value"])), pc_all["z"], -(pc_all["z"]))
                
        else:
            print("pc select bedGraph file cannot be identified.")

    print(timeit(start))

if __name__ == "__main__":
    main()
