#!/usr/bin/env python3
import os 
import subprocess
import argparse 
import pysam 
import pandas as pd 
import bioframe as bf 
import sys
from utility import write_pairs_from_bam,bam2seq 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-bed", "--bed", help="bed file for peaks (read enriched region)")
    parser.add_argument("-bam", "--bam", help = "bam file used to call peaks")
    parser.add_argument("-read_len", "--read_length", type = int, default = 1, help = "read length used to expand its covered region")
    parser.add_argument("-bam_count", "--bam_count", required = False, type = int, help = "read count in the given bam file")
    parser.add_argument("-sample", "--sample", default = 5000000, type = int, help = "number of reads to be sampled (default: 5000000)")
    parser.add_argument("-time", "--time", default = 1, type = int, help = "sample times (default: 1)")
    parser.add_argument("-n", "--threads", default = 1, type = int, help = "number of thread for samtools")
    parser.add_argument("-o", "--outdir", help = "output directory")
    args=parser.parse_args()
    return args

def main():
    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")
    args = args_parser()
    bed = args.bed
    bam = args.bam
    n = args.threads
    # 
    sample = args.sample
    time = args.time 
    outdir = args.outdir
    try:
        os.mkdir(outdir)
    except:
        pass 
    # 
    bed_df = bf.read_table(bed, schema = "bed3")
    bam_handle = pysam.AlignmentFile(bam, "rb")
    chrom_order = bam2seq(bam)
    if args.bam_count == None:
        read_num = sum([bam_handle.count(c) for c in chrom_order])//2
    else:
        read_num = args.bam_count
    score_list = list()
    bam_name = os.path.join(outdir, os.path.basename(bam).split(".bam")[0])
    # 
    subprocess.call("ml samtools", shell = True)
    for i in range(0, time, 1):
        print(f"Run {i+1}/{time}")
        if read_num > sample:
            sub_command = f"samtools view -@ {n} -h -b --subsample-seed {i+1} --subsample {sample/read_num} {bam} -o {bam_name}_sub{i}.bam && samtools index -@ {n} {bam_name}_sub{i}.bam"
            if not os.path.exists(f"{bam}.sub{i}"):
                subprocess.call(sub_command, shell = True)
            open_bam = f"{bam_name}_sub{i}.bam"
            # sample_num = sample
        else:
            open_bam = bam
            # sample_num = read_num
        bam_handle_sample = pysam.AlignmentFile(open_bam, "rb")
        # 
        pair_df = write_pairs_from_bam(open_bam, chrom_order)
        sample_num = len(pair_df)
        # ["readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"]
        r1_df = pair_df[["chr1", "pos1", "readID"]].copy()
        r2_df = pair_df[["chr2", "pos2", "readID"]].copy()
        r1_df["start"] = r1_df["pos1"]-args.read_length; r1_df["end"] = r1_df["pos1"]+args.read_length+1; r1_df = r1_df.rename(columns = {"chr1": "chrom"}).drop("pos1", axis = 1)
        r2_df["start"] = r2_df["pos2"]-args.read_length; r2_df["end"] = r2_df["pos2"]+args.read_length+1; r2_df = r2_df.rename(columns = {"chr2": "chrom"}).drop("pos2", axis = 1)
        # 
        r1_peak = bf.overlap(r1_df, bed_df, how = "inner")
        r2_peak = bf.overlap(r2_df, bed_df, how = "inner")
        # 
        paired_peak = set(r1_peak["readID"]).intersection(set(r2_peak["readID"]))
        # 
        num_read_sample = len(paired_peak)
        spot_score = round(num_read_sample/sample_num, 4)
        score_list.append((bam_name, num_read_sample, sample_num, read_num, spot_score))
    # write score into df
    df = pd.DataFrame(score_list, columns = ["sample", "peak_read_num", "sample_num", "read_num", "spot"])
    df.to_csv(os.path.join(outdir, "spot.txt"), sep = "\t", index = False, header = False, mode = "a")

if __name__ == "__main__":
    main()
    