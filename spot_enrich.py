#!/usr/bin/env python3
import os 
import subprocess
import argparse 
import pysam 
import pandas as pd 
import bioframe as bf 
import sys
import numpy as np 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-bed", "--bed", help="bed file for peaks")
    parser.add_argument("-bam", "--bam", help = "bam file used to call peaks (position sorted and indexed)")
    parser.add_argument("-sample", "--sample", default = 5000000, type = int, help = "number of reads to be sampled (default: 5000000)")
    parser.add_argument("-time", "--time", default = 1, type = int, help = "sample times (default: 1)")
    parser.add_argument("-n", "--threads", default = 1, type = int, help = "number of thread for samtools")
    parser.add_argument("-o", "--outdir", default = ".", help = "output directory")
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
    # 
    bed_df = bf.read_table(bed, schema = "bed3")
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    read_num = bam_handle.mapped//2
    score_list = list()
    bam_name = os.path.basename(bam).split(".bam")[0]
    try:
        os.mkdir(outdir)
    except:
        pass 
    for i in range(0, time, 1):
        if read_num > sample:
            sub_bam = f"{os.path.join(outdir, os.path.basename(bam))}.sub{i}"
            sub_command = f"samtools view -@ {n} -h -b --subsample-seed {i+1} --subsample {sample/read_num} {bam} | samtools sort -@ {n} - -o {sub_bam} && samtools index -@ {n} {sub_bam}"
            if not os.path.exists(sub_bam):
                subprocess.call(sub_command, shell = True)
            open_file = sub_bam
            sample_num = sample
        else:
            open_file = bam
            sample_num = read_num
        bam_handle_sample = pysam.AlignmentFile(open_file, "rb", threads = n)
        read_set = [read.query_name for row in bed_df.itertuples() for read in bam_handle_sample.fetch(str(row.chrom), row.start, row.end)]
        # count read query_name occurance 
        read, read_freq = np.unique(read_set, return_counts = True)
        num_read_sample = len(read[read_freq == 2]) # paired reads reside within enriched loci
        spot_score = round(num_read_sample/sample_num, 4)
        score_list.append((bam_name, sample_num, read_num, spot_score))
    df = pd.DataFrame(score_list, columns = ["sample", "read", "total", "SPOT"])
    if os.path.exists(os.path.join(outdir, "SPOT.txt")):
        df.to_csv(os.path.join(outdir, "SPOT.txt"), sep = "\t", index = False, header = False, mode = "a")
    else:
        df.to_csv(os.path.join(outdir, "SPOT.txt"), sep = "\t", index = False, header = True, mode = "w")

if __name__ == "__main__":
    main()