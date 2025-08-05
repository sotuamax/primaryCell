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
    parser.add_argument("-SE", "--SE", action = "store_true", help = "BAM file is in single end")
    parser.add_argument("-bam", "--bam", help = "bam file used to call peaks (position sorted and indexed)")
    parser.add_argument("-f", "-force", action = "store_true", help = "force regenerate subsampled BAM file")
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
    flagstat_f = bam.replace(".bam", ".flagstat")
    if os.path.exists(flagstat_f) and os.path.getmtime(flagstat_f) > os.path.getmtime(bam):
        pass 
    else:
        subprocess.call(f"samtools flagstat {bam} -@ {n} > {flagstat_f}", shell = True)
    if args.SE:
        from utilities.parse_log import flagstat_parser_SE
        bam_info = flagstat_parser_SE(flagstat_f)
    else:
        from utilities.parse_log import flagstat_parser
        bam_info = flagstat_parser(flagstat_f)
    read_num = bam_info["total"]
    # bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    # read_num = bam_handle.mapped
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
            if not os.path.exists(sub_bam) or args.force:
                subprocess.call(sub_command, shell = True)
            open_file = sub_bam
            flagstat_f0 = sub_bam.replace(".bam", ".flagstat")
            subprocess.call(f"samtools flagstat {sub_bam} -@ {n} > {flagstat_f0}", shell = True)
            #sample_num = sample
            if args.SE: 
                from utilities.parse_log import flagstat_parser_SE
                bam_info = flagstat_parser_SE(flagstat_f0)
            else:
                from utilities.parse_log import flagstat_parser
                bam_info = flagstat_parser(flagstat_f0)
            sample_num = bam_info["total"]
        else:
            open_file = bam
            sample_num = read_num
        bam_handle_sample = pysam.AlignmentFile(open_file, "rb", threads = n)
        a = 0
        for row in bed_df.itertuples(): 
            for _ in bam_handle_sample.fetch(str(row.chrom), row.start, row.end):
                a += 1 
        # count read query_name occurance 
        # read, read_freq = np.unique(read_set, return_counts = True)
        if args.SE:
            num_read_sample = a
        else:
            num_read_sample = a//2
        # num_read_sample = len(read) # paired reads reside within enriched loci
        spot_score = round(num_read_sample/sample_num, 4)
        score_list.append((os.path.basename(bed.replace(".bed", "")), bam_name, sample_num, read_num, num_read_sample, spot_score))
    df = pd.DataFrame(score_list, columns = ["bed", "sample", "read", "total", "overlap_read", "SPOT"])
    if os.path.exists(os.path.join(outdir, "SPOT.txt")):
        df.to_csv(os.path.join(outdir, "SPOT.txt"), sep = "\t", index = False, header = False, mode = "a")
    else:
        df.to_csv(os.path.join(outdir, "SPOT.txt"), sep = "\t", index = False, header = True, mode = "w")

if __name__ == "__main__":
    main()