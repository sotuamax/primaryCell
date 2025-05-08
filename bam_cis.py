#!/usr/bin/env python
import pandas as pd 
import pysam 
import argparse 
import subprocess
import os 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("bam", help = "input bam file (coordinate-sorted and indexed)")
    parser.add_argument("-write", "--write_new", help = "write into new bam file with cis PETs")
    parser.add_argument("-min_quality", "--min_quality", help = "minimum quality of aligned read (default: 0)", default = 0, type = int)
    parser.add_argument("-dist", "--distance", type = int, required = False, help = "minimum distance to filter intrachromosomal PET")
    parser.add_argument("-less_than", action = "store_true", help = "cis PET with distance < cutoff. If not assigned, PET with distance > cutoff")
    parser.add_argument("-t", "--threads", default = 4, type = int, help = "threads to processing bam file")
    parser.add_argument("-o", "--output", required = False, help = "output name")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    bam = args.bam 
    sample = os.path.basename(bam).replace(".bam", "")
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = args.threads)
    # 
    if args.write_new:
        name = args.output
        with pysam.AlignmentFile(name, "wb", template = bam_handle, threads = args.threads) as bam_write:
            if args.distance == None:
                for read in bam_handle.fetch():
                    if read.reference_name == read.next_reference_name:
                        bam_write.write(read)
            else:
                if args.less_than:
                    for read in bam_handle.fetch():
                        if read.reference_name == read.next_reference_name and abs(read.reference_start - read.next_reference_start) < args.distance:
                            bam_write.write(read)
                else:
                    for read in bam_handle.fetch():
                        if read.reference_name == read.next_reference_name and abs(read.reference_start - read.next_reference_start) >= args.distance:
                            bam_write.write(read)
        # techniquely, the bam is still coordinate-sorted
        # subprocess.run(["ml", "samtools"], check = True)
        subprocess.call(f"samtools index -@ {args.threads} {name}", shell = True)
    else:
        total_read = bam_handle.mapped//2
        cis_PETs = {read.query_name for read in bam_handle if read.reference_name == read.next_reference_name}
        bam_handle = pysam.AlignmentFile(bam, "rb", threads = args.threads)
        dist_PETs = {read.query_name for read in bam_handle if read.reference_name == read.next_reference_name and abs(read.reference_start - read.next_reference_start) > 1000}
        cis_num = len(cis_PETs)
        dist_num = len(dist_PETs)
        if not os.path.exists("Cis_ratio.txt"):
            with open("Cis_ratio.txt", "w") as cis_file:
                cis_file.write("sample\ttotal\tcis_PET\tdist_PET\tcis/total\tdist/cis\n")
                cis_file.write(f"{sample}\t{total_read}\t{cis_num}\t{dist_num}\t{round(cis_num/total_read*100, 2)}\t{round(dist_num/cis_num*100, 2)}\n")
        else:
            with open("Cis_ratio.txt", "a") as cis_file:
                cis_file.write(f"{sample}\t{total_read}\t{cis_num}\t{dist_num}\t{round(cis_num/total_read*100, 2)}\t{round(dist_num/cis_num*100, 2)}\n")

if __name__ == "__main__":
    main()

