#!/usr/bin/env python
import pandas as pd 
import pysam 
import argparse 
import subprocess

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("bam", help = "input bam file (coordinate-sorted and indexed)")
    parser.add_argument("-region", type = str, help = "region in format of chr1:12345-13897")
    parser.add_argument("-dist", "--distance", type = int, required = False, help = "minimum distance to filter intrachromosomal PET")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    bam = args.bam 
    # bam = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/merge2/nCAEC.bam"
    bam_handle = pysam.AlignmentFile(bam, "rb")
    
    # name = bam.replace(".bam", "_cis.bam")
    # 
    chr,chr_region = args.region.split(":")
    start,end = chr_region.split("-")
    out = bam.replace(".bam", "_" + "_".join([chr, start, end]) + ".bam")
    # 100 bp as buffer region
    start,end = int(start)-100, int(end)+100
    with pysam.AlignmentFile(out, "wb", template = bam_handle) as bam_write:
        for read in bam_handle.fetch(chr, start, end):
            rseq = read.reference_name 
            mseq = read.next_reference_name 
            if rseq == mseq:
                rstart = read.reference_start 
                mstart = read.next_reference_start
                if rstart > start and rstart < end and mstart > start and mstart < end:
                    if args.distance == None:
                        bam_write.write(read)
                    else:
                        if abs(rstart - mstart) >= args.distance:
                            bam_write.write(read)
    # techniquely, the bam is still coordinate-sorted
    subprocess.call("ml samtools", shell = True)
    subprocess.call(f"samtools index -@ 2 {out}", shell = True)

if __name__ == "__main__":
    main()
