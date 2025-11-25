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
    parser.add_argument("-expand", "--expand", required = False, type = int, help = "expand peak region in bed")
    parser.add_argument("-SE", "--SE", action = "store_true", help = "BAM file is in single end")
    parser.add_argument("-bam", "--bam", help = "bam file used to call peaks (position sorted and indexed)")
    parser.add_argument("-f", "-force", action = "store_true", help = "force regenerate subsampled BAM file")
    parser.add_argument("-barcode", "--barcode", action = "store_true", help = "barcode mode (to count for reads labeled by barcode overlap bed region)")
    parser.add_argument("-sample", "--sample", default = 5000000, type = int, help = "number of reads to be sampled (default: 5000000)")
    # parser.add_argument("-time", "--time", default = 1, type = int, help = "sample times (default: 1)")
    parser.add_argument("-n", "--threads", default = 1, type = int, help = "number of thread for samtools")
    parser.add_argument("-o", "--output", default = ".", help = "output name")
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
    sample = args.sample
    # time = args.time 
    output = args.output
    bed_df = bf.read_table(bed, schema = "bed3")
    if args.expand is not None:
        print(f"Expand bed region by {args.expand} on each side ..")
        bed_df = bf.expand(bed_df, pad = args.expand)
    # for bulk cell level 
    if not args.barcode:
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
        outdir = "."
        try:
            os.mkdir(outdir)
        except:
            pass 
        if read_num > sample:
            sub_bam = f"{os.path.join(outdir, os.path.basename(bam))}.sub"
            sub_command = f"samtools view -@ {n} -h -b --subsample {sample/read_num} {bam} | samtools sort -@ {n} - -o {sub_bam} && samtools index -@ {n} {sub_bam}"
            if not os.path.exists(sub_bam):
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
        if args.SE:
            num_read_sample = a
        else:
            num_read_sample = a//2
        # count read query_name occurance 
        # read, read_freq = np.unique(read_set, return_counts = True)
        # num_read_sample = len(read) # paired reads reside within enriched loci
        spot_score = round(num_read_sample/sample_num, 4)
        score_list.append((os.path.basename(bed.replace(".bed", "")), bam_name, sample_num, read_num, num_read_sample, spot_score))
        df = pd.DataFrame(score_list, columns = ["bed", "sample", "read", "total", "overlap_read", "SPOT"])
        if os.path.exists(os.path.join(outdir, "SPOT.txt")):
            df.to_csv(os.path.join(outdir, "SPOT.txt"), sep = "\t", index = False, header = False, mode = "a")
        else:
            df.to_csv(os.path.join(outdir, "SPOT.txt"), sep = "\t", index = False, header = True, mode = "w")
    # for single cell level (parse each barcode for overlap and total)
    if args.barcode:
        bam_handle_sample = pysam.AlignmentFile(bam, "rb", threads = n)
        barcode_total_dict = dict()
        barcode_overlap_bed_dict = dict()
        for row in bed_df.itertuples():
            for read in bam_handle_sample.fetch(str(row.chrom), row.start, row.end):
                try:
                    barcode_overlap_bed_dict[read.query_name.split(":")[-1]] += 1
                except:
                    barcode_overlap_bed_dict[read.query_name.split(":")[-1]] = 1
        for read in bam_handle_sample.fetch():
            try:
                barcode_total_dict[read.query_name.split(":")[-1]] += 1
            except: 
                barcode_total_dict[read.query_name.split(":")[-1]] = 1 
        for b in barcode_total_dict:
            if b not in barcode_overlap_bed_dict:
                barcode_overlap_bed_dict[b] = 0 
        barcode_total_df = pd.DataFrame.from_dict(barcode_total_dict, orient = "index").rename(columns = {0:"total"})
        barcode_overlap_bed_df = pd.DataFrame.from_dict(barcode_overlap_bed_dict, orient = "index").rename(columns = {0:"overlap"})
        barcode_enrichment = pd.merge(barcode_total_df, barcode_overlap_bed_df, left_index = True, right_index = True)
        barcode_enrichment['overlap_perc'] = barcode_enrichment['overlap']/barcode_enrichment["total"]*100
        if output is not None:
            barcode_enrichment.to_csv(output+".txt", sep = "\t", header = True, index = True)
        else: 
            barcode_enrichment.to_csv("enrich.txt", sep = "\t", header = True, index = True)

if __name__ == "__main__":
    main()