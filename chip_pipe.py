#!/usr/bin/env python3
"""
it execute steps in chip_align.py 

"""
import argparse
import pandas as pd 
import os 
import subprocess

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sample", help = "samples in a xlsx file to process (column: id)")
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help=True)
    sub_parsers = parser2.add_subparsers(dest = "mode", help = "mode to run")
    process = sub_parsers.add_parser("process", help = "process sample data", parents = [parser], add_help = False)
    process.add_argument("-dir", "--directory", help = "directory where to find the fastq data")
    process.add_argument("-outdir", "--outdir", help = "output directory")
    # process.add_argument("-f", "--force", action = "store_true", help = "force regenerate alignment file even exists")
    process.add_argument("-ref", "--reference", required = False, help = "reference genome to use")
    process.add_argument("-n", "--thread", help = "thread number (default: 4)", default = 4, type = int)
    process.add_argument("-m", "--m", choices= ["SE", "PE"], help = "Single-end or Pair-end alignment.", default = "SE")
    summary = sub_parsers.add_parser("summary", help = "generate summary report", parents = [parser], add_help = False)
    summary.add_argument("-outdir", "--outdir", required = False, help = "output directory")
    args = parser2.parse_args()
    return args

def read_df(table):
    # columns: id, sample
    if table.endswith(".xlsx"): 
        sample_df = pd.read_excel(table, header = 0)
        return sample_df 
    else: 
        print("Input sample is not in xlsx format.")
        exit(1)

def main():
    args = args_parser()
    sample_df = read_df(args.sample)
    # process samples 
    if args.mode == "summary": 
        from utilities.parse_log import mark_log_parser, flagstat_tsv_parse 
        import json
        out = args.sample.replace(".xlsx", "_report.xlsx")
        outdir = args.outdir 
        sample_info_list = list()
        for row in sample_df.itertuples():
            stat_file = os.path.join(outdir, f"{row.id}.stat")
            aligned_info = flagstat_tsv_parse(stat_file)
            mark_file = os.path.join(outdir, f"{row.id}.dedup.metric")
            dup_info = mark_log_parser(mark_file)
            all_info = {"sample":row.id} | aligned_info | dup_info
            sample_info_list.append(pd.Series(all_info))
        sample_report = pd.DataFrame(sample_info_list)
        sample_report["aligned_ratio"] = round(sample_report["aligned"].astype(int)/sample_report["total"].astype(int), 2)
        with pd.ExcelWriter(out, mode = "w") as writer:
            sample_report[["sample", "total", "aligned", "aligned_ratio", "dup_rate"]].to_excel(writer, index = False)
    if args.mode == "process":
        dir = args.directory
        outdir = args.outdir 
        n = args.thread
        align_mode = args.m
        for row in sample_df.itertuples():
            print(f"Processing {row.id}")
            read = os.path.join(dir, row.id)
            # alignment step (output .bam and .stat)
            if args.reference == None:
                if align_mode == "SE":
                    align_command = f"chip_align.py align -read {read}_R1.fastq.gz -outdir {outdir} -n {n}"
                if align_mode == "PE":
                    align_command = f"chip_align.py align -read {read}_R1.fastq.gz {read}_R2.fastq.gz -outdir {outdir} -n {n}"
            else:
                print(f"Reference: {args.reference}")
                if align_mode == "SE":
                    align_command = f"chip_align.py align -read {read}_R1.fastq.gz -outdir {outdir} -n {n} -ref {args.reference}"
                if align_mode == "PE":
                    align_command = f"chip_align.py align -read {read}_R1.fastq.gz {read}_R2.fastq.gz -ref {args.reference} -outdir {outdir} -n {n}"
            print(align_command)
            subprocess.call(align_command, shell = True)
            # markdup step 
            align_bam = os.path.join(outdir, row.id + ".bam")
            markdup_command = f"chip_align.py markdup -bam {align_bam} -outdir {outdir}"
            # if args.force:
            #     markdup_command += " -f"
            print(markdup_command)
            subprocess.call(markdup_command, shell = True)
            # markdup 
            mark_bam = os.path.join(outdir, row.id + ".dedup.bam")
            # if args.force:
            #     bw_command += " -f"
            qc_command = f"chip_align.py qc -bam {mark_bam} -n {n} -outdir {outdir}"
            subprocess.call(qc_command, shell = True)
            qc_bam = os.path.join(outdir, row.id + ".dedup.qc.bam")
            bw_command = f"chip_align.py bigwig -bam {qc_bam} -outdir {outdir}"
            print(bw_command)
            subprocess.call(bw_command, shell = True)
            # # log 
            # log_command = f"chip_align.py log -sample {row.id} -mode {align_mode}"
            # print(log_command)
            # subprocess.call(log_command, shell = True)

if __name__ == "__main__":
    main()
