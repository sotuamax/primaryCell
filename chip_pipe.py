#!/usr/bin/env python3
import argparse
import pandas as pd 
import os 
import subprocess

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sample", help = "samples in a xlsx file to process")
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help=True)
    sub_parsers = parser2.add_subparsers(dest = "mode", help = "mode to run")
    process = sub_parsers.add_parser("process", help = "process sample data", parents = [parser], add_help = False)
    process.add_argument("-dir", "--directory", help = "directory where to find the fastq data")
    process.add_argument("-outdir", "--outdir", help = "output directory")
    # process.add_argument("-f", "--force", action = "store_true", help = "force regenerate alignment file even exists")
    process.add_argument("-ref", "--reference", required = False, help = "reference genome to use")
    process.add_argument("-n", "--thread", help = "thread number (default: 4)", default = 4, type = int)
    # process.add_argument("-m", "--m", choices=["SE", "PE"], help = "read alignment mode in SE/PE")
    summary = sub_parsers.add_parser("summary", help = "generate summary report", parents = [parser], add_help = False)
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
    print(args.mode)
    # process samples 
    if args.mode == "summary": 
        import json
        out = args.sample.replace(".xlsx", "_report.xlsx")
        log_dir = "/data/jim4/Seq/CRISPR/log"
        jlist = list()
        for row in sample_df.itertuples():
            logf = os.path.join(log_dir, f"{row.id}.log")
            with open(logf, "r") as logj:
                jdict = json.load(logj)
                jdf = pd.DataFrame.from_dict(jdict, orient = "index")
                jlist.append(jdf)
        sample_report = pd.concat(jlist, axis = 1).transpose()
        # print(sample_report)
        sample_report = pd.merge(sample_df, sample_report, left_on = "id", right_on = "sample").drop("sample", axis = 1)
        with pd.ExcelWriter(out, mode = "w") as writer:
            sample_report.to_excel(writer, index = False)
        
    if args.mode == "process":
        dir = args.directory
        outdir = args.outdir 
        n = args.thread
        # align_mode = args.m
        for row in sample_df.itertuples():
            print(f"Processing {row.id}")
            read = os.path.join(dir, row.id)
            # alignment step (output .bam and .stat)
            if args.reference == None:
                align_command = f"chip_align.py align -read {read}_R1.fastq.gz {read}_R2.fastq.gz -outdir {outdir} -n {n}"
            else:
                print(f"Reference: {args.reference}")
                align_command = f"chip_align.py align -read {read}_R1.fastq.gz {read}_R2.fastq.gz -ref {args.reference} -outdir {outdir} -n {n}"
            # if args.force:
            #     align_command += " -f "
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
