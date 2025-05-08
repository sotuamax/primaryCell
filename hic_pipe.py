#!/usr/bin/env python3
import os 
import pandas as pd 
import argparse 
import sys 
import subprocess

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-sample", "--sample_table", help = "samples in a tsv file to process")
    parser.add_argument("-dir", "--directory", help = "directory where to find the fastq data")
    parser.add_argument("-n", "--thread", help = "thread number (default: 4)", default = 4, type = int)
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    sub_parsers = parser2.add_subparsers(dest = "command", help = "mode to run")
    proc = sub_parsers.add_parser("process", help = "process each library for alignment, markdup, and QC", parents = [parser])
    merge = sub_parsers.add_parser("merge", help = "merge libraries belonging to the same sample", parents = [parser])
    fm = sub_parsers.add_parser("transform", help = "generate pairs file for PET", parents = [parser])
    fm.add_argument("-chrsize", "--chrsize", help = "sorted chrsize file")
    args = parser2.parse_args()
    return args

def read_df(table):
    # columns: id, sample
    if table.endswith(".xlsx"): 
        sample_df = pd.read_excel(table, header = 0)
    return sample_df 

def process_sample(args, sample_df, dir):
    for row in sample_df.itertuples():
        print(f"Processing {row.id}")
        read = os.path.join(dir, row.id)
        # preprocessing step 
        pre_step = f"hic_align.py pre -read {read} -n {args.thread} -o {row.id}"
        subprocess.call(pre_step, shell = True)
        qc_dir = "/data/jim4/Seq/primary_cell_project/fastq/QC/HiTrAC"
        qc_read = os.path.join(qc_dir, row.id)
        # alignment to reference genome
        if row.ref == "human":
            ref = "/data/jim4/Reference/human/GRCh38.p14/fasta/GRCh38.primary_assembly.genome.fa"
        align_step = f"hic_align.py align -read {qc_read} -ref {ref} -n {args.thread} -o {row.id}"
        subprocess.call(align_step, shell = True)
        # markdup alignment file 
        align_dir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/raw/individual"
        align_bam = os.path.join(align_dir, row.id + ".bam")
        mark_step = f"hic_align.py markdup -bam {align_bam} -o {row.id}"
        subprocess.call(mark_step, shell = True)
        # qc alignment file 
        mark_dir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/markdup/individual"
        mark_bam = os.path.join(mark_dir, row.id + ".bam")
        qc_step = f"hic_align.py qc -bam {mark_bam} -n {args.thread} -o {row.id}"
        subprocess.call(qc_step, shell = True)
        # cis bam file 
        qcdir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/individual"
        qc_bam = os.path.join(qcdir, row.id + ".bam")
        cis_step = f"hic_align.py cis -bam {qc_bam} -n {args.thread} -o {row.id}"
        subprocess.call(cis_step, shell = True)

def sample_log(sample_df):
    sample_result = list()
    for row in sample_df.itertuples():
        sample_stat = log_parser(row.id + ".log")
        sample_row = list(row.id) + list(sample_stat)
        sample_result.append(sample_row)
    sample_result_df = pd.DataFrame(sample_result, columns = ["id", "raw", "clean", "aligned", "duplication", "aligned_QC", "cis", "shortPET", "intermediatePET", "longPET"])
    return sample_result_df 

def log_parser(log):
    log_dir = "/data/jim4/Seq/primary_cell_project/log"
    with open(os.path.join(log_dir, log), "r") as log_f:
        for line in log_f.readlines():
            line = line.strip("\n")
            if "#" not in line:
                if "Input" in line.split("\t")[0]:
                    input_total = line.split("\t")[-1]
                if "Clean" in line.split("\t")[0]:
                    clean_total = line.split("\t")[-1]
                if "Aligned" in line.split("\t")[0]:
                    aligned_total = line.split("\t")[-1]
                if "Duplication" in line.split("\t")[0]:
                    duplicate = line.split("\t")[-1]
                if "QC" in line.split("\t")[0]:
                    QC = line.split("\t")[-1]
                if "Cis" in line.split("\t")[0]:
                    cis = line.split("\t")[-1]
                if "short distance" in line.split("\t")[0]:
                    short_dist = line.split("\t")[-1]
                if "intermediate distance" in line.split("\t")[0]:
                    inter_dist = line.split("\t")[-1]
                if "long distance" in line.split("\t")[0]:
                    long_dist = line.split("\t")[-1]
    return (input_total, clean_total, aligned_total, duplicate, QC, cis, short_dist, inter_dist, long_dist)

def bam_merge(args, sample_df, dir):
    for group, group_df in sample_df.groupby("sample"):
        group_id = group_df["id"].to_list()
        RG_label = "\n".join([f'@RG\tID:{id}\tSM:{group}\tLB:{id}\tPL:ILLUMINA' for id in group_id])
        command = f"printf '{RG_label}' > {group}.txt"
        subprocess.call(command, shell = True)
        file_combined = " ".join([os.path.join(dir, g + ".bam") for g in group_id])
        merge_command = f'samtools merge -rh {group}.txt -@ {args.thread} -o - {file_combined} | samtools sort -@ {args.thread} - -o {group}.bam && samtools index -@ {args.thread} {group}.bam && rm {group}.txt'
        if not os.path.exists(group + ".bam"):
            print(merge_command)
            subprocess.call(merge_command, shell = True)

def preserve_bam_header(bam):
    bam_handle = pysam.AlignmentFile(bam, "rb")
    bam_header = bam_handle.header()
    return bam_header 

def main():
    args = args_parser()
    sample_df = read_df(args.sample_table)
    dir = args.directory
    if args.command == "process":
        process_sample(args, sample_df, dir)
    if args.command == "merge":
        bam_merge(args, sample_df, dir)
    if args.command == "transform":
        chrsize = args.chrsize
        for sample, sample_df in sample_df.groupby("sample"):
            sample_bam = os.path.join(dir, sample + ".bam")
            if not os.path.exists(sample + ".pairs"):
                pair_generator = f"hic_format.py pairs {sample_bam} -o {sample} -chrsize {chrsize} -n {args.thread}"
                subprocess.call(pair_generator, shell = True)
            if not os.path.exists(sample + ".cool"):
                cool_generator = f"hic_format.py cool {sample}.pairs -o {sample} -chrsize {chrsize} -n {args.thread}"
                subprocess.call(cool_generator, shell = True)
                blacklist= "/data/jim4/data/blacklist/ENCFF356LFX.bed"
                dchic_command = f"hic_format.py dchic {sample}.mcool -o {sample} -chrsize {chrsize} -n {args.thread} -blacklist {blacklist} -r 25000 50000 100000 250000 500000"
                subprocess.call(dchic_command, shell = True)
            if not os.path.exists(sample + ".hic"):
                hic_generator = f"hic_format.py hic {sample}.pairs -o {sample} -chrsize {chrsize} -n {args.thread}"
                subprocess.call(hic_generator, shell = True)

if __name__ == "__main__":
    main()
