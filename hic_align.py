#!/usr/bin/env python3
"""
To automate the process for Hi-TrAC data processing 

Steps: 
- pre: preprocessing fastq file; 
- align: perform BWA alignment;
- markdup: Markdup using picard; (minimum 20G memory)
- QC: perform QC on bam file (ambiguous alignment removed);

For example: 
================ 
ml fastp 
hic_align.py pre -read <sample> -readdir <directory> -n 4

ml bwa samtools 
hic_align.py align -read <sample> -n 12 

ml picard java 
hic_align.py markdup -bam <sample>.bam

hic_align.py qc -bam <sample>.bam -n 12 
"""
import os
import subprocess
import argparse
import sys
# import logging
# import numpy as np 
import pandas as pd 

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = False)
    #parser.add_argument("-o", "--output", required = True, help = "output name")
    parser.add_argument("-n", "--thread", help = "thread number (default: 4)", default = 4, type = int)
    # subcommand 
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    sub_parsers = parser2.add_subparsers(dest = "command", help = "mode to run")
    # preprocessing step
    pre = sub_parsers.add_parser("pre", help = "preprocessing reads with adpater", parents = [parser])
    pre.add_argument("-read", "--read", help="read name to preprocessing (must be correctly compressed by gzip). ")
    pre.add_argument("-readdir", "--readdir", default = "/data/jim4/Seq/primary_cell_project/fastq/softname/HiTrAC", required = False, help = "read directory")
    pre.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/primary_cell_project/fastq/QC/HiTrAC", required = False, help = "output directory")
    pre.add_argument("-adapter", "--adapter", default = "/data/jim4/Seq/primary_cell_project/data/linker.fa", required = False, help = "adapter fasta file")
    # alignment step
    align = sub_parsers.add_parser("align", help = "align cleaned reads to reference genome", parents = [parser])
    align.add_argument("-ref", "--reference", required = False, default = "/data/jim4/Reference/human/GRCh38.p14/fasta/GRCh38.primary_assembly.genome.fa", help = "reference genome for alignment")
    align.add_argument("-read", "--read", help = "read name for alignment")
    align.add_argument("-readdir", "--readdir", default = "/data/jim4/Seq/primary_cell_project/fastq/QC/HiTrAC", required = False, help = "read directory")
    align.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/raw/individual", required = False, help = "output directory")
    # markdup step
    markdup = sub_parsers.add_parser("markdup", help = "markduplicates of alignment file", parents = [parser])
    markdup.add_argument("-bam", "--bam", help = "alignment bam file")
    markdup.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/markdup/individual", required = False, help = "directory for markdup bam file")
    # qc step 
    qc = sub_parsers.add_parser("qc", help = "quality control of bam file", parents = [parser])
    qc.add_argument("-bam", "--bam", help = "alignment bam file")
    qc.add_argument("-outdir", "--outdir", required = False, default = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/individual", help = "quality control directory for bam file")
    qc.add_argument("-min_quality", "--min_quality", default = 30, help = "minimum quality for aligned read", type = int)
    args = parser2.parse_args()
    return args

def fastp_log_parser(log):
    log_info = list()
    with open(log, "r") as log_open:
        for line in log_open.readlines():
            line = line.strip("\n")
            if line.startswith("total reads"):
                total_read = int(line.split(": ")[-1])
                log_info.append(total_read)
            if line.startswith("reads with adapter trimmed"):
                trimmed_read = int(line.split(": ")[-1])
                log_info.append(trimmed_read)
    return log_info

def mark_log_parser(log):
    i = -2
    with open(log, "r") as log_open: 
        for line in log_open.readlines():
            line = line.strip("\n")
            if i == -1:
                i += 1
            if i == 0:
                pair_read = line.split("\t")[2]
                dup_rate = line.split("\t")[-2]
                break
            if line.startswith("LIBRARY"):
                i += 1
    return pair_read, dup_rate

def main():
    args = args_parser()
    #output = args.output
    # if args.command == "pre":
    #     logging.basicConfig(filename = os.path.join(log_dir, output+".log"), filemode = "w", format = '%(asctime)s %(message)s', datefmt = "%H:%M:%S", level = logging.DEBUG)
    # else:
    #     logging.basicConfig(filename = os.path.join(log_dir, output+".log"), filemode = "a", format = '%(asctime)s %(message)s', datefmt = "%H:%M:%S", level = logging.DEBUG)
    if args.command == "pre": # log file written by fastp 
        print("Preprocessing fastq file ...")
        read = args.read
        readdir = args.readdir
        r1 = os.path.join(readdir, read + "_R1.fastq.gz")
        r2 = os.path.join(readdir, read + "_R2.fastq.gz")
        if os.path.exists(r1) and os.path.exists(r2):
            print("Locate fastq files ... ")
        else: 
            print(f"Cannot find {r1} and {r2}.")
            exit(1)
        outdir = args.outdir
        linker = args.adapter
        n = args.thread
        name = os.path.join(outdir, os.path.basename(read))
        r1_qc = name + "_R1.fastq.gz"
        r2_qc = name + "_R2.fastq.gz"
        print(f"fastp run w/ {n} threads")
        command = f"fastp -i {r1} -I {r2} -o {r1_qc} -O {r2_qc} --detect_adapter_for_pe --adapter_fasta {linker} --thread {n} --json {name}.json --html {name}.html -5 -W 1 &> {name}.log"
        print(command)
        subprocess.call(command, shell = True)
        print("Preprocessing finished!")
        # log_info = fastp_log_parser(name+".log")
        # logging.info(f"Input read\t{log_info[0]}")
        # logging.info(f"QC read\t{log_info[2]}/({round(log_info[2]/log_info[0]*100, 2)}%)")
        # logging.info(f"Adapter read\t{log_info[-1]}/({round(log_info[-1]/log_info[0]*100, 2)}%)")
    if args.command == "align": # no log file 
        print("Perform alignment ...")
        read = args.read
        readdir = args.readdir
        r1 = os.path.join(readdir, read + "_R1.fastq.gz")
        r2 = os.path.join(readdir, read + "_R2.fastq.gz")
        if os.path.exists(r1) and os.path.exists(r2):
            print("Locate fastq files ... ")
        else: 
            print(f"Cannot find {r1} and {r2}.")
            exit(1)
        ref = args.reference 
        n = args.thread
        align_dir = args.outdir
        name = os.path.join(align_dir, os.path.basename(read))
        if not os.path.exists(name + ".bam"):
            print("Start alignment ....")
            print(f"BWA run w/ {n} threads")
            align_command = f"bwa mem -5SP {ref} {r1} {r2} -t {n} | samtools view -@ {n} -Su - | samtools sort -@ {n} - -o {name}.bam && samtools index -@ {n} {name}.bam"
            print(align_command)
            subprocess.call(align_command, shell = True)
        else: 
            print(f"File {name}.bam already existed!")
            exit(1)
        print("Alignment finished!")
    if args.command == "markdup": # log file written by picard ()
        print("Markdup on BAM ...")
        bam = args.bam 
        markdir = args.outdir 
        n = args.thread
        name = os.path.join(os.path.join(markdir, os.path.basename(bam)))
        if not os.path.exists(name):
            picard="/usr/local/apps/picard/3.1.0/picard.jar"
            metric = name.replace(".bam", ".txt")
            markdup_command = f"java -Xmx20g -jar {picard} MarkDuplicates -I {bam} -O {name} -M {metric} --REMOVE_DUPLICATES true && samtools index -@ {n} {name}"
            print(markdup_command)
            subprocess.call(markdup_command, shell = True)
        else:
            print("BAM file already exist! ")
            exit(1)
        print("Markdup finished!")
        # paired_read, dup_rate = mark_log_parser(name.replace(".bam", ".txt"))
        # logging.info(f"Aligned read\t{paired_read}")
        # logging.info("#### Mark duplicates .....")
        # logging.info(f"Duplicate rate\t{dup_rate}")
    if args.command == "qc": # no log file 
        print("Perform QC on BAM file ...")
        bam = args.bam 
        n = args.thread 
        qcdir = args.outdir 
        qc_cutoff = args.min_quality
        name = os.path.join(qcdir, os.path.basename(bam))
        chrsize = pd.read_table("/data/jim4/Reference/human/GRCh38.p14/GRCh38_chrsize.bed", sep = "\t", header = None, names = ["chrom", "size"])
        if not os.path.exists(name):
            print(f"QC w/ {n} threads")
            chroms = " ".join(chrsize["chrom"].tolist())
            qc_command = f"samtools view -@ {n} -h -q {qc_cutoff} -F 2316 -b -o {name} {bam} {chroms} && samtools index -@ {n} {name}"
            print(qc_command)
            subprocess.call(qc_command, shell = True)
        else:
            print("BAM file already exist! ")
            exit(1)
        print("QC finished!")


if __name__ == "__main__":
    main()

