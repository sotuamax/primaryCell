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
ml fastp bwa samtools deeptools picard java 
hic_align.py pre -read <sample_prefix> -n 4
# use fastp to filter reads on adapters 

hic_align.py align -read <sample> -n 12 
# align PETs onto human genome with "-5SP" 

hic_align.py markdup -bam <sample>.bam
# markduplicates and deduplicates for PETs

hic_align.py qc -bam <sample>.bam -n 12
# QC control on bam file MAPQ > 30 and on nuclear chromosome 

hic_align.py nsort -bam <sample>.bam -n 4 
# sort by names and generate bedpe file to count for cis/trans PETs 

hic_align.py log -id <sample> 
# generate log file for each sample for each step processing

"""
import os
import subprocess
import argparse
import sys
# import logging
# import numpy as np 
import pandas as pd 
import glob 
from utilities.misc import ignore_warning
import json

ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = False)
    #parser.add_argument("-o", "--output", required = True, help = "output name")
    parser.add_argument("-n", "--threads", help = "thread number (default: 4)", default = 4, type = int)
    # subcommand 
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    sub_parsers = parser2.add_subparsers(dest = "command", help = "mode to run")
    # preprocessing step
    pre = sub_parsers.add_parser("pre", help = "preprocessing reads with adpater", parents = [parser])
    pre.add_argument("-read", "--read", help="read name to preprocessing (must be correctly compressed by gzip). ")
    #pre.add_argument("-readdir", "--readdir", default = "/data/jim4/Seq/primary_cell_project/fastq/softname/HiTrAC", required = False, help = "read directory")
    pre.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/primary_cell_project/fastq/QC/HiTrAC", required = False, help = "output directory")
    pre.add_argument("-adapter", "--adapter", default = "/data/jim4/Seq/primary_cell_project/data/linker.fa", required = False, help = "adapter fasta file")
    # alignment step
    align = sub_parsers.add_parser("align", help = "align cleaned reads to reference genome", parents = [parser])
    align.add_argument("-ref", "--reference", required = False, default = "/data/jim4/Reference/human/GRCh38.p14/fasta/GRCh38.primary_assembly.genome.fa", help = "reference genome for alignment")
    align.add_argument("-read", "--read", help = "read name for alignment")
    #align.add_argument("-readdir", "--readdir", default = "/data/jim4/Seq/primary_cell_project/fastq/QC/HiTrAC", required = False, help = "read directory")
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
    # name sort 
    nsort = sub_parsers.add_parser("nsort", help = "sort bam file using name (in alpha-numeric ordering)", parents = [parser])
    nsort.add_argument("-bam", "--bam", help = "alignment bam file")
    nsort.add_argument("-outdir", "--outdir", required = False, default = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/sortn", help = "directory for name-sorted bam file")
    # log file
    log = sub_parsers.add_parser("log", help = "parse log file", parents = [parser])
    log.add_argument("-id", "--id", help = "sample ID used to get log")
    args = parser2.parse_args()
    return args

def main():
    args = args_parser()
    n = args.threads
    #output = args.output
    # if args.command == "pre":
    #     logging.basicConfig(filename = os.path.join(log_dir, output+".log"), filemode = "w", format = '%(asctime)s %(message)s', datefmt = "%H:%M:%S", level = logging.DEBUG)
    # else:
    #     logging.basicConfig(filename = os.path.join(log_dir, output+".log"), filemode = "a", format = '%(asctime)s %(message)s', datefmt = "%H:%M:%S", level = logging.DEBUG)
    if args.command == "pre": # log file written by fastp 
        read = args.read
        #readdir = args.readdir
        r1 = glob.glob(f"{read}_*R1.fastq.gz")
        #os.path.join(readdir, read + "_R1.fastq.gz")
        r2 = glob.glob(f"{read}_*R2.fastq.gz") # os.path.join(readdir, read + "_R2.fastq.gz")
        if len(r1) == 1 and len(r2) == 1:
            #print(r1+r2)
            print("Locate fastq files ... ")
        else: 
            print(f"Cannot find paired read for {read}!")
            exit(1)
        outdir = args.outdir
        linker = args.adapter
        name = os.path.join(outdir, os.path.basename(read))
        r1_qc = name + "_R1.fastq.gz"
        r2_qc = name + "_R2.fastq.gz"
        if os.path.exists(r1_qc) and os.path.exists(r2_qc):
            exit(0)
        else:
            print("Preprocessing fastq file ...")
            print(f"fastp run w/ {n} threads")
            command = f"fastp -i {r1[0]} -I {r2[0]} -o {r1_qc} -O {r2_qc} --detect_adapter_for_pe --adapter_fasta {linker} --thread {n} --json {name}.json --html {name}.html -5 -W 1 &> {name}.log"
            print(command)
            subprocess.call(command, shell = True)
            subprocess.call(f"cat {name}.log", shell = True)
            print("Preprocessing finished!")
        # log_info = fastp_log_parser(name+".log")
        # logging.info(f"Input read\t{log_info[0]}")
        # logging.info(f"QC read\t{log_info[2]}/({round(log_info[2]/log_info[0]*100, 2)}%)")
        # logging.info(f"Adapter read\t{log_info[-1]}/({round(log_info[-1]/log_info[0]*100, 2)}%)")
    if args.command == "align": # no log file 
        read = args.read
        #readdir = args.readdir
        r1 = glob.glob(f"{read}_*R1*.fastq.gz")
        r2 = glob.glob(f"{read}_*R2*.fastq.gz")
        if len(r1) == 1 and len(r2) == 1:
            #print(r1+r2)
            print("Locate fastq files ... ")
        else: 
            print(f"Cannot find paired read for {read}!")
            exit(1)
        ref = args.reference 
        align_dir = args.outdir
        name = os.path.join(align_dir, os.path.basename(read))
        if not os.path.exists(name + ".bam.bai"):
            print("Perform alignment ...")
            print(f"BWA run w/ {n} threads")
            align_command = f"bwa mem -5SP {ref} {r1[0]} {r2[0]} -t {n} | samtools view -@ {n} -Su - | samtools sort -@ {n} - -o {name}.bam && samtools index -@ {n} {name}.bam"
            print(align_command)
            subprocess.call(align_command, shell = True)
            print("Alignment finished!")
        align_stat_out = name + ".stat"
        if not os.path.exists(align_stat_out) and os.path.exists(name + ".bam"):
            print("Generate alignment stat ...")
            align_stat = f"samtools flagstat -@ {n} {name}.bam > {align_stat_out}"
            subprocess.call(align_stat, shell = True)
        if os.path.exists(name + ".bam") and os.path.exists(name + ".bam.bai") and os.path.exists(name + ".stat"):
            exit(0)
        else:
            exit(1)
    if args.command == "markdup": # log file written by picard ()
        bam = args.bam 
        markdir = args.outdir 
        name = os.path.join(os.path.join(markdir, os.path.basename(bam)))
        if not os.path.exists(name):
            print("Markdup on BAM ...")
            picard="/usr/local/apps/picard/3.1.0/picard.jar"
            metric = name.replace(".bam", ".txt")
            # picard markdup requires bam sorted by chromosome coordinates 
            markdup_command = f"java -Xmx20g -jar {picard} MarkDuplicates -I {bam} -O {name} -M {metric} --REMOVE_DUPLICATES true && samtools index -@ {n} {name}"
            print(markdup_command)
            subprocess.call(markdup_command, shell = True)
            print("Markdup finished!")
        if os.path.exists(name) and os.path.exists(name.replace(".bam", ".txt")):
            exit(0)
        else:
            exit(1)
        # paired_read, dup_rate = mark_log_parser(name.replace(".bam", ".txt"))
        # logging.info(f"Aligned read\t{paired_read}")
        # logging.info("#### Mark duplicates .....")
        # logging.info(f"Duplicate rate\t{dup_rate}")
    if args.command == "qc": # no log file 
        # QC: mapping quality, flag 2316, chromosome read
        bam = args.bam 
        qcdir = args.outdir 
        qc_cutoff = args.min_quality
        name = os.path.join(qcdir, os.path.basename(bam))
        chrsize = pd.read_table("/data/jim4/Reference/human/GRCh38.p14/GRCh38_chrsize.bed", sep = "\t", header = None, names = ["chrom", "size"])
        if not os.path.exists(name):
            print("Perform QC on BAM file ...")
            chroms = " ".join(chrsize["chrom"].tolist())
            # -q for mapping-quality 
            # -F flag (2316)
            # given chrom to only keep reads on chroms 
            qc_command = f"samtools view -@ {n} -h -q {qc_cutoff} -F 2316 -b -o {name} {bam} {chroms} && samtools index -@ {n} {name}"
            print(qc_command)
            subprocess.call(qc_command, shell = True)
        bw = name.replace(".bam", ".bw")
        if not os.path.exists(bw):
            bw_command = f"bamCoverage -b {name} -o {bw} -of bigwig -p {n} --normalizeUsing CPM"
            subprocess.call(bw_command, shell=True)
        enrich_log = name.replace(".bam", ".stat")
        if not os.path.exists(enrich_log):
            hg38_TSS = "/data/jim4/Reference/human/GRCh38.p14/GRCh38_TSS.bed"
            # bw = f"/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/individual/{id}.bw"
            enrich_score = bw.replace(".bw", "")
            bw_enrich_command = f"bed_profile.py {hg38_TSS} {bw} -o {enrich_score}"
            subprocess.call(bw_enrich_command, shell=True)
        if os.path.exists(name) and os.path.exists(bw) and os.path.exists(enrich_log):
            exit(0)
        else:
            exit(1)
    if args.command == "nsort":
        bam = args.bam
        sortdir = args.outdir 
        name = os.path.join(sortdir, os.path.basename(bam))
        if os.path.exists(name.replace(".bam", ".json")):
            exit(0)
        from utilities.bam_tools import examine_sort
        if examine_sort(bam) == "coordinate" and not os.path.exists(name):
            sort_command = f"samtools sort -@ {n} -n -o {name} {bam}"
            print(sort_command)
            subprocess.call(sort_command, shell = True)
        bed = name.replace(".bam", ".bedpe")
        if examine_sort(name) == "queryname" and not os.path.exists(bed):
            from utilities.bam_tools import bam2bedpe
            bam2bedpe(name, bed)
            # get a small set of data and guess the column info
        qc_json = name.replace(".bam", ".json")
        ## report cis/trans yield 
        import json
        if not os.path.exists(qc_json):
            # final_log
            from utilities.bed_tools import cis_trans_bed
            log_dict = cis_trans_bed(bed)
            with open(qc_json, "w") as ff:
                json.dump(log_dict, ff)
        ## check output
    if args.command == "log": # parse log file 
        id = args.id
        if os.path.exists(f"/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/log/{id}.log"):
            exit(0)
        from utilities.parse_log import fastp_log_parser, flagstat_parser, mark_log_parser
        id_dict = {"ID":id}
        fastp_log = fastp_log_parser(f"/data/jim4/Seq/primary_cell_project/fastq/QC/HiTrAC/{id}.log")
        align_log = flagstat_parser(f"/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/raw/individual/{id}.stat")
        picard_log = mark_log_parser(f"/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/markdup/individual/{id}.txt")
        # bed = f"/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/sortn/{id}.bedpe"
        qc_json = f"/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/sortn/{id}.json"
        ## report cis/trans yield 
        import json
        # if not os.path.exists(qc_json):
        #     # final_log
        #     from utilities.bed_tools import cis_trans_bed
        #     log_dict = cis_trans_bed(bed)
        #     with open(qc_json, "w") as ff:
        #         json.dump(log_dict, ff)
        # else:
        with open(qc_json, "r") as fr:
            log_dict = json.load(fr)
        ## TSS enrichment 
        enrich_log = f"/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/individual/{id}.stat"
        with open(enrich_log, "r") as enf:
            enrich_dict = {line.strip("\n").split("\t")[0]:float(line.strip("\n").split("\t")[-1]) for line in enf.readlines()}
        # collect all log 
        log_all = {**id_dict, **fastp_log, **align_log, **picard_log, **log_dict, **enrich_dict}
        print(log_all)
        with open(f"/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/log/{id}.log", "w") as fo:
            json.dump(log_all, fo)
        # log_df = pd.DataFrame.from_dict(log_dict)

if __name__ == "__main__":
    main()

