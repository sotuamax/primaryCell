#!/usr/bin/env python3
import argparse 
import subprocess
import os 
import glob 


def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # master parameters 
    parser.add_argument("-n", "--threads", help = "threads to use", default = 4, type = int)
    parser.add_argument("-f", "--force", action = "store_true", help = "force regenerate alignment file when it exists")
    parser.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/CRISPR/alignment/individual", help = "output directory for output bam alignment file")
    # separate argumentparser for sub-parsers
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    # initiate sub-parser 
    sub_parsers = parser2.add_subparsers(dest = 'command', help = "mode to run")
    # alignment
    align = sub_parsers.add_parser("align", help = "perform alignment", parents = [parser], add_help = False)
    align.add_argument("-read", "--read", nargs = "+", help="prefix of read file")
    align.add_argument("-ref", "--reference", default = "/data/jim4/Reference/human/GRCh38.p14/fasta/GRCh38.primary_assembly.genome.fa", help = "reference genome for alignment")
    # align.add_argument("-mode", "--mode", choices= ["SE", "PE"], help = "Single-end or Pair-end alignment.", default = "SE")
    # markdup
    markdup = sub_parsers.add_parser("markdup", help = "perform feature count", parents = [parser], add_help = False)
    markdup.add_argument("-bam", "--bam", required = True, help = 'bam alignment file')
    # markdup.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/CRISPR/alignment/markdup", help = "directory for markdup bam file")
    # qc 
    qc = sub_parsers.add_parser("qc", help = "QC on BAM file", parents = [parser], add_help = False)
    qc.add_argument("-bam", "--bam", required = True, help = 'bam alignment file')
    # qc.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/CRISPR/alignment/QC", help = "directory for QC bam file")
    # transform to bigwig format
    bw = sub_parsers.add_parser("bigwig", help = "transform to bigwig format", parents = [parser], add_help = False)
    bw.add_argument("-bam", "--bam", required = True, help = "bam alignment file")
    # bw.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/CRISPR/alignment/bigwig", help = "directory for bigwig file")
    # generate log file 
    # log = sub_parsers.add_parser("log", help = "generate log file for each alignment", parents = [parser], add_help = False)
    # log.add_argument("-sample", "--sample", help = "sample ID")
    # # log.add_argument("-mode", "--mode", choices= ["SE", "PE"], help = "Single-end or Pair-end alignment.", default = "SE")
    # log.add_argument("-outdir", "--outdir", default = "/data/jim4/Seq/CRISPR/log", help = "directory for log file")
    # parse arguments
    args=parser2.parse_args()
    return args

def main():
    args = args_parser()
    n = args.threads
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    if args.command == "align":
        ref = args.reference
        read = args.read 
        r1 = " ".join(read)
        # r1 = glob.glob(f"{read}_*R1.fastq.gz")[0]; r2 = glob.glob(f"{read}_*R2.fastq.gz")[0]
        name = os.path.join(args.outdir, os.path.basename(read[0]).split("_")[0])
        align_command = f"bwa mem {ref} {r1} -t {n} | samtools view -@ {n} -Su - | samtools sort -@ {n} - -o {name}.bam && samtools index -@ {n} {name}.bam"
        if not os.path.exists(name + ".bam"):
            print("Perform alignment ...")
            subprocess.call(align_command, shell = True)
        
        align_stat = f"samtools flagstat -@ {n} {name}.bam > {name}.stat"
        if not os.path.exists(name + ".stat"):
            subprocess.call(align_stat, shell=True)
        if args.force:
            print("Enforce alignment ...")
            subprocess.call(align_command, shell = True)
            subprocess.call(align_stat, shell = True)
    if args.command == "markdup":
        bam = args.bam 
        picard="/usr/local/apps/picard/3.1.0/picard.jar"
        name = os.path.join(args.outdir, os.path.basename(bam).replace(".bam", ".dedup.bam"))
        metric = name.replace(".bam", ".metric")
        # picard markdup requires bam sorted by chromosome coordinates 
        markdup_command = f"java -Xmx50g -jar {picard} MarkDuplicates -I {bam} -O {name} -M {metric} --REMOVE_DUPLICATES true && samtools index -@ {n} {name}"
        print(markdup_command)
        if not os.path.exists(metric):
            print("Mark duplicates in BAM ...")
            subprocess.call(markdup_command, shell = True)
        else:
            print("deduplicate BAM exists!")
        if args.force:
            print("Enforce markdup ...")
            subprocess.call(markdup_command, shell = True)
    if args.command == "qc":
        bam = args.bam; outdir = args.outdir
        from utilities.bam_tools import bamfilter
        name = os.path.join(outdir, os.path.basename(bam).replace(".bam", ".qc.bam"))
        if not os.path.exists(name) or args.force:
            print("QC on BAM ...")
            bamfilter(bam, name, clip_check=False, threads = n, chrom = True)
        else:
            print("QC BAM exists.")
    if args.command == "bigwig":
        bam = args.bam
        name = os.path.join(args.outdir, os.path.basename(bam).replace(".bam", ".bw"))
        bw_command = f"bamCoverage -b {bam} -o {name} -of bigwig -p {n} --normalizeUsing CPM"
        # bw None normalization ()
        if not os.path.exists(name):
            print("Generate bigwig file ...")
            subprocess.call(bw_command, shell=True)
        else:
            print("Bigwig exists!")
        if args.force:
            print("Enforce bigwig ...")
            subprocess.call(bw_command, shell=True)
    
    # if args.command == "log":
    #     from utilities.parse_log import flagstat_parser_SE, mark_log_parser, flagstat_parser
    #     sample = args.sample 
    #     id_dict = {"sample":sample}
    #     if args.mode == "SE":
    #         align_log = flagstat_parser_SE(f"/data/jim4/Seq/CRISPR/alignment/individual/{sample}.stat")
    #     else:
    #         align_log = flagstat_parser(f"/data/jim4/Seq/CRISPR/alignment/individual/{sample}.stat")
    #     picard_log = mark_log_parser(f"/data/jim4/Seq/CRISPR/alignment/markdup/{sample}.metric")
    #     del picard_log["pass_markdup_reads"] #.remove("pass_markdup_reads")
    #     log_all = {**id_dict, **align_log, **picard_log}
    #     import json 
    #     print("Write into log ...")
    #     with open(os.path.join(args.outdir, f"{sample}.log"), "w") as fo:
    #         json.dump(log_all, fo)

if __name__ == "__main__":
    main()
