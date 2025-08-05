#!/usr/bin/env python3
"""
Take sample excel sheet, and extract sample to run. 

"""
import pandas as pd 
import numpy as np 
import argparse 
import subprocess 
import os 
import sys
from datetime import date

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # 
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help = True)
    sub_parsers = parser2.add_subparsers(dest = 'command', help = "mode to run")
    run = sub_parsers.add_parser("run", help = "run standard steps to identify tagged genes", parents = [parser], add_help = False)
    run.add_argument("sample", help = "samples to process in a excel file")
    run.add_argument("-dir", "--directory", help = "directory where to find the fastq data")
    run.add_argument("-outdir", "--outdir", help = "directory where to write output")
    run.add_argument("-primer")
    run.add_argument("-gRNA", "--gRNA", required = False, help = "gRNA table")
    run.add_argument("-n", "--threads", default = 8, type = int, help = "threads for alignment process")
    # 
    summary = sub_parsers.add_parser("summary", help = "summarize all samples report", parents = [parser], add_help = False)
    summary.add_argument("result_folder", help = "result folder to summarize sample information from")
    summary.add_argument("-o", "--output", required = False, help = "summary output prefix name")
    args = parser2.parse_args()
    return args

def gene_parse(gene_file):
    """
    Parse tagged genes dataframe 
    """
    gene_df = pd.read_table(gene_file, sep = "\t", header = 0)
    total_gene_n = len(set(gene_df["symbol"]))
    return (gene_df, total_gene_n)

def uni_barcode_gene_num(gene_file, min_cnt = 0):
    """
    Parse gene dataframe on conditon of unique barcode
    """
    gene_df = pd.read_table(gene_file, sep = "\t", header = 0)
    gene_df_filtered = gene_df.query("read_count > @min_cnt").copy()
    uni_barcode = (gene_df_filtered.groupby("barcode")["symbol"].nunique() == 1).index
    unibarcode_gene_df = gene_df_filtered.query("barcode in @uni_barcode")
    return len(set(unibarcode_gene_df["symbol"]))

def collect_sample_info(result_dir):
    """
    Given directory of experiment results folder, read the "log" and "gene" file, to report sample summary.
    Returns: 
    master_df (dataframe): contains read-statistics and #genes for each sample 
    inframe_df (dataframe): inframe genes detected for all samples are pooled together
    gene_df_dict (dictionary): sample name as key, and its genes of all frame as item. 
    """
    master_tab = list()
    gene_df_dict = dict()
    inframe_gene_dict = dict()
    for dirpath,dirname,filename in os.walk(result_dir):
        for d in dirname:
            sample = d.split(".")[0]
            print(sample)
            # per d is a per sample folder
            for dp,dn,fn in os.walk(os.path.join(result_dir, d)):
                for f in fn: 
                    if not f.startswith("."):
                        if f.endswith(".log"):
                            log_file = os.path.join(dp, f)
                        if f.endswith("_gene.txt"):
                            gene_file = os.path.join(dp, f)
                        else:
                            gene_df, gene_n = 0, 0
                        if f.endswith("_inframe.txt"):
                            tag_gene_file = os.path.join(dp, f)
                        else:
                            tag_gene_df, tag_gene_n = 0, 0
            tot_r, cln_r, ali_r, fet_r, frame,frame0,frame0i, frame1,frame1i,frame2,frame2i = log_parse(log_file)
            gene_df, gene_n = gene_parse(gene_file)
            tag_gene_df, tag_gene_n = gene_parse(tag_gene_file)
            if "barcode" in tag_gene_df.columns:
                sample_df = pd.DataFrame({"sample":[sample], "total_read":[int(tot_r)], "clean_read":[cln_r], "align_read":[ali_r], "CDS_align":[fet_r], "frame":frame, "F0 gene":frame0, "F0 read":frame0i, "F1 gene":frame1, "F1 read":frame1i, "F2 gene":frame2, "F2 read":frame2i, "total_gene":[gene_n], "gene_inframe":[tag_gene_n], "barcode_single_target(read_count>0)":uni_barcode_gene_num(tag_gene_file, 0), "barcode_single_target(read_count>10)":uni_barcode_gene_num(tag_gene_file, 10), "barcode_single_target(read_count>20)":uni_barcode_gene_num(tag_gene_file, 20), "barcode_single_target(read_count>50)":uni_barcode_gene_num(tag_gene_file, 50)})
            else:
                sample_df = pd.DataFrame({"sample":[sample], "total_read":[int(tot_r)], "clean_read":[cln_r], "align_read":[ali_r], "CDS_align":[fet_r], "frame":frame,"F0 gene":frame0, "F0 read":frame0i, "F1 gene":frame1, "F1 read":frame1i, "F2 gene":frame2, "F2 read":frame2i, "total_gene":[gene_n], "gene_inframe":[tag_gene_n]})
            master_tab.append(sample_df)
            gene_df_dict[sample] = gene_df #.to_excel(writer, sheet_name = sample, index = False)
            inframe_gene_dict[sample] = tag_gene_df # .to_excel(writer, sheet_name = sample+"_inframe")
    #### concatenate df in list
    master_df = pd.concat(master_tab, axis = 0).sort_values(by = ["sample"], ignore_index = True)
    inframe_list = list()
    for s in inframe_gene_dict:
        sample_i_df = inframe_gene_dict[s].copy()
        sample_i_df["sample"] = s 
        inframe_list.append(sample_i_df)
    inframe_df = pd.concat(inframe_list, axis = 0, ignore_index = True)
    symbol, symbol_freq = np.unique(inframe_df["symbol"], return_counts = True)
    symbol_freq_df = pd.DataFrame.from_dict(dict(zip(symbol,symbol_freq)), orient = "index", columns = ["freq"])
    inframe_df = pd.merge(inframe_df, symbol_freq_df, left_on = "symbol", right_index = True, how = "left")
    return (master_df, inframe_df, gene_df_dict)

def filter_true(inframe_df):
    if "gRNA_target" in inframe_df.columns:
        inframe_df = inframe_df[inframe_df["gRNA_target"] == True]
    if "pass_donor_filter" in inframe_df.columns:
        inframe_df = inframe_df[inframe_df["pass_donor_filter"] == True]
    return inframe_df

def run_process(args, sample_df, current_dir):
    """
    Run each step for gene tag detection. 
    Step 1: pre
    """
    for row in sample_df.itertuples():
        # create the directory of experiment as container of sample's results 
        sample = row.sample
        GFP_sample = sample + "GFP"
        ChiC_sample = sample + "0"
        primer = row.primer
        frame = row.frame
        ref = row.reference 
        n = args.core
        # step 1: preprocessing step
        pre_step = f"mpiexec -n {n} tag_preprocessing.py -sample {GFP_sample} -dir {dir} -outdir {}"
        # genome matched gtf and star_index 
        if ref == "mouse":
            gtf = "/data/jim4/Reference/mouse/GRCm39/GTF/gencode.vM34.primary_assembly.annotation.gtf"
            star_ref = "/data/jim4/Reference/mouse/GRCm39/STAR_overhang50"
            bwa_ref = "/data/jim4/Reference/mouse/GRCm39/fasta/GRCm39.primary_assembly.genome.fa"
        if ref == "human":
            # update on 05/06/2025 (GTF file for human protein coding genes only)
            gtf = "/data/jim4/Reference/human/GRCh38.p14/GTF/gencode.v45.primary_assembly.annotation.protein_coding.gtf"
            star_ref = "/data/jim4/Reference/human/GRCh38.p14/STAR_overhang_protein"
            bwa_ref = "/data/jim4/Reference/human/GRCh38.p14/fasta/GRCh38.primary_assembly.genome.fa"
        # step 2: alignment 
        sample_dir = sample + date.today().strftime(".%m-%d-%Y")
        if dtype == "cDNA":
            align_step = f"g_target.py align -dtype {dtype} -fastq {sample_dir}/{sample}_R1.fq {sample_dir}/{sample}_R2.fq -ref {star_ref} -gtf {gtf} -o {sample} -n {args.align_threads}"
        if dtype == "DNA":
            align_step = f"g_target.py align -dtype {dtype} -fastq {sample_dir}/{sample}_R1.fq {sample_dir}/{sample}_R2.fq -ref {bwa_ref} -gtf {gtf} -o {sample} -n {args.align_threads}"
        # set up alignment folder 
        if args.align_folder != None:
            align_folder = args.align_folder
            align_step = align_step + f" -dir {align_folder}"
        else:
            align_folder = "/data/jim4/Seq/CRISPR/alignment"
        # set up single cell mode 
        if args.sc:
            barcode = "/data/jim4/Seq/CRISPR/data/multiplex.xlsx"
            align_step = align_step + f" -multiplex {barcode}"
        r1_len = fq2len(f"{sample_dir}/{sample}_R1.fq") # get read length for alignment 
        print("Alignment step ....")
        print(align_step)
        subprocess.call(align_step, shell = True)
        # store all alignment bam file here: /data/jim4/Seq/CRISPR/alignment
        if dtype == "cDNA":
            bam_f = os.path.join(align_folder, sample+f"L{r1_len}_Aligned.sortedByCoord.out.bam")
        if dtype == "DNA":
            bam_f = os.path.join(align_folder, sample+f"L{r1_len}.bam")
        bed_f = bam_f.replace(".bam", ".bed")
        # step 3: screen alignment file for genes intag
        if os.path.exists(bam_f) and os.path.exists(bed_f):
            screen_step = f"g_target.py screen -dtype {dtype} -bam {bam_f} -gtf {gtf} -frame {frame[-1]} -o {sample}"
            if args.sc:
                barcode = "/data/jim4/Seq/CRISPR/data/multiplex.xlsx"
                screen_step += f" -multiplex {barcode}"
            if donor != "-":
                screen_step = screen_step + f" -donor {donor}"
            if assay.upper() == "SA":
                screen_step = screen_step + f" -side_end"
            if not args.gRNA is None:
                screen_step = screen_step + f" -gRNA {args.gRNA}"
            print("Screen step ....")
            print(screen_step)
            subprocess.call(screen_step, shell = True)
        else:
            print("bam/bed file does not exist!")
            exit(1)
        os.chdir(current_dir)

def main():
    args = args_parser()
    # preprocessing 
    sample = args.sample
    dir = args.directory 
    SA = args.primer
    output_dir = args.outdir 
    n = args.threads
    GFP_sample = sample + "GFP"
    ChiC_sample = sample 
    # preprocessing 
    pre_GFP = f"mpiexec -n {n} tag_preprocessing.py -sample {GFP_sample} -dir {dir} -outdir {output_dir} -primer {SA} -barcode /data/jim4/Seq/ChIP_seq/data/barcode.txt"
    pre_ChiC = f"mpiexec -n {n} tag_preprocessing.py -sample {sample} -dir {dir} -outdir {output_dir} -barcode /data/jim4/Seq/ChIP_seq/data/barcode.txt"
    # align 
    hg38_gtf_pick="/data/jim4/Reference/human/GRCh38.p14/GTF/gencode.v45.primary_assembly.annotation.protein_coding.gtf.gz"
    align_GFP = f"tag_align.py -sample {GFP_sample} -seq_type cDNA -outdir {output_dir} -ref $hg38_star -gtf {hg38_gtf_pick} -n {n}"
    align_ChiC = f"tag_align.py -sample {sample} -seq_type DNA -outdir {output_dir} -ref $hg38 -n {n} -qc -report"
    # generate flag summary for GFP
    tag_GFP = f"tag_summary.py -sample {GFP_sample} -outdir {output_dir}"
    # generate html report for GFP and ChiC 
    html = f"tag_html.py -sample {sample} -outdir {output_dir}"
    # 

if __name__ == "__main__":
    main()
