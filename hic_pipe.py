#!/usr/bin/env python3
"""
This is a pipeline executing the following steps in hic_align.py. 
1) trim read ; 
2) read aligment; 
3) markduplicate; 
4) bam QC; 

"""
import os 
import pandas as pd 
import argparse 
import sys 
import subprocess
from utilities.misc import ignore_warning
import glob 
import os 

ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter, 
                                     usage = "")
    parser.add_argument("sample", help = "samples in a xlsx file to process")
    parser2 = argparse.ArgumentParser(prog = "PROG", add_help=True)
    sub_parsers = parser2.add_subparsers(dest = "mode", help = "mode to run")
    process = sub_parsers.add_parser("process", help = "process sample data", parents = [parser], add_help = False)
    process.add_argument("-dir", "--directory", help = "directory where to find the fastq data")
    process.add_argument("-n", "--thread", help = "thread number (default: 4)", default = 4, type = int)
    process.add_argument("-f", "--force", action = "store_true", help = "force to rerun if when output file already exists!")
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

def format_excel(out, sample_report):
    """
    Write sample report in style
    out: str. output excel name
    sample_report: df. dataframe to be styled. 
    """
    import string
    letters = list(string.ascii_uppercase)
    # 
    yield_col = [col for col in sample_report.columns if col.endswith("yield")]
    cis_yield_col = [c for c in sample_report.columns if "cis_yield" in c][0]
    cis_valid_yield_col = [c for c in sample_report.columns if "cis_valid_yield" in c][0]
    sample_report.rename(columns = {cis_yield_col:"cis_yield", cis_valid_yield_col:"cis_valid_yield"}, inplace = True)
    col_list = sample_report.columns.tolist()
    # print(sample_report)
    with pd.ExcelWriter(out, mode = "w") as writer:
        sample_report.to_excel(writer, index = False, sheet_name = "Sheet1")
        workbook = writer.book 
        worksheet = writer.sheets["Sheet1"]
        #int_format = workbook.add_format({"num_format": "#,##0"})
        bg_format = workbook.add_format({"bg_color":"#fff7c7"})
        bg_format_low = workbook.add_format({"bg_color":"#898481"})
        #tot_col = col_list.index("total_reads")
        #worksheet.set_column(f"{tot_col}:{tot_col}", None, int_format)
        # add background color to yield columns 
        for y in yield_col:
            format_col = letters[col_list.index(y)]
            worksheet.set_column(f"{format_col}:{format_col}", None, bg_format)
        for row in sample_report.itertuples():
            if row.trim_yield < 0.8:
                worksheet.write(row.Index+1, col_list.index("trim_yield"), sample_report.iloc[row.Index, col_list.index("trim_yield")], bg_format_low)
            if row.align_ratio < 0.8:
                "low quality aligned"
                worksheet.write(row.Index + 1, col_list.index("align_ratio"), sample_report.iloc[row.Index, col_list.index("align_ratio")], bg_format_low)
            if row.dup_rate > 0.6: 
                "high duplication"
                worksheet.write(row.Index + 1, col_list.index("dup_rate"), sample_report.iloc[row.Index, col_list.index("dup_rate")], bg_format_low)
            if row.cis_yield < 0.5:
                "low cis"
                worksheet.write(row.Index + 1, col_list.index("cis_yield"), sample_report.iloc[row.Index, col_list.index("cis_yield")], bg_format_low)
            # if row.cis_valid_yield < 0.15:
            #     "low cis with long distance"
            #     worksheet.write(row.Index + 1, col_list.index("cis_valid_yield"), sample_report.iloc[row.Index, col_list.index("cis_valid_yield")], bg_format_low)
            # if row.enrich < 1.25:
            #     worksheet.write(row.Index + 1, col_list.index("enrich"), sample_report.iloc[row.Index, col_list.index("enrich")], bg_format_low)

def main():
    args = args_parser()
    sample_df = read_df(args.sample)
    # process samples 
    if args.mode == "summary": 
        import json
        log_dir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/log"
        jlist = list()
        for row in sample_df.itertuples():
            logf = os.path.join(log_dir, f"{row.id}.log")
            try:
                with open(logf, "r") as logj:
                    jdict = json.load(logj)
                    jdf = pd.DataFrame.from_dict(jdict, orient = "index")
                    jlist.append(jdf)
            except Exception as e:
                print(e)
        sample_report = pd.concat(jlist, axis = 1).transpose()
        sample_report = pd.merge(sample_df, sample_report, left_on = "id", right_on = "ID").drop("ID", axis = 1)
        short_col = [c for c in sample_report.columns if "short_distance" in c][0]
        sample_report["final_usable_PETs"] = sample_report["cis"] - sample_report["short_distance (<1 kb)"]
        sample_report["final_yield%(long_distance_PETs/total_reads)"] = (sample_report["cis"]-sample_report[short_col])/sample_report["total_reads"]*100
        out = args.sample.replace(".xlsx", "_report.xlsx")
        sample_report.drop(["total", "pass_markdup_reads", "short_distance (<1 kb)"], axis = 1, inplace = True)
        # with pd.ExcelWriter(out, mode = "w") as writer:
        #     sample_report.to_excel(writer, index = False, sheet_name = "Sheet1")
        format_excel(out, sample_report)
    if args.mode == "process":
        dir = args.directory
        n = args.thread
        for row in sample_df.itertuples():
            print(f"Processing {row.id}")
            read = os.path.join(dir, row.id)
            print(read)
            # preprocessing step 
            pre_step = f"hic_align.py pre -read {read} -n {n}"
            # subprocess.call(pre_step, shell = True)
            qc_dir = "/data/jim4/Seq/primary_cell_project/fastq/QC/HiTrAC"
            qc_read = os.path.join(qc_dir, row.id)
            # alignment to reference genome
            # if row.ref == "human":
            #     ref = "/data/jim4/Reference/human/GRCh38.p14/fasta/GRCh38.primary_assembly.genome.fa"
            align_step = f"hic_align.py align -read {qc_read} -n {n}"
            # subprocess.call(align_step, shell = True)
            # markdup alignment file 
            align_dir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/raw/individual"
            align_bam = os.path.join(align_dir, row.id + ".bam")
            mark_step = f"hic_align.py markdup -bam {align_bam}"
            # subprocess.call(mark_step, shell = True)
            # qc alignment file 
            mark_dir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/markdup/individual"
            mark_bam = os.path.join(mark_dir, row.id + ".bam")
            qc_step = f"hic_align.py qc -bam {mark_bam} -n {n}"
            # subprocess.call(qc_step, shell = True)
            qc_dir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/individual"
            qc_bam = os.path.join(qc_dir, row.id + '.bam')
            name_sort = f"hic_align.py nsort -bam {qc_bam} -n {n}"
            # subprocess.call(name_sort, shell = True)
            log_command = f"hic_align.py log -id {row.id}"
            # subprocess.call(log_command, shell = True)
            if args.force or (not os.path.exists(qc_read + "_R1.fastq.gz") and not os.path.exists(qc_read + "_R2.fastq.gz")):
                print(pre_step)
                subprocess.call(pre_step, shell = True)
            if args.force or  (not os.path.exists(align_bam)):
                print(align_step)
                subprocess.call(align_step, shell = True)
            if args.force or (not os.path.exists(mark_bam)):
                print(mark_step)
                subprocess.call(mark_step, shell = True)
            if args.force or (not os.path.exists(qc_bam) or not os.path.exists(qc_bam.replace(".bam", ".stat"))):
                print(qc_step)
                subprocess.call(qc_step, shell = True)
            subprocess.call(name_sort, shell = True)
            subprocess.call(log_command, shell = True)
            # if 
            #if not os.path.exists()
            #consecutive_steps = " && ".join([pre_step, align_step, mark_step, qc_step, name_sort, log_command])
            #subprocess.call(consecutive_steps, shell = True)
        # # cis bam file 
        # qcdir = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/individual"
        # qc_bam = os.path.join(qcdir, row.id + ".bam")
        # cis_step = f"hic_align.py cis -bam {qc_bam} -n {args.thread} -o {row.id}"
        # subprocess.call(cis_step, shell = True)
    # if args.command == "process":
    #     process_sample(args, sample_df, dir)
    # if args.command == "merge":
    #     bam_merge(args, sample_df, dir)
    # if args.command == "transform":
    #     chrsize = args.chrsize
    #     for sample, sample_df in sample_df.groupby("sample"):
    #         sample_bam = os.path.join(dir, sample + ".bam")
    #         if not os.path.exists(sample + ".pairs"):
    #             pair_generator = f"hic_format.py pairs {sample_bam} -o {sample} -chrsize {chrsize} -n {args.thread}"
    #             subprocess.call(pair_generator, shell = True)
    #         if not os.path.exists(sample + ".cool"):
    #             cool_generator = f"hic_format.py cool {sample}.pairs -o {sample} -chrsize {chrsize} -n {args.thread}"
    #             subprocess.call(cool_generator, shell = True)
    #             blacklist= "/data/jim4/data/blacklist/ENCFF356LFX.bed"
    #             dchic_command = f"hic_format.py dchic {sample}.mcool -o {sample} -chrsize {chrsize} -n {args.thread} -blacklist {blacklist} -r 25000 50000 100000 250000 500000"
    #             subprocess.call(dchic_command, shell = True)
    #         if not os.path.exists(sample + ".hic"):
    #             hic_generator = f"hic_format.py hic {sample}.pairs -o {sample} -chrsize {chrsize} -n {args.thread}"
    #             subprocess.call(hic_generator, shell = True)

if __name__ == "__main__":
    main()
