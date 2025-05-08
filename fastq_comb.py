#!/usr/bin/env python3

"""Give the directory of fastq file, identical GC number files are merged
Example:
fastq_comb.py -dir run0001/ -dest merged_dir/
"""

import subprocess 
import argparse
import pandas as pd 
import os 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-dir", "--directory", nargs = "+", help = "master directory sequencing data")
    parser.add_argument("-dest", "--destination", default = ".", help = "destination to store the merged files")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    dirs = args.directory
    dest = args.destination
    for dir in dirs:
        file_folder = [(dirpath, f) for dirpath, dirname, filename in os.walk(dir) for f in filename]
        file_folder_df = pd.DataFrame(file_folder, columns = ["folder", "file"])
    # extract GC and read (R1/R2) label
    file_folder_df["GC"] = file_folder_df["file"].str.split("_", expand = True)[0]
    file_folder_df["R"] = file_folder_df["file"].str.split("_", expand = True)[3]
    #
    for group in file_folder_df.groupby(by = ["GC", "R"]):
        group_files = [os.path.join(row.folder, row.file) for row in group[1].itertuples()]
        all_file = " ".join(sorted(group_files))
        new_file = os.path.join(dest, '_'.join(group[0]))
        if len(group_files) > 1:
            command = f"cat {all_file} > {new_file}.fastq.gz"
        else:
            command = f"mv {all_file} {new_file}.fastq.gz"
        if not os.path.exists(new_file+".fastq.gz"):
            print(command)
            subprocess.call(command, shell = True)

if __name__ == "__main__":
    main()


