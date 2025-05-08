#!/usr/bin/env python3
import os 
import json 
import subprocess
import argparse
import pandas as pd 

"""
Given the file path in the source location, use rsync to transfer indicated file (GC) to the destination.

example: 
rsync -C -P -e ssh jim4@137.187.135.165:/mnt/usbdisk/fastq_backup/file_*.json . # get json file 
rsync_proc.py file_*.json GC.txt -d . > rsyn.sh # generate bash file for data transfering 
sh ./rsyn.sh # run rsync 
"""

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="\n \
            Given the file path in the source location, use rsync to transfer indicated files (GC) to the destination. \n \
            example: \n \
            rsync -C -P -e ssh jim4@137.187.135.165:/mnt/usbdisk/fastq_backup/file_*.json . \n \
            rsync_proc.py file_*.json GC.txt > rsyn.sh \n \
            sh ./rsyn.sh ")
    parser.add_argument("json_path", help="json file contain all file path")
    parser.add_argument("GC_file", help = "GC number of interest for transfering/downloading")
    # parser.add_argument("-dir", "--directory", help = "master directory for all files", required = False)
    parser.add_argument("-d", "--destination", default = ".", help = "local destination to save the data")
    args=parser.parse_args()
    return args

def GC_list(GC_file):
    gc_file = pd.read_table(GC_file, header = None, sep = "\t", names = ["GC"])
    GC_list = gc_file["GC"].tolist()
    return GC_list

def remote_transfer(json_path, GC_list, destination):
    with open(json_path, "r") as file:
        all_files = json.load(file)
    for dir in all_files:
        file_list = [os.path.join(dir, f) for f in all_files[dir] if f.split("_")[0] in GC_list]
        if len(file_list) > 0:
            all_file4transfer = " :".join(file_list)
            run_number = os.path.basename(dir).split("_")[0]
            if not os.path.exists(os.path.join(destination, run_number)): 
                os.mkdir(os.path.join(destination, run_number))
            output_folder = os.path.join(destination, run_number)
            command = f"rsync -C -P -e ssh jim4@137.187.135.165:{all_file4transfer} {output_folder}"
            print(command)
            # subprocess.call(command, shell = True)

def main():
    args = args_parser()
    if args.json_path != None: 
        remote_transfer(args.json_path, GC_list(args.GC_file), args.destination)

if __name__ == "__main__":
    main()
