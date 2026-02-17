#!/usr/bin/env python3
import subprocess 
import argparse 
import os 
import time 
from utilities.misc import timeit

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = True, 
                                     usage="hic_format.py -o <out> inputfile -n 10")
    parser.add_argument("-i", "--input", nargs="+", help = "multiple input files used to merge.")
    parser.add_argument("-r", "--resolution", nargs = "*", default = [1000], type = int, required = False, help = "output hic/cool file in the requested resolutions. ")
    parser.add_argument("-n", "--thread", default = 2, type = int, help = "number of processes")
    parser.add_argument("-o", "--output", required = True, help = "output file name (prefix)")
    args = parser.parse_args()
    return args

def main():
    args = args_parser()
    input_files = args.input; files = " ".join(input_files)
    n = args.thread 
    output = args.output
    maxt = 0 
    for f in input_files:
        ftime = os.path.getmtime(f)
        maxt = maxt if ftime < maxt else ftime 
    if not os.path.exists(output + ".cool") or os.path.getmtime(output + ".cool") < maxt:
        print("To merge cool files into one ...")
        command = f"cooler merge {output}.cool {files}"
        print(command)
        subprocess.call(command, shell = True)
    else:
        print(f"{output}.cool already exists!")
    if len(args.resolution) > 1:
        # when assigned resolution more than 1, after merge, also generate a multiple resolutions cooler file (mcool)
        all_resolutions = ",".join([str(r) for r in args.resolution])
        mcool_command = f"cooler zoomify -n {n} -r {all_resolutions} -o {output}.mcool {output}.cool"
        if not os.path.exists(output + ".mcool") or os.path.getmtime(output + ".mcool") < os.path.getmtime(output + ".cool"):
            print(mcool_command)
            print(f"Run cooler w/ {n} threads ...")
            subprocess.call(mcool_command, shell = True)
        else:
            print(f"{output}.mcool already exists!")

if __name__ == "__main__":
    main()