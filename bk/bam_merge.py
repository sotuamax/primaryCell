#!/usr/bin/env python3
import os 
import subprocess
import argparse 

"""merge multiple bam files, and generate index file for merged bam."""

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-bam", "--bam", nargs = "+", help="bams to merge")
    parser.add_argument("-o", "--output", help = "output bam name in prefix", required = False)
    parser.add_argument("-r", "--header", help = "read group header (e.g., 'SM:samplename')", type = str, required = False)
    parser.add_argument("-n", "--name", help = "sort by name", action = "store_true")
    parser.add_argument("-t", "--thread", help = "thread number (default: 12)", default = 24)
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    bam = args.bam
    bam_dir = os.path.dirname(os.path.abspath(bam[0]))
    o = args.output
    if o == None: 
        o = "_".join([b.split(".bam")[0] for b in bam])
    n = args.thread
    h = args.header
    # 
    prefix = os.path.basename(o).split(".bam")[0]
    subprocess.call("ml samtools ", shell = True)
    file_combined = " ".join(bam)
    if h != None:
        l_dict = {l.split(":")[0]:l.split(":")[-1] for l in h.split(";")}
        RG_label = "".join([f'@RG\tID:{f}\tSM:{l_dict["SM"]}\tLB:{f}\tPL:ILLUMINA\n' for f in bam])
        command = f"printf '{RG_label}' > {bam_dir}/{prefix}.txt"
        subprocess.call(command, shell = True)
        if args.name:
            merge_command = f'samtools merge -rh {bam_dir}/{prefix}.txt -@ {n} -o - {file_combined} | samtools sort -@ {n} - -n -o {o}.bam && rm {bam_dir}/{prefix}.txt'
        else:
            merge_command = f'samtools merge -rh {bam_dir}/{prefix}.txt -@ {n} -o {o}.bam {file_combined} && samtools index -@ {n} {o}.bam && rm {bam_dir}/{prefix}.txt'
    else:
        if args.name:
            merge_command = f'samtools merge -@ {n} -o - {file_combined} | samtools sort -@ {n} - -n -o {o}.bam'
        else:
            merge_command = f'samtools merge -@ {n} -o {o}.bam {file_combined} && samtools index -@ {n} {o}.bam'

    print(merge_command)
    subprocess.call(merge_command, shell = True)

if __name__ == "__main__":
    main()

