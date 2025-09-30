#!/usr/bin/env python3
import os 
import subprocess 
import argparse 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-read", "--read_name", help="read name to used for alignment. ")
    parser.add_argument("-ref", "--reference", help = "assign reference used for the alignment (human/mouse)", default = "human")
    parser.add_argument("-outdir", "--outdir", help = "output directory")
    parser.add_argument("-n", "--thread", help = "thread number (default: 12)", default = 12)
    parser.add_argument("-mode", "--mode", help = "input fastq format in PE or SE mode", default = "PE")
    args=parser.parse_args()
    return args

def retrieve_fastq(read, mode = "PE"):
    if mode == "PE":
        r1 = read + "_R1.fastq.gz"
        r2 = read + "_R2.fastq.gz"
        return (r1, r2)

def align(args):
    read = args.read_name
    outdir = args.outdir
    n = args.thread
    r1, r2 = retrieve_fastq(read, mode = args.mode)
    # 
    bam_out = os.path.basename(read)
    bam_out = os.path.join(outdir, bam_out)
    # 
    if args.reference == "human":
        print("Run for human ....... ")
        if os.path.exists(r1) and os.path.exists(r2):
            if not os.path.exists(bam_out + ".bam"):
                align_command = f"bwa mem -t {n} $hg38 {r1} {r2} | samtools view -@ {n} -Su - | samtools sort -@ {n} - -o {bam_out}.bam && samtools index -@ {n} {bam_out}.bam"
                print(align_command)
                subprocess.call(align_command, shell = True)
        else:
            print(r1 + " and " + r2 + " not found")
    elif args.reference == "mouse": 
        print("Run for mouse ....... ")
        align_command = f"bwa mem -t {n} $mm39 {r1} {r2} | samtools view -@ {n} -Su - | samtools sort -@ {n} - -o {bam_out}.bam && samtools index -@ {n} {bam_out}.bam"
        subprocess.call(align_command, shell = True)
    else:
        print("The reference is not supported! Please provide valid reference for alignment. ")
        exit(0)
        if len(read) == 1:
            print("SE bwa alignment is not supported yet.")
            exit(0)

def main():
    args = args_parser()
    align(args)

if __name__ == "__main__":
    main()
