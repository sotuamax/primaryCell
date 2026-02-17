#!/usr/bin/env python3
"""
update(01/29/2026.12PM): support intra-chromosomal reads output in pairs file 

Transform BAM input file into pairs/cool/hic.

Input: 
    required: 
    - bam: bam file (recommend in name sorted)
    - ouput: output name 
    optional:
    - assembly: genome assembly name 
    - threads: threads to use 
    - chrsize: chromosome size file 

===============
Example: 
source myconda 
mamba activate bio 
ml pairtools juicer

hic_format.py -n 12 -o <sample> inputfile 


# for BAM to pairs 
hic_format.py -n 12 file.bam -of pairs -o sample 
# from pairs to cool 
hic_format.py -n 12 file.pairs.gz -of cool -o sample 

"""
import subprocess 
import argparse 
import os 
import time 
from utilities.misc import timeit

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = True, 
                                     usage="hic_format.py -o <out> inputfile -n 10")
    parser.add_argument("input", help = "input file (either BAM file coordinate or name sorted, provide name sorted file to save time; or pairs file sorted and compressed for UU reads)")
    parser.add_argument("-of", "--output_format", help = "output file format (output cool is not multiple-threaded)", choices = ["pairs", "cool", "hic", "transform"])
    parser.add_argument("-cis", "--cis", action = "store_true", help = "use cis (intra-chromosomal) data only.")
    parser.add_argument("-r", "--resolution", nargs = "*", default = [1000], type = int, required = False, help = "output hic/cool file in the requested resolutions. ")
    parser.add_argument("-o", "--output", required = True, help = "output file name (prefix)")
    parser.add_argument("-n", "--thread", default = 2, type = int, help = "thread to process bam")
    parser.add_argument("-chrsize", "--chrsize", required = False, default = "/data/jim4/Reference/human/GRCh38.p14/GRCh38_chrsize.bed", help = "file for ordered chromosome and size")
    parser.add_argument("-assembly", "--assembly", default = "hg38", help = "genome assembly name")
    args = parser.parse_args()
    return args

def main():
    start = time.time()
    args = args_parser()
    input_file = args.input
    n = args.thread
    out = args.output
    chrsize = args.chrsize 
    assembly = args.assembly
    if args.output_format == "pairs":
        if not input_file.endswith(".bam"):
            raise ValueError("Input file is not in BAM format!")
        bam = input_file 
        from utilities.bam_tools import examine_sort
        try:
            b_sort = examine_sort(bam)
            print("BAM is sorted by", b_sort)
        except Exception as e: 
            print(e)
            exit(1)
        if b_sort == "coordinate" and not args.cis:
            # only keep bwa header line in the pairs file
            pair_command = f"""samtools sort -@ {n} -n {bam} | samtools view -@ {n} -h | pairtools parse -c {chrsize} --assembly {assembly} --nproc-in {n} --nproc-out {n} --drop-sam | pairtools select --nproc-in {n} --nproc-out {n} '(pair_type==\"UU\")' | pairtools sort --memory 50G --nproc-in {n} --nproc-out {n} -o {out}.pairs.gz && pairix -p pairs {out}.pairs.gz"""
        if b_sort == "coordinate" and args.cis:
            pair_command = f"""samtools sort -@ {n} -n {bam} | samtools view -@ {n} -h | pairtools parse -c {chrsize} --assembly {assembly} --nproc-in {n} --nproc-out {n} --drop-sam | pairtools select --nproc-in {n} --nproc-out {n} '(pair_type=="UU") and (chrom1==chrom2)' | pairtools sort --memory 50G --nproc-in {n} --nproc-out {n} -o {out}.pairs.gz && pairix -p pairs {out}.pairs.gz"""
        if b_sort == "queryname" and not args.cis:
            pair_command = f"""samtools view -@ {n} -h {bam} | awk '$1!="@PG" || ($0 ~ /ID:bwa/)' | pairtools parse -c {chrsize} --assembly {assembly} --nproc-in {n} --nproc-out {n} --drop-sam | pairtools select --nproc-in {n} --nproc-out {n} '(pair_type==\"UU\")' | pairtools sort --memory 50G --nproc-in {n} --nproc-out {n} -o {out}.pairs.gz && pairix -p pairs {out}.pairs.gz"""
        if b_sort == "queryname" and args.cis:
            pair_command = f"""pairtools parse --cmd-in 'samtools view -@ {n} -h {bam}' -c {chrsize} --assembly {assembly} --nproc-in {n} --nproc-out {n} --drop-sam | pairtools select --nproc-in {n} --nproc-out {n} '(pair_type==\"UU\") and (chrom1==chrom2)' | pairtools sort --nproc-in {n} --nproc-out {n} --memory 50G -o {out}.pairs.gz && pairix -p pairs {out}.pairs.gz"""
        if not os.path.exists(out + ".pairs.gz") or os.path.getmtime(out + ".pairs.gz") < os.path.getmtime(input_file):
            # use standard pairtools to generate pairs file 
            print("Generate pairs ...")
            print(pair_command)
            subprocess.call(pair_command, shell = True)
        else:
            print("Pairs file already exists!")
            exit(0)
    if args.output_format == "cool":
        if not input_file.endswith(".pairs.gz"):
            raise ValueError("Input file is not in pairs format!")
        if not os.path.exists(out + ".cool") or os.path.getmtime(out + ".cool") < os.path.getmtime(input_file):
            print("Load pairs to generate cool ...")
            # https://cooler.readthedocs.io/en/latest/cli.html#cooler-cload-pairs
            # pairs do NOT to be sorted. Accept compressed file. 
            if len(args.resolution) == 1:
                cool_command = f"cooler cload pairs --assembly {assembly} {chrsize}:{args.resolution[0]} {input_file} {out}.cool -c1 2 -p1 3 -c2 4 -p2 5"
                print(cool_command)
                subprocess.call(cool_command, shell = True)
            else:
                min_resolution = min(args.resolution)
                print(f"Generate cool at {min_resolution} resolutions ...")
                cool_command = f"cooler cload pairs --assembly {assembly} {chrsize}:{min_resolution} {input_file} {out}.cool -c1 2 -p1 3 -c2 4 -p2 5"
                print(cool_command)
                if not os.path.exists(f"{out}.cool"):
                    subprocess.call(cool_command, shell = True)
                all_resolutions = ",".join([str(r) for r in args.resolution])
                print(f"Generate mcool at {all_resolutions} resolutions ...")
                mcool_command = f"cooler zoomify {out}.cool -n {n} -r {all_resolutions} -o {out}.mcool"
                print(mcool_command)
                print(f"Run cooler w/ {n} threads ...")
                subprocess.call(mcool_command, shell = True)
        else:
            print(f"{out}.cool already exists!")
    # generate hic
    if args.output_format == "hic":
        if not input_file.endswith(".pairs.gz"):
            raise ValueError("Input file is not in pairs format!")
        if not os.path.exists(out + ".hic") or os.path.getmtime(out + ".hic") < os.path.getmtime(input_file):
            # although cool no need to sort pairs, juicer must sort pairs by chromosome 
            print("Load pairs to generate HiC ... ")
            juicertools="/usr/local/apps/juicer/juicer-1.6/scripts/juicer_tools.jar"
            all_resolutions = ",".join([str(r) for r in args.resolution])
            print(f"Load pairs to generate HiC at {all_resolutions} resolutions ... ")
            juicer_command = f"java -Xmx48g -jar {juicertools} pre -r {all_resolutions} -k KR,ICE,VC {input_file} {out}.hic --threads {n} {assembly}"
            print(juicer_command)
            print(f"Run juicer w/ {n} threads ...")
            subprocess.call(juicer_command, shell = True)
    if args.output_format == "transform":
        if input_file.endswith("cool"):
            out_end = "hic"
        if input_file.endswith("hic"):
            out_end = "mcool"
        convert_command = f"hictk convert {input_file} {out}.{out_end}"
        if not os.path.exists(f"{out}.{out_end}") or os.path.getmtime(f"{out}.{out_end}") < os.path.getmtime(input_file):
            print(convert_command)
            subprocess.call(convert_command, shell = True)
        else:
            print(f"{out}.{out_end} already exists!")
    print(timeit(start))

if __name__ == "__main__":
    main()
