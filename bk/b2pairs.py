#!/usr/bin/env python3
"""
Transform bedpe or bam file to pairs format. 

Update history:
07 Oct 2024 - support bam input.

"""
import pandas as pd 
import bioframe as bf 
import argparse 
import os 
import numpy as np 
import subprocess 
import pysam  

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("b", help="input bam/bed file")
    parser.add_argument("-chrsize", "--chrsize", required = False, help = "file for ordered chromosome and size")
    parser.add_argument("-assembly", "--assembly", default = "hg38", help = "genome assembly name")
    parser.add_argument("-f", "--format", required = False, choices = ["BAM", "BED"], help = "input b file format (BAM/BED)")
    parser.add_argument("-d", "--distance", required = False, type = int, help = "PET distance used for filtering")
    parser.add_argument("-less", "--less_than", action = "store_true", help = "filter PET distance < d, if not assigned, PET distance > d")
    parser.add_argument("-map_quality", "--map_quality", type = int, required = False, help = "minimum mapping quality required")
    parser.add_argument("-o", "--out", required = False, help = "output name (prefix)")
    args=parser.parse_args()
    return args

def pairs_header(chr_df, assembly):
    chromsize_list = [f"#chromsize: {row.chrom} {row.size}" for row in chr_df.itertuples()]
    chromsize_total = "\n".join(chromsize_list)
    # note that empty line is not allowed for pairs file
    header = f"## pairs format v1.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n#genome_assembly: {assembly}\n{chromsize_total}\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
    return header 

def chunk_bed(bed, chunksize):
    bed_df = pd.read_table(bed, chunksize = chunksize, names = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"])
    for bed_chunk in bed_df:
        yield bed_chunk 

def write_pairs_from_bed(bed, chrom_order):
    # read bedpe file
    bed_df = bf.read_table(bed, schema = "bedpe")
    bed_df = bed_df[(bed_df["chrom1"].isin(chrom_order)) & (bed_df["chrom2"].isin(chrom_order))].copy()
    # find the 5' position of reads 
    bed_df["pos1"] = np.where(bed_df["strand1"] == "+", bed_df["start1"], bed_df["end1"])
    bed_df["pos2"] = np.where(bed_df["strand2"] == "+", bed_df["start2"], bed_df["end2"])
    # create container for filling elements 
    df_container = list()
    for row in bed_df.itertuples():
        if row.chrom1 == row.chrom2 and row.pos1 < row.pos2:
            df_container.append([row.name, row.chrom1, row.pos1, row.chrom2, row.pos2, row.strand1, row.strand2])
        if row.chrom1 == row.chrom2 and row.pos1 > row.pos2:
            df_container.append([row.name, row.chrom2, row.pos2, row.chrom1, row.pos1, row.strand2, row.strand1])
        if row.chrom1 != row.chrom2:
            i1 = chrom_order.index(row.chrom1)
            i2 = chrom_order.index(row.chrom2)
            if i1 < i2:
                df_container.append([row.name, row.chrom1, row.pos1, row.chrom2, row.pos2, row.strand1, row.strand2])
            if i1 > i2:
                df_container.append([row.name, row.chrom2, row.pos2, row.chrom1, row.pos1, row.strand2, row.strand1])
    df_pairs = pd.DataFrame(df_container, columns = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"])
    # give the chromorder for dataframe ordering
    sorterIndex = dict(zip(chrom_order, range(len(chrom_order))))
    df_pairs['chrom1_rank'] = df_pairs['chrom1'].map(sorterIndex)
    df_pairs["chrom2_rank"] = df_pairs["chrom2"].map(sorterIndex)
    df_pairs.sort_values(by = ["chrom1_rank", "chrom2_rank", "pos1", "pos2"], inplace = True, ignore_index = True)
    df_pairs.drop(["chrom1_rank", "chrom2_rank"], axis = 1, inplace = True)
    return df_pairs

def write_pairs_from_bam(bam, chrom_order, out, map_quality):
    # read bam file
    strand_dict = {True:"+", False:"-"}
    bam_handle = pysam.AlignmentFile(bam, "rb")
    fw = open(out + ".tmp", "w")
    # parse read by read on coordinate-sorted file
    if map_quality != None:
        for read in bam_handle.fetch():
            if read.mapping_quality > map_quality:
                rname = read.query_name
                # get the coordinates of read alignment 
                rseq = read.reference_name 
                mseq = read.next_reference_name 
                if rseq in chrom_order and mseq in chrom_order:
                    rstart = read.reference_start
                    mstart = read.next_reference_start
                    # check read is aligned before its mate
                    rindex = chrom_order.index(rseq)
                    mindex = chrom_order.index(mseq)
                    # correction: when index == 0 (in case of chr1), no cis read can be written. 0 cannot < 0
                    if rindex <= mindex and rstart*rindex < mstart*mindex:
                        fw.write(f"{rname}\t{rseq}\t{rstart}\t{mseq}\t{mstart}\t{strand_dict[read.is_forward]}\t{strand_dict[read.mate_is_forward]}\n")
    else:
        rname = read.query_name
        # get the coordinates of read alignment 
        rseq = read.reference_name 
        mseq = read.next_reference_name 
        if rseq in chrom_order and mseq in chrom_order:
            rstart = read.reference_start
            mstart = read.next_reference_start
            # check read is aligned before its mate
            rindex = chrom_order.index(rseq)
            mindex = chrom_order.index(mseq)
            if rindex <= mindex and rstart*rindex < mstart*mindex:
                fw.write(f"{rname}\t{rseq}\t{rstart}\t{mseq}\t{mstart}\t{strand_dict[read.is_forward]}\t{strand_dict[read.mate_is_forward]}\n")
    fw.close()

def compress_index(out):
    # compress pairs text file using bgzip 
    subprocess.call(f"bgzip -c {out}.pairs > {out}.pairs.gz", shell = True)
    # generate index for compressed pairs file  (index file with px2 suffix)
    subprocess.call(f"pairix -p pairs {out}.pairs.gz", shell = True)

def main():
    args = args_parser()
    b = args.b
    chrsize = args.chrsize 
    assembly = args.assembly 
    if chrsize != None:
        chr_df = pd.read_table(chrsize, sep = "\t", header = None, names = ["chrom", "size"])
    # 
    if args.format == None:
        formatt = b.split(".")[-1].upper()
    else:
        formatt = args.format.upper()
    # 
    if formatt == "BAM"
        args.format.upper() == "BAM":
        bam_handle = pysam.AlignmentFile(b, "rb")
        chr_list = [line.split("\t")[1:] for line in str(bam_handle.header).split("\n") if line.startswith("@SQ")]
        chr_df = pd.DataFrame(chr_list, columns = ["chrom", "size"])
        chr_df["chrom"] = chr_df["chrom"].str.split(":", expand = True)[1]
        chr_df["size"] = chr_df["size"].str.split(":", expand = True)[1]
    if chr_df not in locals():
        print("chrsize file is required. ")
        exit(1)
    chrom_order = chr_df["chrom"].tolist()
    header = pairs_header(chr_df, assembly)
    if args.out == None:
        out = os.path.basename(b).split(".")[0]
    else:
        out = args.out
    # 
    if formatt == "BED":
        df_pairs = write_pairs_from_bed(b, chrom_order)
        df_pairs.to_csv(out + ".tmp", sep = "\t", header = False, index = False)
    elif formatt == "BAM":
        write_pairs_from_bam(b, chrom_order, out, args.map_quality)
    else:
        print("format is not supported")
        exit(1)
    # write header (parsed from input chrsize file)
    with open(out + ".header", "w") as fw:
        fw.write(header)
    try: 
        subprocess.call(f"cat {out}.header {out}.tmp > {out}.pairs", shell = True)
        os.remove(f"{out}.header")
        os.remove(f"{out}.tmp")
    except: 
        print("header not added to pairs dataframe")
        exit(1)
    try: 
        compress_index(out)
    except:
        print("pairs file is not indexed successfully using pairix")
        pass 
        
if __name__ == "__main__":
    main()