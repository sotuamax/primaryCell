# 
import pandas as pd 
import bioframe as bf 
import argparse 
import os 
import numpy as np 
import subprocess 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("-bed", "--bed", help="bed file")
    parser.add_argument("-chrsize", "--chrsize", help = "file for ordered chromosome and size")
    parser.add_argument("-assembly", "--assembly", default = "hg38", help = "genome assembly name")
    parser.add_argument("-o", "--out", required = False, help = "output name")
    args=parser.parse_args()
    return args

def get_header(chrsize, assembly):
    chr_df = pd.read_table(chrsize, sep = "\t", header = None, names = ["chrom", "size"])
    chromsize_list = [f"#chromsize: {row.chrom} {row.size}" for row in chr_df.itertuples()]
    chromsize_total = "\n".join(chromsize_list)
    # note that empty line is not allowed for pairs file
    header = f"## pairs format v1.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n#genome_assembly: {assembly}\n{chromsize_total}\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
    return header 

def chunk_bed(bed, chunksize):
    bed_df = pd.read_table(bed, chunksize = chunksize, names = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"])
    for bed_chunk in bed_df:
        yield bed_chunk 

def write_pairs(bed, chrsize):
    chr_df = pd.read_table(chrsize, sep = "\t", header = None, names = ["chrom", "size"])
    chrom_order = chr_df["chrom"].tolist()
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

def compress_index(out):
    # compress pairs text file using bgzip 
    subprocess.call(f"bgzip -c {out}.pairs > {out}.pairs.gz", shell = True)
    # generate index for compressed pairs file  (index file with px2 suffix)
    subprocess.call(f"pairix -p pairs {out}.pairs.gz", shell = True)

def main():
    args = args_parser()
    bed = args.bed 
    chrsize = args.chrsize 
    assembly = args.assembly 
    header = get_header(chrsize, assembly)
    if args.out == None:
        out = os.path.basename(bed).split(".bed")[0]
    else:
        out = args.out
    df_pairs = write_pairs(bed, chrsize)
    df_pairs.to_csv(out + ".tmp", sep = "\t", header = False, index = False)
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
        print("pairs file is not compressed and indexed successfully using pairix")
        pass 
        
if __name__ == "__main__":
    main()