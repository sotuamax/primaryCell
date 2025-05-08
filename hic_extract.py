#!/usr/bin/env python3
import hicstraw
import argparse 
import numpy as np 
import os 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("hic", help="hic file. ")
    parser.add_argument("-b", "--binsize", type = int, default = 25000, help = "binsize for the resolution of hic (default: 25000)")
    parser.add_argument("-chr", "--chrom", type = str, help = "contact frequency of chromosome to retrieve from", required = False)
    args=parser.parse_args()
    return args

def chrsize(hic):
    chrlen = {chr.name:chr.length for chr in hic.getChromosomes()}
    del chrlen["M"]
    del chrlen["ALL"]
    del chrlen["Y"]
    return chrlen 

def count_retrieve(hic_handle, chr, chrlen, binsize):
    cs = chrlen[chr]
    matrix_chr = hic_handle.getMatrixZoomData(chr, chr, "observed", "KR", "BP", binsize)
    if cs/binsize < 5000:
        matrix_count = matrix_chr.getRecordsAsMatrix(0, cs, 0, cs)
        yield (0, cs, matrix_count)
    else: 
        rs = [r for r in range(0, cs, 5000*binsize)] + [cs]
        for i in range(len(rs)-1):
            matrix_count = matrix_chr.getRecordsAsMatrix(rs[i], rs[i+1]-1, rs[i], rs[i+1]-1)
            yield (rs[i], rs[i+1], matrix_count)

def main():
    args = args_parser()
    hic = args.hic 
    binsize = args.binsize 
    chr = args.chrom
    # 
    hic_handle = hicstraw.HiCFile(hic)
    chrlen = chrsize(hic_handle)
    resolution = hic_handle.getResolutions() 
    # 
    chrmax = max(chrlen.values())
    if binsize in resolution: 
        pass 
    else: 
        print("the binsize (resolution) is not available in current hic file")
        exit(1)
    # 
    try:
        hic_path = hic.split(".hic")[0]
        os.mkdir(hic_path)
    except: 
        pass
    if chr != None:
        for start,end,cnt in count_retrieve(hic_handle, chr, chrlen, binsize):
            # chrom, start, binsize, cnt_shape
            filename = f"chr{chr}_pos{int(start/1000)}k_res{int(binsize/1000)}k_bin{cnt.shape[0]}.txt"
            np.savetxt(os.path.join(hic_path, filename), cnt, delimiter = "\t")
    else:
        for chr in chrlen:
            for start,end,cnt in count_retrieve(hic_handle, chr, chrlen, binsize):
                # chrom, start, binsize, cnt_shape
                filename = f"chr{chr}_pos{int(start/1000)}k_res{int(binsize/1000)}k_bin{cnt.shape[0]}.txt"
                np.savetxt(os.path.join(hic_path, filename), cnt, delimiter = "\t")

if __name__ == "__main__":
    main()
