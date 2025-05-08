#!/usr/bin/env python3
"""
This script is for collecting scores at specified bed surrounding region (flanking region). 
Input:
    - BED: regions to collect score; 
    - score: bigwig/bam format; 
    Note: for bigwig value, no normalization applied, thus, recommend input bigwig is already normalized; for bam, coverage values are calculated, and normalization will apply on the count value. 

Output:
    - score matrix (each column match to )
"""
import pandas as pd 
import numpy as np 
import argparse 
import os 
from utilities.misc import timeit
import time 
from utilities.misc import ignore_warning

ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("bed", help = "bed file (for intervals to examine)")
    parser.add_argument("score", help = "provide scores on genomic position \n bigwig, or bam to caculate coverage from")
    parser.add_argument("-scale", action = "store_true", help = "scale bed file to use midpoint as reference")
    parser.add_argument("-flank", "--flank", default = 1000, type = int, help = "flanking side size centered at bed region pospoint (default: 1000)")
    parser.add_argument("-filter", "--filter", action = "store_true", help = "when assigned, filter bed input for chromosome")
    parser.add_argument("-step", "--step", default = 50, type = int, help = "sample step size within the 2xflanking region")
    parser.add_argument("-min", "--min_coverage", default = 0, type = int, help = "minimum read for a region to be included in the average caculation")
    parser.add_argument("-plot", "--plot", action = "store_true", help = "when assigned, plot the screened region")
    parser.add_argument("-p", "--bam_processor", default = 4, type = int, help = "number of processors to use for processing bam file")
    parser.add_argument("-TSS", "--TSS", action = "store_true", help = "perform TSS enrichment call. ")
    parser.add_argument("-o", "--output", help = "prefix output name")
    args=parser.parse_args()
    return args

def m_filter(matrix, min_cnt):
    m_filtered = matrix[~np.all(matrix <= min_cnt, axis = 1)]
    return m_filtered 

def scale_mm(matrix_filtered):
    from sklearn.preprocessing import minmax_scale 
    matrix_scaled = minmax_scale(matrix_filtered, axis = 1)
    return matrix_scaled

def encode_normalize(matrix_filtered, step):
    """"
    For normalization idea, borrow from: 
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_tss_enrich.py 
    """
    endflank = 100//step + 1
    step_mean = matrix_filtered.mean(axis = 0)
    avg_noise = (step_mean[:endflank].sum() + step_mean[-endflank:].sum())/(2*endflank)
    matrix_normalized = matrix_filtered/avg_noise
    return matrix_normalized

def main():
    # all ranks run for argument parsing 
    start = time.time()
    args = args_parser()
    bed = args.bed
    score = args.score
    output = args.output
    flank = args.flank
    step = args.step 
    min_cnt = args.min_coverage
    #bw = args.bigwig 
    if score.endswith("bigwig") or score.endswith(".bw"):
        import pyBigWig
        score_handle = pyBigWig.open(score)
        if score_handle.isBigWig():
            print("Input score is in bigwig format")
            mode = "bw"
    if score.endswith(".bam"):
        import pysam 
        score_handle = pysam.AlignmentFile(score, "rb", threads = args.bam_processor)
        print("Input score is in BAM format")
        mode = "bam" 
    import bioframe as bf 
    bed_df = bf.read_table(bed, schema = "bed4")
    if args.scale: 
        bed_df = bf.expand(bed_df, scale = 0)
        bed_df["end"] = bed_df["end"] + 1
    score_list = list()
    if mode == "bw":
        print("Collecting scores at bed surrounding region ...")
        for i in range(-flank, flank, step):
            if i == 0:
                score_array = np.array([score_handle.stats(row.chrom, row.start, row.end) for row in bed_df.itertuples()]).astype(float)
            if i < 0: 
                score_array = np.array([score_handle.stats(row.chrom, row.start+i, row.start+i+step) for row in bed_df.itertuples()]).astype(float)
            if i > 0: 
                score_array = np.array([score_handle.stats(row.chrom, row.end+i, row.end+i+step) for row in bed_df.itertuples()]).astype(float)
            score_list.append(score_array)
        # score_list contains for each list as a relative position of all bed intervals (e.g., all values at 1000 bp downstream of bed intervals)
        print("All positions collected!")
        score_stack = np.hstack(score_list)
        # stack in vertical format (relative position on each row)
        print('Write scores into csv!')
        np.savetxt(output + ".csv", score_stack, delimiter = ",")
        score_stack = score_stack[~np.isnan(score_stack).any(axis = 1), :]
        # get the average score at each relative position 
        score_average = np.mean(score_stack, axis = 0).T
    if mode == "bam":
        print("bam format is under development")
        exit(1)
        for i in range(-flank, flank+1, step):
            if i == 0:
                score_array = np.array([score_handle.count(row.chrom, row.start, row.end) for row in bed_df.itertuples()])
            if i < 0: 
                score_array = np.array([score_handle.stats(row.chrom, row.start+i, row.start+i+step, type = "max") for row in bed_df.itertuples()])
            if i > 0:
                score_array = np.array([score_handle.stats(row.chrom, row.end+i, row.end+i+step, type = "max") for row in bed_df.itertuples()])
            score_list.append(score_array)
        matrix_filtered = m_filter(bed_matrix, min_cnt)
        pd.DataFrame(matrix_filtered).to_csv(output+".mat", sep = "\t", header = False, index = False)
        scaled_matrix = scale_mm(matrix_filtered)
    # 
    if args.TSS:
        normalized_matrix = encode_normalize(matrix_filtered, step)
        pos_mean_encode = normalized_matrix.mean(axis = 0)
    # non_scale_matrix = z_scale(bed_matrix, 3)
    # pos_mean_raw = matrix_filtered.mean(axis = 0)
    # pos_mean_mm = scaled_matrix.mean(axis = 0)
    # #pos_median_no_scale = np.median(non_scale_matrix, axis = 0)
    # score_df = pd.DataFrame()
    # score_df["pos"] = screen_relative_pos
    # score_df["minmax_mean"] = pos_mean_mm
    # if args.TSS:
    #     score_df["norm_mean"] = pos_mean_encode
    # score_df["mean"] = pos_mean_raw
    # score_df.to_csv(f"{output}.txt", sep = "\t", index = False, header = True)
    # if args.TSS:
    #     sample_score = score_df.query("pos == 0").copy()
    #     sample_score["sample"] = os.path.basename(bam).split(".bam")[0]
    #     if not os.path.exists("TSS_enrich.txt"):
    #         sample_score.to_csv("TSS_enrich.txt", sep = "\t", header = True, index = False, mode = "w")
    #     else:
    #         sample_score.to_csv("TSS_enrich.txt", sep = "\t", header = False, index = False, mode = "a")
    if args.plot:
        print("Plot scores!")
        score_df = pd.DataFrame(score_average, index = range(-flank, flank, step), columns = ["score"])
        import matplotlib.pyplot as plt 
        fig, axis = plt.subplots(1,1)
        axis.plot(score_df.index, score_df["score"])
        # axis.axvline(x = 0, ymin = score_df["scaled_coverage_mean"].min(), ymax = score_df["scaled_coverage_mean"].max(), linestype = "--", size = 0.5, )
        tick_break = int(200/step)
        axis.set_xticks(np.array(score_df.index[::tick_break]))
        axis.set_xlabel("distance from peak interval")
        axis.set_ylabel("score")
        plt.savefig(f"{output}.pdf", format = "pdf")
        plt.savefig(f"{output}.png", dpi = 1000)
        ### 
# fig, (ax_line, ax_heatmap) = plt.subplots(2, 1, figsize=(10, 6),sharex=True, gridspec_kw={'height_ratios': [1, 5]})

# indices = np.random.choice(score_stack.shape[0], size=10000, replace=False)
# # Extract those rows
# random_rows = score_stack[indices]
# random_rows = random_rows[np.argsort(random_rows[:, random_rows.shape[1]//2])[::1]]
# # --- Line Plot (for selected row) ---
# ax_line.plot(np.arrange(score_stack.shape[1]), score_stack.mean(axis = 0), color='red')
# ax_line.set_ylabel("average score")
# # --- Heatmap ---
# im = ax_heatmap.imshow(random_rows, aspect='auto', cmap='viridis', extent = [0, random_rows.shape[1], 0, random_rows.shape[0]])
# ax_heatmap.set_ylabel("score")
# #ax_heatmap.set_yticks(range(len(row_labels)))
# #ax_heatmap.set_yticklabels(row_labels)
# #ax_heatmap.set_title("Heatmap + Line Plot")
# plt.colorbar(im, ax=ax_heatmap, orientation='vertical')
# plt.savefig(f"{output}.png", dpi = 1000)
    # finish and report process time
    print(timeit(start))

if __name__ == "__main__":
    main()
