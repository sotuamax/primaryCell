#!/usr/bin/env python3
"""
This script is for collecting scores at specified bed surrounding region (flanking region). 
Input:
    - BED: regions to collect score; 
    - score: bigwig/bam format; 
    Note: for bigwig value, no normalization applied, thus, recommend input bigwig is already normalized; 
    # for bam, coverage values are calculated, and normalization will apply on the count value. (deleted this function)

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
import bioframe as bf 

ignore_warning()

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="bed_profile.py <bed> <bw> ")
    parser.add_argument("bed", help = "bed file (for intervals to examine)")
    parser.add_argument("bw", help = "provide scores on genomic position in bigwig format")
    parser.add_argument("-scale", choices = ["mm", "bg"], required = False, help = "methods used to scale score ")
    parser.add_argument("-flank", "--flank", default = 1000, type = int, help = "flanking side size centered at bed region pospoint (default: 1000)")
    #parser.add_argument("-filter", "--filter", action = "store_true", help = "when assigned, filter bed input for chromosome")
    parser.add_argument("-step", "--step", default = 50, type = int, help = "sample step size within the 2xflanking region")
    parser.add_argument("-min", "--min_coverage", default = 0, type = int, help = "minimum read for a region to be included in the average caculation")
    parser.add_argument("-plot", "--plot", action = "store_true", help = "when assigned, plot the screened region")
    # parser.add_argument("-p", "--bam_processor", default = 4, type = int, help = "number of processors to use for processing bam file")
    # parser.add_argument("-TSS", "--TSS", action = "store_true", help = "perform TSS enrichment call. ")
    parser.add_argument("-o", "--output", help = "output prefix")
    args=parser.parse_args()
    return args

def scale_mm(mtx):
    """"
    Scale score on the min and max value to 0-1 range
    """
    from sklearn.preprocessing import minmax_scale 
    matrix_scaled = minmax_scale(mtx, axis = 1)
    return matrix_scaled

def scale_bg(mtx):
    """"
    Scale score based on background signal (that calculated from all intervals side ends)

    For normalization idea, borrow from: 
    https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_tss_enrich.py 
    """
    #endflank = 100//step + 1
    score_mean = mtx.mean(axis = 0)
    bg_score = (score_mean[0] + score_mean[-1])/2
    matrix_scaled = mtx/bg_score
    return matrix_scaled

def main():
    # all ranks run for argument parsing 
    start = time.time()
    args = args_parser()
    bed = args.bed
    bw = args.bw
    output = args.output
    flank = args.flank
    step = args.step
    min_cnt = args.min_coverage
    #bw = args.bigwig 
    bed_df = bf.read_table(bed, schema = "bed4")
    refseq = bf.assembly_info("hg38")
    standard_seq = refseq.seqinfo["name"].tolist()
    bed_df = bed_df[(bed_df["chrom"].isin(standard_seq)) & (bed_df["chrom"] != "chrM") & (bed_df["chrom"] != "chrY")].copy()
    bed_df.to_csv(bed.replace(".bed", ".tmp.bed"), sep = "\t", header = False, index = False)
    import pyBigWig
    score_handle = pyBigWig.open(bw)
    if score_handle.isBigWig():
        #score_list = list()
        print("Collecting scores at bed surrounding region ...")
        score_list = [np.array([score_handle.stats(row.chrom, i, i+50)[0] for i in range((row.start+row.end)//2-flank, (row.start+row.end)//2+flank+1, step)]) for row in bed_df.itertuples()]
        # for i in range(-flank, flank, step):
        #     if i == 0:
        #         score_array = np.array([score_handle.stats(row.chrom, row.start, row.end) for row in bed_df.itertuples()]).astype(float)
        #     if i < 0: 
        #         score_array = np.array([score_handle.stats(row.chrom, row.start+i, row.start+i+step) for row in bed_df.itertuples()]).astype(float)
        #     if i > 0: 
        #         score_array = np.array([score_handle.stats(row.chrom, row.end+i, row.end+i+step) for row in bed_df.itertuples()]).astype(float)
            #score_list.append(score_array)
        # score_list contains for each list as a relative position of all bed intervals (e.g., all values at 1000 bp downstream of bed intervals)
        print("All positions collected!")
        score_stack = np.vstack(score_list) # stack in vertical format (relative position on each row)
        if score_stack.shape[0] != len(bed_df):
            raise ValueError("score dimension is not match to bed file")
        np.savetxt(output + ".mtx", score_stack, delimiter = "\t", fmt = "%d")
        # filter score, bed row score w/ all zeros removed
        score_filtered = score_stack[~np.all(score_stack <= min_cnt, axis = 1)]
        if args.scale == "mm":
            print("Perform minmax scale on signal ...")
            score_filtered = scale_mm(score_filtered)
        if args.scale == "bg":
            print("Perform background sclae on signal ...")
            score_filtered = scale_bg(score_filtered)
        score_mean = pd.DataFrame(score_filtered.mean(axis = 0), columns = ["signal"])
        score_mean["pos"] = range(-flank, flank+1, step)
        score_mean[['pos', "signal"]].to_csv(output + "_mean.txt", sep = "\t", header = True, index = False)
        bg_score = (score_filtered[:, 0:10].mean() + score_filtered[:, -10:-1].mean())/2
        center_score = score_filtered[:, score_filtered.shape[-1]//2].mean()
        enrich_value = center_score/bg_score
        print(f"Enrichment score: {enrich_value}")
        with open(output+".stat", "w") as fw:
            fw.write(f"enrich\t{score_filtered.shape[0]}\t{enrich_value}\n")
        #print('Write scores into csv!')
        # np.savetxt(output + ".csv", score_stack, delimiter = ",")
        # score_stack = score_stack[~np.isnan(score_stack).any(axis = 1), :]
        # get the average score at each relative position 
        # score_mean = np.mean(score_filtered, axis = 0)
        # score_normalized = encode_normalize(score_filtered, step)
        #pos_mean_encode = normalized_matrix.mean(axis = 0)
    # if mode == "bam":
    #     print("bam format is under development")
    #     exit(1)
    #     for i in range(-flank, flank+1, step):
    #         if i == 0:
    #             score_array = np.array([score_handle.count(row.chrom, row.start, row.end) for row in bed_df.itertuples()])
    #         if i < 0: 
    #             score_array = np.array([score_handle.stats(row.chrom, row.start+i, row.start+i+step, type = "max") for row in bed_df.itertuples()])
    #         if i > 0:
    #             score_array = np.array([score_handle.stats(row.chrom, row.end+i, row.end+i+step, type = "max") for row in bed_df.itertuples()])
    #         score_list.append(score_array)
    #     matrix_filtered = m_filter(bed_matrix, min_cnt)
    #     pd.DataFrame(matrix_filtered).to_csv(output+".mat", sep = "\t", header = False, index = False)
    #     scaled_matrix = scale_mm(matrix_filtered)
    # 
    #if args.TSS:
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
        import matplotlib.pyplot as plt 
        #fig, axis = plt.subplots(1,1)
        plt.plot(score_mean["pos"], score_mean["signal"])
        plt.scatter(score_mean["pos"], score_mean["signal"])
        # axis.axvline(x = 0, ymin = score_df["scaled_coverage_mean"].min(), ymax = score_df["scaled_coverage_mean"].max(), linestype = "--", size = 0.5, )
        # tick_break = int(200/step)
        # axis.set_xticks(np.array(score_df.index[::tick_break]))
        # axis.set_xlabel("distance from peak interval")
        # axis.set_ylabel("score")
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
