#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program is to clean sequencing data based on the primer given (R1, enriched), 
and assign cell identify based on cell barcode given in R2.
Output file is clean and cell barcode labeled R1 fastq file 
1). start with primer; 
2). valid cell barcode. 

g_target.py align -ref $hg38_star -gtf $gtf -sample GC760804GFP -seq_type cDNA -outdir output/GC760804 -n 24
g_target.py align -ref $hg38 -sample GC760804 -seq_type DNA -outdir output/GC760804 -n 24


"""

import os 
import sys
import pandas as pd
import numpy as np
import argparse
import pysam 
import subprocess
import matplotlib.pyplot as plt
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def args_parser():
    '''parser the argument from terminal command'''
    parser = argparse.ArgumentParser(prog = "PROG", add_help = True, formatter_class = argparse.RawDescriptionHelpFormatter)
    # important parameters 
    parser.add_argument("-sample", "--sample", help="sample prefix used to find fastq file")
    parser.add_argument("-outdir", "--outdir", default = ".", help="output directory name (also for sample fastq location)")
    parser.add_argument("-n", "--threads", help = "number of threads for alignment", default = 1, type = int)
    parser.add_argument("-seq_type", "--seq_type", help = "sequencing data type", choices = ["cDNA", "DNA"])
    parser.add_argument("-ref", "--reference", help = "reference genome index used for alignment")
    parser.add_argument("-qc", action = "store_true", help = "perform QC on BAM alignment file")
    parser.add_argument("-report", action = "store_true", help = "generate report for read count on cell barcode scale")
    parser.add_argument("-f", "--force", action = "store_true", help = "force alignment even when bam file exists")
    parser.add_argument("-gtf", "--gtf", required = False, help = "reference gtf file used for alignment")
    args=parser.parse_args()
    return args

def main():
    args = args_parser()
    sample = args.sample
    output_dir = args.outdir
    bam_file = os.path.join(output_dir, sample + ".bam")
    n= args.threads
    try:
        os.mkdir(output_dir)
        print(f"Create {output_dir}")
    except:
        pass
    align_r1 = os.path.join(output_dir, sample + "_R1.fastq.gz")
    if not os.path.exists(align_r1):
        raise IOError("Cannot find R1 for alignment.")
    print(f"Locate {sample} file ...")
    ### align for DNA seq
    if args.seq_type == "DNA":
        from utilities.align_tools import BWA_SE
        if not os.path.exists(bam_file) or args.force:
            print("Perform alignment ...")
            BWA_SE(align_r1, args.reference, os.path.join(output_dir, args.sample), n)
            flagstat = bam_file.replace(".bam", ".flagstat")
            subprocess.call(f"samtools flagstat {bam_file} -@ {n} > {flagstat}", shell = True)
        if not os.path.exists(bam_file.replace(".bam", ".dedup.bam")):
            # enable single-cell deduplicate process
            from utilities.bam_tools import sc_markdup
            print("Deduplicate on single cell level ...")
            print(f"Process {bam_file} ... ")
            # markdup(bam_file, os.path.join(output_dir, args.sample) + ".dedup")
            # make sure that single cell deduplication is on cell barcode and coordinates
            sc_markdup(bam_file, bam_file.replace(".bam", ".dedup.bam"), n)
        if args.qc:
            # check soft-clip side at 5' end for QC 
            from utilities.bam_tools import bamfilter 
            if os.path.exists(bam_file.replace(".bam", ".dedup.bam")):
                bam_file = bam_file.replace(".bam", ".dedup.bam")
            qc_bam = bam_file.replace(".bam", ".qc.bam")
            if not os.path.exists(qc_bam):
                print("Perform QC ...")
                print(f"Process {bam_file} ...")
                bamfilter(bam_file, qc_bam, clip_check=True, threads = n)
            from utilities.bam_tools import bam2bw
            bw = qc_bam.replace(".bam", ".bw")
            if not os.path.exists(bw):
                bam2bw(qc_bam, bw, n)
        if args.report:
            bam_file = os.path.join(output_dir, sample + ".dedup.qc.bam")
            bam_handle = pysam.AlignmentFile(bam_file, "rb", threads = n)
            if not os.path.exists(os.path.join(output_dir, sample + "_barcode.stat")):
                print(f"Generate report for cell barcode on {sample}...")
                print(f"Process {os.path.basename(bam_file)} ... ")
                read_barcode_list = [(read.query_name).split(":")[-1] for read in bam_handle.fetch()]
                barcode, barcode_freq = np.unique(read_barcode_list, return_counts = True)
                barcode_freq_df = pd.DataFrame({"barcode":barcode, "freq":barcode_freq})
                barcode_freq_df.to_csv(os.path.join(output_dir, sample + "_barcode.stat"), sep = "\t", header = True, index = False)
            else:
                barcode_freq_df = pd.read_table(os.path.join(output_dir, sample + "_barcode.stat"), sep = "\t", header = 0)
            barcode_freq_df.sort_values(by = "freq", ascending=False, inplace = True, ignore_index = True)
            # make plots for ChiC read per cell
            fig, axes = plt.subplots(1, 2, figsize = (10,4))
            axes[0].plot(range(len(barcode_freq_df)), barcode_freq_df["freq"], "o", markersize = 0.2)
            axes[0].set_xscale('log')
            axes[0].set_yscale('log')
            axes[0].set_xlabel("cell barcode ranked by read count")
            axes[0].set_ylabel("read count per cell barcode")
            axes[1].hist(np.log10(barcode_freq_df["freq"]), bins = 50, orientation='horizontal')
            axes[1].set_ylabel("log(read count per cell barcode)")
            axes[1].set_xlabel("freq")
            plt.savefig(os.path.join(output_dir, sample + "_barcode.png"), dpi = 300)
            # combine with GFP (when GFP data available)
            # ChiC read per cell for cells with valid GFP
            if os.path.exists(os.path.join(output_dir, sample + "GFP.xlsx")):
                print("Generate report for cell barcode w/ GFP filtering ...")
                tag_df = pd.read_excel(os.path.join(output_dir, sample + "GFP.xlsx"), sheet_name  = "flag_gene")
                if "library" not in tag_df.columns:
                    library_df = pd.read_excel(os.path.join(output_dir, sample + "GFP.xlsx"), sheet_name  = "summary")
                    tag_df["index"] = tag_df["barcode"].str.slice(0,15)
                    tag_df = pd.merge(library_df[["library", "index"]].drop_duplicates(keep = "first"), tag_df, on = "index", how = "right")
                tag_df.columns = [c.strip("_") for c in tag_df.columns]
                # ChiC_barcode_freq = pd.read_table(os.path.join(output_dir, sample + "_barcode.stat"), sep = "\t", header = 0)
                # ChiC_barcode_freq.sort_values(by = ['freq'], ascending=False, inplace = True, ignore_index = True)
                gfp_tag_copy = tag_df[["library", "index", "barcode", "gene", "symbol", "gene_frame_read"]].copy()
                gfp_tag_copy["i7"] = gfp_tag_copy['barcode'].str.split("+", expand = True)[0]
                barcode_freq_df["i7"] = barcode_freq_df["barcode"].str.split("+", expand = True)[0]
                gfp_tag_copy["b"] = gfp_tag_copy["barcode"].str.split("+", expand = True)[2]
                barcode_freq_df["b"] = barcode_freq_df["barcode"].str.split("+", expand = True)[2]
                GFP_ChiC_df = pd.merge(gfp_tag_copy.rename(columns = {"barcode":"GFP_barcode"}), barcode_freq_df.rename(columns = {"barcode":"ChiC_barcode"}), on = ["b", "i7"], how = "inner")
                GFP_ChiC_df.drop(["i7", "b"], axis = 1, inplace = True)
                GFP_ChiC_df.rename(columns = {"gene_frame_read":"GFP_read_count", "freq":"ChiC_read_count"}, inplace = True)
                GFP_ChiC_df.sort_values(["ChiC_read_count"], ascending=[False], inplace = True, ignore_index=True)
                GFP_ChiC_df.to_csv(os.path.join(output_dir, sample + "_clean_barcode.stat"), sep = "\t", header = True, index = False)
                # enrich 
                # enrich_df = pd.read_table(os.path.join(output_dir, sample + ".enrich.txt"), sep = "\t", header = 0)
                # enrich_sum_df = pd.merge(enrich_df.groupby("barcode")["total"].sum().reset_index(), enrich_df.groupby("barcode")["overlap_peak"].sum().reset_index(), on ="barcode")
                # chrM_count = enrich_df.query("chrom == 'chrM'")[["barcode", "total"]].rename(columns = {"total":"chrM"})
                # enrich_sum_df = pd.merge(enrich_sum_df, chrM_count, on = "barcode")
                # enrich_sum_df.rename(columns = {"barcode":"ChiC_barcode"}, inplace = True)
                # enrich_sum_df["enrichment"] = enrich_sum_df["overlap_peak"]/enrich_sum_df["total"]
                # enrich_sum_df["chrM_ratio"] = enrich_sum_df["chrM"]/enrich_sum_df["total"]
                # 
                # GFP_ChiC_df = pd.merge(GFP_ChiC_df, enrich_sum_df, on = "ChiC_barcode")
                GFP_ChiC_df["ChiC_read_count_log"] = np.log10(GFP_ChiC_df["ChiC_read_count"])
                GFP_ChiC_df["GFP_read_count_log"] = np.log10(GFP_ChiC_df["GFP_read_count"])
                # GFP_ChiC = GFP_ChiC_df[["ChiC_read_count_log", "GFP_read_count_log", "enrichment", "chrM_ratio"]]
                # import seaborn as sns
                # pair = sns.pairplot(GFP_ChiC, plot_kws={'s': 5}, markers = ".")
                # pair.figure.savefig(os.path.join(output_dir, sample + "pair.png"), dpi=300, bbox_inches='tight')
                ### using GM (Gaussian Mixture model to separate cells and find the real single cell)
                 # plot the relationship between GFP and ChiC read count 
                from sklearn.mixture import GaussianMixture
                gmm = GaussianMixture(n_components=2, random_state=0)
                data = np.array(GFP_ChiC_df["ChiC_read_count_log"]).reshape(-1,1)
                gmm.fit(data)
                # gmm.predict(data)
                # i = np.argmax(gmm.means_) # get the index for max-mean 
                GFP_ChiC_df["label"] = -1
                gmm_parameters = sorted([(gmm.means_[i][0],np.sqrt(gmm.covariances_[i][0][0])) for i in range(2)])
                # default label == -1, group with smaller mean label == 0, group with larger mean label == 1 (1 is what we can use for further analysis)
                for index, (m,sd) in enumerate(gmm_parameters):
                    min_i = m-sd*1.645; max_i = m+sd*1.645
                    GFP_ChiC_df["label"] = np.where((GFP_ChiC_df["ChiC_read_count_log"] > min_i) & (GFP_ChiC_df["ChiC_read_count_log"] < max_i), index, GFP_ChiC_df["label"])
                
                # GFP_ChiC_df["label"] = np.where(GFP_ChiC_df["chrM_ratio"] < 0.2, GFP_ChiC_df["label"], -1)
                ### make plots of the selected cell
                fig, axes = plt.subplots(1, 2, figsize = (10, 4))
                color_map = {0: 'red', 1:"red", -1: 'grey'}
                colors = GFP_ChiC_df['label'].map(color_map)
                axes[0].scatter(range(len(GFP_ChiC_df)), GFP_ChiC_df["ChiC_read_count"], c=colors, s = 0.2)
                # .plot(, "o", markersize = 0.2, )
                axes[0].set_xscale('log')
                axes[0].set_yscale('log')
                axes[0].set_xlabel("cell ranked by read count")
                axes[0].set_ylabel("ChiC read count per cell")
                axes[1].hist(GFP_ChiC_df["ChiC_read_count_log"], bins = 50, orientation='horizontal')
                for i in range(2):
                    m = gmm.means_[i][0]
                    sd = np.sqrt(gmm.covariances_[i][0][0])
                    min_i = m-sd*1.645; max_i = m+sd*1.645
                    axes[1].axhline(y=min_i, color='black', linestyle='--', linewidth=1)
                    axes[1].axhline(y=max_i, color='black', linestyle='--', linewidth=1)
                axes[1].set_ylabel("log10(ChiC read count per cell)")
                axes[1].set_xlabel("freq")
                # axes[2].scatter(GFP_ChiC_df["GFP_read_count"], GFP_ChiC_df["ChiC_read_count"], c = colors, s = GFP_ChiC_df["enrichment"]*10)
                # axes[2].set_xscale("log")
                # axes[2].set_yscale("log")
                # axes[2].set_xlabel("GFP read count per cell")
                # axes[2].set_ylabel("ChiC read count per cell")
                # # 
                # axes[4].scatter(GFP_ChiC_df["enrichment"], GFP_ChiC_df["ChiC_read_count"], c = colors, s = 0.2)
                # axes[4].set_xlabel("Enrichment")
                # axes[4].set_yscale("log")
                # axes[4].set_ylabel("ChiC read count per cell")
                # # 
                # axes[5].scatter(GFP_ChiC_df["chrM_ratio"], GFP_ChiC_df["ChiC_read_count"], c = colors, s = 0.2)
                # axes[5].set_xlabel("chrM ratio")
                # axes[5].set_yscale("log")
                # axes[5].set_ylabel("ChiC read count per cell")
                # 
                plt.savefig(os.path.join(output_dir, sample + "_clean_barcode.png"), dpi = 300)
                GFP_ChiC_df.drop(["index", "ChiC_read_count_log", "GFP_read_count_log"], axis = 1).to_csv(os.path.join(output_dir, sample + "_clean_barcode.stat"), sep = "\t", header = True, index = False)
    ### align for RNA seq 
    if args.seq_type == "cDNA":
        from utilities.align_tools import STAR_SE
        #from utilities.parse_log import STAR_log
        gtf = args.gtf
        if gtf.endswith(".gz"):
            gtf = gtf.split(".gz")[0]
        ### 
        if not os.path.exists(bam_file) or args.force:
            print("Perform alignment ...")
            STAR_SE(align_r1, args.reference, gtf, os.path.join(output_dir, sample), n = n, enable_novel = False)
        from utilities.gtf_tools import parse_gtf
        gtf_gz = gtf + ".gz"
        print("Parse CDS region ...")
        feature_df = parse_gtf(gtf_gz, "CDS").drop("transcript", axis = 1).drop_duplicates(keep = "first", ignore_index=True)
        # annotated frame is for 5' end 
        # update frame on 3' end (where donor seq to add)
        feature_df["frame"] = (feature_df["end"] - feature_df["frame"] - feature_df["start"])%3
        #
        from utilities.bam_tools import bam2bed
        import bioframe as bf 
        #filtered_list = list()
        print("Join alignment and CDS ...")
        #tmp_bam = bam_file.replace(".bam", ".tmp.bam")
        tmp_bed = bam_file.replace(".bam", ".bed")
        # delete all previous writing if exists
        if os.path.exists(tmp_bed):
            os.remove(tmp_bed)
        # write chrom by chrom
        i = 0
        for c,feature_c in feature_df.groupby("chrom"):
            c_read = bam2bed(bam_file, chrom = c)
            if c_read is not None:
                read_feature_overlap = bf.overlap(c_read, feature_c, how = "inner")
                # read_feature_overlap.drop_duplicates(keep = "first", inplace = True, ignore_index = True)
                # read_feature_overlap = read_feature_overlap[~read_feature_overlap["cigar"].str.contains("N")]
                ### orientation & end overlap
                scenario1 = (read_feature_overlap["strand"] == "+") & (read_feature_overlap['start'] == read_feature_overlap["start_"]) & (read_feature_overlap["strand_"] == "-")
                scenario2 = (read_feature_overlap["strand"] == "-") & (read_feature_overlap['end'] == read_feature_overlap["end_"]) & (read_feature_overlap["strand_"] == "+")
                read_feature_overlap_filtered = read_feature_overlap[scenario1 | scenario2]
                ### per read match to single frame 
                single_frame = read_feature_overlap_filtered.groupby("name")["frame_"].transform("nunique") == 1
                read_feature_overlap_filtered = read_feature_overlap_filtered[single_frame]
                ### per read (labeled by its name) match to single target
                single_target = read_feature_overlap_filtered.groupby("name")["gene_"].transform("nunique") == 1
                read_feature_overlap_filtered = read_feature_overlap_filtered[single_target].copy()
                read_feature_overlap_filtered.drop(["score", "chrom_", "start_", "end_", "strand_"], axis = 1, inplace = True)
                read_feature_overlap_filtered.drop_duplicates(keep = "first", ignore_index = True, inplace = True)
                if i == 0:
                    read_feature_overlap_filtered.to_csv(tmp_bed, mode = "a", sep = "\t", header = True, index = False)
                else:
                    read_feature_overlap_filtered.to_csv(tmp_bed, mode = "a", sep = "\t", header = False, index = False)
                i += 1
                #filtered_list.append(read_feature_overlap_filtered)
        ### write qualified reads 
        # final_filtered = pd.concat(filtered_list, axis = 0)
        # final_filtered.to_csv(tmp_bed, sep = "\t", header = True, index = False)
        # print("Collect read tag info ...")
        ### read target gene 
        #read_gene_dict = final_filtered.set_index("name")[["gene_", "symbol_"]].to_dict()
        # read_gene = read_gene_dict["gene_"]
        # read_symbol = read_gene_dict["symbol_"]
        # ### read target gene frame 
        # read_frame = final_filtered.set_index("name")[["frame_"]].to_dict()["frame_"]
        # ### read cluster 
        # print("Generate genomic block ...")
        # d = 100000
        # read_cluster = bf.cluster(final_filtered[["chrom", "start", "end", "name"]].drop_duplicates(keep = "first"), min_dist = d)[["chrom", "cluster_start", "cluster_end", "name"]].drop_duplicates(keep = "first").rename(columns = {"cluster_start":"start", "cluster_end":"end"})
        # while len(set(read_cluster["name"])) != len(read_cluster):
        #     d = d*2
        #     read_cluster = bf.cluster(final_filtered[["chrom", "start", "end", "name"]].drop_duplicates(keep = "first"), min_dist = d)[["chrom", "cluster_start", "cluster_end", "name"]].drop_duplicates(keep = "first").rename(columns = {"cluster_start":"start", "cluster_end":"end"})
        ### 
        # seqnum = len(set(read_cluster["name"]))
        # logging.info(f"Reads flag CDS\t{seqnum}")
        # fm,fm_cn = np.unique(list(read_frame.values()), return_counts = True)
        # for f in enumerate(fm):
        #     ff = fm_cn[f[0]]
        #     logging.info(f"Frame {f[-1]}\t{ff}")
        #     logging.info(f"Frame {f[-1]}%\t{round(ff/seqnum*100, 1)}")
        #### 
        # print("Write read tag info ...")
        # bam_handle = pysam.AlignmentFile(bam_file, "rb", threads = args.threads//2)
        # with pysam.AlignmentFile(tmp_bam, "wb", template = bam_handle, threads = args.threads//2) as newbam:
        #     for (c,cstart,cend),c_df in read_cluster.groupby(["chrom", "start", "end"]):
        #         for read in bam_handle.fetch(c, cstart, cend):
        #             if read.query_name in c_df["name"].tolist():
        #                 read.set_tag("TG", read_gene[read.query_name], value_type = "Z") # target gene 
        #                 read.set_tag("SY", read_symbol[read.query_name], value_type = "Z") # gene symbol
        #                 read.set_tag("FR", read_frame[read.query_name], value_type = "i") # frame 
        #                 newbam.write(read)
        # new_bam = bam_file.replace(".bam", ".clean.bam")
        # subprocess.call(f"samtools sort {tmp_bam} -o {new_bam} -@ {args.threads} && samtools index -@ {args.threads} {new_bam}", shell = True)
        # try: 
        #     os.remove(tmp_bam)
        # except:
        #     pass
        # #### 
        # print("Get flag genes ...")
        # new_bam_handle = pysam.AlignmentFile(new_bam, "rb", threads = args.threads)
        # read_target_df = pd.DataFrame([[read.query_name, read.get_tag("TG"), read.get_tag("SY")] for read in new_bam_handle.fetch() if read.get_tag("FR") == args.frame], columns = ["name", "gene", "symbol"])
        # read_target_df["library"] = read_target_df["name"].str.split(":", expand = True).iloc[:, -1]
        # read_target_df.to_csv(os.path.join(output_dir, sample+".tmp"), sep = "\t", header = True, index = False)
        # #### 
        # read_target_df_count = pd.merge(read_target_df.groupby(["library", "gene", "symbol"]).size().reset_index(name = "count"), read_target_df.groupby("library").size().reset_index(name = "total"), on = "library")
        # read_target_df_count["ratio"] = round(read_target_df_count["count"]/read_target_df_count["total"]*100)
        # read_target_df_count.to_csv(os.path.join(output_dir, sample+".txt"), sep = "\t", header = True, index = False)
        # read_target_df_count.query("count >= 3 and ratio >= 80").to_csv(os.path.join(output_dir, sample+".clean.txt"), sep = "\t", header = True, index = False)
        

if __name__ == "__main__":
    main()
