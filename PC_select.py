#!/usr/bin/env python3
import pandas as pd 
import argparse 
import glob 
import bioframe as bf 
import os 
import numpy as np 
from itertools import product
from sklearn.preprocessing import StandardScaler

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, usage="")
    parser.add_argument("-dir", "--directory", help = "dchic output directory") 
    parser.add_argument("-exp", "--expression", help = "expression data", required = False)
    parser.add_argument("-dhs", "--dhs", help = "DHS signal file in bedGraph format", required = False)
    parser.add_argument("-tss", "--tss", help = "TSS file in bedGraph format", default = "/data/jim4/Seq/primary_cell_project/analysis/Compartment/data/hg38_goldenpathData/hg38_TSScount_40k.bedGraph")
    parser.add_argument("-gc", "--gc", help = "GC ratio in bedGraph format", default = "/data/jim4/Seq/primary_cell_project/analysis/Compartment/data/hg38_goldenpathData/hg38_GCratio_40k.bedGraph")
    parser.add_argument("-blacklist", "--blacklist", default = "/data/jim4/data/blacklist/hg38_blacklist_jointed.bed", help = "blacklist bed file")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args


def value_scale(value):
    value = np.array(value)
    scaled_value = StandardScaler().fit_transform(value.reshape(-1,1))
    return scaled_value.reshape(-1)

def main():
    args = args_parser()
    dchic_folder = args.directory
    tss = args.tss 
    gc = args.gc 
    dhs = args.dhs 
    exp = args.expression 
    output = args.output 
    print("Read PC values ...")
    all_pc_df = pd.concat([bf.read_table(file, sep = "\t", header = 0).drop("index", axis = 1) for file in glob.glob(os.path.join(dchic_folder, "*pc.txt"))], axis = 0)
    all_pc_df.rename(columns = {"chr":"chrom"}, inplace = True)
    print(f"PC bins: {len(all_pc_df)}")
    print("Remove blacklist ...")
    from utilities.bed_tools import rmblacklist
    blacklist_df = bf.read_table(args.blacklist, schema = "bed3")
    all_pc_df = rmblacklist(all_pc_df, blacklist_df)
    print(f"PC bins: {len(all_pc_df)}")
    # read TSS / GC 
    tss_bg = bf.read_table(tss, schema = "bedGraph"); tss_bg.rename(columns = {"value":"TSS"}, inplace = True)
    gc_bg = bf.read_table(gc, schema = "bedGraph"); gc_bg.rename(columns = {"value":"GC"}, inplace = True)
    #cols = ["chrom", "PC", "GC_cor", "TSS_cor"]
    print("Add GC ratio ...")
    all_pc_df = pd.merge(all_pc_df, gc_bg, on = ["chrom", "start", "end"], how = "inner")
    print(f"PC bins: {len(all_pc_df)}")
    print("Add TSS density ...")
    all_pc_df = pd.merge(all_pc_df, tss_bg, on = ["chrom", "start", "end"], how = "inner")
    print(f"PC bins: {len(all_pc_df)}")
    if dhs is not None:
        dhs_bg = bf.read_table(dhs, schema = "bedGraph"); dhs_bg.rename(columns = {"value":"DHS"}, inplace = True)
        #cols.append("DHS_cor")
        print("Add DHS score ...")
        all_pc_df = pd.merge(all_pc_df, dhs_bg, on = ["chrom", "start", "end"], how = "inner")
        print(f"PC bins: {len(all_pc_df)}")
    PCs = [col for col in all_pc_df.columns if col.startswith("PC")]
    PC_covar = [col for col in all_pc_df.columns if not col.startswith("PC") and col not in ["chrom", "start", "end", "index"]]
    PC_comb = list(product(PCs, PC_covar))
    PC_predf = pd.DataFrame(PC_comb, columns = ["PC", "var"])
    all_chr_PC = list()
    for chr, chr_df in all_pc_df.groupby("chrom"):
        cor_list = [chr_df[pc].corr(chr_df[var]) for pc, var in PC_comb]
        PC_predf["cor"] = cor_list
        PC_cor = PC_predf.pivot(index="PC", columns="var", values="cor")
        PC_cor.reset_index(inplace = True)
        PC_cor["chrom"] = chr
        all_chr_PC.append(PC_cor)
    pc_pearson_corr_df = pd.concat(all_chr_PC, axis = 0)
    # sort by chromosome and PC
    pc_pearson_corr_df["chrom"] = pd.Categorical(pc_pearson_corr_df["chrom"], categories=bf.fetch_chromsizes("hg38").index.tolist(), ordered=True)
    pc_pearson_corr_df.sort_values(["chrom", "PC"], inplace = True, ignore_index=True)
    if dhs is None:
        pc_pearson_corr_df["DHS"] = 0
    pc_pearson_corr_df["sign"] = np.where(pc_pearson_corr_df["GC"] < 0, -1, 1)
    for row in pc_pearson_corr_df.itertuples():
        if row.sign == -1:
            mask = all_pc_df["chrom"] == row.chrom
            all_pc_df.loc[mask, row.PC] *= -1
    for col in PC_covar:
        pc_pearson_corr_df[col] = pc_pearson_corr_df[col]*pc_pearson_corr_df["sign"]
    if exp is not None:
        print('Add expression report ..')
        exp_list = list()
        exp_df = bf.read_table(exp, schema = "bedGraph"); exp_df.rename(columns = {"value":"exp"}, inplace = True)
        pc_exp = bf.overlap(all_pc_df, exp_df, how = "inner")
        for chr, chr_sub in pc_exp.groupby("chrom"):
            for pc in PCs:
                active_A = len(chr_sub[(chr_sub[pc] > 0) & (chr_sub["exp_"] > 1)])
                active_B = len(chr_sub[(chr_sub[pc] < 0) & (chr_sub["exp_"] > 1)])
                exp_list.append((chr, pc, active_A, active_B)) 
        active_df = pd.DataFrame(exp_list, columns = ["chrom", "PC", "active_A", "active_B"])
        active_df["A%"] = round(active_df["active_A"]*100/(active_df["active_A"]+active_df["active_B"]), 2)
        pc_pearson_corr_df = pd.merge(pc_pearson_corr_df, active_df, on = ["chrom", "PC"])
        # active expression occupy A minimum 50%
        # pc_pearson_corr_df = pc_pearson_corr_df[pc_pearson_corr_df["A%"] > 0.5]
    else:
        pc_pearson_corr_df["active_A"] = 0; pc_pearson_corr_df["active_B"] = 0; pc_pearson_corr_df["A%"] = 0
    pc_pearson_corr_df.drop(["active_A", "active_B"], axis = 1, inplace = True)
    
    print("Write correlation results ...")
    # sort columns 
    pc_pearson_corr_df = pc_pearson_corr_df[["chrom", "PC", "sign"] + sorted([col for col in pc_pearson_corr_df.columns if col not in ["chrom", "PC", "sign"]])].copy()
    # all GC / TSS (DHS) show positive correlation
    #pc_pearson_corr_df["filter"] = np.where((pc_pearson_corr_df[[col for col in pc_pearson_corr_df.columns if col not in ["chrom", "PC", "sign"]]] >= 0).mean(axis = 1) == 1, 1, -1)
    print("select PC ...")
    if exp is not None and dhs is not None:
        pc_pearson_corr_df["select"] = np.where(((pc_pearson_corr_df["GC"] > 0.1) | (pc_pearson_corr_df["TSS"] > 0.1)) & (pc_pearson_corr_df["A%"] > 60) & (pc_pearson_corr_df["DHS"] > 0.1), 1, 0)
    if exp is not None and dhs is None:
        pc_pearson_corr_df["select"] = np.where(((pc_pearson_corr_df["GC"] > 0.1) | (pc_pearson_corr_df["TSS"] > 0.1)) & (pc_pearson_corr_df["A%"] > 60), 1, 0)
    if exp is None and dhs is not None:
        pc_pearson_corr_df["select"] = np.where(((pc_pearson_corr_df["GC"] > 0.1) | (pc_pearson_corr_df["TSS"] > 0.1)) & (pc_pearson_corr_df["DHS"] > 0.1), 1, 0)
    pc_pearson_corr_df.to_csv(f"{output}_PC.txt", sep = "\t", header = True, index = False)
    try:
        os.remove(f"{output}.bedGraph")
    except:
        pass 
    pc_select = pc_pearson_corr_df.query("select == 1");pc_select.index = range(len(pc_select))
    # when multiple PC remained, select top PC 
    if len(pc_select) != pc_select["chrom"].nunique(): 
        print("Selected PC is not unique per chromosome, select the top PC!")
        if exp:
            pc_select = pc_select[pc_select["A%"] == pc_select.groupby("chrom")["A%"].transform("max")]
        if not exp and dhs:
            pc_select = pc_select[pc_select["DHS"] == pc_select.groupby("chrom")["DHS"].transform("max")]
    if len(pc_select) == 23:
        print('Write into bedGraph ...')
        pc_select["PC#"] = pc_select["PC"].str.split("PC", expand = True)[1]
        print(pc_select)
        pc_select.to_csv(f"{output}_PC_select.txt", sep = "\t", header = True, index = False)
        for row in pc_select.itertuples():
            all_pc_df.query("chrom == @row.chrom")[["chrom", "start", "end", row.PC]].to_csv(f"{output}.bedGraph", sep = "\t", header = False, index = False, mode = "a")
        pc_value_df = bf.read_table(f"{output}.bedGraph", schema = "bedGraph")
        pc_value_df["z"] = pc_value_df.groupby("chrom")["value"].transform(value_scale)
        if np.equal(np.sign(pc_value_df["z"]), np.sign(pc_value_df["value"])).sum() == len(pc_value_df):
            print("After normalization (z-score), sign was not changed.")
        else:
            print("After normalization (z-score), sign changed. ")
            print("Correcting sign ...")
            pc_value_df["z"] = np.where(np.equal(np.sign(pc_value_df["z"]), np.sign(pc_value_df["value"])), pc_value_df["z"], -(pc_value_df["z"]))
        pc_value_df[["chrom", "start", "end", "z"]].to_csv(f"{output}_z.bedGraph", sep = "\t", header = False, index = False)

    else:
        print(f"### chromosome PC not unique for {dchic_folder}, reselect on lower score ...")
        exit(1)
        pc_pearson_corr_df["srank"] = pc_pearson_corr_df.groupby("chrom")["score"].rank(ascending=False)
        pc_pearson_corr_df["select"] = np.where(((pc_pearson_corr_df["GC_select"] == "pass") | (pc_pearson_corr_df["TSS_select"] == "pass")) & (pc_pearson_corr_df["srank"] == 1), 1, np.where(((pc_pearson_corr_df["GC_select"] == "pass") | (pc_pearson_corr_df["TSS_select"] == "pass")) & (pc_pearson_corr_df["srank"] == 2), 1, 0))
        pc_pearson_corr_df.drop(["srank"], axis = 1).to_csv(f"{output}.txt", sep = "\t", header = True, index = False)
        pc_select = pc_pearson_corr_df.query("select == 1");pc_select.index = range(len(pc_select))
        pc_select = pc_select[pc_select["PC"] == pc_select.groupby("chrom")["PC"].transform("min")]
        if len(pc_select) == 23:
            print(pc_select)
            print('Write into bedGraph ...')
            for row in pc_select.itertuples():
                all_pc_df.query("chrom == @row.chrom")[["chrom", "start", "end", row.PC]].to_csv(f"{output}.bedGraph", sep = "\t", header = False, index = False, mode = "a")

if __name__ == "__main__":
    main()
