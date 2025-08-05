#!/usr/bin/env python3

"""
stratified linkage disequilibrium score regression (S-LDSC)
A method to test if specific genomic region contributes significantly to a phenotypic inheritance. 

ldscore regression weights:
for each SNPs, weights measure its confidence values (higher weights, better confidence)
SNPs with high LD and low MAF generate more noises in statistic test, so lower weights 
downweight low conficence SNPs to improve estimate of h2 
for each chromosome, ldscore, M, M_5_50, log (ldscore w/ header : CHR, SNP, BP, L2)

reference SNP: bim/bed/fam from plink (based on 1000G) 
observation on ~1000 chrom (2 per individual), and ~10 Million SNPs across chr1-22 
for each chromosome, bed (genotype), bim, fam, frq 
bim wo/ header : SNP; CHR, SNP, cM, BP, A1, A2 (genome version matters)
fam wo/ header : family; FID (familyID), IID (individual ID), PID (Paternal ID, 0 for missing), MID (maternal ID), sex (0 for unknown), phenotype (-9/0 missing) 489 
frq w/ header : SNP frequency; CHR, SNP, A1, A2, MAF, NCHROBS (number of chrom observed, 2*number of individuals)

baseline model:
for each chromosome, (annot.gz), l2.ldscore.gz, M, M_5_50, log
annot.gz w/ header: CHR, BP, SNP, CM, base (all 1s), annot1, annot2, annot3, annot4, annot5, .... (match to *bim for each chromosome)
l2.ldscore w/ header : CHR, SNP, BP, baseL2, annot1, annot2, annot3, annot4, annot5, ... (hm3 SNP only to reduce heritabilty calculation time)

"""
import argparse 
import pandas as pd 
import glob 
import subprocess
import os 
import bioframe as bf

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", 
                                   formatter_class = argparse.RawDescriptionHelpFormatter, 
                                   description="Given bed file of interest genomic region, to test if regions significantly contributed to a specific phenotype (contain higher inheritance compared to random region)", 
                                   help = "mpiexec -n 12 h_enrich.py -bed region.bed  -gwas /data/jim4/Seq/primary_cell_project/data/GWAS/sumstats_63/*.gz -outdir output")
    parser.add_argument("-gwas", "--gwas", nargs = "+", help = "summary statistics for GWAS")
    parser.add_argument("-bed", "--bed", nargs = "+", help = "bed file for regions of interest to test for partitioned inheritance")
    parser.add_argument("-ref_SNP", "--ref_SNP", help = "folder for reference SNP used to annotated bed file (bim/bed/frq/fam)", default = "/data/jim4/Seq/primary_cell_project/data/GWAS/GRCh38/plink_files/")
    parser.add_argument("-baseline", "--baseline", help = "folder for baseline model used as conditioning", default = "/data/jim4/Seq/primary_cell_project/data/GWAS/GRCh38/baselineLD_v2.2")
    parser.add_argument("-w", "--weights", help = "folder for regression weights", default = "/data/jim4/Seq/primary_cell_project/data/GWAS/GRCh38/weights")
    parser.add_argument("-outdir", "--outdir", help = "output directory", default = ".")
    args=parser.parse_args()
    return args

def bed_examine(bed):
    """
    bed examination 
    """
    bed_df = bf.read_table(bed, schema = "bed3")
    bed_cluster = bf.cluster(bed_df, min_dist = 1)
    bed_cluster = bed_cluster[["chrom", "cluster_start", "cluster_end"]].drop_duplicates(keep = "first")
    if len(bed_df) != len(bed_cluster):
        return bed_cluster 
    else:
        return None

def main():
    from mpi4py import MPI
    from mpi4py.util import pkl5
    comm = pkl5.Intracomm(MPI.COMM_WORLD) # to overcome the overflow error when comm data > 2 GB
    rank = comm.Get_rank()
    size = comm.Get_size()
    # 
    args = args_parser()
    gwas = args.gwas # summary statistics of GWAS 
    bed_list = args.bed # specific region of interest for test in inheritance 
    # prepare annotation file (per chromosome)
    # output: SNPs labeled a 0/1 as overlapping/non-overlapping with region in bed (it should match to SNP in *bim)
    if rank == 0:
        for bed in bed_list:
            annot_dir = bed.replace(".bed", "_annot")
            try:
                os.mkdir(annot_dir)
            except:
                pass 
            annot_list = list(); annot_files = list()
            for bim_file in glob.glob(os.path.join(args.ref_SNP, "*.bim")):
                chr = bim_file.split(".")[-2]
                out_name = os.path.basename(bed).replace("bed", chr+".annot.gz")
                out_annot = os.path.join(annot_dir, out_name)
                annot_run = f"make_annot.py --bed-file {bed} --bimfile {bim_file} --annot-file {out_annot}"
                annot_files.append(out_annot)
                if not os.path.exists(out_annot):
                    annot_list.append(annot_run)
        w_idx = 0
        worker_tasks1 = {w:[] for w in range(size)}
        for ann in annot_list:
            worker_tasks1[w_idx].append(ann)
            w_idx = (w_idx + 1) % size
    else:
        worker_tasks1 = None
        annot_files = None
    worker_tasks1 = comm.bcast(worker_tasks1, root = 0)
    annot_files = comm.bcast(annot_files, root = 0)
    for t in worker_tasks1[rank]:
        subprocess.call(t, shell = True)
    # calculate partitioned LD score (per chromosome)
    # for each chromosome, each matched snp, bim/bed/fam and annot file are given 
    # to build on top of baseline mode (SNPs ldscore should match as in 1000G and Hapmap3)
    if rank == 0:
        l2_list = list()
        for a in annot_files:
            chr = a.split(".")[-3] # *.annot.gz 
            l2_file = a.replace(".annot.gz", "")
            com_n = os.path.commonprefix([f for f in os.listdir(args.ref_SNP)])
            bfile = os.path.join(args.ref_SNP, f"{com_n}{chr}") # prefix for bim file
            snp = os.path.join(args.baseline, f"baselineLD.{chr}.snp")
            l2_command = f"ldsc.py --l2 --bfile {bfile} --ld-wind-cm 1 --annot {a} --thin-annot --out {l2_file} --print-snps {snp}" # --thin-annot (input annot file will take bim file for its row name)
            if not os.path.exists(l2_file + ".l2.ldscore.gz") or not os.path.exists(l2_file + ".l2.M") or not os.path.exists(l2_file + ".l2.M_5_50"):
                l2_list.append(l2_command)
                # subprocess.call(l2_command, shell = True)
        w_idx = 0
        worker_tasks2 = {w:[] for w in range(size)}
        for l2 in l2_list:
            worker_tasks2[w_idx].append(l2)
            w_idx = (w_idx+1) % size 
    else:
        worker_tasks2 = None 
    worker_tasks2 = comm.bcast(worker_tasks2, root = 0)
    for t in worker_tasks2[rank]:
        subprocess.call(t, shell = True)
    # partitioned heritability (conditioned on baseline model)
    if rank == 0:
        test = 1
    else:
        test = None 
    test = comm.bcast(test, root = 0)
    # estimate h2 score for gwas summary statistics score (z > 6 pass to downstream analysis)
    h2_estimate = f"ldsc.py --h2 {g} --ref-ld-chr {eur_ldscore} --w-ld-chr {weight_ldscore} --frqfile-chr {frq_file} --out {out}"

    if rank == 0:
        try:
            os.mkdir(args.outdir)
        except:
            pass
        baseline_prefix = os.path.join(args.baseline, os.path.commonprefix([f for f in os.listdir(args.baseline)]))
        frq_prefix = os.path.join(args.ref_SNP, os.path.commonprefix([f for f in os.listdir(args.ref_SNP)]))
        weight_prefix = os.path.join(args.weights, os.path.commonprefix([f for f in os.listdir(args.weights)]))
        # combine all annot
        annot_prefix_list = list()
        for bed in bed_list:
            annot_dir = bed.replace(".bed", "_annot")
            annot_prefix = os.path.join(annot_dir, os.path.commonprefix([f for f in os.listdir(annot_dir)]))
            annot_prefix_list.append(annot_prefix)
        annot_prefix_collect = ",".join(annot_prefix_list)
        h2_list = list()
        for g in gwas:
            out = os.path.join(args.outdir, os.path.basename(g).replace(".sumstats.gz", ""))
            h2_command = f"ldsc.py --h2 {g} --overlap-annot --ref-ld-chr {annot_prefix_collect},{baseline_prefix} --frqfile-chr {frq_prefix} --w-ld-chr {weight_prefix} --out {out} --print-coefficients"
            if not os.path.exists(out + ".results"):
                # print(h2_command)
                h2_list.append(h2_command)
                # subprocess.call(h2_command, shell = True)
        w_idx = 0
        worker_tasks3 = {w:[] for w in range(size)}
        for h2 in h2_list:
            worker_tasks3[w_idx].append(h2)
            w_idx = (w_idx+1) % size 
    else:
        worker_tasks3 = None
        annot_prefix_collect = None
    worker_tasks3 = comm.bcast(worker_tasks3, root = 0)
    annot_prefix_collect = comm.bcast(annot_prefix_collect, root = 0)
    for t in worker_tasks3[rank]:
        subprocess.call(t, shell = True)
        # parse results file 
        out = t.split("--out ")[-1].split(" ")[0]
        df = pd.read_table(out + ".results", sep = "\t", header = 0)
        # rename L2 to annot_file_applied
        for aa in enumerate(annot_prefix_collect.split(",")):
            df.iloc[aa[0], 0] = os.path.basename(aa[-1])
        df.to_csv(out+".results", sep = "\t", header = True, index = False)

if __name__ == "__main__":
    main()