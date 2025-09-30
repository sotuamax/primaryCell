import os 
from bam_tools import bam_count 
import pandas as pd 
import pysam 


def bam_count_loop(folder):
    bam_count_list = list()
    for bam in os.listdir(folder):
        if bam.endswith(".bam"):
            bam_file = os.path.join(folder, bam)
            if os.path.exists(bam_file + ".bai"):
                bam_count_list.append((bam_file.split(".bam")[0], bam_count(bam_file, n = 36)))
            else:
                pysam.index(bam_file)
    bam_count_df = pd.DataFrame(bam_count_list, columns = ["sample", "count"])
    return bam_count_df 


raw_bam = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/raw/individual"
qc_bam = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/individual"

raw_bam_count_df = bam_count_loop(raw_bam)
qc_bam_count_df = bam_count_loop(qc_bam)

bam_count_all = pd.concat([raw_bam_count_df, qc_bam_count_df], axis = 0)

bam_count_all.to_csv("bam_count_hitrac.txt", sep = "\t", header = True, index = False)
