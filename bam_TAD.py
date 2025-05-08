"""filter bam based on TAD: PET must within TAD region"""
import pysam 
import pandas as pd 
import bioframe as bf 
from joblib import Parallel, delayed 
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

domain = "/data/jim4/Seq/primary_cell_project/analysis/hi-trac/domain_nCAEC.txt" # only for chrom1
bam = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/QC/merge2/nCAEC.bam"
pairs = "/data/jim4/Seq/primary_cell_project/alignment/HiTrAC/format/nCAEC.pairs.gz"


domain_df = pd.read_table(domain, sep = "\t", header = 0)
domain_df = domain_df[domain_df["tag"] == "domain"]
domain_df.columns = ["chr", "id1", "start", "id2", "end", "tag", "size"]

bam_handle = pysam.AlignmentFile(bam, "rb")
pair_df = pd.read_table(pairs, sep = "\t", header = None, comment = "#", names = ["readID","chr1","pos1", "chr2", "pos2", "strand1", "strand2"], chunksize = 1e6)

chunk_df = list()
for chunk in pair_df:
    if chunk["chr1"].iloc[0] == "chr1":
        chunk1 = chunk[chunk["chr1"] == chunk["chr2"]].copy()
        chunk1 = chunk1[["chr1", "pos1", "pos2", "readID"]]
        for row in domain_df.itertuples():
            region = (row.chr, row.start, row.end)
            chunk2 = chunk1[(chunk1["pos1"] >= region[1]) & (chunk["pos2"] <= region[2])]
            chunk_df += chunk2["readID"].tolist()
    else:
        break 


with pysam.Alignmentfile("new.bam", "wb", template = header) as newbam:
    for read in bam_handle.fetch("chr1"):
        if read in chunk_df:
            newbam.write(read)
chunk_selected = pd.concat(pd.DataFrame(chunk_df), axis = 0)

