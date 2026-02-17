#!/usr/bin/env python3
"""
Given single end BAM file, to generate pseudo-bedpe file by replicate single-end read twice. 

"""
import pysam 
import sys 
import pandas as pd 

bam = sys.argv[-1]
strand={True:"+", False:"-"}

print("Read BAM file ...")
bam_handle = pysam.AlignmentFile(bam, "rb", threads = 8)
bedpe_list=list()
print('Remove any reads on chrY and chrM ...')
print('Remove any reads on small scaffold ...')
with open(bam.replace(".bam", ".bedpe"), "w") as fw:
    for read in bam_handle.fetch():
        if (read.reference_name.startswith("chr")) and (read.reference_name not in ["chrY", "chrM"]):
            chr,start,end = (read.reference_name, read.reference_start, read.reference_end) 
            fw.write(f"{chr}\t{start}\t{end}\t{chr}\t{start}\t{end}\t{read.query_name}\t{read.mapping_quality}\t{strand[read.is_forward]}\t{strand[read.is_reverse]}\n") #.append((chr,start,end,chr,start,end,read.query_name,read.mapping_quality,strand[read.is_forward],strand[read.is_reverse]))

# bedpe_df = pd.DataFrame(bedpe_list)
# print("Write into bedpe file ...")
# bedpe_df.to_csv(bam.replace(".bam", ".bedpe"), sep = "\t", header = False, index = False)
