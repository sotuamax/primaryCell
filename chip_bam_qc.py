#!/usr/bin/env python
import os 
import sys 
import pysam 

bam=sys.argv[-2]
mq=sys.argv[-3]
new_bam=sys.argv[-1]

bam_handle = pysam.AlignmentFile(bam, "rb", threads = 24)
header = str(bam_handle.header)
new_header = list()
for line in header.split("\n"):
    if line.startswith("@SQ"):
        if line.split("\t")[1].split(":")[-1].startswith("chr"):
            new_header.append(line)
    else:
        new_header.append(line)
new_header = "\n".join(new_header)
with pysam.AlignmentFile(new_bam, "wb", text = new_header, threads = 24) as newbam:
    for read in bam_handle.fetch():
        if (not read.is_secondary) and (not read.is_supplementary) and (read.mapping_quality > int(mq)) and (read.reference_name.startswith("chr")) and (read.reference_name not in ["chrM"]) and (read.is_read1):
            newbam.write(read)
