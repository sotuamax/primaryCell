from utilities.bam_tools import bamfilter
import glob 
import os 
import subprocess
import pysam 

for bam in glob.glob("/data/jim4/Seq/primary_cell_project/alignment/DNase/markdup/individual/*bam"):
    print(bam)
    output_bam = os.path.join("/data/jim4/Seq/primary_cell_project/alignment/DNase/QC/individual", os.path.basename(bam))
    bamfilter(bam, output_bam, two_end_clip_check=True, threads = 32)


