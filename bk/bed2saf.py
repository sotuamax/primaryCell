#!/usr/bin/env python3
import pandas as pd 
import bioframe as bf 
import sys 
import os 

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
    
bed = sys.argv[-1]
name = os.path.basename(bed).replace("bed", "saf")

bed_df = bf.read_table(bed, schema = "bed5")
saf = bed_df[["name", "chrom", "start", "end"]]
saf["strand"] = "."

saf.to_csv(name, sep = "\t", header = False, index = False)
