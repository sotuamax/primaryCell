import pandas as pd 

def featurecount_format(f):
    import os 
    df = pd.read_table(f, sep = '\t', header = 0, comment="#")
    df.columns = [os.path.basename(c).replace(".bam", "") for c in df.columns]
    bed_df = df[["Chr", "Start", "End", "Geneid"]]
    df.drop(["Chr", "Start", "End", "Strand", "Length"], axis = 1, inplace = True)
    return bed_df, df 

