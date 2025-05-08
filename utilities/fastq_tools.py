import pysam 

def read_count(fastq):
    """
    Parse fastq for read number
    Input: 
    fastq: fastq file
    Output:
    read count
    """
    q_handle = pysam.FastxFile(fastq)
    # read_df = pd.DataFrame([[entry.name, entry.comment, entry.sequence, entry.quality] for entry in q_handle]) 
    i = 0
    for entry in q_handle:
        i += 1
    return i

def read_len(fastq):
    """
    Load fastq file for read length 
    Input: 
    fastq: fastq file
    Output:
    read length
    """
    q_handle = pysam.FastxFile(fastq)
    for entry in q_handle:
        entry_len = len(entry.sequence)
        return entry_len

