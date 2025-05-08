import pysam 
import subprocess
import os 
import pandas as pd 

def examine_sort(bam):
    """
    Examine BAM as name or coordinates sorted
    """
    bamfile = pysam.AlignmentFile(bam, "rb")
    # Access the header dictionary
    header = bamfile.header
    # Check sorting order from @HD line
    sort_order = header.get('HD').get('SO')
    # standard output: queryname (sort by name); coordinate (sort by coordinate)
    return sort_order 

def nsort_clean(bam:str, new_bam:str, n = 1):
    """
    When QC applied on each read, read and its mate may not be retained at the same time. This func is to clean BAM file for reads that forms in pair. 
    Input: 
    bam: name sorted bam file 
    new_bam: new bam file name to write into
    n: threads to use
    Return: 
    write into new bam file with read and its mate sorted by name. 
    """
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    for read in bam_handle:
        read0 = read 
        break
    with pysam.AlignmentFile(new_bam, "wb", template = bam_handle, threads = n) as newbam:
        for read in bam_handle:
            if read.query_name == read0.query_name:
                newbam.write(read0)
                newbam.write(read)
            read0 = read 

def bam2bed(bam, label = "cDNA", read_seq = True):
    """
    Transform bam into bed format
    Input: 
    bam: bam file 
    label: bam file type (cDNA/DNA)
    Output: 
    read bed file (chrom, start, end, name, score, strand, cigar, read-seq as optional)
    """
    bam_handle = pysam.AlignmentFile(bam, "rb")
    strandness = {True:"+", False:"-"}
    # chrom, start, end, name, score, strandness, cigar
    # read_df = pd.DataFrame(read_df, columns = ["read", "seq", 'strand', "cigar"])
    if label.upper() in ["CDNA", "RNA"]:
        # for RNA-seq alignment, ("NH", 1) to get unique alignment (no quality filter, MAPQ see score)
        if read_seq == True:
            read_bed = pd.DataFrame([[read.reference_name, read.reference_start,read.reference_end, read.query_name, read.mapping_quality, strandness[read.is_forward],read.cigarstring, read.query_sequence] for read in bam_handle.fetch() if ("NH", 1) in read.tags])
            read_bed.columns = ["chrom", "start", "end", "name", "score", "strand", "cigar", "seq"]
        else:
            read_bed = pd.DataFrame([[read.reference_name, read.reference_start,read.reference_end, read.query_name, read.mapping_quality, strandness[read.is_forward],read.cigarstring] for read in bam_handle.fetch() if ("NH", 1) in read.tags])
            read_bed.columns = ["chrom", "start", "end", "name", "score", "strand", "cigar"]
    if label.upper() == "DNA":
        # for DNA-seq alignment, filtering on primary alignment (no quality filter, MAPQ see score)
        if read_seq == True:
            read_bed = pd.DataFrame([[read.reference_name,read.reference_start,read.reference_end,read.query_name, read.mapping_quality, strandness[read.is_forward], read.cigarstring, read.query_sequence] for read in bam_handle.fetch() if (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped)])
            read_bed.columns = ["chrom", "start", "end", "name", "score", "strand", "cigar", "seq"]
        else:
            read_bed = pd.DataFrame([[read.reference_name,read.reference_start,read.reference_end,read.query_name, read.mapping_quality, strandness[read.is_forward], read.cigarstring] for read in bam_handle.fetch() if (not read.is_secondary) and (not read.is_supplementary) and (not read.is_unmapped)])
            read_bed.columns = ["chrom", "start", "end", "name", "score", "strand", "cigar"]
    return read_bed

def bam_count(bam, chr = None, start = None, end = None, n = 1):
    """
    Count number of reads within specified region of given bam file
    Input: 
    bam: bam file 
    chr: chromosome region 
    start: chromosome start region 
    end: chromosome end region 
    n: thread to process bam file 
    Output: 
    number of reads within the specified region of bam file
    """
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    read_set = {read.query_name for read in bam_handle.fetch(chr, start, end)}
    # for read in bam_handle.fetch()
    return len(read_set)

def bam_count_fast(bam, chr = None, start = None, end = None, n = 1):
    """
    For PE reads, count total number of reads and then divide it by 2 to get estimated read number
    # Note: it is not accurate, as some paired-reads only have 1-end of it pass quality control. 
    """
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    read_count = bam_handle.count(chr, start,end)
    return read_count//2

def bam_merge(bams:list, name, thread, labels = None):
    file_combined = " ".join(bams)
    if labels != None:
        RG_label = "\n".join([f'@RG\tID:{id}\tSM:{group}\tLB:{id}\tPL:ILLUMINA' for id in group_id])
        command = f"printf '{RG_label}' > {group}.txt"
        subprocess.call(command, shell = True)
        merge_command = f'samtools merge -rh {name}.txt -@ {thread} -o - {file_combined} | samtools sort -@ {thread} - -o {name}.bam && samtools index -@ {thread} {name}.bam && rm {name}.txt'
    else:
        merge_command = f"samtools merge -@ {thread} -o - {file_combined} | samtools sort -@ {thread} - -o {name}.bam && samtools index -@ {thread} {name}.bam"
    if not os.path.exists(name + ".bam"):
        print(merge_command)
        subprocess.call(merge_command, shell = True)

def bed2saf(bed):
    """
    Transfer bed file input saf format (see featureCount)
    Input: 
    bed: bed file 
    Output:
    None
    """
    import bioframe as bf 
    name = os.path.basename(bed).replace("bed", "saf")
    bed_df = bf.read_table(bed, schema = "bed5")
    saf = bed_df[["name", "chrom", "start", "end"]]
    saf["strand"] = "."
    saf.to_csv(name, sep = "\t", header = False, index = False)

def bamfilter(bam, select_read, output):
    """
    Rewrite new bam file by given selected reads name"""
    bam_handle = pysam.AlignmentFile(bam, "rb")
    new_bam_name = os.path.join(output, output + "_clean.bam")
    with pysam.AlignmentFile(new_bam_name, "wb", template = bam_handle) as newbam:
        for read in bam_handle.fetch():
            if read.query_name in select_read:
                newbam.write(read)
    try:
        subprocess.call(f"samtools index -@ 4 {new_bam_name}", shell = True)
    except Exception as e:
        print(e)
