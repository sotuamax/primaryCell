import pysam 
import pandas as pd 
import bioframe as bf 
import numpy as np 

def parse_gtf(gtf, feature = "gene", ENCODE = True):
    """
    Get the feature coordinates information from GTF file
    Input: 
    gtf: GTF file 
    feature: cds/gene/exon
    Output: 
    feature coordinates information in the given GTF file
    """
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    if ENCODE:
        # collection all CDS (its coordinates, strand, and frame)
        if feature.upper() == "EXON":
            gtf_list = [(g.contig, g.start, g.end, g.frame, g.transcript_id, g.gene_id, g.gene_name, g.strand) for g in gtf_handle.fetch() if g.feature.upper() == feature.upper()]
            gtf_df = pd.DataFrame(gtf_list, columns = ["chrom", "start", "end", "frame", "transcript", "gene", "symbol", "strand"])
        if feature.upper() == "CDS":
            gtf_list = [(g.contig, g.start, g.end, g.frame, g.transcript_id, g.gene_id, g.gene_name, g.strand) for g in gtf_handle.fetch() if g.feature.upper() == feature.upper()]
            gtf_df = pd.DataFrame(gtf_list, columns = ["chrom", "start", "end", "frame", "transcript", "gene", "symbol", "strand"])
        if feature.upper() == "GENE":
            gtf_list = {(g.contig, g.start, g.end, g.gene_id, g.gene_name, g.strand) for g in gtf_handle.fetch() if g.feature == feature}
            gtf_df = pd.DataFrame(gtf_list, columns = ["chrom", "start", "end", "gene", "symbol", "strand"])
        if feature.upper() == "TRANSCRIPT":
            gtf_list = [(g.contig, g.start, g.end, g.transcript_id, g.gene_type, g.gene_id, g.gene_name, g.strand) for g in gtf_handle.fetch() if g.feature.upper() == feature.upper()]
            gtf_df = pd.DataFrame(gtf_list, columns = ["chrom", "start", "end", "transcript", "type", "gene", "symbol", "strand"])
        if feature.upper() == "STOP":
            pass
    else:
        if feature.upper() == "TRANSCRIPT":
            gtf_list = [(g.contig, g.start, g.end, g.transcript_id, g.gene_id, g.gbkey, g.strand, g.gene, g.product) for g in gtf_handle.fetch() if g.feature.upper() == feature.upper()]
            gtf_df = pd.DataFrame(gtf_list, columns = ["chrom", "start", "end", "transcript", "gene", "type", "strand", "symbol", "product"])
    return gtf_df

def protein_coding(gtf, out):
    """
    Filter GTF for protein coding genes only 
    """
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF())
    gtf_header = "\n".join(gtf_handle.header)
    with open(out, "w") as fw:
        fw.write(gtf_header+ "\n") # write header (add \n after header written)
        for g in gtf_handle.fetch(): 
            if g.gene_type == "protein_coding":
                # wrote protein coding features 
                fw.write(str(g) + "\n")


def parse_TSS(gtf, feature = "gene"):
    """
    Find Transcript Start Site (TSS) in GTF file
    Input: 
    gtf: gtf file 
    Output: 
    TSS positon of gene
    """
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF(), threads = 2)
    if feature.lower() == "transcript":
        gene_TSS = [(row.contig, row.start, row.end, row.gene_id, row.transcript_id, row.gene_name, row.strand) for row in gtf_handle.fetch() if row.feature == "transcript"]
        TSS_df = pd.DataFrame(gene_TSS, columns = ["chrom", "start", "end", "gene", "transcript", "symbol", "strand"])
        TSS_df["pos"] = np.where(TSS_df["strand"] == "+", TSS_df["start"], TSS_df["end"])
        TSS_df = TSS_df[["chrom", "pos", "gene", "transcript", "symbol", "strand"]].copy()
    if feature.lower() == "gene":
        gene_TSS = [(row.contig, row.start, row.end, row.gene_id, row.gene_name, row.strand) for row in gtf_handle.fetch() if row.feature == "gene"]
        TSS_df = pd.DataFrame(gene_TSS, columns = ["chrom", "start", "end", "gene", "symbol", "strand"])
        TSS_df["pos"] = np.where(TSS_df["strand"] == "+", TSS_df["start"], TSS_df["end"])
        TSS_df = TSS_df[["chrom", "pos", "gene", "symbol", "strand"]].copy()
    return TSS_df

def parse_intron(gtf):
    """"
    Retrieve intron region in GTF file 
    Input: 
    gtf: gtf file 
    Output: 
    intron coordinates 
    """
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF(), threads = 2)
    exon_bed = pd.DataFrame([(row.contig, row.start, row.end, row.transcript_id) for row in gtf_handle.fetch() if row.feature == "exon"], columns = ["chrom", "start", "end", "transcript", "strand"])
    intron_df_list = list()
    for label, transcript_group in exon_bed.groupby(["chrom", "transcript", "strand"]):
        chrom = label[0]
        transcript_s = label[1]
        strand_s = label[2]
        intron_list = sorted(transcript_group["start"].tolist()+transcript_group["end"].tolist())[1:-1]
        intron_start = intron_list[::2];intron_end = intron_list[1::2]
        intron_df = pd.DataFrame({"chrom" :chrom, "start":pd.Series(intron_start, dtype = "int"), "end": pd.Series(intron_end, dtype = "int"), "transcript":transcript_s, "strand":strand_s})
        intron_df_list.append(intron_df)
    intron_df_all = pd.concat(intron_df_list, axis = 0)
    return intron_df_all 

def noncoding(gtf, chr):
    """
    Retrieve non-coding exon region
    Input: 
    gtf: GTF file"""
    gtf_handle = pysam.TabixFile(gtf, parser = pysam.asGTF(), threads = 2)
    cds_bed = pd.DataFrame([(row.contig, row.start, row.end, row.transcript_id) for row in gtf_handle.fetch(chr) if row.feature == "CDS"], columns = ["chrom", "start", "end", "transcript"])
    exon_bed = pd.DataFrame([(row.contig, row.start, row.end, row.transcript_id) for row in gtf_handle.fetch(chr) if row.feature == "exon"], columns = ["chrom", "start", "end", "transcript"])
    cds_bed["label"] = 0
    exon_bed["label"] = range(1, len(exon_bed)+1)
    cds_exon_combined = pd.concat([cds_bed, exon_bed], axis = 0)
    nocoding_list = list()
    for transcript, transcript_group in cds_exon_combined.groupby("transcript"):
        nocoding_transcript = bf.complement(transcript_group[transcript_group["label"].to_numpy() == 0], transcript_group[transcript_group["label"].to_numpy() != 0], view_name_col = "label").drop("view_region", axis = 1)
        nocoding_transcript["transcript"] = transcript
        nocoding_list.append(nocoding_transcript)
    # combine all list 
    nocoding_df = pd.concat(nocoding_list, axis = 0)
    # nocoding_df.drop_duplicates(keep = "first", inplace = True, ignore_index = True)
    # nocoding_df["label"] = "ncoding"
    return nocoding_df
