import bioframe as bf
import pandas as pd 

def chrom_size(ref = "hg38", filter = True):
    """
    Docstring for chrom_size
    To get chromosome size for a reference genome.
    :param ref: reference genome 
    :param filter: if remove chrY and chrM
    """
    ref_chrsize = bf.fetch_chromsizes(ref)
    ref_chrsize = pd.DataFrame.from_dict({"chrom":ref_chrsize.index, "start":([0]*len(ref_chrsize)), "end":ref_chrsize.values})
    if filter:
        ref_chrsize = ref_chrsize.query("chrom != 'chrM' and chrom != 'chrY'").copy()
    return ref_chrsize

def chrom_arms(ref = "hg38", filter = True):
    """
    Retrive chromosome arms for a given reference genome
    Input: 
    ref(str): reference genome to retrieve info from (default: hg38)
    Return:
    dataframe: chromosome arms of reference genome
    Note: 
    chrY & chrM are removed.
    """
    ref_chrsize = bf.fetch_chromsizes(ref)
    ref_chrsize = pd.DataFrame.from_dict({"chrom":ref_chrsize.index, "start":([0]*len(ref_chrsize)), "end":ref_chrsize.values})
    ref_cens = bf.fetch_centromeres(ref) 
    ref_arms = bf.subtract(ref_chrsize, ref_cens)
    if filter:
        ref_arms = ref_arms.query("chrom != 'chrM' and chrom != 'chrY'").copy()
    ref_arms.index = range(len(ref_arms))
    return ref_arms

