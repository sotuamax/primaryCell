import pandas as pd 
import pysam 

def pair_stat(pair, chrs):
    pair_chunk = pd.read_table(pair, sep = "\t", comment = "#", header = None, chunksize = 1000000, names = ["id", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2", "pairtype"])
    pair_count = {pair:0 for pair in combinations_with_replacement(chrs, 2) if chrs.index(pair[0]) <= chrs.index(pair[-1])}
    for chunk in pair_chunk:
        for g, chunk_g in chunk.groupby(["chr1", "chr2"]):
            pair_count[g] += len(chunk_g)
    count_num = pd.DataFrame.from_dict(pair_count, orient = "index", columns = ["count"])
    return count_num


def compress_index(out):
    # compress pairs text file using bgzip 
    subprocess.call(f"bgzip -c {out}.pairs > {out}.pairs.gz", shell = True)
    # generate index for compressed pairs file  (index file with px2 suffix)
    subprocess.call(f"pairix -p pairs {out}.pairs.gz", shell = True)


def pairs_header(chr_df, assembly):
    """"
    Generate pairs header from chrsize df
    """
    chromsize_list = [f"#chromsize: {row.Index} {row.size}" for row in chr_df.itertuples()]
    chromsize_total = "\n".join(chromsize_list)
    # note that empty line is not allowed for pairs file
    header = f"## pairs format v1.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n#genome_assembly: {assembly}\n{chromsize_total}\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
    return header 

def write_pairs_from_bam(bam, chrom_order, out, n):
    """"
    Self-defined function to write paris format from bam input 
    deprecated: use following instead
    pairtools parse --assembly hg38 -c chrsize input.bam -o output.pairs
    """
    # strand dictionary 
    strand_dict = {True:"+", False:"-"}
    # parse read by read on coordinate-sorted file
    # chr_pair = [c for c in combinations_with_replacement(chrom_order, 2)]
    bam_handle = pysam.AlignmentFile(bam, "rb", threads = n)
    all_read = {read.query_name for read in bam_handle.fetch()}
    print("Retrieve all reads ....")
    with open(out + ".tmp", "w") as fw:
        for chr in chrom_order:
            for read in bam_handle.fetch(chr):
                try:
                    all_read.remove(read.query_name)
                    # sort it and then write into tmp
                    if read.reference_name == read.next_reference_name: 
                        if read.reference_start < read.next_reference_start:
                            # chr_list.append()
                            fw.write(f"{read.query_name}\t{read.reference_name}\t{read.reference_start}\t{read.next_reference_name}\t{read.next_reference_start}\t{strand_dict[read.is_forward]}\t{strand_dict[read.mate_is_forward]}\tUU\n")
                        else:
                            fw.write(f"{read.query_name}\t{read.next_reference_name}\t{read.next_reference_start}\t{read.reference_name}\t{read.reference_start}\t{strand_dict[read.mate_is_forward]}\t{strand_dict[read.is_forward]}\tUU\n")
                    else:
                        if read.next_reference_name in chrom_order:
                            if chrom_order.index(read.reference_name) < chrom_order.index(read.next_reference_name):
                                fw.write(f"{read.query_name}\t{read.reference_name}\t{read.reference_start}\t{read.next_reference_name}\t{read.next_reference_start}\t{strand_dict[read.is_forward]}\t{strand_dict[read.mate_is_forward]}\tUU\n")
                            else:
                                fw.write(f"{read.query_name}\t{read.next_reference_name}\t{read.next_reference_start}\t{read.reference_name}\t{read.reference_start}\t{strand_dict[read.mate_is_forward]}\t{strand_dict[read.is_forward]}\tUU\n")
                except KeyError:
                    pass 

