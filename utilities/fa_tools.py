import pysam 
from Bio import SeqIO


def ref_chrom(fa, out):
    """
    Process reference genome for chromosome only
    """
    seq = SeqIO.parse(fa, "fasta")
    with open(out, "w") as fan:
        for s in seq:
            if s.name.startswith("chr") and s != "chrM":
                fan.write(">"+s.description + "\n")
                fan.write(str(s.seq) + "\n")



