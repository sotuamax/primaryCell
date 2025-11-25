
import subprocess 
import os 
import pandas as pd

def BT_SE(r1, ref, output, n):
    """
    BOWTIE aligner
    """
    bowtie_command = f"bowtie2 -x {ref} -U {r1} -S {output} "
    return (bowtie_command)

def is_gzip_file(file):
    with open(file, 'rb') as f:
        magic = f.read(2)
        return magic == b'\x1f\x8b'

def flagstat(bam, n, out):
    command = f"samtools flagstat {bam} -@ {n} > {out}"
    subprocess.call(command, shell = True)

def STAR_SE(r1, index_dir, gtf, output, n, enable_novel = True):
    """"
    Perform single-end (SE) read alignment using STAR
    Input: 
    r1: fastq file for alignment
    index_dir: STAR index directory
    gtf: GTF file 
    output: output name prefix 
    n: threads 
    Output: 
    None
    """
    command = f"STAR --runThreadN {n} --genomeDir {index_dir} --sjdbGTFfile {gtf} --readFilesIn {r1} --outFileNamePrefix {output} --outSAMtype BAM SortedByCoordinate"
    if is_gzip_file(r1):
        command += " --readFilesCommand zcat " # disable twopassMode: --twopassMode Basic 
    if not enable_novel:
        command += " --alignSJDBoverhangMin 1 \
                    --alignSJoverhangMin 1000000 \
                    --alignIntronMin 1000000 \
                    --alignIntronMax 1 \
                    --outFilterIntronMotifs RemoveNoncanonical "
    # generate index for BAM 
    command += f" && mv {output}Aligned.sortedByCoord.out.bam {output}.bam && samtools index -@ {n} {output}.bam"
    print(command)
    subprocess.call(command, shell = True)
    # remove temporatory directories when exists
    import shutil 
    try:
        shutil.rmtree(f"{output}_STARgenome")
        shutil.rmtree(f"{output}_STARpass1")
    except Exception as e:
        pass

def STAR_PE(r1, r2, star_index, gtf, name, thread):
    """
    Perform paired-end (PE) read alignment using STAR
    Input: 
    r1: fastq file of read1
    r2: fastq file of read2
    star_index: STAR index directory 
    gtf: GTF file 
    name: output name prefix 
    thread: threads 
    Output:
    None 
    """
    STAR_process = f"STAR --runThreadN {thread} --genomeDir {star_index} --sjdbGTFfile {gtf} --readFilesIn {r1} {r2} --twopassMode Basic --sjdbOverhang 50 --outFileNamePrefix {name} --outSAMtype BAM SortedByCoordinate "
    if r1.endswith(".gz"):
        STAR_process += " --readFilesCommand zcat"
    STAR_process += f" && samtools index -@ {thread} {name}Aligned.sortedByCoord.out.bam"
    subprocess.call(STAR_process, shell = True)
    # remove temporatory directories when exists
    import shutil 
    try:
        shutil.rmtree(f"{name}_STARgenome")
        shutil.rmtree(f"{name}_STARpass1")
    except Exception as e:
        pass

def BWA_SE(r1, ref, output, n):
    """
    Perform single-end (SE) read alignment using BWA
    Input: 
    r1: fastq file for alignment
    ref: reference genome
    output: output name prefix 
    n: threads 
    Output: 
    None
    """
    command = f"bwa mem -t {n} {ref} {r1} | samtools view -@ {n} -Su - | samtools sort -@ {n} - -o {output}.bam && samtools index -@ {n} {output}.bam"
    # logging.info(f"Command\t{command}")
    print(command)
    subprocess.call(command, shell = True)

def BWA_PE(r1, r2, ref, output, n):
    """
    Perform paired-end (PE) read alignment using BWA
    Input: 
    r1: read1 fastq file 
    r2: read2 fastq file
    ref: reference genome
    output: output name prefix 
    n: threads 
    Output: 
    None
    """
    align_command = f"bwa mem -t {n} {ref} {r1} {r2} | samtools view -@ {n} -Su - | samtools sort -@ {n} - -o {output}.bam && samtools index -@ {n} {output}.bam"
    subprocess.call(align_command, shell = True)

def markdup(bam, output):
    """
    Markdup bam file (applies to bulk cell data)
    Input:
    bam: bam file 
    output: output prefix name 
    Output: 
    None
    """
    # 
    picard="/usr/local/apps/picard/3.1.0/picard.jar"
    markdup_command = f"java -jar {picard} MarkDuplicates -I {bam} -O {output}.bam -M {output}.txt --REMOVE_DUPLICATES true"
    subprocess.call(markdup_command, shell = True)

def featurecount(bam, gtf, name, thread):
    """
    Perform featureCounts on the bam input (PE mode)
    Input: 
    bam: bam file 
    gtf: GTF file 
    name: output file prefix 
    thread: thread
    Output: 
    None
    """
    count_command = f"featureCounts -p --countReadPairs -C -B -a {gtf} -F GTF --extraAttributes gene_name -s 0 -T {thread} {bam} -o {name}"
    print(count_command)
    subprocess.call(count_command, shell = True)

def parse_featurecount(name):
    """
    clean the format of featureCount output
    Input: 
    name: file name to open with 
    Output:
    pd.DataFrame: dataframe of clean featureCount"""
    feature_df = pd.read_table(name, sep = "\t", header = 0, comment = "#")
    new_columns = list()
    for c in feature_df.columns:
        if c.endswith(".bam"):
            new_columns.append(os.path.basename(c))
        else:
            new_columns.append(c)
    feature_df.columns = new_columns 
    try:
        feature_df.drop(["Chr", "Start", "End", "Strand", "Length"], axis = 1, inplace = True)
        return feature_df
    except:
        return feature_df

def rsem_quantify(r1, r2, ref_dir, name, thread):
    """
    Perform feature quantification using RSEM (PE mode)
    Input: 
    r1: read1 fastq file
    r2: read2 fastq file
    ref: reference genome
    name: output name 
    thread: thread to use
    Output:
    None
    """
    # module = " ".join(modules)
    # p = os.system(f"module load {module} | echo $PATH", shell = True)
    # p.communicate()
    # current_env = module_load(modules)
    rsem_command = f"rsem-calculate-expression --paired-end --star --strandedness none -p {thread} "
    if r1.endswith(".gz"):
        rsem_command += " --star-gzipped-read-file"
    rsem_command += f" {r1} {r2} {ref_dir} {name}"
    subprocess.call(rsem_command, shell = True)

def parse_rsem(name):
    """
    Parse output of rsem
    Input: 
    name: output name assigned for rsem
    Output:
    cleaned results of rsem
    """
    gene_results = pd.read_table(name + ".genes.results", sep = "\t", header = 0)
    gene_results = gene_results[["gene_id", "TPM", "FPKM"]]
    return gene_results
