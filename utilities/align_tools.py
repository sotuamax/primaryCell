
import subprocess 
import os 
import logging
import pandas as pd
from utility import module_load 


def STAR_SE(r1, index_dir, gtf, output, n):
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
    command = f"STAR --runThreadN {n} --genomeDir {index_dir} --sjdbGTFfile {gtf} --readFilesIn {r1} --outFileNamePrefix {output} --outSAMtype BAM SortedByCoordinate && samtools index -@ {n} {output}Aligned.sortedByCoord.out.bam" # disable twopassMode: --twopassMode Basic 
    if r1.endswith(".gz"):
        commmand += " --readFilesCommand zcat"
    modules = ["STAR", "samtools"]
    my_env = module_load(modules)
    # 
    subprocess.call(command, shell = True, env = my_env)
    # remove temporatory directories when exists
    import shutil 
    if os.path.exists(f"{output}_STARgenome"):
        try:
            shutil.rmtree(f"{output}_STARgenome")
        except Exception as e:
            print(e)
    if os.path.exists(f"{output}_STARpass1"):
        try:
            shutil.rmtree(f"{output}_STARpass1")
        except Exception as e:
            print(e)

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
    STAR_process = f"STAR --runThreadN {thread} --genomeDir {star_index} --sjdbGTFfile {gtf} --readFilesIn {r1} {r2} --twopassMode Basic --sjdbOverhang 50 --outFileNamePrefix {name} --outSAMtype BAM SortedByCoordinate && samtools index -@ {thread} {name}Aligned.sortedByCoord.out.bam"
    if r1.endswith(".gz"):
        STAR_process += " --readFilesCommand zcat"
    modules = ["STAR", "samtools"]
    my_env = module_load(modules)
    subprocess.call(STAR_process, shell = True, env = my_env)
    # remove temporatory directories when exists
    import shutil 
    if os.path.exists(f"{name}_STARgenome"):
        try:
            shutil.rmtree(f"{name}_STARgenome")
        except Exception as e:
            print(e)
    if os.path.exists(f"{name}_STARpass1"):
        try:
            shutil.rmtree(f"{name}_STARpass1")
        except Exception as e:
            print(e)

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
    command = f"bwa mem -t {n} {ref} {r1} -L 0,0 | samtools view -@ {n} -Su - | samtools sort -@ {n} - -o {output}.bam && samtools index -@ {n} {output}.bam"
    # logging.info(f"Command\t{command}")
    # print(command)
    modules = ["bwa", "samtools"]
    my_env = module_load(modules)
    subprocess.call(command, shell = True, env = my_env)

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
    modules = ["bwa", "samtools"]
    my_env = module_load(modules)
    subprocess.call(align_command, shell = True, env = my_env)

def STAR_log(output):
    with open(output, "r") as fr:
        for line in fr.readlines():
            line = line.strip()
            if "Number of input reads |" in line:
                input_r = line.split("Number of input reads |\t")[-1]
            if "Uniquely mapped reads number |" in line or "Uniquely mapped reads % |" in line:
                if "Uniquely mapped reads number |" in line:
                    unique_r = line.split("Uniquely mapped reads number |\t")[-1]
                if "Uniquely mapped reads % |" in line:
                    unique_p = line.split("Uniquely mapped reads % |\t")[-1].strip("%")
    logging.info(f"Align input\t{input_r}")
    logging.info(f"Align unique\t{unique_r}({unique_p}%)")

def BWA_log(bam):
    """
    Parse bam file for alignment report 
    Input:
    bam: bam file 
    Output: 
    alignment stat
    """
    from bam_tools import bam2bed
    read_bed = bam2bed(bam, label = "DNA", read_seq = False)
    read_num = len(set(read_bed["name"]))
    logging.info(f"Align unique\t{read_num}")

def markdup(bam, output):
    """
    Markdup bam file 
    Input:
    bam: bam file 
    output: output prefix name 
    Output: 
    None
    """
    modules = ["picard", "java"]
    my_env = module_load(modules)
    # 
    picard="/usr/local/apps/picard/3.1.0/picard.jar"
    markdup_command = f"java -jar {picard} MarkDuplicates -I {bam} -O {output}.bam -M {output}.txt --REMOVE_DUPLICATES true"
    subprocess.call(markdup_command, shell = True, env = my_env)

def markdup_log(file):
    """
    Parse picard markdup output file 
    Input: 
    file: file to open 
    Output:
    pd.DataFrame of the log stats
    """
    log_list = list()
    with open(file, "r") as log_open:
        lines = log_open.readlines(); content = "".join(lines)
        stat_block = content.split("\n\n")[1]
        for line in stat_block.split("\n"):
            if not line.startswith("#"):
                log_list.append(line.split("\t"))
    # 
    df = pd.DataFrame(log_list[1], index = log_list[0], columns = ["stat"])
    return df

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
    modules = ["subread"]
    my_env = module_load(modules)
    count_command = f"featureCounts -p --countReadPairs -C -B -a {gtf} -F GTF --extraAttributes gene_name -s 0 -T {thread} {bam} -o {name}"
    print(count_command)
    subprocess.call(count_command, shell = True, env = my_env)

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
    modules = ["rsem/1.3.3", "STAR", "samtools"]
    my_env = module_load(modules)
    # module = " ".join(modules)
    # p = os.system(f"module load {module} | echo $PATH", shell = True)
    # p.communicate()
    # current_env = module_load(modules)
    rsem_command = f"rsem-calculate-expression --paired-end --star --strandedness none -p {thread} "
    if r1.endswith(".gz"):
        rsem_command += " --star-gzipped-read-file"
    rsem_command += f" {r1} {r2} {ref_dir} {name}"
    subprocess.call(rsem_command, shell = True, env = my_env)

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
