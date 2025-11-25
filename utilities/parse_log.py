import pandas as pd 

def log_parse(log_file):
    """
    Parse gene_tag.py log file
    it is strictly match to log format (when changes in log occurs, double check its applicability).
    """
    print(log_file)
    total_read, clean_read, align_unique, feature_align, frame = 0,0,0,0,0
    frame0_g, frame1_g, frame2_g = 0, 0, 0
    frame0_r, frame1_r, frame2_r = 0, 0, 0
    with open(log_file, "r") as fr:
        lines = fr.readlines()
        for line in lines:
            line = line.strip("\n")
            if "\t" in line:
                obj, value = line.split("\t")
                if "Total read" in obj:
                    total_read = value
                    continue
                if "Clean Read" in obj:
                    clean_read = value
                    continue
                if "Align unique" in obj:
                    align_unique = value 
                    continue
                if "CDS alignment" in obj:
                    feature_align = value 
                    continue
                if "Target frame" in obj:
                    frame = value
                    continue
                if "Frame0 gene" in obj:
                    frame0_g = value 
                    continue
                if "Frame1 gene" in obj:
                    frame1_g = value 
                    continue
                if "Frame2 gene" in obj:
                    frame2_g = value
                    continue
                if "Frame0 Read" in obj:
                    frame0_r = value 
                    continue
                if "Frame1 Read" in obj:
                    frame1_r = value 
                    continue
                if "Frame2 Read" in obj:
                    frame2_r = value
                    continue
    return (total_read, clean_read, align_unique, feature_align, frame, frame0_g, frame0_r, frame1_g, frame1_r, frame2_g, frame2_r)

def fastp_log_parser(log):
    """
    parse fastp log 

    Input:
    log: log file to read (from fastp)
    Return:
    log info: dictionary (total_reads, trimmed_reads)
    """
    log_info = dict()
    with open(log, "r") as log_open:
        for line in log_open.readlines():
            line = line.strip("\n")
            if line.startswith("total reads"):
                total_read = int(line.split(": ")[-1])
                log_info["total_reads"] = total_read
                break
    with open(log, "r") as log_open:
        for line in log_open.readlines():
            line = line.strip("\n")
            if line.startswith("total reads"):
                total_read = int(line.split(": ")[-1])
                log_info["total_reads(after filtering)"] = total_read
            # if line.startswith("reads with adapter trimmed"):
            #     trimmed_read = int(line.split(": ")[-1])
            #     log_info["trimmed_reads"] = trimmed_read
    log_info["trim_yield"] = log_info["total_reads(after filtering)"]/log_info["total_reads"]
    return log_info

def mark_log_parser(log):
    """
    parse picard markdup metric 

    Input: 
    log: log file to read
    Return: 
    log info: dictionary (pass_markdup_reads, dup_rate)
    """
    i = -2
    log_info = dict()
    with open(log, "r") as log_open: 
        for line in log_open.readlines():
            line = line.strip("\n")
            if i == -1:
                i += 1
            if i == 0:
                pair_read = line.split("\t")[2]
                dup_rate = line.split("\t")[-2]
                log_info["pass_markdup_reads"] = int(pair_read)
                log_info["dup_rate"] = float(dup_rate)
                # log_info["dup_yield"] = 1-log_info["dup_rate"]
                break
            if line.startswith("LIBRARY"):
                i += 1
    return log_info

def markdup_log_parser(file):
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
    df = pd.DataFrame(log_list[1], index = log_list[0], columns = ["stat"])
    return df

def flagstat_tsv_parse(f):
    map_info = dict()
    df = pd.read_table(f, sep = "\t", header = None, names = ["pass", "fail", "category"])
    total_mapped = df.query("category == 'primary mapped'")["pass"].item()
    map_info["aligned"] = total_mapped
    total = df.query("category == 'total (QC-passed reads + QC-failed reads)'")["pass"].item()
    map_info["total"] = total
    return map_info 
    # r1 = df.query("category == 'read1'")["pass"].values[0]
    # r2 = df.query("category == 'read2'")["pass"].values[0]
    # trans = df.query("category == 'with mate mapped to a different chr'")["pass"].values[0]
    # cis = (total_mapped-trans)//2*0.95
    # return cis 

def flagstat_parser(log):
    """
    parse samtools flagstat output (PE alignment)
    """
    log_info = dict()
    with open(log, 'r') as log_open:
        for line in log_open.readlines():
            line = line.strip("\n")
            if "read1" in line:
                input_read = int(line.split(" +")[0])
                log_info["total"] = input_read
            if "with itself and mate mapped" in line:
                align_pair = int(line.split(" +")[0])//2
                log_info["aligned"] = align_pair 
    log_info["align_ratio"] = log_info["aligned"]/log_info["total"]
    return log_info

def flagstat_parser_SE(log):
    """
    Parse samtools flagstat output (SE alignment)
    """
    log_info = dict()
    with open(log, "r") as log_open:
        for line in log_open.readlines():
            line = line.strip("\n")
            if "in total" in line:
                total = int(line.split(" +")[0])
                log_info["total"] = total
            if "primary mapped" in line:
                aligned = int(line.split(" +")[0])
                log_info["aligned"] = aligned
        log_info["align_ratio"] = log_info["aligned"]/log_info["total"]
    return log_info 

def h2_log_parse(log):
    """
    Parse ldsc --h2 (heritability) output 
    """
    log_info = dict()
    with open(log, "r") as log_open:
        for line in log_open.readlines():
            line = line.strip("\n")
            if "" in line:
                h2 = line.split(":")
                h2_se = line.split(":")
    log_info["h2"] = h2 
    log_info["h2_se"] = h2_se
    return log_info 

def STAR_log(log):
    """
    Parse log from STAR aligner 

    Input: 
    log file from STAR aligner 

    Return: 
    a dictionary with total_input_reads, unique_align, unique_ratio
    """
    log_info = dict()
    with open(log, "r") as fr:
        for line in fr.readlines():
            line = line.strip()
            if "Number of input reads |" in line:
                input_r = line.split("Number of input reads |\t")[-1]
            if "Uniquely mapped reads number |" in line:
                unique_r = line.split("Uniquely mapped reads number |\t")[-1]
            if "Uniquely mapped reads % |" in line:
                unique_p = line.split("Uniquely mapped reads % |\t")[-1].strip("%")
    log_info["Total reads"] = input_r 
    log_info["Aligned reads"] = unique_r
    log_info["Aligned %"] = unique_p
    return log_info

