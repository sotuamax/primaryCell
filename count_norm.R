#!/usr/bin/env Rscript
suppressMessages(library("tidyverse"))
suppressMessages(library("argparse"))
suppressMessages(library("DESeq2"))

parse_arg <- function() {
    parser <- ArgumentParser(description = "Usage: Rscript count_norm.R")
    parser$add_argument("-count", "--count", nargs = "+", help = "feature count file ")
    parser$add_argument("-norm", "--normalization", default = "DESeq2", help = "normalization methods", choices = c("RPM", "DESeq2"))
    parser$add_argument("-updatecol", "--updatecol", action = "store_true", help = "update column names")
    parser$add_argument("-o", "--output", help = "output name")
    parser$parse_args()
}

args <- parse_arg()
count <- args$count 

print("Read count file ...")
count_list <- list()
for (c in count) {
    feature_count <- read.table(c, sep = "\t", header = T, row.names = 1, comment.char = "#", check.names = F) 
    if ("Chr" %in% colnames(feature_count)) {
        feature_mat <- feature_count |> select(-c(Chr,Start,End,Strand,Length))
    } else {feature_mat <- feature_count}
    count_list[[c]] <- feature_mat
}

names(count_list) <- NULL
feature_mat <- do.call(cbind, count_list)

# feature_size <- feature_count |> select(c("Length"))


# remove gene_name if ever exists 
if ("gene_name" %in% colnames(feature_mat)) {
    rownames(feature_mat) <- paste0(rownames(feature_mat), "|", feature_mat$gene_name)
    feature_mat <- feature_mat |> select(-c(gene_name))
}

if (args$updatecol) {
    print("Update column values ... ")
    old_col <- colnames(feature_mat)
    new_col <- tools::file_path_sans_ext(basename(old_col))
    colnames(feature_mat) <- new_col
}

# get RPM function 
RPM <- function(count_df) {
  count_total <- apply(count_df, 2, sum)
  count_norm <- sweep(count_df, 2, count_total, FUN = "/")*1000000
  return (count_norm)
}


if (args$normalization == "RPM") {
    colnames(feature_mat) <- basename(colnames(feature_mat))
    colnames(feature_mat) <- gsub(".bam", "", colnames(feature_mat))
    # feature_mat <- feature_mat[rowSums(feature_mat < 8) < ncol(feature_mat), ]
    # 
    # feature_mat |> write.table(paste0(args$output, ".txt"), sep = "\t", quote = F, row.names = T)
    feature_rpm <- RPM(feature_mat)
    print("Write normalized count ...")
    feature_rpm |> write.table(paste0(args$output, "_RPM.txt"), sep = "\t", quote = F, row.names = T)
}

# if (args$normalization == "RPKM") {
#     print("RPKM is not developed yet. ")
    # exit(1)
    # sample_df <- data.frame(sample = colnames(feature_count)) 
    # dds <- DESeqDataSetFromMatrix(countData = feature_count, colData = sample_df, design = ~1)
    # dds <- estimateSizeFactors(dds)
    # # get normalized count (normalized by sequencing depth)
    # feature_normalized <- counts(dds, normalized = T)
    # write.table(feature_normalized, file = gsub(".txt", "_norm.txt", count), sep = "\t", quote = F, row.names = T)
#}

# c1 <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(feature_count))
# c2 <- gsub(".*\\.", "", c1)
# reformat feature count data and write into a matrix format 
# colnames(feature_count) <- c1
# write.table(feature_count, file = gsub(".txt", "_mat.txt", count), sep = "\t", quote = F, row.names = T)

if (args$normalization == "DESeq2") {
    print("Perform normalization using DESeq2 ...")
    print('Write raw count ...')
    feature_mat |> write.table(paste0(args$output, "_mat.txt"), sep = "\t", quote = F, row.names = T)
    colname <- colnames(feature_mat)
    col_sample <- data.frame(sample = colname)
    print("Write samples ...")
    col_sample |> write.table(paste0(args$output, "_sample.txt"), sep = "\t", quote = F, row.names = F)
    dds <- DESeqDataSetFromMatrix(countData = feature_mat, colData = col_sample, design = ~ 1)
    dds <- estimateSizeFactors(dds)
    norm_cts <- counts(dds, normalized = T)
    print("Write normalized count ...")
    norm_cts |> write.table(paste0(args$output, "_DESeq2.txt"), sep = "\t", quote = F, row.names = T)
}

# normalize by size (this is specifically for peak, not gene, gene size is not easily scaled)
# sizefactor <- feature_size$Length/100 
# feature_normby_size <- sweep(feature_normalized, 1, sizefactor, "/")
# write.table(feature_normby_size, file = gsub(".txt", "_norm_size.txt", count), sep = "\t", quote = F, row.names = T)

