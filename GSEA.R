#!/usr/bin/env Rscript
suppressMessages(library("tidyverse"))
suppressMessages(library("argparse"))
suppressMessages(library("clusterProfiler"))


label_description <- "
    H: Hallmark (well-defined biological states or processes and display coherent expression);
    CP: canonical pathway; 
    CGP: chemical and genetic pertubations;
    MIR: microRNA target;
    TFT: transcription factor target;
    BP: GO biological process ontology;
    CC: GO cellular component ontology
    MF: GO molecular funciton ontology
    PO: phenotype ontology;
    CT: cell types from sc-RNA;
    OS: oncogenic signature (human only);
    IMMU: chemical and genetic perturbation of the immune system (human only);
    CA: cancer-cell atlas (human only);
    ALL: for all gene sets (not recommend)
    "


parse_arg <- function() {
    parser <- ArgumentParser(description = "Usage: Rscript GSEA.R -h")
    # parser$add_argument("-label", "--gene_set_label", choices = c("H", "CP", "CGP", "MIR", "TFT", "BP", "CC", "MF", "PO", "CT", "OS", "IMMU", "CA"), default = "CP", help = paste0('gene set label for enrichment analysis.', label_description))
    parser$add_argument("-gmt", "--gmt", help = "gmt for gene annotation terms ")
    parser$add_argument("-gene", "--gene", help = "file with genes of interest")
    parser$add_argument("-symbol", "--symbol", required = F, help = "if given gene set is not symbol, symbol file should be provided.")
    parser$add_argument("-bg", "--background", required = F, help = "gene set background")
    parser$add_argument("-padj", "--padj", default = 0.01, help = "p cutoff value (default: 0.01)")
    parser$add_argument("-O", "--output", help = "output file prefix. ")
    parser$parse_args()
}

# parse arguments 
args <- parse_arg()
# label <- args$gene_set_label
gmt <- args$gmt
gene <- args$gene
padj <- as.numeric(args$padj)
out <- args$output

# read gene set 
gene_df <- read.table(gene, sep = "\t", header = T) |> rename_with(~c("gene"), 1)
# print(paste0("Input unique genes ", length(unique(gene_df$gene))))

if (!is.null(args$symbol)) {
    symbol_df <- read.table(args$symbol, sep = "\t", header = F, col.names = c("gene", "symbol"))
    gene_df <- inner_join(gene_df, symbol_df, by = "gene")
} else {
    gene_df <- gene_df |> rename_with(~c("symbol"), 1)
}
print(paste0("Input unique genes ", length(unique(gene_df$symbol))))

# geneList <- gene_df |> select(gene_symbol, ranking) |> arrange(ranking) |> deframe()

# get gene set file by assigning key work for specific gene set
#anno_df <- read.table(file.path(anno_dir, "annotation.txt"), sep = "\t", header = F, col.names = c("file", "description", "label"))
# file_pattern <- anno_df[anno_df$label == label, "file"] 
# key_f <- list.files(anno_dir, pattern = file_pattern)
gmt_df <- read.gmt(gmt)

# clean gene set based on background (if given)
if (!is.null(args$background)) {
    background_df <- read.table(args$background, sep = "\t", header = T)
    gmt_df <- gmt_df[gmt_df$gene %in% background_df$symbol, ]
    print(paste0("Background genes ", length(unique(background_df$gene)), ", ", length(unique(gmt_df$gene)), " with annotation"))
    # gmt_bg <- inner_join(gmt_df, background_df, by = "gene")
    # print(paste0("Actual background genes (with term) ", length(unique(gmt_df$gene))))
} else {
    print(paste0("Annotated genes ", length(unique(gmt_df$gene))))
}

gene_anno <- gene_df[gene_df$symbol %in% gmt_df$gene, ]
print(paste0("Input with annotation ", length(unique(gene_anno$symbol))))

# term2gid and term2nm are dataframe, 
# term2gid with annotation term and gID (the database TERM2GENE table will be used as background)
# term2nm with annotation term and annotation name (optional)
enrich <- enricher(gene_df$symbol, TERM2GENE = gmt_df, pvalueCutoff = padj, pAdjustMethod = "BH")

# num_term <- length(unique(gmt_df$term))
# num_gene <- length(unique(gmt_df$gene))
# print(paste0(num_term, " gene sets and ", num_gene, " genes in the current selection"))
# overlap_prerank <- intersect(unique(gmt_df$gene), unique(geneList$gene_symbol))
# print(paste0(length(overlap_prerank), " (", round(length(overlap_prerank)/num_gene*100, 2), "%) overlapped between prerank genes and selected gene set"))
# enrich_result <- GSEA(geneList, TERM2GENE = gmt_df, minGSSize = 10, maxGSSize = 800, pvalueCutoff = padj, pAdjustMethod = "BH", nPerm = 5000)

enrich_sig <- enrich@result |> arrange(p.adjust) |> filter(p.adjust < padj) |> select(-(ID))

print(paste0(nrow(enrich_sig), " significant enrichment with padj ", padj))
enrich_sig$Description <- tolower(enrich_sig$Description)
write.table(enrich_sig, file = paste0(out, ".txt"), sep = "\t", quote = F, row.names = F)
