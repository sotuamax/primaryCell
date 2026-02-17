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
    parser <- ArgumentParser(description = "Usage: Rscript GSEA.R -gmt ann.gmt -gene target_gene.txt -symbol gene_symbol.txt -bg background_gene.txt -O gene_enrich_res")
    # parser$add_argument("-label", "--gene_set_label", choices = c("H", "CP", "CGP", "MIR", "TFT", "BP", "CC", "MF", "PO", "CT", "OS", "IMMU", "CA"), default = "CP", help = paste0('gene set label for enrichment analysis.', label_description))
    parser$add_argument("-gmt", "--gmt", default = "/data/jim4/data/Gene_set/human/h.all.v2025.1.Hs.symbols.gmt", help = "gmt for gene annotation terms ")
    parser$add_argument("-gene", "--gene", help = "file with genes of interest (w/ header and rowname as gene list)")
    parser$add_argument("-bg", "--background", required = F, help = "gene set background (w/ header and rowname as gene list)")
    parser$add_argument("-symbol", "--symbol", required = F, help = "if given gene set is not symbol, symbol file should be provided (gene-symbol no header).")
    parser$add_argument("-padj", "--padj", default = 0.01, help = "p cutoff value (default: 0.01)")
    parser$add_argument("-o", "--output", help = "output file prefix. ")
    parser$parse_args()
}

# parse arguments 
args <- parse_arg()
# label <- args$gene_set_label
gmt <- args$gmt
gene <- args$gene
bg <- args$background
symbol_file <- args$symbol
padj <- as.numeric(args$padj)
out <- args$output

cat("read gmt ... \n")
gmt_df <- read.gmt(gmt)
cat("read target gene set ...\n")
gene_df <- read.table(gene, sep = "\t", header = T, row.names = 1) 
gene_list <- row.names(gene_df)
# clean gene set based on background (if given)
if (!is.null(bg)) {
    cat("read background gene set ...\n")
    background_df <- read.table(bg, sep = "\t", header = T, row.names = 1)
    background_gene_list <- row.names(background_df)
    #head(background_gene_list)
}
# test if any gene overlap gmt_gene 
gmt_filtered <- gmt_df[gmt_df$gene %in% gene_list, ]
if (nrow(gmt_filtered) == 0) {
    cat("gene symbol required!\n")
    if (!is.null(symbol_file)) {
        cat("Read gene symbol ...\n")
        symbol_df <- read.table(symbol_file, sep = "\t", header = F, col.names = c("gene", "symbol"))
        gene_list <- symbol_df[symbol_df$gene %in% gene_list, ]$symbol
        if (!is.null(bg)) background_gene_list <- symbol_df[symbol_df$gene %in% background_gene_list, ]$symbol
    } else {
        cat("Provide gene symbol!")
        quit()
    }
} 

gmt_df <- gmt_df[gmt_df$gene %in% background_gene_list, ]
cat("Perform enrichment analysis ...\n")
cat(paste0("target gene: ", length(gene_list), "\n"))
cat(paste0("total background gene: ", nrow(gmt_df), "\n"))
cat(paste0("p.adjust:", padj, "\n"))
enrich <- enricher(gene_list, TERM2GENE = gmt_df, pvalueCutoff = padj, pAdjustMethod = "BH")

cat("Write significant enrichment ...\n")
enrich_sig <- enrich@result |> arrange(p.adjust) |> filter(p.adjust < padj) |> select(-(ID))
enrich_sig$Description <- tolower(enrich_sig$Description)
cat(paste0("significant enrichment:", nrow(enrich_sig), "\n"))
write.table(enrich_sig, file = paste0(out, ".tsv"), sep = "\t", quote = F, row.names = F)

#print(paste0("Background genes ", length(unique(background_df$gene)), ", ", length(unique(gmt_df$gene)), " with annotation"))

# print(paste0("Input unique genes ", length(unique(gene_df$symbol))))

# geneList <- gene_df |> select(gene_symbol, ranking) |> arrange(ranking) |> deframe()

# get gene set file by assigning key work for specific gene set
#anno_df <- read.table(file.path(anno_dir, "annotation.txt"), sep = "\t", header = F, col.names = c("file", "description", "label"))
# file_pattern <- anno_df[anno_df$label == label, "file"] 
# key_f <- list.files(anno_dir, pattern = file_pattern)


#gene_anno <- gene_df[gene_df$symbol %in% gmt_df$gene, ]
#print(paste0("Input with annotation ", length(unique(gene_anno$symbol))))

# term2gid and term2nm are dataframe, 
# term2gid with annotation term and gID (the database TERM2GENE table will be used as background)
# term2nm with annotation term and annotation name (optional)


# num_term <- length(unique(gmt_df$term))
# num_gene <- length(unique(gmt_df$gene))
# print(paste0(num_term, " gene sets and ", num_gene, " genes in the current selection"))
# overlap_prerank <- intersect(unique(gmt_df$gene), unique(geneList$gene_symbol))
# print(paste0(length(overlap_prerank), " (", round(length(overlap_prerank)/num_gene*100, 2), "%) overlapped between prerank genes and selected gene set"))
# enrich_result <- GSEA(geneList, TERM2GENE = gmt_df, minGSSize = 10, maxGSSize = 800, pvalueCutoff = padj, pAdjustMethod = "BH", nPerm = 5000)

