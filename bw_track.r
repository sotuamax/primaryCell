#!/usr/bin/env Rscript
suppressMessages(library("tidyverse"))
suppressMessages(library("argparse"))
suppressMessages(library("plotgardener"))
suppressMessages(library("RColorBrewer"))

parse_arg <- function() {
    parser <- ArgumentParser(description = "Usage: Rscript GSEA.R -h")
    # parser$add_argument("-label", "--gene_set_label", choices = c("H", "CP", "CGP", "MIR", "TFT", "BP", "CC", "MF", "PO", "CT", "OS", "IMMU", "CA"), default = "CP", help = paste0('gene set label for enrichment analysis.', label_description))
    parser$add_argument("-track", "--track", nargs = "+", help = "track file in bigwig format")
    parser$add_argument("-track_col", "--track_color", nargs = "+", help = "assign track color with number to indicate group of samples")
    parser$add_argument("-chrom", "--chrom", help = "chromosome to plot")
    parser$add_argument("-start", "--start", help = "start position to plot")
    parser$add_argument("-end", "--end", help = "end position to plot")
    parser$add_argument("-extend", "--extend", help = "extend range around the start and end")
    parser$add_argument("--assembly", "--assembly", help = "genome assembly to use for gene model", choices = c("hg38", "mm39"))
    parser$add_argument("-o", "--output", help = "output name for pdf plot")
    parser$add_argument("-width", "--width", help = "page width")
    parser$add_argument("-height", "--height", help = "page height")
    parser$parse_args()
}

# parse arguments 
args <- parse_arg()
track <- args$track
track_color <- args$track_color
assembly <- args$assembly 
chrom <- args$chrom
start <- as.integer(args$start)
end <- as.integer(args$end)
output <- args$output

# prepare track color
uniq_col <- unique(track_color)
pick_color <- brewer.pal(length(uniq_col), "Dark2")
sample_df <- data.frame(track = track, label = track_color)
col_df <- data.frame(label = uniq_col, col = pick_color)
sample_col <- sample_df |> left_join(col_df, by = "label")

# prepare assembly 
if (assembly == "mm39") {
    assembly <- assembly("mm39", TxDb = "TxDb.Mmusculus.UCSC.mm39.refGene", OrgDb = "org.Mm.eg.db")
}

# determine signal range 
score_max <- c()
for (f in track) {
    bw_read <- readBigwig(f, params = pm)
    score_max <- c(score_max, bw_read$score |> max())
}
score_max <- score_max |> max()

# create page 
w <- as.integer(args$width)
h <- as.integer(args$height)
pdf(paste0(output, ".pdf"), width = w, height = h)
pm <- pgParams(assembly = assembly, chrom = chrom, chromstart = start, chromend = end)
pageCreate(
    width = w, height = h, default.units = "inches",
    showGuides = F, xgrid = 0, ygrid = 0,
    params = pm
)
# plot region
side_empty <- 0.15
plot_w <- w-2*side_empty
plotGenomeLabel(params = pm, x = side_empty, y = 0.01, length = plot_w, scale = "Kb")
plotGenes(params = pm, x = side_empty, y = 0.05, width = plot_w, height = 0.7, fontcolor = c("#669fd9", "black"), fill = c("#669fd9", "black"))

j <- 0
for (f in track) {
    bw_read <- readBigwig(f, params = pm)
    panel_h <- round((h-2)/length(track), 2)
    fill_c <- sample_col |> filter(track == f) |> pull(col)
    plotSignal(bw_read, range = c(0,round(score_max * 0.6)), scale = T, linecolor = NA, fill = fill_c, label = strsplit(f, ".bw") |> unlist(), params = pm, x = side_empty, y = 1+(panel_h+0.005)*j, width = plot_w, height = panel_h)
    plotRect(x = side_empty, 1+(panel_h+0.005)*j, width = plot_w, height = panel_h)
    j <- j+1
}


dev.off()