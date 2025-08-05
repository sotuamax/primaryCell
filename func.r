
# draw loops 
loop_arc_fun <- function(loops, start, end, row = F, col = F, score = F) {
  # input loops as df w/ columns : anchor1 & anchor2 & score (optional)
    # Map loop start/end to [0,1] relative positions along the annotation axis (absolute start/end)
    loops$a1 <- (loops$anchor1 - start)/(end-start)
    loops$a2 <- (loops$anchor2 - start)/(end-start)
    loops$c1 <- (loops$a2 - loops$a1)/3 + loops$a1
    loops$c2 <- (loops$a2 - loops$a1)/3*2 + loops$a1
    # Arc height (adjust as needed)
    h <- max(loops$a2 - loops$a1)
    loops$h <- (loops$a2 - loops$a1)/h
    loops$id <- seq(nrow(loops))
    
    if (row) {
      # Draw bezier curve as arc annotated on row (flip x/y and reverse axis coordinates)
      if (score) {
      for (i in seq(nrow(loops))) {
      grid.bezier(
        # start, control1, control2, end 
        y = unit(1-c(loops[loops$id == i, ][, c("a1", "c1", "c2", "a2")] |> unlist()), "npc"), 
        x = unit(1-c(0, loops[loops$id == i, "h"], loops[loops$id == i, "h"], 0), "npc"), 
        gp = gpar(lwd = loops[loops$id == i, "score"])
      )
      }} else {
      for (i in seq(nrow(loops))) {
      grid.bezier(
        y = unit(1-c(loops[loops$id == i, ][, c("a1", "c1", "c2", "a2")] |> unlist()), "npc"), 
        x = unit(1-c(0, loops[loops$id == i, "h"], loops[loops$id == i, "h"], 0), "npc")
      )
    }
      }
        } 
    if (col) {
      # Draw bezier curve as arc annotated on col (retain its own axis)
      if (score) {
      for (i in seq(nrow(loops))) {
      grid.bezier(
        x = unit(c(loops[loops$id == i, ][, c("a1", "c1", "c2", "a2")] |> unlist()), "npc"), 
        y = unit(c(0, loops[loops$id == i, "h"], loops[loops$id == i, "h"], 0), "npc"), 
        gp = gpar(lwd = loops[loops$id == i, "score"])
      )
      }} else {
      for (i in seq(nrow(loops))) {
      grid.bezier(
        x = unit(c(loops[loops$id == i, ][, c("a1", "c1", "c2", "a2")] |> unlist()), "npc"), 
        y = unit(c(0, loops[loops$id == i, "h"], loops[loops$id == i, "h"], 0), "npc")
      )
    }
      }
        }
}


# draw contact data as heatmap 
hic_pixel_heatmap <- function(mat, sample, chr, start, end, r, anno = NULL, upper = F, lower = F, anno = "row") {
  vec <- as.vector(mat[upper.tri(mat1, diag = F)])
  # get the maximum value
  if (lower) {
    mat[upper.tri(mat, diag = T)] <- NA  # remove upper triangle (show lower triangle)
  }
  if (upper) {
    mat[lower.tri(mat, diag = T)] <- NA # remove lower triangle (show upper triangle)
  }
  # color code 
  # "#F7D340","#C43C4E", "#FBB91F",
  # make sure (25% - 75% data show up w/ code #F17020)
  col <- c("grey98", "#FBB91F", "#F17020", "#F17020", "#F17020", "#D64B40", "#D64B40")
  col_value <- c(0, quantile(log2(vec+1), c(0.05, 0.25, 0.5, 0.75, 0.9, 0.95)))
  col_fun <- colorRamp2(col_value, col)
  # title 
  title <- paste0(sample, "", chr, ":", start/1e6, "-", end/1e6, "Mb")
  # add annotation 
  if (!is.null(anno)) {
    compartment_score <- anno[["compartment"]]
    if (anno == "row") {
      # for row annotation, the axis orientation is not match to matrix, reverse score here
      row_ha <- rowAnnotation(compartment = anno_barplot(-compartment_score, 
                                                  gp = gpar(fill = ifelse(-compartment_score < 0, "firebrick", "steelblue"), col = NA), 
                                                border = T, 
                                                baseline = 0, 
                                                add_numbers = F, 
                                                axis = F 
                                                ),
                              DHS = anno_barplot(-anno[["DHS"]], 
                                                gp = gpar(fill = "skyblue", col = NA), 
                                                border = T, 
                                                add_numbers = F, 
                                                axis = F
                                                
                              ))
    } else row_ha <- NULL
    if (anno == "col") {
      col_ha <- HeatmapAnnotation(compartment = anno_barplot(compartment_score, 
                                  gp = gpar(fill = ifelse(compartment_score > 0, "firebrick", "steelblue"), col = NA), 
                                  border = T, 
                                  baseline = 0, 
                                  add_numbers = F, 
                                  axis = F, 
                                  which = "column"),
                                DHS = anno_barplot(anno[["DHS"]], 
                                                   gp = gpar(fill = "skyblue", col = NA), 
                                                   border = T, 
                                                   add_numbers = F, 
                                                   axis = F
                                                   ))
    } else col_ha <- NULL
  }
  h <- ComplexHeatmap::Heatmap(
    log2(mat+1), 
    cluster_rows = F, 
    cluster_columns = F, 
    show_row_dend = F, 
    show_column_dend = F, 
    show_column_names = F, 
    show_row_names = F, 
    border = T, 
    col = col_fun, 
    na_col = "transparent",
    show_heatmap_legend = T,
    heatmap_legend_param = list(title = "log2(contactFreq+1)", 
                                title_position = "topcenter", 
                                direction = "horizontal"
                                ),
    width = unit(0.9, "snpc"), 
    height = unit(0.9, "snpc"), 
    row_title = title,
    left_annotation = row_ha,
    top_annotation = col_ha
    )
    return (h)
}

# retrieve HiC pixel data for specific region (chr, start, end)
retrieve_hic_pixel <- function(hic_file, chr, start, end, r, norm = "NONE") {
  # parameters 
  pm <- pgParams(
    assembly = "hg38", 
    chrom = chr, 
    chromstart = start, 
    chromend = end
  )
  # read hic into df (sparse data, only non-zero write into df)
  hic_df <- readHic(
  hic_file, 
  params = pm, 
  resolution = r, 
  norm = norm
  )
  colnames(hic_df) <- c("A", "B", "counts")
  # from df to matrix 
  # to change absolution position in the Matrix into relative position
  # only upper matrix data retained, to make it a matrix, to add upper + lower matrix
  # Do NOT add diagonal matrix data twice (for lower matrix, diagonal as 0)
  hic_mat_upper <- Matrix::sparseMatrix(
    i = (hic_df$A-start)/r + 1, 
    j = (hic_df$B-start)/r + 1, 
    x = hic_df$counts) |> as.matrix()
  hic_mat_lower <- Matrix::sparseMatrix(
    j = (hic_df$A-start)/r + 1,
    i = (hic_df$B-start)/r + 1, 
    x = hic_df$counts
  ) |> as.matrix()
  diag(hic_mat_lower) <- 0 # to prevent diagonal value get doubled 
  return (hic_mat_upper + hic_mat_lower)
}


# perform hypergeometric test 
hyper_enrich <- function(sample_success, sample_total, bg_success, bg_total) {
  x <- sample_success
  m <- bg_success
  n <- bg_total - bg_success 
  k <- sample_total 
  phyper(x, m, n, k, lower.tail = F)
  # Where,
  # x: The number of successes in the sample.
  # m: The number of items of interest in the population.
  # n: The number of items not of interest in the population.
  # k: The number of items drawn in the sample.
  # lower.tail = F (probability of X > x)
}

# GSEA result plot
GSEA_plot <- function(df, plot = "bar") {
  # input df have columns: name (names of term) & padj (significance value)
  # plot GSEA result using different styles 
  if (plot == "bar") {
    gsea_p <- ggplot(df) + geom_bar(aes(x = reorder(name, -log10(padj)), y = -log10(padj))) + theme_minimal() + coord_flip() + xlab(NULL)
  }
  if (plot == "circle") {
    # extra column: count (gene count in the term) & total (total gene for testing)
    df <- df |> mutate(size = Count/total)
    gsea_p <- ggplot(df) + geom_point(aes(x = reorder(name, -log10(padj)), y = -log10(padj), size = Count/total*100, color = -log10(padj))) + theme_minimal() + coord_flip() + xlab(NULL) + scale_color_viridis_c(option = "F", direction = -1)
  }
  return (gsea_p) 
}

barplot <- function(df) {
  # input df with columns in group + freq
  # add label of its proportion
  df$pos <- df$freq/2
  df$prop <- round(df$freq/sum(df$freq)*100, 2)
  df$lab <- paste0(df$group, "/", df$prop, "%")
  p <- ggplot(df) + geom_bar(aes(x = group, y = freq)) + geom_text(aes(x=group, y = pos, label = lab)) + theme_minimal()
  return (p)
}

# correlation analysis on sample replicates 
corr_func <- function(df, filter = NULL, method = "pearson") {
  cor_res <- list()
  if (method == "pearson") {
    df <- log2(df+1)
  } 
  if (!is.null(filter)) df <- df[rowSums(df > 1) > 0, ]
  # 
  cor_mat <- cor(df, method = method) # given a matrix/df, the coefficient value for every pair of variables is calculated
  cor_res$correlation.mat <- cor_mat
  cor_p <- pheatmap(cor_mat, width = 8, height = 8, cluster_rows = T, cluster_cols = T, fontsize = 5)
  cor_res$correlation.heat <- cor_p 
  # 
  cor_mat[lower.tri(cor_mat, diag = T)] <- NA
  cor_df <- reshape2::melt(cor_mat) |> filter(!is.na(value))
  cor_df <- cor_df |> dplyr::rename(c("sample1" = "Var1",  "sample2" = "Var2")) |> filter(sample1 != sample2)
  cor_res$correlation.df <- cor_df 
  cor_df$cell1 <- gsub("_.*", "", cor_df$sample1)
  cor_df$cell2 <- gsub("_.*", "", cor_df$sample2)
  replicates_cor <- cor_df |> filter(cell1 == cell2)
  cor_res$correlation.replicate <- replicates_cor
  return (cor_res)
}

# fisher's exact test for enrihment 
fisher_enrich <- function(df, padjust = T) {
  # df should have two feature columns 
  # each row is a variable (e.g., gene), which can be assigned into types in the two features' classification
  # recommend: in the given table, the two features are given as factors 
  df <- as.data.frame(df)
  tbl_freq <- table(df)
  # run chisqure test to see if correlation between the two features classification (chi-square test of independence)
  cq <- chisq.test(tbl_freq)
  cq_pvalue <- cq$p.value
  if (cq_pvalue < 0.01) {
    print("Two variables are dependent")
    # get the residual value (the differences between obs - exp)
    residualTbl <- cq$residuals
    observeTbl <- cq$observed 
    expectTbl <- cq$expected
  } else {
    print("Two variables are independent (H0)")
    return (NULL)
  }
  if (cq_pvalue < 0.01) {
    pTbl <- tbl_freq # p.value table 
    orTbl <- tbl_freq # odds ratio table 
    lowerci <- tbl_freq # lower CI 
    upperci <- tbl_freq # upper CI 
    rowCats <- levels(factor(as.vector(df[, 1])))
    colCats <- levels(factor(df[, 2]))
    # for each combination of features, use two sided fisher's exact test for enrichment test
    for (rc in rowCats) 
      for (cc in colCats) {
        rclass <- (df[, 1] == rc)
        cclass <- (df[, 2] == cc)
        # when fisher's test p.value significant, true odds ratio != 1 (there is differences in )
        fc_tmp <- fisher.test(rclass, cclass, alternative = "two.sided")
        CI95 <- fc_tmp$conf.int |> as.numeric()
        CI95_lower <- CI95[1]
        CI95_upper <- CI95[2]
        # 
        pTbl[rc, cc] <- fc_tmp$p.value
        orTbl[rc, cc] <- fc_tmp$estimate
        lowerci[rc, cc] <- CI95_lower
        upperci[rc, cc] <- CI95_upper
      }
  }
  enrich_list <- list()
  enrich_list[["observed"]] <- observeTbl
  enrich_list[["expect"]] <- expectTbl
  enrich_list[["odds.ratio"]] <- orTbl
  enrich_list[["residual"]] <- residualTbl
  enrich_list[["p"]] <- pTbl
  enrich_list[["CI.lower"]] <- lowerci 
  enrich_list[["CI.upper"]] <- upperci
  if (padjust) {
    padjTbl <- pTbl
    padj <- p.adjust(pTbl, method = "BH")
    i <- 0
    for (cc in colCats) 
      for (rc in rowCats) {
        i <- i+1
        padjTbl[rc, cc] <- padj[i]
      }
    enrich_list[["padj"]] <- padjTbl
    }
  return (enrich_list)
  # visualize enrichment as a heatmap
  # corrplot::corrplot(enrich_list$odds.ratio, p.mat=enrich_list$padj, insig="label_sig", sig.level=c(0.001,0.01,0.05), is.corr=F, method="color",pch.col="white",tl.col="black")
  # visualize enrichment as a barplot 
  # or_df <- d_res2$odds.ratio |> as.data.frame() |> filter(disease_gene == "D")
  # or_df$CI_L <- d_res2$CI.lower |> as.data.frame() |> filter(disease_gene == "D") |> pull(Freq)
  # or_df$CI_U <- d_res2$CI.upper |> as.data.frame() |> filter(disease_gene == "D") |> pull(Freq)
  # ggplot(or_df) + geom_bar(stat = "identity", aes(x = loop_g, y = Freq), fill = "gold") + xlab("Loop number") + ylab("Fold enrichment of disease genes") + theme_minimal() + geom_errorbar(aes(ymin = CI_L, ymax = CI_U, x = loop_g), width = 0.3)
}

fisher_pipe <- function(tf_df) {
  tf_query <- tf_df |> select(-c(chrom, start, end,name))
  n <- length(colnames(tf_query)[-1])
  tf_enrich <- data.frame(obs = rep(0, n), exp=rep(0, n), oddsr = rep(0, n), p = rep(1, n), row.names = colnames(tf_query)[-1])
  for (tf in colnames(tf_query)[2:n]) {
    fish_res <- fisher_enrich(tf_query |> select(all_of(c("query", tf))), padjust = F)
    if (!is.null(fish_res)) {
      tf_enrich[tf, "obs"] <- fish_res$observed[2,2]
      tf_enrich[tf, "exp"] <- fish_res$expect[2,2]
      tf_enrich[tf, "oddsr"] <- fish_res$odds.ratio[2,2]
      tf_enrich[tf, "p"] <- fish_res$p[2,2]
    }
  }
  return (tf_enrich)
}

# t-statistic for the identification of sample-specific features (when sample size is large, and follow normality assumption)
t_stat_comp <- function(count_mat, group_sample, group_name) {
  # input count_mat, column for sample, row for peak/gene/loop
  # t: (X – μ) / [ s/√(n) ]
  # X: sample mean 
  # μ: population mean 
  # s: sample sd 
  # n: sample size 
  sample_mat <- count_mat |> select(all_of(c(group_sample, "mean"))) 
  print(head(sample_mat))
  t_return <- apply(sample_mat, 1, function(x) {
    t_res <- t.test(x[-length(x)], mu = x[length(x)], alternative = "greater")
    return(t_res$p.value)
  })
  count_mat[, paste0("p.", group_name)] <- t_return
  return (count_mat)
  # df_mean <- apply(df, 1, mean)
  # df_ssd <- (apply(sweep(df, 1, df_mean, "-")^2, 1, sum)/(ncol(df)-1))^0.5 # sample SD 
  # df_SE <- df_ssd/(ncol(df)^0.5) # standard error 
  # df_t <- sweep(sweep(df, 1, df_mean, "-"), 1, df_SE, "/")
  # return (df_t)
}

# welch two-sample t.test (compare means of two independent groups when variance are unequal and sample sizes may differ)
welch_stat_comp <- function(count_mat, group_sample, compare_sample, group_name) {
  sample_bool <- names(count_mat) %in% group_sample # bool vector 
  compare_bool <- names(count_mat) %in% compare_sample # bool vector 
  wilcox_return <- apply(count_mat, 1, function(x) {
    wilcox_res <- suppressWarnings(t.test(x[sample_bool], x[compare_bool], alternative = "g"))
    return(wilcox_res$p.value)
  })
  count_mat[, paste0("p.", group_name)] <- wilcox_return
  return (count_mat)
}

# Wilcoxon test for sample-specific features (when sample size is small and t test assumptions violated)
wilcox_stat_comp <- function(count_mat, group_sample, compare_sample, group_name) {
  # wilcox is a ranked test, and when ties, it cannot give exact p-value 
  # when value differences for the top ranks, it cannot differentiate (only ranks considered)
  sample_bool <- names(count_mat) %in% group_sample # bool vector 
  compare_bool <- names(count_mat) %in% compare_sample # bool vector 
  wilcox_return <- apply(count_mat, 1, function(x) {
    wilcox_res <- suppressWarnings(wilcox.test(x[sample_bool], x[compare_bool], alternative = "g"))
    return(c(p=wilcox_res$p.value, h1=mean(x[sample_bool]), h0=mean(x[compare_bool])))
  })
  count_mat[, paste0("p.", group_name)] <- wilcox_return["p", ]
  count_mat[, paste0("mean.", group_name)] <- wilcox_return["h1", ]
  return (count_mat)
}

# perform padjust on each column 
p_adjust_col <- function(df) {
  p_col <- grep("^p.", colnames(df), value = T)
  for (col in p_col) {
    df[, gsub("^p.", "padj.", col)] <- p.adjust(df |> pull(col), method = "BH")
  }
  return (df)
}

# pick columns where its p-value is significant and return in a new column 
pick_group <- function(padj_df) {
  # get padj columns 
  padj_col <- grep("^padj", colnames(padj_df), value = T) 
  padj_sub <- padj_df |> select(all_of(padj_col)) 
  # log transform padj
  padj_log <- -log10(padj_sub)
  # ratio of log(padj) 
  padj_ratio <- apply(padj_log, 1, function(x) x/sum(x)) |> t() |> as.data.frame()
  colnames(padj_ratio) <- gsub("padj.", "", colnames(padj_ratio))
  # dominant padj as the assigned cluster
  padj_ratio$cluster <- apply(padj_ratio, 1, function(x) {
    matched <- unlist(names(x)[x > 0.9][1])
    if (length(matched) == 0) "." else matched })
  # how many clusters significant here
  padj_ratio$n <- rowSums(padj_sub < 0.1)
  # get peak and its best assigned cluster
  padj_cluster <- padj_ratio |> filter(n > 0) |> filter(cluster != ".") |> rownames_to_column("peak")
  padj_cluster <- padj_df |> rownames_to_column("peak") |> left_join(padj_cluster[, c("peak", "cluster")], by = "peak")
  return (padj_cluster)
}

# based on p-value, to get the top # peaks for specific cluster 
group_specific <- function(group_df, group, top) {
  cluster_df <- group_df[group_df$cluster == group, ]
  # score is weighted by -log(padj) times signal of mean in the group
  cluster_df$score <- -log10(cluster_df[[paste0("padj.", group)]])*cluster_df[[paste0("mean.", group)]]
  top_df <- cluster_df |> arrange(desc(score)) |> head(top)
  # return the top # features 
  return (top_df)
}


# pca function, return pca summary 
pca_fun <- function(df, filter = FALSE, top = FALSE, scale = FALSE, log = FALSE) {
  if (filter) {
    print("Filtering ...")
    df <- df[rowSums(df > filter) > 0, ]} # at least one feature have value > 1
  if (top) {
    # select top n most variable features 
    df_var <- rowSds(as.matrix(df)) |> as.data.frame() |> rename_with(~c("sd"), 1) |> arrange(desc(sd))
    df_var <- df_var[1:top, , drop = F]
    df <- df[rownames(df_var), ]
  }
  if (scale) {
    # scale features to the same scale
    print("Scaling ...")
    res <- prcomp(t(df), center = T, scale. = T)
  } 
  if (log) {
    # or instead, use log2 to scale the value
    df_scaled <- t(log2(df+1)) # apply(t(htseq_exp), 2, scale)
    res <- prcomp(df_scaled)
  } else {
    res <- prcomp(df)
  }
  summ <- summary(res)
  return (summ)
}

# pca distance (weighted by variance explained per component), return pca distance df 
pca_dist <- function(pca_res, pc_variance = 0) {
  prop_var <- t(pca_res$importance) |> as.data.frame()
  colnames(prop_var) <- c("sd", "var", "cumulate")
  pc_var <- prop_var[prop_var$var > pc_variance, ]
  # get pc-sample matrix into a dataframe
  sample_pc_pos <- t(pca_res$x) |> as.data.frame()
  # 
  pca_prop_val <- sample_pc_pos[rownames(pc_var), ] * pc_var[, "var"]
  pca_dist <- dist(t(pca_prop_val), upper = T, diag = T)
  return (pca_dist)
}

# pca plot on pca summary result
pca_plot <- function(pca_res, style = "point") {
  if (style == "point") {
    p <- ggplot(data = pca_res$x, aes(x = PC1, y = PC2)) + geom_point(shape = 1, size = 3) + ggrepel::geom_text_repel(aes(label = rownames(pca_res$x)), size = 4) + theme(legend.position = "none", panel.background = element_rect(fill = "transparent", color = "grey"), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) + xlab(paste0("PC1: ", round(t(pca_res$importance)["PC1", 2]*100, 1), "%")) + ylab(paste0("PC2: ", round(t(pca_res$importance)["PC2", 2]*100, 1), "%"))
    return (p)
  } 
  if (style == "heatmap") {
    prop_var <- t(pca_res$importance)[, 2]
    pca_prop_val <- t(pca_res$x) * prop_var
    pca_dist <- dist(t(pca_prop_val), upper = T, diag = T)
    color_scale <- circlize::colorRamp2(seq(0, max(pca_dist), by = max(pca_dist)%/%5), 
                                    c("#4575B4","#B8D9E9","#FFFFBF","#FDB674","#D73027", "red"))
    p <- Heatmap(as.matrix(pca_dist), show_column_names = F, name = "PCA distance", heatmap_legend_param = list(title_gp = gpar(rot = 90)), row_title = NULL, column_title = NULL, col = color_scale)
    return (p)
  }
}

# RPM normalization (read per million)
RPM <- function(count_df) {
  count_total <- apply(count_df, 2, sum)
  count_norm <- sweep(count_df, 2, count_total, FUN = "/")*1000000
  return (count_norm)
}

# z-score based specific signal 
z_specificity <- function(count_norm, p_cutoff) {
  count_z <- t(apply(count_norm, 1, scale))
  colnames(count_z) <- colnames(count_norm)
  # get p-value for z-score
  z_p <- apply(count_z, 1, function(x) pnorm(x, lower.tail = F)) |> t()
  z_sig_p <- z_p[rowSums(z_p < p_cutoff) > 0, ]
  z_sig_p[z_sig_p > p_cutoff + 0.1] <- 1
  z_sig_p_log <- -log2(z_sig_p) 
  z_p_sum <- apply(z_sig_p_log, 1, sum) 
  # get weighted-p strength 
  z_sig_ratio <- sweep(z_sig_p_log, 1, z_p_sum, "/")
  z_specific <- z_sig_ratio[rowSums(z_sig_ratio > 0.9) == 1, ]
  # get the colnames of the maximum weighted-p value
  z_specific_col <- apply(z_specific, 1, function(x) colnames(z_specific)[x > 0.9]) 
  z_specific_df <- data.frame(s = I(z_specific_col))
  return (z_specific_df)
}


