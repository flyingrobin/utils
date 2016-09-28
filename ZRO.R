library(taigr)
library(dplyr)
library(stringr)
library(devtools)
source_url("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")

setwd("/Users/wangli/Documents/Broad/")
source("Utils/fisherZ_transformation.R")

# This function import PRISM or any drug screen data matrix and CCLE gene expression matrix, and compute:
# (1) fisher Z transformed correlation of compound killing profile (z socre, viability or AUC) vs gene expression,
#     over the panel of cell lines
# (2) For a given compound, the rank of every correlation entry from (1), among all cell lines
# (3) For a given compound, the outliers of Z score, classified as outlier ( > 1.5IQR) and extreme outlier (>3 IQR)

# outlier detection
outlier_IQR <- function(x){
  # Args:
  #   x: one row/column of data matrix
  list[lower, upper] <- quantile(x, c(0.25,0.75), na.rm = T)
  iqr <- upper - lower
  breaks <- c(lower - 3* iqr, lower - 1.5 * iqr, upper + 1.5 * iqr, upper + 3 * iqr)
  # define -1 as lower end outliers and 1 as upper end outliers
  # define -2 as extreme lower end outliuer and 2 as extreme upper
  outlier <- ifelse(x < breaks[1], -2, ifelse(x < breaks[2], -1, 0)) +
    ifelse(x > breaks[4], 2, ifelse(x > breaks[3], 1, 0))
  outlier
}

ZRO <- function(mx, ccle.expr){
  # shared cell lines and excluede suspension cell lines
  ccl.common <- intersect(colnames(mx), colnames(ccle.expr))

  # fisher Z transformation of corr. coeff.
  z <- fisherZ(mx[, ccl.common] %>% t(), ccle.expr[, ccl.common] %>% t())

  # rank of each correlation coeff. for each compound
  rank <- apply(z, 1, function(x) rank(-abs(x))) %>% t()
  colnames(rank) <- colnames(z)

  # find outliers for each compound
  outlier <- apply(z, 1, outlier_IQR) %>% t()

  # return Z, R and O matrix
  list(z, rank, outlier)
}

# Z-Rank plot for a list of genes
#' @param gene.list input vector of genes, order of the gene.list represent the rank of genes. genes are column names of Z, R and O matrix 
#' @param z matrix of z 
#' @param rank matrix of rank
#' @param outlier matrix of outlier
#' @param gene.rank: ranking of the gene. if not provided the gene.list will be rank, starting with rank 1
#' @param cpd.group 
#' @param gene.lookup named vector to switch geneID with gene symbol back and forth
#' @param return ggplot of ZRank and legend separately

plot_ZRO <- function(gene.list, z, rank, outlier, 
                     gene.rank = gene.list, cpd.group, 
                     gene.lookup = NULL, repurposing = FALSE){
  source("~/Documents/Broad/Utils/grid_arrange_shared_legend.R")
  # aesthetic settings
  labels <- c('-2' = 'Extreme Lower', 
              '-1' = 'Lower', 
              '0' = 'None', 
              '1' =  'Upper', 
              '2'  = 'Extreme upper')
  colors <-  c('-2' = 'orange', 
               '-1' = 'blue', 
               '0' = 'black', 
               '1' =  'red', 
               '2'  = 'green')
  sizes <-  c('TRUE' = 2.5, 'FALSE' = 1)
  shapes <- c('TRUE' = 17, 'FALSE' = 19)
  
  # ZRank plot list
  g.ZRank.list <- lapply(gene.list, function(gene){
    # generate a legend with all outlier category
    df <- data.frame(z = z[, gene],
                     rank = rank[, gene],
                     outlier = outlier[, gene] %>% as.factor(),
                     TLCD1.cpds = rownames(rank) %in% cpd.group)
    # plot title: row1 gene/cpd name; row2 outlier ?out of 8, rank
    if(repurposing){
      title.row1 <- get_ZRankplot_title_REP(gene, gene.lookup)
    }else{
      title.row1 <- paste(gene, gene.lookup[gene])
    }
    
    title.row2 <- paste(sum(df[df$outlier != '0', 'TLCD1.cpds']), '/8, rank = ', which(gene.rank == gene), sep = '')
    ggplot(df, aes(x = z, y = rank)) +
      geom_point(aes(color = outlier,  size = TLCD1.cpds, shape = TLCD1.cpds)) + 
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )  + 
      scale_color_manual(labels= labels, values = colors) +
      scale_size_manual(values = sizes) +
      scale_shape_manual(values = shapes) +
      ggtitle(paste(title.row1, title.row2, sep = '\n')) + 
      theme(legend.position="none")
  })
  
  # generate legend
  g.legend <- ggplot(data.frame(x = rnorm(10), y = rnorm(10), outlier = as.factor(c(-2,-1,0,1,2)) ,TLCD1.cpds = c(TRUE, FALSE)), aes(x, y)) + geom_point(aes(color = outlier, size = TLCD1.cpds, shape = TLCD1.cpds)) + scale_color_manual(labels= labels, values = colors) +
    scale_size_manual(values = sizes) +
    scale_shape_manual(values = shapes) 
  legend <- gg_legend_extract(g.legend, legend.panel = F)
  list(g.ZRank.list, legend)
}
