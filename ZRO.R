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

# plotting functions of ZRO
# beeswarm <- beeswarm(boxplot.stats(z.auc2expr[1, ])$out)
# ggplot(data = data.frame(z = z.auc2expr[1, ]), aes(x = 1, y = z)) +
#   stat_boxplot(geom ='errorbar', width = 0.5) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(data = beeswarm, aes(x = x, y = y)) +
#   ggtitle(rownames(z.auc2expr)[1]) + xlab("")
