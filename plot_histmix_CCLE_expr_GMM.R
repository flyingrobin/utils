library(ggplot2)
library(reshape2)
library(taigr)
library(dplyr)
library(readr)
library(useful)
library(tidyr)
library(gridExtra)
source("/Users/wangli/Documents/Broad/Utils/CCLE_expr_parser.R")

ccle.rpkm.log2  <- load.from.taiga(data.name = 'ccle-rnaseq-gene-expression-rpkm', data.version = 5)
list[ccle.rpkm.log2, gene.lookup] <- ccle_expr_cleaning_fromTaiga(ccle.rpkm.log2, log = T)

load("/Users/wangli/Documents/Broad/PRISM/RNA-seq/data/pcg_rmNA.GMM.RData")

plot_histmix <- function(gene, cell.line.group = NULL){
  # require ccle.rpkm.log2, med.gmm.pred and gene.lookup
  # Args:
  #   gene: gene ID
  
  # creae histogram mix plot
  df <- data.frame(value = ccle.rpkm.log2[gene, ], 
                   Gaussian.group = mx.G.pred[gene, 1:1019] %>% as.factor(),
                   Cell.line.group = colnames(ccle.rpkm.log2) %in% cell.line.group)
  
  ggplot(df, aes(x = value)) + 
    geom_histogram(aes(y = ..count.. *8), bins = 100) + 
    geom_density(aes(y = ..count.., fill = Gaussian.group, color = Gaussian.group), alpha = 0.4) + 
    xlab("log2 RPKM") + ylab("scaled count") + ggtitle(gene.lookup[gene]) +
    geom_rug() + 
    theme(legend.position = "bottom") + guides(fill = guide_legend(nrow=1,byrow=TRUE))

}
