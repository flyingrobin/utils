library(taigr)
library(dplyr)
library(useful)
library(stringr)
library(devtools)
source_url("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")


ccle.expr_log_transform <- function(mx.ccle, min.rpkm = -7){
  mx.ccle <- log2(mx.ccle)
  mx.ccle[is.infinite(mx.ccle)] <- NA
  mx.ccle[mx.ccle < min.rpkm] <- min.rpkm
  mx.ccle
}

# parse ccle expression mx imported from Taiga
ccle_expr_cleaning_fromTaiga <- function(mx.ccle, log = FALSE, min.rpkm = -7){
  # Args:
  #   mx.ccle: ccle gene expression matrix
  #   log: boolean indicating whether to log2 transform the expression matrix
  #   min.rpkm: if log is TRUE, convert -Inf entries to NA and cap min of log2(rpkm) as -7

  colnames(mx.ccle) <- gsub(')', '', colnames(mx.ccle))
  gene.lookup <- str_split_fixed(colnames(mx.ccle), ' \\(', 2) %>% as.data.frame()
  gene.ids <- gsub('\\.[0-9]+$', '', gene.lookup$V2)

  gene.lookup <- gene.lookup$V1 %>% as.character()
  names(gene.lookup) <- gene.ids # gene lookup table as a vector with names
  
  # switch gene ID and gene symbol
  geneID.lookup <- names(gene.lookup)
  names(geneID.lookup) <- gene.lookup
  
  # combine the two
  gene.lookup <- c(gene.lookup, geneID.lookup)
  
  colnames(mx.ccle) <- gene.ids

  if(log){
    mx.ccle <- log2(mx.ccle)
    mx.ccle[is.infinite(mx.ccle)] <- NA
    mx.ccle[mx.ccle < min.rpkm] <- min.rpkm
  }

  # return transpose matrix (row as genes and column as cell lines) and gene lookup table
  list(t(mx.ccle), gene.lookup)
}
prism_expr_cleaning_fromTaiga <- function(mx.prism, log = FALSE, min.rpkm = -7){
  # Args:
  #   mx.ccle: ccle gene expression matrix
  #   log: boolean indicating whether to log2 transform the expression matrix
  #   min.rpkm: if log is TRUE, convert -Inf entries to NA and cap min of log2(rpkm) as -7
  
  mx.prism <- t(mx.prism)
  if(log){
    mx.prism <- log2(mx.prism)
    mx.prism[is.infinite(mx.prism)] <- NA
    mx.prism[mx.prism < min.rpkm] <- min.rpkm
  }
  mx.prism
}

ccle_expr_cleaning <- function(df.ccle){

  # gene id-symbol lookup table
  gene.lookup <- df.ccle$Description # save gene id/symbol columns to a new variable
  names(gene.lookup) <- gsub('\\.[0-9]+$', '', df.ccle$Name)

  # simplify ccle names
  colnames(df.ccle) <- gsub('fh_|\\-Tumor$','',colnames(df.ccle))

  df.ccle <- df.ccle[, c(-1, -2)] %>% as.matrix() # remove gene name columns and convert to numerical matrix
  rownames(df.ccle) <- names(gene.lookup)

  list(ccle = df.ccle, lookup = gene.lookup)
}

ccle_expr_toTaiga <- function(df.ccle){
  colnames(df.ccle) <- gsub('fh_|\\-Tumor$','',colnames(df.ccle))
  mx.ccle <- df.ccle[, c(-1, -2)] %>% t()
  colnames(mx.ccle) <- apply(df.ccle[, c(1,2)], 1, function(x) paste(x[[2]], ' (',  x[[1]], ')', sep = ''))
  mx.ccle
}

ccle_prism_data_cleaning <- function(df.ccle, df.prism, prism.cellinfo = 'RNAseq_sample_CCLE_name.csv'){
  # Args:
  #   df.ccle: CCLE data table (count, rpkm etc.)
  #   df.prism: PRISM data table (count, rpkm etc.)
  #   prism.cellinfo: a file containing names of prism cell lines and the matching CCLE cell line

  # simplify ccle matrix
  list[df.ccle, gene.lookup] <- ccle_expr_cleaning(df.ccle)

  df.prism <- df.prism[, c(-1, -2)] %>% as.matrix() # remove gene name columns and convert to numerical matrix
  rownames(df.prism) <- names(gene.lookup)

  # replace prism cell line names with standard CCLE names
  cellinfo <- read.csv(prism.cellinfo)
  rownames(cellinfo) <- cellinfo$CCLE.Name
  nm <- gsub('[\\. \\-]','_', cellinfo$Cell.Line)
  ccle_name_prism25 <- cellinfo$CCLE.Name %>% as.character()
  names(ccle_name_prism25)  <- nm

  nm <- gsub('[\\.-]','_', colnames(df.prism))
  colnames(df.prism) <- ccle_name_prism25[nm] %>% as.character()

  # returns
  #   ccle: cleaned ccle data matrix
  #   prism: cleaned prism data matrix
  #   lookup: gene lookup table
  list(ccle = df.ccle, prism = df.prism, lookup = gene.lookup, cell.info = cellinfo) # return a list and a data frame

}

log2_rpkm <- function(mx, min.rpkm = -7){
  mx <- log2(mx)
  mx[is.infinite(mx)] <- NA
  mx[mx < min.rpkm] <- min.rpkm
  mx
}
