library(dplyr)
library(readr)
library(tidyr)
library(Matrix)
library(devtools)  

#' @description Fisher extact test (implemented as hypergeomatric test) of a list of elements over a list of binary feature matrices
#' @param input.list a list of input elements, could be a list of genes for GO analysis
#' @param feat.list a list of feature binary matrix. The name of each element in the list is the name for the feature block
#' @param A scalar number, min.input.size minimum size of input to run fisher exact test
#' @param min.feat.size minimum size of feature to run fisher exact test, must be a vector with the length equals the length of the feature list 
#' @param cutoff: must be a vector with the length equals the length of the res list 
#' @return a data frame containing feature blocks, row and colmumn number and names and the value

batch_fisher_test <- function(input.list, feat.list, min.input.size = 10, min.feat.size = 2, cutoff = 1e-10){
  if(length(min.feat.size) == 1){
    min.feat.size <- rep(min.feat.size, length(feat.list))
  }
  res.list <- lapply(1:length(feat.list), function(i){
    feat_mx <- feat.list[[i]]
    feat_mx <- feat_mx[, colSums(feat_mx) >= min.feat.size[i]] %>% as.matrix() # filter out features based on min.feat.size
    
    if(is.null(colnames(feat_mx))){ # if none of features pass min.feature.size, return NULL
      return(NULL)
    }
    # Output feature size and nodes number
    print(paste('Calculating enrichment of', ncol(feat_mx), 'features for', length(input.list),'nodes', sep = ' '))
    
    # In the hypergeometric setting
    #q number of white balls drawn without replacement 
    #  from an urn which contains both black and white balls (# of overlap genes)
    #m the number of white balls in the urn. (size of a gene set)
    #n the number of black balls in the urn. (total number of the genes - size of a gene set)
    #k the number of balls drawn from the urn. (size of the input gene list)
    
    q <- sapply(input.list, function(x) {colSums(feat_mx[rownames(feat_mx) %in% x, ])}) %>% t()
    m <- matrix(colSums(feat_mx), nrow(q), ncol(q), byrow = T)
    n <- matrix(nrow(feat_mx), nrow(q), ncol(q)) - m # n is the number of genes that not in the geneset
    k <- matrix(sapply(input.list, length), nrow(q), ncol(q))
    
    # hypergeometric test for over-representation
    p.val <- phyper(q-1, m, n, k, lower.tail = F) # the -1 is to get the analytically correct solution
    q.val <- apply(p.val, 1, function(x) p.adjust(x, method = 'BH')) %>% t() # FDR
    list(q.val = q.val, query.size = k, geneset.size = m, overlap = q) 
  })
  names(res.list) <- names(feat.list)
  
  # get sparse results
  if(length(cutoff) == 1){
    cutoff <- rep(cutoff, length(res.list))
  }
  res.sparse <- data.frame(geme.list = character(),
                           feat.block = character(),
                           feat.name = character(),
                           FDR = numeric())
  for(i in 1:length(res.list)){
    q.val <- res.list[[i]]$q.val
    query.size = res.list[[i]]$query.size
    geneset.size = res.list[[i]]$geneset.size
    overlap = res.list[[i]]$overlap
    
    q.val[q.val > cutoff[i]] <- 0
    q.val[is.na(q.val)] <- 0
    sparse.mx <- Matrix(as.matrix(q.val), sparse = TRUE)

    # convert dgCMatrix object to a easily readible data frame
    sparse.mx <- summary(sparse.mx) %>% as.data.frame()
    if(nrow(sparse.mx) == 0){
      sparse.df <- data.frame(query = NA, query.size = NA, 
                              feat.block =  names(res.list)[i],
                              geneset.name = NA, geneset.size = NA, overlap = NA,
                              FDR = NA)
    }else{
      sparse.df <- data.frame(query = rownames(q.val)[sparse.mx$i],
                              query.size = apply(sparse.mx, 1, function(x) query.size[x[[1]], x[[2]]]),
                              feat.block =  names(res.list)[i],
                              geneset.name = colnames(q.val)[sparse.mx$j] %>% as.character(),
                              geneset.size = apply(sparse.mx, 1, function(x) geneset.size[x[[1]], x[[2]]]),
                              overlap =apply(sparse.mx, 1, function(x) overlap[x[[1]], x[[2]]]),
                              FDR = sparse.mx$x)
    }

    res.sparse <- rbind(res.sparse, sparse.df)
  }
  # res.sparse %>%
  #   separate('query' ,c("cell.line", 'UpDn', 'filter'), sep = '\\.') %>%
  #   transform(geneset.name = as.character(geneset.name))
  res.sparse
}
