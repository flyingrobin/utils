require(dplyr)

source('~/Documents/Broad/Utils/batch_fisher_test.R')
source('~/Documents/Broad/Utils/df_to_lookup.R')

# load precomputed MSigDB binary matrix
if(!exists('list.msigdb')){
  load("~/Documents/Broad/Utils/data/msigdb_binary_mx.RData")
}

if(!exists('list_rep_anno')){
  list_rep_anno <- readRDS("~/Documents/Broad/Utils/data/repurposing_binary_mx.rds")
}


# convert mx.msigdb to list of feat.blocks ("GO", "Hallmk", "Pertb", "Pthway")
msigdb_to_list <- function(mx.msigdb){
  feat.block <-  levels(gsub('[[:digit:]]*', '', colnames(mx.msigdb)) %>% as.factor())
  list.msigdb <- lapply(feat.block, function(c){
    mx.msigdb[, grep(c, colnames(mx.msigdb))]
  })
  names(list.msigdb) <- feat.block
  list.msigdb
}
#list.msigdb <- msigdb_to_list(mx.msigdb)

#' @description Run Fisher extact test (implemented as hypergeomatric test) of MSigDB binary feature matrices
#' @param input.list a list of input elements, could be a list of genes for GO analysis
#' @param feat a list of feature binary matrix. The name of each element in the list is the name for the feature block
#' @param A scalar number, min.input.size minimum size of input to run fisher exact test
#' @param min.feat.size minimum size of feature to run fisher exact test, must be a vector with the length equals the length of the feature list 
#' @param cutoff: must be a vector with the length equals the length of the res list 
#' @return a data frame containing feature blocks, row and colmumn number and names and the value

overlapMsigDB <- function(input.list, min.input.size = 10, min.feat.size = 2, cutoff = 1e-10){
  batch_fisher_test(input.list, feat.list = list.msigdb, 
                    min.input.size = min.input.size, 
                    min.feat.size = min.feat.size, 
                    cutoff = cutoff) %>%
    mutate(geneset.name = genesets.lookup[geneset.name])
}


overlapREP <- function(input.list, min.input.size = 10, min.feat.size = 2, cutoff = 1e-10){
  batch_fisher_test(input.list, feat.list = list_rep_anno, 
                    min.input.size = min.input.size, 
                    min.feat.size = min.feat.size, 
                    cutoff = cutoff)
}

