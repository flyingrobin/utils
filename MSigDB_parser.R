library(taigr)
library(dplyr)
library(readr)
library(useful)
library(tidyr)
library(magrittr)
library(rvest)
library(stringr)



##############################################################
# code for generating the msigdb binary matrix
# gene alias lookup table
source("~/Documents/Broad/Utils/alias2geneSymbol.R")
alias.lookup <- alias2geneSymbol()

# replace gene alias with official gene symbol. 
# if gene symbol is not among anno.gene$GeneSymbol, remove this row
anno.gene <- load.from.taiga(data.name = 'ensembl-gene-type-annotation', data.version = 1)

correct_geneSymbol <- function(x){
  ifelse(x %in% anno.gene$GeneSymbol, x, 
         ifelse(x %in% names(alias.lookup), alias.lookup[x], NA)) %>%
    na.omit()
}

gmt.path <- "/Users/wangli/Documents/Broad/Utils/data/"
gmt.files <- list.files(path = gmt.path, pattern = ".*\\.gmt")

gmt.list <- lapply(gmt.files, function(x){
  GSA::GSA.read.gmt(paste(gmt.path, x, sep = '/'))
}) 
names(gmt.list) <- c('GO', 'Hallmk', 'Pthwy', 'Pertb')

# convert alias to offical gene symbol
gmt.list.geneset.corrected <-
  lapply(gmt.list, function(x){
  parallel::mclapply(x$genesets, correct_geneSymbol)
})

gmt.list$GO$genesets <- gmt.list.geneset.corrected$GO
gmt.list$Hallmark$genesets <- gmt.list.geneset.corrected$Hallmark
gmt.list$Pathway$genesets <- gmt.list.geneset.corrected$Pathway
gmt.list$Perturb$genesets <- gmt.list.geneset.corrected$Perturb

# convert gene set list to binary matrix
gmt.list.flat <- unlist(gmt.list, recursive = F)

gene.union <- unlist(gmt.list.flat[c(1,4,7,10)]) %>% unique()
names(gmt.list.flat)[c(1,4,7,10)] <- c('GO', 'Hallmk', 'Pthwy', 'Pertb')
genesets.union <- unlist(gmt.list.flat[c(1,4,7,10)], recursive = F)

# lookup table for geneset names -vs- short names
genesets.lookup <- data.frame(
  geneset.name = unlist(gmt.list.flat[c(2,5,8,11)]),
  geneset = names(unlist(gmt.list.flat[c(2,5,8,11)]))
) %>% mutate(geneset = gsub('\\.geneset.names', '', geneset))


# binary matrix is n (# of genes) x m (# of genesets)
mx.msigdb <- matrix(0, nrow = length(gene.union), ncol = length(genesets.union))
rownames(mx.msigdb) <- gene.union %>% sort()
colnames(mx.msigdb) <- names(genesets.union)
for(i in 1:length(genesets.union)){
  rows <- mx.msigdb[genesets.union[[i]], names(genesets.union)[i]] <- 1
}


# filter out row names that not among anno.gene$GeneSymbol
mx.msigdb <- mx.msigdb[rownames(mx.msigdb) %in% anno.gene$GeneSymbol, ]

# convert to sparsematrix
mx.msigdb.sparse <- Matrix(mx.msigdb, sparse = TRUE)

list_msigdb <- 
  lapply(c('GO', 'Hallmk', 'Pthwy', 'Pertb'), function(i){
    mx.msigdb.sparse[, grepl(i, colnames(mx.msigdb.sparse))]
}) %>% set_names(c('GO', 'Hallmk', 'Pthwy', 'Pertb'))

# Gene set description, scraping form GSEA webside using rvest
genesets.dscrp <- lapply(c(3,6,9,12), function(i){
  sapply(gmt.list.flat[[i]], function(url){
    sample.page <- read_html(url)
    dscrp.brief <- sample.page %>% html_nodes("tr:nth-child(3) td") %>% html_text()
    dscrp.brief[1]
  }, USE.NAMES = F)
})
genesets.dscrp <- unlist(genesets.dscrp)
names(genesets.dscrp) <- unlist(gmt.list.flat[c(2,5,8,11)])

save(list_msigdb, genesets.lookup, genesets.dscrp, file = "~/Documents/Broad/Utils/data/msigdb_binary_mx.RData")
saveRDS(list_msigdb, file = "~/Documents/Broad/Utils/data/msigdb_binary_mx.RDS")
