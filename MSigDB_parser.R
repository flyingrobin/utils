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


gmt.path <- "/Users/wangli/Documents/Broad/PRISM/RNA-seq/data"
gmt.files <- list.files(path = gmt.path, pattern = "MSigDB.*symbol\\.gmt")

gmt.list <- lapply(gmt.files, function(x){
  GSA::GSA.read.gmt(paste(gmt.path, x, sep = '/'))
}) 
names(gmt.list) <- c('GO', 'Pathway', 'Perturb', 'Hallmark')

# convert gene set list to binary matrix
gmt.list.flat <- unlist(gmt.list, recursive = F)

gene.union <- unlist(gmt.list.flat[c(1,4,7,10)]) %>% unique()
names(gmt.list.flat)[c(1,4,7,10)] <- c('GO', 'Pthway', 'Pertb', 'Hallmk')
genesets.union <- unlist(gmt.list.flat[c(1,4,7,10)], recursive = F)

# lookup table for geneset names -vs- short names
genesets.lookup <- names(genesets.union)
names(genesets.lookup) <- unlist(gmt.list.flat[c(2,5,8,11)])
genesets.lookup2 <- names(genesets.lookup)
names(genesets.lookup2) <- genesets.lookup
genesets.lookup <- c(genesets.lookup, genesets.lookup2)

# binary matrix is n (# of genes) x m (# of genesets)
mx.msigdb <- matrix(0, nrow = length(gene.union), ncol = length(genesets.union))
rownames(mx.msigdb) <- gene.union %>% sort()
colnames(mx.msigdb) <- names(genesets.union)
for(i in 1:length(genesets.union)){
  rows <- 
  mx.msigdb[genesets.union[[i]], names(genesets.union)[i]] <- 1
}

# replace gene alias with official gene symbol. 
# If an alias does not have a match, remove this row
# if gene symbol is not among anno.gene$GeneSymbol, remove this row
anno.gene <- load.from.taiga(data.name = 'ensembl-gene-type-annotation', data.version = 1)

# alias names that not offical symbol and exists in the lookup table
alias.rownms <- rownames(mx.msigdb)[rownames(mx.msigdb) %in% names(alias.lookup) &
                                      !rownames(mx.msigdb) %in% anno.gene$GeneSymbol]
# convert alias to offical gene symbol
rownames(mx.msigdb)[rownames(mx.msigdb) %in% alias.rownms] <- alias.lookup[alias.rownms]

# filter out row names that not among anno.gene$GeneSymbol
mx.msigdb <- mx.msigdb[rownames(mx.msigdb) %in% anno.gene$GeneSymbol, ]


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

save(mx.msigdb, genesets.lookup, genesets.dscrp, file = "Documents/Broad/Utils/data/msigdb_binary_mx.RData")
