library(dplyr)
library(useful)
library(magrittr)
library(Matrix)
library(doMC)
doMC::registerDoMC(cores=36) # or however many cores you have access to

linlog_hybrid <- function(x, c = 20, base = 10){
  ifelse(x < c, x, c*log(x/c, base = base) + c)
}

ks_test_msigdb <- function(cormat, mx.msigdb){
  df_cormat <- as.data.frame(cormat) %>% mutate(geneSymbol = gsub('Exp_', '', rownames(.)))
  common_genes <- intersect(df_cormat$geneSymbol, rownames(mx.msigdb))
  
  # reorder df_cormat according to the order of rownames of mx.msigdb
  df_cormat %<>% filter(geneSymbol %in% common_genes) %>%
    mutate(geneSymbol = factor(geneSymbol, levels = common_genes))
  mx.msigdb <- mx.msigdb[common_genes, ]
  
  pbmcapply::pbmclapply(1:(ncol(df_cormat)-1), function(i){
    plyr::aaply(mx.msigdb, 2, function(x){
      within.set <- x > 0
      # up = ks.test(df_cormat[!within.set, i], df_cormat[within.set, i], alternative = 'greater')$p.value
      # down = ks.test(df_cormat[!within.set, i], df_cormat[within.set, i], alternative = 'less')$p.value
      # ifelse(up < down, up, -down)
      ks.test(df_cormat[!within.set, i], df_cormat[within.set, i])$p.value
    }, .parallel = T) 
  }) %>% do.call(cbind, .) %>% set_colnames(colnames(cormat))
}

list_msigdb <- readRDS("../../Utils/data/msigdb_binary_mx.RDS")
cor_vp_ge <- readRDS("../cormat_Vp_GE.rds")

x <- ks_test_msigdb(cor_vp_ge, list_msigdb$Hallmark)
# 
# kstest.deseq.reslist <- function(deseq.reslist, mx.msigdb, genesets, genesets.lookup, gene.lookup){
#   mx.kstest <- lapply(genesets, function(gs.nm){
#     gene.nms <- rownames(mx.msigdb)[mx.msigdb[, genesets.lookup[gs.nm]] > 0]
#     
#     # iterate each cell line
#     sapply(deseq.reslist, function(df) {
#       df$in.geneset <- rownames(df) %in% gene.lookup[gene.nms]
#       c(
#         # down = ks.test(filter(df, !in.geneset)$log2FoldChange, 
#         #                filter(df, in.geneset)$log2FoldChange, alternative = "less")$p.value,
#         # up = ks.test(filter(df, !in.geneset)$log2FoldChange, 
#         #              filter(df, in.geneset)$log2FoldChange, alternative = "greater")$p.value
#         down = ksplus::ks.stat(filter(df, in.geneset)$log2FoldChange, filter(df, !in.geneset)$log2FoldChange, 
#                                alternative = "lower", cutoff = -1, p_method = 'empirical')$pvalue,
#         up = ksplus::ks.stat(filter(df, in.geneset)$log2FoldChange, filter(df, !in.geneset)$log2FoldChange, 
#                              alternative = "upper", cutoff = 1, p_method = 'empirical')$pvalue
#       )
#     })
#   })
#   names(mx.kstest) <- genesets
#   
#   # log10 transform and then lin-log transform at threshold - 10
#   mx.down <- sapply(mx.kstest, function(i) i['down', ]) %>% log10()*-1 
#   mx.up <-  sapply(mx.kstest, function(i) i['up', ]) %>% log10()*-1
#   ifelse(mx.up > mx.down, linlog_hybrid(mx.up), -1*linlog_hybrid(mx.down))
# }
# kstest.pergs.percl.permu <- kstest.deseq.reslist(deseq.reslist = single.cl.de.list,
#                                                  mx.msigdb = mx.msigdb, 
#                                                  genesets = overlap.grouped.geneset$geneset.name, 
#                                                  genesets.lookup = genesets.lookup, 
#                                                  gene.lookup = gene.lookup2)
# rownames(kstest.pergs.percl.permu) <- stringr::str_split(rownames(kstest.pergs.percl.permu), '_', simplify = T)[, 1]
# write.csv(as.data.frame(kstest.pergs.percl.permu), "kstest.pergs.percl.permu.csv", row.names = T)
