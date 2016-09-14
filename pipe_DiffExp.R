# Differential expression pipeline
###############################################
library(edgeR)
library(DESeq2)
library(BiocParallel)
library(devtools)
source_url("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")

run_edgeR <- function(count.data, group, design, test = 'exact', multi.factor = FALSE){
  print("Start running edgeR...")
  cds <- DGEList(counts=count.data, group = group)
  
  # We need to filter out low count reads 
  # since it would be impossible to detect differential expression. 
  # The method used in the edgeR vignette is to keep only those genes that 
  # have at least 1 read per million in at least 3 samples. 
  #cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  print("Normalize for RNA composition...")
  cds <- calcNormFactors(cds) # normalize for RNA composition
  
  # Estimating Dispersions
  print("Estimating dispersions...")
  cds <- estimateDisp(cds, design = design, robust=TRUE)
  
  # if multi.factor is true, return cds object and for specific tests
  if(multi.factor){
    return(cds)
  }
  print("Detecting differentially expressed genes...")
  stopifnot(test %in% c('exact', 'likelihood', 'quasi'))
  if(test == 'exact'){
    de <- exactTest(cds, pair = unique(group) ) # if group is [WT, WT, KO, KO] then we compare KO/WT
    res_edger <-  topTags(de, n = Inf) %>% as.data.frame()
  
  }else if(test == 'quasi'){
    # Quasi-likelihood F-tests
    fit <- glmQLFit(cds, design)
    qlf <- glmQLFTest(fit, coef= 3:26)
    res_edger <-  topTags(qlf, n = Inf) %>% as.data.frame()
  
  }else if(test == 'likelihood'){
    # Quasi-likelihood F-tests
    fit <- glmFit(cds, design)
    lrt <- glmLRT(fit)
    res_edger <-  topTags(lrt, n = Inf) %>% as.data.frame()
  }
  print("Done!")
  return(list(DGEList = cds, results = res_edger))
}

###############################################################
# count.data <- as.matrix(df.merged[rownames(df.merged) %in% pcg$ID,])
# col.data <- data.frame(condition = c(rep('CCLE',25),rep('PRISM',25)), row.names = names(df.merged))
# dds2 <- DESeqDataSetFromMatrix(countData = count.data2,
#                                colData = col.data,
#                                design = ~ condition)
# 
# register(MulticoreParam(4))
# dds2 <- DESeq(dds2, parallel = T)
# res2 <- results(dds2, parallel = T) %>% as.data.frame()
# res2 <- cbind(geneSymbol = gene.lookup[rownames(res2)], res2)
# 
# # significant level = 0.05
# res.hit2 <- res2[order(res2$padj),]  %>% subset(padj < 0.05)
# 
# vsd <- vst(dds2, blind = FALSE)
# vsd.mx <- assay(vsd)
# 
# rlog <- rlog(dds2, blind = FALSE)
# rlog.mx <- assay(rlog)
# 
# res.up10fold <- res.hit %>% subset(log2FoldChange > log2(10) & padj < 0.01)
# res.dn10fold <- res.hit %>% subset(log2FoldChange < log2(0.1) & padj < 0.01)
# 
# 
# pdf('../report/rlogExpr_fc10fold_PRISM-v-CCLE.pdf', w = 8.1, h = 11)
# par(mfrow = c(4,4), mar = c(2,2,4,2))
# s <- sapply(rownames(res.up10fold), function(x) gene_boxplot(x, count = combined_mx))
# s <- sapply(rownames(res.dn10fold), function(x) gene_boxplot(x, count = combined_mx))
# 
# # explore the relationship of lineage specificity and differential expression
# ccle.info <- load.from.taiga(data.name="ccle-type-refined-lineage-table", data.version=1)
# 
# 
# rld <- rlog(dds2, blind=FALSE)
# vsd <- varianceStabilizingTransformation(dds2, blind=FALSE)
# vsd.fast <- vst(dds2, blind=FALSE)
# 
# library(vsn)
# notAllZero <- (rowMedians(counts(dds2))>10)
# meanSdPlot(log2(counts(dds2,normalized=TRUE)[notAllZero,] + 1))
# meanSdPlot(assay(rld[notAllZero,]))
# meanSdPlot(assay(vsd[notAllZero,]))
# meanSdPlot(assay(vsd.fast[notAllZero,]))
# 
# 
# ##################################################
# #
# # Differential expression analysis using limma
# library(limma)
# dge <- DGEList(counts=count.data )
# dge <- calcNormFactors(dge)
# 
# condition <- c(rep('CCLE',25),rep('PRISM',25))
# design <- model.matrix(~ condition)
# v <- voom(cds,design,plot=TRUE)
# 
# fit <- lmFit(v,design)
# fit <- eBayes(fit)
# res.limma <- topTable(fit, coef=ncol(design), number="Inf", sort.by="P")
# res.hits.limma <- topTable(fit,coef=ncol(design), p = 0.05, n = 5000)
# 
# # remove high variance genes
# genesd <- apply(count.data, 1, sd)
# count.data2 <- count.data[row.names(count.data) %in% names(genesd[genesd < 10000]), ]
# 
# 
# # qval vs mean