# Differential expression pipeline
###############################################
library(edgeR)
library(DESeq2)
library(BiocParallel)
library(devtools)
source_url("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")

fit_edgeR <- function(count.data, group, design, fit = 'glm'){
  print("Start running edgeR...")
  cds <- DGEList(counts=count.data, group = group)
  
  #cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  print("Normalize for RNA composition...")
  cds <- calcNormFactors(cds) # normalize for RNA composition
  
  # Estimating Dispersions
  print("Estimating dispersions...")
  cds <- estimateDisp(cds, design = design, robust = T)
  
  stopifnot(fit %in% c('glm', 'quasi-likelihood'))
  if(fit == 'glm'){
    # Quasi-likelihood F-tests
    fit <- glmFit(cds, design)
  }else if(fit == 'quasi'){
    fit <- glmQLFit(cds, design)
  }
  fit
}

run_edgeR <- function(count.data, group, design, test = 'exact'){
  print("Start running edgeR...")
  cds <- DGEList(counts=count.data, group = group)
  
  #cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  print("Normalize for RNA composition...")
  cds <- calcNormFactors(cds) # normalize for RNA composition
  
  # Estimating Dispersions
  print("Estimating dispersions...")
  cds <- estimateDisp(cds, design = design, robust = T)
  
  print("Detecting differentially expressed genes...")
  stopifnot(test %in% c('exact', 'likelihood', 'quasi'))
  if(test == 'exact'){
    de <- exactTest(cds, pair = unique(group) ) # if group is [WT, WT, KO, KO] then we compare KO/WT
    res_edger <-  topTags(de, n = Inf) %>% as.data.frame()
  
  }else if(test == 'quasi'){
    # Quasi-likelihood F-tests
    fit <- glmQLFit(cds, design)
    qlf <- glmQLFTest(fit, coef= 2)
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
# DESeq2 for batch factor test
run_DESeq2 <- function(counts, col.data, full.model, reduced.model){
  library(DESeq2)
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = col.data,
                                design = full.model)
  
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  register(MulticoreParam(4))
  dds <- DESeq(dds, test = "LRT", reduced = reduced.model, parallel = T)
  
  res <- results(dds, parallel = T) %>% as.data.frame()
  res[order(res$padj), ]
}

# ##################################################
# #
# # Differential expression analysis using limma
library(limma)
fit_limma <- function(count.data, group, design){
  print("Start running edgeR...")
  cds <- DGEList(counts=count.data, group = group)
  
  #cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  print("Normalize for RNA composition...")
  cds <- calcNormFactors(cds) # normalize for RNA composition
  
  # Estimating Dispersions
  print("Estimating dispersions...")
  v <- voom(cds, design, plot = F)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  fit
}

# dge <- DGEList(counts=count.data )
# dge <- calcNormFactors(dge)
# 
# condition <- c(rep('CCLE',25),rep('PRISM',25))
# design <- model.matrix(~ condition)
# v <- voom(cds,design,plot=TRUE)
# 
# fit <- lmFit(v,design)
# fit <- ?(fit)
# res.limma <- topTable(fit, coef=ncol(design), number="Inf", sort.by="P")
# res.hits.limma <- topTable(fit,coef=ncol(design), p = 0.05, n = 5000)
# 
# # remove high variance genes
# genesd <- apply(count.data, 1, sd)
# count.data2 <- count.data[row.names(count.data) %in% names(genesd[genesd < 10000]), ]
# 
# 
# # qval vs mean