library(GSEABase) #1.32.0
library(GSVA) #1.18.0

#' Transformation of Dependency datasets using ssGSEA method
#'
#' The transformation results in a per cell line value for each geneset. The ssGSEA method is the 
#' original Barbie et al. (2009) method implemented by the GSVA package. If lower values 
#' represent stronger dependencies in the input dataset and the alternative is set to "lower", the 
#' resulting matrix will also have lower values for more dependent genesets. 
#' 
#' Row names in the dependency dataset should be symbols/names that match symbol identifiers in the gmt
#' 
#' @param dep_dataset gene (row) by cell line (col) matrix of dependency values
#' @param geneset path to geneset gmt file 
#' @param alternative indicates whether more negative (lower) or more positive (upper) dependency 
#'        values should get higher scores
#' @return matrix genesets (row) by cell line (col)
ssGSEA_dependency_transformation <- function(dep_dataset,
                                             geneset,
                                             alternative=c("upper","lower","either")){
  alternative <- match.arg(alternative)
  
  gsc <- getGmt(geneset,geneIdType=SymbolIdentifier())
  
  absolute <- ifelse(alternative == "either", T, F)
  if(identical(alternative,"lower")){dep_dataset <- -dep_dataset}
 
  result <- gsva(dep_dataset,gsc,method="ssgsea",abs.ranking=absolute)
  
  #keep more dependent genesets as negative values
  result <- if(identical(alternative,"lower")){result <- -result}
  
  return(result)
}



