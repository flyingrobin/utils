library(dplyr)

# Count complete pairs of each pairwise correlation
count_complete_pairs <-function(x, y = NULL){
  # binerize the matrix that all NAs are converted to 0 and otherwise 1
  x.bin <- ifelse(!is.na(x), 1, 0)
  
  if(is.null(y)){
    # Do matrix multiplication of transposed A and A, 
    # the off-diagonal entries are the number of complete pairs of corresponding pair of vectors
    return(t(x.bin) %*% x.bin)
  }else{
    y.bin <- ifelse(!is.na(y), 1, 0)
    
    # the off-diagonal entries are the number of complete pairs of corresponding pair of vectors
    return(t(x.bin) %*% y.bin)
  }
  
  
}

# Fisher transformation of correlation coefficient
fisherZ <- function(x, y = NULL){
  # Args:
  #   x: input matrix of p x m. if y is null, m == p 
  #   y: input matrix of p x n
  # Return: 
  #   Z-transformed correlation coeff. of column vectors of x and y
  
  print("Computing correlation...")
  if(is.null(y)){
    mx_cor <- cor(x, use = "pairwise.complete.obs")
    cor_complete_pair_count <- count_complete_pairs(x)
  }else{
    mx_cor <- cor(x, y, use = "pairwise.complete.obs")
    cor_complete_pair_count <- count_complete_pairs(x, y)
    
  }
  
  print("Computing Fisher transformation...")
  z <- 0.5 * (log(1+mx_cor) - log(1-mx_cor)) * sqrt(cor_complete_pair_count-3)
  z
}