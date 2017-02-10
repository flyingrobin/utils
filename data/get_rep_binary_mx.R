library(taigr)
library(stringr)

####################################
# helper functions
my_summary <- function(v){ # A summary function that suppress reporting of the number of NAs, if any
  if(!any(is.na(v))){
    res <- summary(v, maxsum = Inf)
  } else{
    summary_v <- summary(v, maxsum = Inf)
    res <- summary_v[-length(summary_v)]
  }
  return(res)
}

get_cpd_featblock <- function(onecol, row.names, min.occur = 2, sep = ','){

  # creat a list in which each element is one entry of the column vector, as a form of char. vector
  tolist <- str_split(onecol, pattern = sep)
  names(tolist) <- paste(seq(1,length(onecol)), '.',sep = '') # names of the list is number + '.'
  # remove empty item in tolist
  tolist <- tolist[sapply(tolist, function(x) max(str_length(x))) > 0]
  
  flat<- unlist(tolist) %>% as.factor()
  mx <- matrix(0, length(onecol), length(levels(flat))) # creat an empty matrix of nrow x no. of levels
  colnames(mx) <- levels(flat) # column names of the matrix the the level name
  rownames(mx) <- seq(1,length(onecol)) # row names is the index (row name) of the input column

  row.i <- names(flat) %>% as.numeric() %>% floor() # names(flat) returns character as '1.5', then converted to number and take the floor
  col.j <- iconv(flat)
  mx[cbind(row.i, col.j)] <- 1

  rownames(mx) <- row.names
  return(mx[, my_summary(flat) > min.occur])
}

####################################
# clean annotations
mx_rep_anno <- load.from.taiga(data.name="prism-repurposing-compounds-annotation", data.version = 1)
mx_rep_anno <- unique(mx_rep_anno) # remove duplicate lines

# Binarize compound features

cpd.feat.target <-
  get_cpd_featblock(mx_rep_anno$target, row.names = mx_rep_anno$Broad_ID, min.occur = 2, sep = ',') %>%
  set_colnames(paste('target', colnames(.), sep = '_'))

cpd.feat.moa <-
  get_cpd_featblock(mx_rep_anno$moa, row.names = mx_rep_anno$Broad_ID, min.occur = 2, sep = ',')

cpd.feat.indications <-
  get_cpd_featblock(mx_rep_anno$indications, row.names = mx_rep_anno$Broad_ID, min.occur = 2, sep = '\\|')

list_rep_anno <- list(TargetGene = cpd.feat.target, MOA = cpd.feat.moa, Indication = cpd.feat.indications)
saveRDS(list_rep_anno, file = '../../Utils/data/repurposing_binary_mx.rds')
