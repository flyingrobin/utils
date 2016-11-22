
#' @title Convert a data frame to a named vector as a lookup table
#' @description  Given a data frame, in which there is a key column and a value column, convert it to a vector whose value is the value column and key is the key columns. The key column must have all unique entries.
#' @param df input data frame
#' @param key.col a string or column number of the key column 
#' @param value.col a string or column number of the value column 
#' @param both.way = F logical of whether returning a lookup table of both key-to-value and value-to-key. If ture, require 1-to-1 relationship of key-value pairs
df_to_lookup <- function(df, key.col, value.col, both.way = F){
  df <- as.data.frame(df)
  stopifnot(nrow(df) == nrow(df[, key.col] %>% unique())) # key column must have all unique entries
  values <- df[, value.col]
  names(values) <- df[, key.col]
  
  if(both.way){
    stopifnot(nrow(df) != length(df[, value.col] %>% unique()))
    values2 <- df[, key.col]
    names(values2) <- df[, value.col]
    values <- c(values, values2)
  }
  values
}

# switch key and value for a lookup vector
rev_lookup <- function(lookup){
  stopifnot(length(lookup) == length(unique(lookup))) 
  values <- names(lookup)
  names(values) <- values
  values
}