
# convert gene alias to official gene symbol 
# based on database of HGNC and NCBI

alias2geneSymbol <- function(){
  # get all gene synonyms (alias)
  # HGNC
  hgnc <- read.delim(url("http://www.genenames.org/cgi-bin/hgnc_downloads.cgi?title=HGNC+output+data&hgnc_dbtag=on&col=gd_app_sym&col=gd_aliases&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit&.cgifields=&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag"))
  names(hgnc)[1] <- 'Symbol'
  hgnc.alias <- hgnc %>% filter(!grepl('withdrawn', Symbol), str_length(Synonyms) > 0) # remove wuithdrawn gene symbols
  
  # NCBI
  ncbi <- read_delim("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", delim = '\t')
  ncbi.alias <- ncbi[!is.na(ncbi$Symbol), c(2,3,5)]
  names(ncbi.alias)[1] <- 'EntrezGeneID'
  
  # merge HGNC and NCBI
  #' @param symbol: gene symbols
  #' @param alias: synonyms or aliases that match gene symbols
  #' @param sep: separation character   
  alias_parser <- function(symbol, alias, sep){
    # symbols and alias has to match
    stopifnot(length(symbol) == length(alias))
    parsed <- str_split(alias, pattern = sep)
    names(parsed) <- symbol
    
    flatten.key <- unlist(parsed, use.names = F)
    flatten.count <- sapply(parsed, length)
    flatten.value <- rep(names(flatten.count), flatten.count)
    data.frame(Synonym = flatten.key, Symbol = flatten.value)
  }
  
  
  hgnc.alias.parsed <- alias_parser(hgnc.alias$Symbol, hgnc.alias$Synonyms, sep = ', ')
  ncbi.alias.parsed <- alias_parser(ncbi.alias$Symbol, ncbi.alias$Synonyms, sep = '\\|')
  alias.merged <- rbind(hgnc.alias.parsed, ncbi.alias.parsed) %>% unique()
  alias.merged <- with(alias.merged, transform(alias.merged, Synonym = as.character(Synonym), Symbol = as.character(Symbol)))
  # remove one-alias-matches-all cases
  alias.merged.unique <- alias.merged %>% group_by(Synonym) %>% summarise(count = n()) %>% filter(count == 1)
  
  alias.lookup <- alias.merged[alias.merged$Synonym %in% alias.merged.unique$Synonym, 'Symbol']
  names(alias.lookup) <-  alias.merged[alias.merged$Synonym %in% alias.merged.unique$Synonym, 'Synonym']
  alias.lookup
}