##### Generic Genomic functions


### From a tibble of gene identifieres remove everything after a .
# Use this to remove information from genomic build (ie 
# SolycXXXX.3.2,SolycYYY.2.5 -> SolycXXXX,SolycYYYY )

removeDotGene <- function(dotGeneName,ColumnName="",verbose=F){
  if (verbose) cat ("removeDotGene function called\n")
  tmp <- gsub("\\.[0-9].*$","",unlist(dotGeneName))
  return(tmp)
}
