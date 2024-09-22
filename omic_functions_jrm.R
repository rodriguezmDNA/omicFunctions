##### reMake functions using tidyverse

## Original functions are in the metaFunctions_forNetworkAnalysis.R \
# and work with base R functions.
# The new functions work similar but are written using tidyverse and intended to work with tibbles and similar data types
###########################################################################
### Test, works on a for loop but not in a lapply
# each <- DEListTibble[[1]]
# binaryTable
# binaryTable <- binaryTable %>%
#   left_join( each %>%
#                filter( .[grep("adj.P.Val",colnames(each)) ] < 0.05) %>%
#                select( "Genes",grep("logFC",colnames(each))) %>%
#                rename(!!str_replace(colnames(each)[grep("logFC",colnames(each))], pattern = ".logFC$", replacement = "") := !!colnames(each)[grep("logFC",colnames(each))] )
#              
#   ) 

## Not working, unkown problem // Used a for loop instead. Check later
# lapply(DEListTibble, function(each){
#   
#     binaryTable <- binaryTable %>%
#     left_join( each %>%
#                  filter( .[grep("adj.P.Val",colnames(each)) ] < 0.05) %>%
#                  select( "Genes",grep("logFC",colnames(each))) %>%
#                  rename(!!str_replace(colnames(each)[grep("logFC",colnames(each))], pattern = ".logFC$", replacement = "") := !!colnames(each)[grep("logFC",colnames(each))] )
#                
#     ) 
# })


get_significant_logFC <- function(listDFs) {
  logFCvalues <- Reduce(union,lapply(listDFs, "[","Genes")) ## Get all genes in list of DEG tables
  #columns <- names(DEListTibble) # Get all the names of the tables
  #for (each in columns){ binaryTable <- add_column(binaryTable,!!(each):=NA)} # Add empty columns
  str(logFCvalues)
  
  ## For each 
  for (each in seq_along(listDFs) ){   
    each <- listDFs[[each]]
    # print (binaryTable)
    logFCvalues  <- logFCvalues %>%
      left_join( each %>%
                   filter( .[grep("adj.P.Val",colnames(each)) ] < 0.05) %>%
                   select( "Genes",grep("logFC",colnames(each))) %>%
                   rename(!!str_replace(colnames(each)[grep("logFC",colnames(each))], pattern = ".logFC$", replacement = "") 
                          := !!colnames(each)[grep("logFC",colnames(each))] )
                 
      )
    # print (binaryTable) 
  }
  #logFCvalues
  return(logFCvalues)
}

make_binary_table <- function(logFCvalues) {
  binaryTable <- mutate_at(logFCvalues, vars( -1), funs( ifelse(is.na(.),0,1)))
  return(binaryTable)
}


collapse_columns_tibbles <- function(listDFs,byColumn="logFC"){
  collapseColumns <- Reduce(union,lapply(listDFs, "[","Genes")) ## Get all genes in list of DEG tables
  #columns <- names(DEListTibble) # Get all the names of the tables
  #for (each in columns){ binaryTable <- add_column(binaryTable,!!(each):=NA)} # Add empty columns
  str(collapseColumns)
  collapseColumns
  ## For each 
  for (each in seq_along(listDFs) ){   
    each <- listDFs[[each]]
    # print (binaryTable)
    
    collapseColumns <- collapseColumns %>%
      left_join( each %>%
                   #filter( .[grep("adj.P.Val",colnames(each)) ] < 0.05) %>%
                   select( "Genes",grep(byColumn,colnames(each))) %>%
                   rename(!!str_replace(colnames(each)[grep(byColumn,colnames(each))], pattern = paste0(".", byColumn), replacement = "") 
                          := !!colnames(each)[grep(byColumn,colnames(each))])
                 
      )
    # each["Solyc06g074010.3.1",]
  }
  #logFCvalues
  return(collapseColumns)
}



## Match DE genes in the list with the empty table. Add true/1 if they're DE 





