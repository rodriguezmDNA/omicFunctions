
### Obtained from gramene/biomart 
require(tidyverse)
cat("Reading table of AGI - ITAG orthologues (Downloaded from Ensenmbl/gramene using biomart)\n
* Ath - TAIR10
* Sly - ITAG2.5\n\n")
orthologues <- read_tsv("~/Desktop/OrthoSly/biomart_gramene/grameneBioMart-Sly2.5Genes_AthTAIR10-GeneStableID_GeneName.txt")
colnames(orthologues) <- c("Sly2.5_ITAG","TAIR10_ID","TAIR10_Symbol")
orthologues



######## Convert AGI to Sly
## Takes a tibble with a column with AGI identifiers and returns a table with matched Sly orthologues (0, 1 or many) and the associated TAIR10_Symbol.
findOrtho_AGI2Sly <- function(DataWithAGI,columnName="TAIR10_ID"){
  tmpMatched <- DataWithAGI %>%
    rename("TAIR10_ID":=!!columnName) %>%
    left_join( orthologues %>% 
                 mutate(TAIR10_ID=removeDotGene(TAIR10_ID)), "TAIR10_ID") %>%
    rename(!!columnName:="TAIR10_ID")
  return(tmpMatched)
}

## Reciprocal
findOrtho_Sly2AGI <- function(DataWithITAG,columnName=colnames(DataWithITAG)){
  tmpMatched <- DataWithITAG %>%
    #rename("Sly2.5_ITAG":=!!columnName) %>% 
    dplyr::rename_(. ,.dots=setNames(list(columnName),"Sly2.5_ITAG")) %>%
    mutate(Sly2.5_ITAG=removeDotGene(Sly2.5_ITAG)) %>%
    left_join( orthologues %>% 
                 mutate(Sly2.5_ITAG=removeDotGene(Sly2.5_ITAG)), by = c("Sly2.5_ITAG"="Sly2.5_ITAG")) %>%
    dplyr::rename_(. ,.dots=setNames(list("Sly2.5_ITAG"),columnName))
  return(tmpMatched)
}









