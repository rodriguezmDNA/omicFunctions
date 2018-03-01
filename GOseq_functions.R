## To install goseq:
# Make sure a SQL client is installed.
# for mac use `brew install mariadb-connector-c``
### 
##source("https://bioconductor.org/biocLite.R")
# install.packages('RMySQL', type='source')
#biocLite(c("goseq","GenomicRanges","Rsamtools","rtracklayer"))
library(goseq)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)

##### About the readGO function #####
### The file format that the readGO function processes is as follows:
# First column are the gene identifiers
# Second column the GO terms associated with the gene, separated by a coma and might contain a space or not

#                 ITAG                        GO
# 1 Solyc00g005000.2.1 0019538, 0009056, 0016787
# 1 Solyc00g005000.2.1 0006508,0004190

############ !!!!!! ############
### If the file you're trying to read has a different format, modify the regex into an appropiate form to match your format. - DON'T PANIC.


##### About the call_goseq function #####
# The removeDots is a logical argument (TRUE/FALSE). It can be turned on if the GFF or GO annotation is from a different genome release than the used in the experiments 
# (ie, when GO annotation is from a previous release -gene id ID#####.2.1- and RNASeq was done in a release v3.1 -gene id ID#####.3.1-)
# The parameter is intended to 'normalize' the annotations to avoid resulting in few overlaps 
####


######################## Here be functions ########################

######## Read GFF and get gene lengths ########

## This function takes in the name of a GFF file to read and process the lengths of the genes to be used by goseq 
# the name must be quoted and must be inside the GOdata/ directory 
# GFFfile <- "Name_of_file.gff"

getGeneLengthfromGFF <- function(gffFileName) {
  ## Requires rtracklayer
  require(rtracklayer)
  ###
  GFFPath <- paste0("GOdata/",gffFileName)
  GFF <- import.gff(GFFPath,version="3",feature.type="gene")
  grl <- GenomicRanges::reduce(split(GFF, mcols(GFF)$Name)) #Conflict with reduce between GRanges and tidyverse
  reducedGTF <- unlist(grl, use.names=T)
  mcols(reducedGTF)$Name <- rep(names(grl), elementNROWS(grl))
  reducedGTF
  AllLengths <- width(reducedGTF)
  names(AllLengths) <- mcols(reducedGTF)$Name
  return(AllLengths)
}



######## Read GO file and process it ########
## This function takes in the name of a GO file to read and process it to a format that goseq likes
# the name must be quoted and must be inside the "GOdata/" directory 
# goFileName <- "Name_of_file.gff"

readGO <- function(goFileName){
  
  GOitag <- read.table( paste0("GOdata/",goFileName) ,
                        stringsAsFactors = F,header = T)
  
  GOitagSplit <- split(GOitag,GOitag[,1])
  tmp <- GOitagSplit[[1]]
  tst <- lapply(GOitagSplit,function(tmp){cbind( 
    tmp[,1],unlist(     # Split by regex:
      strsplit(tmp[,2],",[[:space:]]{0,}"))) #Split by a coma and any or no space after. 
    # This makes it work for both files I got from Ted's data, but just in case any other GO file might be using tabs or any type of space.
    # If
  })
  GOitagNew <- data.frame(do.call(rbind,tst),stringsAsFactors = F)
  colnames(GOitagNew) <- c("ITAG","GOID")
  head(GOitagNew)
  GOitagNew[,2] <- paste0("GO:",GOitagNew[,2])
  return(GOitagNew)
}




######## Call goseq functions and perform enrichment testing ########

## Calls go seq functions
call_goseq <- function(genesToAnalyze,assayed.genes,AllLengths,go.goseq,removeDots=TRUE,plotFit=FALSE){
  ## Intended to normalize names if using info from different builds. Doesn't hurt if the genes names are the same across different data sets
  if (removeDots){
    cat ("--- \n Removing genome version info from gene names (GeneId1234.genomeVersion3 --> GeneId1234)\n---\n")
    
    assayed.genes <- removeDotGene(assayed.genes)
    genesToAnalyze <- removeDotGene(genesToAnalyze)
    names(AllLengths) <- removeDotGene ( names(AllLengths) )
    go.goseq$ITAG <- removeDotGene (go.goseq$ITAG)
  } else cat ("Leaving gene ID's as is")
  
  
  cat("Total number of genes analyzed (universe):", length(assayed.genes),"\n")
  cat("Number of genes used for testing", length(genesToAnalyze),"\n")
  
  #
  gene.vector=as.integer(assayed.genes%in%genesToAnalyze)
  names(gene.vector)=assayed.genes
  #head(gene.vector)
  #table(gene.vector)
  
  ### Filter length info from all the genes in the genome to those analyzed (present in the universe)
  GeneLengths <- AllLengths[names(gene.vector)]
  
  ### Calls go seq function
  pwf=nullp(gene.vector, bias.data=GeneLengths,plot.fit = plotFit) #By default doens't plot the fit 
  go.all <- goseq(pwf, gene2cat=go.goseq)
  
  # Return results
  return(go.all)
}


######## Find which genes from each ########

## signGOList is a list of all the GO terms that were significant when testing genes. 
## GenePatterns is a matrix where each column is a list of genes of interest that were used for the GO enrichment testing

get_Genes_Per_Category <- function(signGOList,GenePatterns) {
  getGenesPerCategory <- lapply(seq_along(signGOList),function(DP){
    require(tidyverse)
    cat("Genes per category in pattern", names(signGOList)[DP],"\n")
    
    ## Get DE genes
    de.genes <- GenePatterns[GenePatterns[,DP]!="",DP]
    
    ### Get categories
    categories <- signGOList[[DP]][,"category"]
    
    ## From the full GO seq list, get the genes that overlap with the DE
    DEgenesInCategory <- go.goseq %>%
      filter(.$GOID %in% categories) %>% group_by(GOID) %>% filter(removeDotGene(ITAG) %in% removeDotGene(de.genes)) %>% nest()
    ### Convert the table with significant GO categories to a tibble and join the list-columns table
    test <- as.tibble(signGOList[[DP]])
    test <- test %>% 
      left_join(DEgenesInCategory,by=c("category"="GOID")) %>% mutate("Set"=names(signGOList)[DP]) #%>% rename("DEGenesinCategory"="data") 
    return(test)
  })
  names(getGenesPerCategory) <- names(signGOList)
  return(getGenesPerCategory)
}


######## Find which genes from each ########
## Calls the findOrtho_Sly2AGI from an accompanying script
# Need to add > source(orthologueFindingFunctions.R) to the script
findOrthologues_GOgenes <- function(significantGOtibble){
  significantGOtibble <- significantGOtibble %>%  # To this table
    left_join( # Join the following transformation:
      significantGOtibble %>% select(uniqueID,data) %>%
        # Group by category (maybe this isn't necessary)
        mutate(                  # Create a new column that is
          "Sly2AGI_orthologues" = map(data, .f= findOrtho_Sly2AGI)) %>% 
        select(uniqueID,3) 
      , by = c("uniqueID" = "uniqueID")
    ) %>% rename("DE_ITAGsinCats"="data")
  return(significantGOtibble)
}


gg_makeTileMap <- function(significantGOtibble,extraName="",keyword=NULL,filterbyOnto=c("MF","BP","CC")){
  require(gtools)
  plotTitle <- "Heatmap"
  if (extraName != "") plotTitle <- paste(plotTitle,extraName,sep=" - ")
  if (!is.null(keyword)) {
    searchThis <- paste(keyword,collapse = "|") 
    plotTitle <-paste(plotTitle, paste("keywords: ", paste(keyword,collapse=",") ,sep = ""), sep="\n")
  } else searchThis <- "*"
  ## Heatmap
  plotHeat <- significantGOtibble %>%  filter(grepl(searchThis,term,ignore.case = T)) %>%
    ## Do some filtering
    select(Set,term, category,ontology,2) %>% 
    mutate(Set=factor(Set,levels = mixedsort(unique(Set),decreasing = F))) %>%
    filter(ontology %in% filterbyOnto) %>% #If only one ontology term is wanted
    #### Dealing with NA's
    #mutate("term" = ifelse(is.na(term),category,term)) %>% # Either assign the GOcategory to a term that is NA,
    filter(!is.na(term)) %>%# or filter them out
    
    ## Group
    group_by(term,ontology,Set) %>%
    #gather(., key=contrast,value = pVal,over_represented_pvalue)
    ## In case we want significance points
    mutate("signif"=ifelse(over_represented_pvalue < 0.05,"*","")) %>% 
    ## Transform 
    mutate("log10pVal"=-log10(over_represented_pvalue)) %>%
    ### Make the heatmap
    ggplot(., aes(Set, term)) + 
    geom_tile(aes(fill=log10pVal),colour="white",width=1, height=1) + 
    #scale_fill_distiller(palette = "BuGn",direction = 1) +
    scale_fill_distiller(palette = "Blues",direction = 1,na.value="white",name="-log10(P)")+
    facet_grid(ontology~.,scales="free") +
    
    ### Add text
    #geom_text(aes(label=signif),size=2.5, hjust = 0.5) + 
    ## Names of axes and other graph aesthetics
    ggtitle(plotTitle) + xlab("Set") + ylab("GO terms") + 
    theme_minimal() + 
    theme(#axis.title.x = element_blank(),
      axis.ticks = element_blank(), 
      panel.background=element_rect(fill="white", colour="lightgray"),
      panel.grid.major.x =  element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color="black", size=8),
      axis.text.x = element_text(angle = 0, hjust = 1))
    #scale_alpha(guide = 'none') + 
    #guides(fill=guide_legend(title="-log10(pVal)"))
  return(plotHeat)
}

######## Make a wrapper to deal with the results ########

wrapSignificantGOterms <- function(GOresults,GeneList,saveToFile=F,resultsPath="GenesPerCategory",extraName="",keyword=NULL,filterbyOnto=c("MF","BP","CC"),minFreq=1,findOrhto=T) {
  require(purrrlyr)
  ##### Filter by p-value of overrepresented (enriched) categories
  listSignificantGOterms <- lapply(GOresults, function(x){  x[x$over_represented_pvalue < 0.05,] })
  
  
  ### Calls a function to get the genes that were present in each category
  getGenesPerCategory <- get_Genes_Per_Category(listSignificantGOterms,GeneList)
  
  # Make a table out of the lists
  significantGOtibble <- bind_rows(getGenesPerCategory) %>%  # Bind lists into a single table
    mutate(uniqueID = row_number())
  
  ### Get the orthologues of the genes
  if (findOrhto) {
  ### Find the Ath orthologues for each Sly gene
  significantGOtibble <- findOrthologues_GOgenes(significantGOtibble) 
  }
  significantGOtibble
  
  
  ####### Make pretty graphs ########
  plotTitle <- "Frequency of Categories"
  if (extraName != "") plotTitle <- paste(plotTitle,extraName,sep=" - ")
  ## GO term frequency
  
  tmpCounts <- significantGOtibble %>% filter(!is.na(term)) %>%
    group_by(ontology,term) %>% summarise("Count"=n()) %>% arrange(ontology,desc(Count)) 
  ## 
  if (max(tmpCounts$Count)>=minFreq) tmpCounts <- tmpCounts %>% filter(Count >= minFreq) 
  ## Make plot
  plotByOnthology <- tmpCounts %>% 
    ggplot(.,aes(x=reorder(term,Count), # Order X axis labels by Count
                 y=Count, #Counts go in the Y axis
                 fill=ontology)) + #Color by ontology
    geom_bar(stat="identity") + coord_flip() + # Do a bar plot 
    xlab("GO category") + xlab("Frequency") + labs(title=plotTitle)  + 
    facet_grid(ontology~.,scales="free") #Separate by ontology
  
  plotTitle <- "Number of enriched GO terms per set"
  if (extraName != "") plotTitle <- paste(plotTitle,extraName,sep="\n")
  ## Terms per Set
  plotOfCatNumber <- significantGOtibble %>% select(Set) %>% group_by(Set) %>% count() %>%
    ggplot(.,aes(x=as.factor(Set),y=n)) + labs(title=plotTitle)  + 
    geom_bar(stat="identity") + #,width = 10/length(unique(significantGOtibble$Set)) ) + # Do a bar plot 
    xlab("Set") + ylab("# of significant GO terms") + theme_light()
  
  
  ## Heatmap
  plotHeat <- gg_makeTileMap(significantGOtibble,extraName,keyword,filterbyOnto)
  
  ######## Save to files (or not) ########
  if (saveToFile) { 
    cat("---\nSaving to file\n---\n")
    ##### Create directories to save outputs
    if (extraName != "") resultsPath <- paste(resultsPath,extraName,sep = "-")
    
    itagPath <- paste0(resultsPath,"/ITAG_byCategory/")
    orthoPath <- paste0(resultsPath,"/AthOrthologues_byCategory/")
    dir.create(itagPath,recursive = T) 
    dir.create(orthoPath,recursive = T)
    
    unite_(significantGOtibble,"tmpName",c("term","ontology","Set"),sep="_",remove=F) %>%  mutate("tmpName"=gsub("[\\>'()\\[\\\\/]","",tmpName)) %>%
      ## Create new names for the files
      mutate("itagFile"=paste0(itagPath,"/",gsub("/|[[:space:]]","-",tmpName),".txt")) %>%
      by_row(~write.table(.$DE_ITAGsinCats, file = .$itagFile,sep = "\t",row.names = F,quote = F)) %>%
      mutate("orthoFile"=paste0(orthoPath,"/",gsub("/|[[:space:]]","-",tmpName),".txt")) %>%
      by_row(~write.table(.$Sly2AGI_orthologues, file = .$orthoFile,sep = "\t",row.names = F,quote = F))
  } else cat("---\nNot saving to file\n---\n")
  
  return(list("significantGOtibble"=significantGOtibble,
              "plotByOnthology"=plotByOnthology,
              "plotOfCatNumber"=plotOfCatNumber,
              "plotHeatmap"=plotHeat))
}

