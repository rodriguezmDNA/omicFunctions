########################
##### 

plotAllPoints <- function(tmp_logFCPair,Xname,Yname){
  ###### Make plot
  tmpTitle <- paste(Xname, "vs", Yname,sep=" ")
  tmpFullPlot <- tmp_logFCPair %>%
    ggplot(.,mapping = aes_string(x=Xname, y=Yname)) +
    geom_point( shape=20,
                alpha=ifelse(nrow(tmp_logFCPair)>700,0.15,1),
                size=ifelse(nrow(tmp_logFCPair)>900,0.5,0.7)) + 
    geom_hline(yintercept = 0, linetype="dashed",size=0.25,color="darkblue") +
    geom_vline(xintercept = 0, linetype="dashed",size=0.25,color="darkblue") +
    theme_gray() + ggtitle(tmpTitle) + xlab("logFC") + ylab("logFC") #+
    #theme(
      #panel.grid = element_blank()
    #  panel.background = element_rect(fill = "white", colour = "grey50")
    #)
  return(tmpFullPlot)
}


##########
subset_logFC_direction <- function(tmp_logFCPair,Xname,Yname,Xind,Yind,addFit,directionX=1,directionY=1,maxX=NULL,maxY=NULL){
  ## Subsets a tibble with pairs of logFC values for two contrasts by direction (+/+, +/-,etc...) 
  # This function is intended to be called by another (see below) function
  
  ###### Set a default value of max X and Y
  
  if (is.null(maxX) && is.null(maxY)){
    cat ("Setting max limits to default \n")
    maxX <- max(tmp_logFCPair[2])
    maxY <- max(tmp_logFCPair[3])
  } 
  
  ###### Subset by direction 
  subsetByDirection <- tmp_logFCPair %>%
    filter(sign(.[2])==directionX , sign(.[3])==directionY )
  
  ###### Make plot
  tmpTitle <- paste(Xname, "vs", Yname,"\nSignificant - by Quadrant",sep=" ")
  tmpPlot <- subsetByDirection %>%
    ggplot(.,mapping = aes_string(x=Xname, y=Yname)) +
    coord_cartesian(xlim = c(0,maxX), ylim = c(0,maxY), expand = TRUE) +
    geom_point( shape=20,
                alpha=ifelse(nrow(subsetByDirection)>500,1,1),
                size=ifelse(nrow(subsetByDirection)>500,0.7,0.7)) + 
    theme_gray() + #ggtitle(tmpTitle) + 
     xlab("logFC") + ylab("logFC") #+
    # theme(
    #   #panel.grid = element_blank()
    #   panel.background = element_rect(fill = "white", colour = "grey50")
    #   # panel.grid.minor.y =  element_blank(),
    #   # panel.grid.minor.x =  element_blank()
    #   #panel.grid.major.y =  element_blank(),
    #   #panel.grid.major.x =  element_blank()
    # )
  if (addFit) tmpPlot <- tmpPlot + geom_smooth(method='lm',formula=y~x,se=F,size=0.25)
  ####
  return(list("subset"=subsetByDirection,"plot"=tmpPlot))
}
########################################################################################################################
########################################################################################################################


########################
##### 
##########
get_significant_logFCRelationship <- function(logFCTibble,Xpattern,Ypattern,pValueTable,addFit=T){ 
  # Needs a pattern to filter exactly one contrast for X axis and another to exactly match one contrast for the Y axis
  # Direction 1/0 to extract genes with a positive (1) or negative (-1) sign
  
  ### Process names
  # Get indices
  Xind <- grep(Xpattern,names(logFCTibble)) 
  Yind <- grep(Ypattern,names(logFCTibble))
  # Get names
  Xname <- names(logFCTibble)[Xind]
  Yname <- names(logFCTibble)[Yind]
  ### Verify
  if (length(Xind)+length(Yind) != 2){
    cat ("Make each pattern matches one single contrast \n")
    cat ("X pattern matches:", Xname, "\n")
    cat ("Y pattern matches:", Yname, "\n")
  } else {
    cat ("Processing",Xname, " vs ", Yname, "\n")
  }
  
  ### Subset table into columns wanted
  tmp_logFC <- logFCTibble %>% 
    select(Genes,Xind,Yind)
  
  ### Divide into sectors, by direction of change of the X and Y axes. 4 possible combinations
  maxX <- ceiling(max(tmp_logFC[2]) * 1.15)
  maxY <- ceiling(max(tmp_logFC[3]) * 1.15)
  
  ###
  cat ("Full plot of pairwise logFC relationship \n")
  tmp_logFCFullPlot <- plotAllPoints(tmp_logFC,Xname,Yname)
  
  cat ("Getting plots \n")
  plotFull <- plot_grid(
    tmp_logFCFullPlot,
    labels = c("A"), align = "h")
  
  
  ################################################
  ### Now, filter by p Value and subset plots.
  ################################################
  cat ("Filtering significant values\n")
  
  significant <- pValueTable %>% 
    select(Genes,Xind,Yind) %>%
    filter(.[2]<=0.05 , .[3] <= 0.05 ) #Get genes that are significant for both contrasts
  
  tmp_significant_logFC <- tmp_logFC %>% 
    filter(Genes %in% significant$Genes)
  
  ### Divide into sectors, by direction of change of the X and Y axes. 4 possible combinations
  cat ("Plotting quadrants by significant values \n")
  signif_logFC_IQ <-   subset_logFC_direction(tmp_significant_logFC,Xname,Yname,Xind,Yind,addFit,1,1,maxX,maxY)
  signif_logFC_IIQ <- subset_logFC_direction(tmp_significant_logFC,Xname,Yname,Xind,Yind,addFit,-1,1,-maxX,maxY)
  signif_logFC_IIIQ <- subset_logFC_direction(tmp_significant_logFC,Xname,Yname,Xind,Yind,addFit,-1,-1,-maxX,-maxY)
  signif_logFC_IVQ <- subset_logFC_direction(tmp_significant_logFC,Xname,Yname,Xind,Yind,addFit,1,-1,maxX,-maxY)
  
  signifTitle <- ggdraw() + draw_label( paste0(Xname," vs ",Yname), fontface='bold')
  
  plotQuadrantSignif <- plot_grid(
    signif_logFC_IIQ[["plot"]],
    signif_logFC_IQ[["plot"]],
    signif_logFC_IIIQ[["plot"]],
    signif_logFC_IVQ[["plot"]],
    labels = c("II", "I","III","IV"), align = "h")
  
  plotQuadrantSignif <- plot_grid(signifTitle, plotQuadrantSignif, ncol=1, rel_heights=c(0.1, 1)) 
  # Test for truth
  nrow(signif_logFC_IQ[["subset"]]) +
    nrow(signif_logFC_IIQ[["subset"]]) +
    nrow(signif_logFC_IIIQ[["subset"]]) +
    nrow(signif_logFC_IVQ[["subset"]]) == nrow(tmp_significant_logFC)
  
  significantGenes <- list(
    "QI"=signif_logFC_IQ[["subset"]],
    "QII"=signif_logFC_IIQ[["subset"]],
    "QIII"=signif_logFC_IIIQ[["subset"]],
    "QIV"=signif_logFC_IVQ[["subset"]])
  
  ####
  return(list("significantGenes"=significantGenes,"plotFull"=plotFull,"plotQuadrantSignif"=plotQuadrantSignif))
}
########################################################################################################################
########################################################################################################################


########################
##### 
##########
get_gene_relativeComplement <- function(logFCTibble,Xpattern,Ypattern,pValueTable){ 
  ## The relative complement of A with respect to a set B, also termed the difference of sets A and B, 
  # written B ∖ A, is the set of elements in B but not in A.
  
  ##### Start function here
  ############################################################
  
  ## Check we're dealing with only two elements
  Xind <- grep(Xpattern,names(logFCTibble)) 
  Yind <- grep(Ypattern,names(logFCTibble))
  # Get names
  Xname <- names(logFCTibble)[Xind]
  Yname <- names(logFCTibble)[Yind]
  ### Verify
  if (length(Xind)+length(Yind) != 2){
    cat ("Make each pattern matches one single contrast \n")
    cat ("X pattern matches:", Xname, "\n")
    cat ("Y pattern matches:", Yname, "\n")
  } else {
    cat ("Processing",Xname, "vs", Yname, "\n")
  }
  
  ### Subset table into columns wanted
  ################################
  
  ##### Filter by pValue
  onlyTotal <- pValueTable %>% 
    select(Genes,Xind,Yind) %>%
    filter(.[2]<=0.05 , .[3] > 0.05) %>%
    select(Genes)
  
  onlyTRAP <- pValueTable %>% 
    select(Genes,Xind,Yind) %>%
    filter(.[2] > .05 , .[3] <= 0.05) %>%
    select(Genes)
  
  onlyBoth <- pValueTable %>% 
    select(Genes,Xind,Yind) %>%
    filter(.[2]<=0.05 , .[3] <= 0.05) %>%
    select(Genes)
  
  ## Filter logFC table by genes significantly DE in either Total, TRAP or both
  
  genes_onlyTotal <- logFCTibble %>% 
    select(Genes,Xind,Yind) %>%
    filter (Genes %in% onlyTotal$Genes)
  
  genes_onlyTRAP <- logFCTibble %>% 
    select(Genes,Xind,Yind) %>%
    filter (Genes %in% onlyTRAP$Genes)
  
  genes_onlyBoth <- logFCTibble %>% 
    select(Genes,Xind,Yind) %>%
    filter (Genes %in% onlyBoth$Genes)
  
  ######
  complementResult <- list("genes_onlyTotal"=genes_onlyTotal,
                           "genes_onlyTRAP"=genes_onlyTRAP,
                           "genes_onlyBoth"=genes_onlyBoth)
  return(complementResult)
  
}
################################################################################################################################################################################################################################################