
#### This piece of code is just to generate a random set of data, to show how the code works
pValues <- data.frame(replicate(10,runif(50,0,1))) # Creates a matrix of random p Values (between 0 and 1)

## Add names to the dataframe
GOterms <- replicate(50,paste0(sample(letters,3),collapse = "")) # these would be the GO terms
samples <- replicate(10,paste0(sample(letters,10),collapse = "")) # sample names
rownames(pValues) <- GOterms
colnames(pValues) <- samples

## With this block I make sure to have some cells that will be below a threshold
rowIdx <- sample(seq(1,nrow(pValues)),10)
colIdx <- sample(seq(1,ncol(pValues)),4)
pValues[rowIdx,colIdx] <- pValues[rowIdx,colIdx] * .01



### Define colors
cols <-  colorRampPalette(RColorBrewer::brewer.pal(name = "Blues",9))
## Add asterisks to the plot if the p-value is below a threshold
tmp_signif <- ifelse(pValues < 0.05,"*","")

## Transform data
hmData <- -log10(pValues)

## Set the cell size and font
cell = 8#nrow(dat)*.45
fsize = cell*.5

## Load library
library(pheatmap)

pheatmap(as.matrix(hmData),#cluster_rows = F,cluster_cols = F,
         cellwidth = cell,cellheight = cell,
         display_numbers = tmp_signif,na_col = "gray",
         fontsize_row = fsize,fontsize_col = fsize, number_color = "black", fontsize_number = fsize,
         color = cols(600)) 
