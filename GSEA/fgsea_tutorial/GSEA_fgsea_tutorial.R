#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("fgsea")
library(fgsea)
library(tidyverse)
## R vignette
# http://www.bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

### Example run
fgseaResexample <- fgsea(pathways = examplePathways, 
                  stats = exampleRanks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)

head(fgseaResexample[order(pval), ])

plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")



##### With real data
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

load("Annotated_Results_LvV.rdata")
load("mouse_H_v5.rdata")

head(annotLvV)
head(Mm.H)

gseaDat <- filter(shrinkLvV, !is.na(Entrez)) #Remove genes without an Entrez 
ranks <- gseaDat$logFC
names(ranks) <- gseaDat$Entrez
head(ranks)

barplot(sort(ranks, decreasing = T))

## Get the pathways to test for enrichment
pathwaysH <- Mm.H

## Perform GSEA
fgseaRes <- fgsea(pathwaysH, 
                  ranks, # The function automatically sorts the list
                  minSize=15, maxSize = 500, nperm=1000)


head(fgseaRes[order(padj, -abs(NES)), ], n=10)

plotEnrichment(pathwaysH[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
               ranks) + labs(title="HALLMARK_OXIDATIVE_PHOSPHORYLATION")

plotEnrichment(pathwaysH[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]], 
               ranks) + labs(title="HALLMARK_ESTROGEN_RESPONSE_EARLY")

fgseaRes[fgseaRes$pathway == "HALLMARK_ESTROGEN_RESPONSE_EARLY",,drop=F]


topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)