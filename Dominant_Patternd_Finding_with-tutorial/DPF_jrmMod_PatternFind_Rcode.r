resultsFolder <- "DPF_Results/"
dir.create(resultsFolder,showWarnings = T)
##-------------- Parameters for first step - fuzzy K-means clustering
##--------------------------------------------------------------------

inputFile <- "Data/MeansNormalized_MAexpression_Hu_2007_ATH1_chl1_Significant.csv"
#inputFile <- "Normalized_MAexpression_Hu_2007_ATH1_chl1_Full.csv"
minExpFilter <- 5 #FALSE OR 1
minVarFilter <- 25 #FALSE OR 50
kChoice <- 5 #20
fuzzyKmemb <- 1.04 #1.05
alreadyLog2 <- FALSE #TRUE
methodResultFile <- paste(resultsFolder,"DPF_Hu2007_chl1Vwt_patternIdent_result.rDump",sep="")
diagnosticFile <- paste(resultsFolder,"DPF_Hu2007_chl1Vwt_clustering_diagnostic.txt",sep="")


##---------This controls how much data is in the final R object
##--------------------------------------------------------------------
saveRawClusterObject <- FALSE #FALSE, if TRUE the file will be large
saveClusterDistanceMat <- FALSE #FALSE, if TRUE the file will be large
saveCompleteExpData <- TRUE #FALSE
saveFunctions <- TRUE #FALSE
saveParamVariables <- TRUE #FALSE

##-------------- Parameters for 2nd step - building patterns
##--------------------------------------------------------------------

autoSelectClusterCutoff <- TRUE #TRUE
clusterCutoff <- 0.4 #0.4
patternSimilarityCutoff <- 0.9 #0.9
pearsonCutoff <- 0.85 #0.85

userPatternInputFile <- FALSE #FALSE
patternOutputFile <- "DPF_Hu2007_chl1Vwt_patternOutput.txt"
groupOutputFile <- "DPF_Hu2007_chl1Vwt_genesToPatterns.txt"
diagnosticFile <- "DPF_Hu2007_chl1Vwt_patternIdent_diagnostic.txt"

patternOutputFile <- paste(resultsFolder,patternOutputFile,sep="")
groupOutputFile <- paste(resultsFolder,groupOutputFile,sep="")
diagnosticFile <- paste(resultsFolder,diagnosticFile,sep="")

#-------------------------------------------------
#-------------------------------------------------

library(cluster)
removeLowE <- 0
removeLowV <- 0
dateRun<- date()
cat("Starting Fuzzy K-Means clustering on ",dateRun,"\n",file=diagnosticFile)
expressionData <-read.csv(file=inputFile,header=TRUE,row.names = 1,stringsAsFactors = F,as.is = T)
expressionData <- as.matrix(expressionData)
cat("Expression data read from ",inputFile,"\n",file=diagnosticFile,append=TRUE)
cat("\t",nrow(expressionData)," genes with ",ncol(expressionData)," observations.\n",file=diagnosticFile,append=TRUE)
expDataFiltered <- expressionData
filterSettings <- paste("Filter Settings:\n\tInput File = ",inputFile,sep="")

#-------------------------------------------------
#-------------------------------------------------

cat("Removing low expressed genes - ",file=diagnosticFile,append=TRUE)
filterSettings <- paste(filterSettings,"\n\tLow Expression Filter = ",sep="")
if(is.numeric(minExpFilter)){
	keep <- c(1:nrow(expDataFiltered))[apply(expDataFiltered,1,max)>=minExpFilter]
	removeLowE <- nrow(expDataFiltered)-length(keep)
	expDataFiltered <- expDataFiltered[keep,]
	cat("Done\n",file=diagnosticFile,append=TRUE)
	filterSettings <- paste(filterSettings,minExpFilter," (",removeLowE," genes removed)\n",sep="")
}else{
	cat("Skipped\n",file=diagnosticFile,append=TRUE)
	filterSettings <- paste(filterSettings,"FALSE\n",sep="")
}

#-------------------------------------------------
#-------------------------------------------------

cat("Removing low varying genes - ",file=diagnosticFile,append=TRUE)
filterSettings <- paste(filterSettings,"\tLow Variance Filter = ",sep="")
if(is.numeric(minVarFilter)){
	lowVarCut <- quantile(apply(expDataFiltered,1,var),probs=(minVarFilter/100))
	keep <- c(1:nrow(expDataFiltered))[apply(expDataFiltered,1,var)>=lowVarCut]
	removeLowV <- nrow(expDataFiltered)-length(keep)
	expDataFiltered <- expDataFiltered[keep,]
	cat("Done\n",file=diagnosticFile,append=TRUE)
	filterSettings <- paste(filterSettings,minVarFilter,"% (",removeLowV," genes removed)\n",sep="")
}else{
	cat("Skipped\n",file=diagnosticFile,append=TRUE)
	filterSettings <- paste(filterSettings,"FALSE\n",sep="")
}
cat(filterSettings,file=diagnosticFile,append=TRUE)

#-------------------------------------------------
#-------------------------------------------------

cat("Log2 Transforming expression data - ",file=diagnosticFile,append=TRUE)
if(!alreadyLog2){
    expDataFiltered[expDataFiltered==0] <- 1e-10
    expDataFiltered <- log2(expDataFiltered/apply(expDataFiltered,1,mean))
    expressionData[expressionData==0] <- 1e-10
    expressionData <- log2(expressionData/apply(expressionData,1,mean))
    cat("Done\n",file=diagnosticFile,append=TRUE)
}else{
    cat("Skipped\n",file=diagnosticFile,append=TRUE)
}

#-------------------------------------------------
#-------------------------------------------------

cat("Building distance matrix for clustering - ",file=diagnosticFile,append=TRUE)
distanceMat<-as.dist((1-cor(t(expDataFiltered),method="pearson"))/2)
cat(" Done\n",file=diagnosticFile,append=TRUE)

fuzzyKSettings <- paste("Fuzzy K-Means settings:\n\tK = ",kChoice,"\n\tmemb.exp = ",fuzzyKmemb,sep="")
fuzzyKSettings <- paste(fuzzyKSettings,"\n\tmaxit = ",(nrow(expDataFiltered)*4),"\n",sep="")

cat("Running fuzzy K-means to produce ",kChoice," clusters from ",nrow(expDataFiltered)," genes.\n",file=diagnosticFile,append=TRUE)
cat(fuzzyKSettings,file=diagnosticFile,append=TRUE)
cat("\n!!! This may take a long time to complete. !!!\n",file=diagnosticFile,append=TRUE)

initClust <- fanny(distanceMat, kChoice, diss = TRUE, memb.exp= fuzzyKmemb,maxit=nrow(expDataFiltered)*4)
cat("Done clustering!\nSaving results to ",methodResultFile,"\n",file=diagnosticFile,append=TRUE)

saveList <- c("expDataFiltered","filterSettings","fuzzyKSettings","dateRun")
if(saveRawClusterObject){
	saveList <- c("initClust",saveList)
}
if(saveClusterDistanceMat){
	saveList <- c("distanceMat",saveList)
}else{
	rm("distanceMat")
}
if(saveCompleteExpData){
	saveList <- c("expressionData",saveList)
}else{
	rm("expressionData")
}

save(file=methodResultFile,list=saveList)
cat("Done. Initial clustering completed.\n\n",file=diagnosticFile,append=TRUE)


#-----------------Functions used by pattern collapse and assignment code -------------
#-------------------------------------------------------------------------------------
getClustProfiles <- function(expDataIn,clustIn,clusterCutoff){
    numClusters <- ncol(clustIn$membership)
    clustProfiles <- matrix(data=0,ncol=ncol(expDataIn),nrow=numClusters)
    for(i in 1:numClusters){
        genesInClust <- c(1:nrow(expDataIn))[clustIn$membership[,i]>=clusterCutoff]
        if(length(genesInClust)>0){
                    if(length(genesInClust)>1){
                             clustProfiles[i,]<-as.numeric(apply(expDataIn[genesInClust,],2,median))
                    }else{
                                clustProfiles[i,]<-as.numeric(expDataIn[genesInClust,])
                    }
        }
    }
    whichCluster <- c(1:numClusters)[apply(clustProfiles,1,var)>0]
    clustProfiles <- clustProfiles[apply(clustProfiles,1,var)>0,]
    return(list(profiles=clustProfiles,whichClust=whichCluster))
}


makeCollapsedProfiles <- function(expDataIn,clustProfsIn, clusterCutoff,grouping,origClust,clustIn){
    collapsedProfiles <- matrix(data=0,ncol=ncol(expDataIn),nrow=max(grouping))
    for(i in 1:max(grouping)){
        inClust <- c(1:length(grouping))[grouping==i]
        if(length(inClust)==1){
            collapsedProfiles[i,]<-clustProfsIn[inClust,]
            genesInClust <- c(1:nrow(expDataIn))[clustIn$membership[,origClust[inClust]]>=0.4]
        }else{
            cat("\tCollapsing profiles ",paste(inClust,collapse=","),"into new profile ",i,"\n")
            clustLookAt <- origClust[inClust]
            genesInClust <- c(1:nrow(expDataIn))[apply(clustIn$membership[,clustLookAt],1,max)>= clusterCutoff]
            collapsedProfiles[i,]<-apply(expDataIn[genesInClust,],2,median)
        }
    }
return(collapsedProfiles)
}

correlateGenesToProfiles <- function(expDataIn,profilesIn){
    profileCor <- matrix(data=0,ncol=nrow(profilesIn),nrow=nrow(expDataIn))
    for(i in 1:nrow(profilesIn)){
        profileCor[,i] <- apply(expDataIn,1,cor,y=profilesIn[i,])
        cat("\tCorrelating genes to profile ",i," out of ",nrow(profilesIn),"\n")
    }
    rownames(profileCor)<-rownames(expDataIn)
    colnames(profileCor)<-rownames(profilesIn)
    return(profileCor)
}


writeMemberList <- function(patternCor,minDist,fileOut){
    maxSize <- max(apply(patternCor>=minDist,2,sum))
    tempMat <- matrix(data="",ncol=ncol(patternCor),nrow=maxSize)
    colnames(tempMat)<-colnames(patternCor)
    for(i in 1:ncol(tempMat)){
        geneList <- rownames(patternCor)[patternCor[,i]>=minDist]
        if(length(geneList)>0){
                tempMat[1:length(geneList),i]<-geneList
        }
    }
    write.table(file=fileOut,tempMat,sep="\t",quote=FALSE,row.names=FALSE,col.names=colnames(tempMat))
}

#-----------------------------------------------------------------------------------------------------
#-----------------Start of executed code-------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

cat("Starting Pattern Identification & Gene Assignment\n",file=diagnosticFile)
load(methodResultFile)
if(autoSelectClusterCutoff){
    clusterCutoff <- initClust$membership[sort(initClust$membership,index.return=TRUE,decreasing=TRUE)$ix[nrow(initClust$membership)]]
}
cantAssignGenes<-(apply(initClust$membership,1,max)<clusterCutoff)
goodGenes<-(apply(initClust$membership,1,max)>=clusterCutoff)
initClust$cluster[cantAssignGenes]<-0
initClust$membership[cantAssignGenes,]<-0
cat("Number of genes assigned to a cluster above ",clusterCutoff," =",file=diagnosticFile,append=TRUE)
cat(sum(goodGenes),"/",length(initClust$cluster),"\n",file=diagnosticFile,append=TRUE)
cat("Number of genes not assigned to a cluster above ",clusterCutoff," =",file=diagnosticFile,append=TRUE)
cat(sum(cantAssignGenes),"/",length(initClust$cluster),"\n",file=diagnosticFile,append=TRUE)
cat("Generating inital set of patterns.\n",file=diagnosticFile,append=TRUE)
clustProfiles <- getClustProfiles(expDataFiltered,initClust,clusterCutoff)
whichClusters <- clustProfiles$whichClust
clustProfiles <- clustProfiles$profiles

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

cat("Calculating distances between each pattern.\n",file=diagnosticFile,append=TRUE)
profileCor <- cor(t(clustProfiles))
diag(profileCor)<-NA
hClust <- hclust(as.dist(1-profileCor),method="single")
collapseGrp <- cutree(hClust,h=(1-patternSimilarityCutoff))

cat("Collapsing similar patterns.\n",file=diagnosticFile,append=TRUE)
finalProfiles<-makeCollapsedProfiles(expDataFiltered,clustProfiles, clusterCutoff ,collapseGrp,whichClusters,initClust)
rownames(finalProfiles)<- paste("Pattern_",1:nrow(finalProfiles),sep="")
colnames(finalProfiles)<-colnames(expDataFiltered)
colnames(finalProfiles)[1] <- paste("\t",colnames(finalProfiles)[1],sep="")
if(is.character(userPatternInputFile)){
	cat("Adding user defined patterns from \"",userPatternInputFile,"\".\n",file=diagnosticFile,append=TRUE)
	userPatterns <- read.delim(file=userPatternInputFile,sep="\t",header=TRUE)
	tempFinal <- matrix(data=0,ncol=ncol(finalProfiles),nrow=nrow(finalProfiles)+nrow(userPatterns))
	tempFinal[1:nrow(finalProfiles),]<-finalProfiles
	tempFinal[(nrow(finalProfiles)+1):nrow(tempFinal),]<-as.matrix(userPatterns[,2:ncol(userPatterns)])
	rownames(tempFinal) <- c(rownames(finalProfiles),userPatterns[,1])
	colnames(tempFinal) <- colnames(finalProfiles)
	finalProfiles <- tempFinal
}
cat("Writing patterns to output file.\n",file=diagnosticFile,append=TRUE)
write.table(file=patternOutputFile,finalProfiles,quote=FALSE,sep="\t")

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

cat("Assigning genes to patterns.\n",file=diagnosticFile,append=TRUE)
geneToPatternCor <- correlateGenesToProfiles(expDataFiltered,finalProfiles)

cat("Writing gene to pattern output file.\n",file=diagnosticFile,append=TRUE)
writeMemberList(geneToPatternCor,pearsonCutoff,groupOutputFile)
cat("Pattern Identification and Gene Assignment completed.\n",file=diagnosticFile,append=TRUE)



if(!saveRawClusterObject){
	rm("initClust")
}
removeVarList <- c("cantAssignGenes","goodGenes","profileCor","hClust","collapseGrp","clustProfiles","removeVarList","removeFunctionList","inputFile","minExpFilter","minVarFilter","fuzzyKmemb","alreadyLog2","diagnosticFile","saveRawClusterObject","saveClusterDistanceMat","saveCompleteExpData","saveFunctions","saveParamVariables","autoSelectClusterCutoff","userPatternInputFile","patternOutputFile","groupOutputFile","removeLowE","removeLowV","kChoice","whichClusters")
removeFunctionList <- c("getClustProfiles","makeCollapsedProfiles","writeMemberList","correlateGenesToProfiles")
if(!saveParamVariables){
	rm(list=removeFunctionList)
}
if(!saveFunctions){
	rm(list=removeVarList)
}

save(file=methodResultFile,list=ls())

# --------- jrm
# I commented out this snippet because breaks R
##This little snippet just checks to see if its possible to draw images on the current machine
##this is only really a problem on cluster machines
#gsexe <- Sys.getenv("R_GSCMD")
#if(is.null(gsexe) || !nchar(gsexe)) {
#        gsexe <- "gs"
#        rc <- system(paste(gsexe, "-help > /dev/null"))
#        if(rc != 0){q("no")} #We cant make images. quit and end program
#}
# ---------





#--------------------------------------------------------------------------------------------
#----This code snippet will draw the heatmaps and line graphs for the patterns & genes-------
#--------------------------------------------------------------------------------------------
#--Change the numbers inside par(mar=c(x,x,x,x)) to adjust margin size on bottom,left,top,right---
#---in case the axis labels do not print, just make the appropriate axis margin larger ---

drawGeneGroupCurves <- function(inputData,title="",xAxisNames=colnames(inputData),centerCurve=ifelse(nrow(inputData)>1,apply(inputData,2,median),inputData)){
	yRng <- max(abs(c(inputData,centerCurve)))
	yRng <- c(-yRng,yRng)
	matplot(1:ncol(inputData),t(inputData),ylim=yRng,xaxt="n",ylab="Expression",main=title,xlab="",type="l",lwd=0.5,col=8,lty=1)
	lines(1:ncol(inputData),centerCurve,col=1,lwd=2)
	if(!is.na(xAxisNames[1])){
		axis(1,1:ncol(inputData),xAxisNames,las=2)
	}
}

drawExpressionHeatmap <- function(inputData,title="",yAxisNames=rownames(inputData),xAxisNames=colnames(inputData),autoOrder=TRUE,maxExp=2,minExp=-2,standAlone=TRUE){
	cols <- c(rgb(0,49:0/49,49:0/49),rgb(1:50/50,1:50/50,0))
	if(autoOrder & nrow(inputData)>2){
		newOrd <- hclust(dist(inputData))$order
		inputData <- inputData[newOrd,]
		yAxisNames <- yAxisNames[newOrd]
	}
	brks <- seq(minExp,maxExp,length=101)
	if(standAlone){
		nf <- layout(matrix(data=c(1,1,1,1,2),nrow=1))
	}
	inputData[inputData>maxExp]<-maxExp
	inputData[inputData<minExp]<-minExp
	image(1:ncol(inputData),1:nrow(inputData),t(inputData),xaxt="n",yaxt="n",xlab="",ylab="",main=title,breaks=brks,col=cols)
	if(!is.na(xAxisNames[1])){
		axis(1,1:ncol(inputData),xAxisNames,las=2)
	}
	if(!is.na(yAxisNames[1])){
		axis(2,1:nrow(inputData),yAxisNames,las=2)
	}
	if(standAlone){
		tempMat <- matrix(data=seq(minExp,maxExp,length=100),ncol=2,nrow=100)
		image(1:2,seq(minExp,maxExp,length=100),t(tempMat),breaks=brks,col=cols,xaxt="n",xlab="",ylab="",main="Key")
	}
}


## Changed bitmap to PDF

load(file=paste(resultsFolder,"DPF_Hu2007_chl1Vwt_patternIdent_result.rDump",sep=""))

pdf(file=paste(resultsFolder,"DPF_Hu2007_chl1Vwt_profileHeatmap.pdf",sep=""),height=4,width=6)
par(mar=c(11,7,2,1))
drawExpressionHeatmap(finalProfiles,title="Patterns")
dev.off()

for(i in 1:nrow(finalProfiles)){
	if(sum(geneToPatternCor[,i]>=pearsonCutoff)>0){
		drawData  <- as.matrix(expDataFiltered[geneToPatternCor[,i]>=pearsonCutoff,])
		if(nrow(drawData)==length(drawData)){ # only one gene
			drawData <- t(drawData)		#make the matrix 1 row instead of 1 col. just a hack
		}
		
		pdf(file=paste(resultsFolder,"DPF_Hu2007_chl1Vwt_pattern_",i,"_geneCurves.pdf",sep=""),height=4,width=6)
		par(mfrow=c(1,1),mar=c(11,5,2,1))
		drawGeneGroupCurves(drawData,centerCurve=finalProfiles[i,],title=paste("Pattern ",i," Genes",sep=""))
		dev.off()

		pdf(file=paste(resultsFolder,"DPF_Hu2007_chl1Vwt_pattern_",i,"_geneHeatmap.pdf",sep=""),height=4,width=6)
		par(mfrow=c(1,1),mar=c(11,4,2,1))
		drawExpressionHeatmap(drawData,yAxisNames=NA,title=paste("Pattern ",i," Genes",sep=""))
		dev.off()
	}else{
		cat("No genes assigned to pattern ",i," --- Skipping drawing\n")
	}
}




