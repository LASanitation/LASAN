#Evaluate the nLAEI.
rm(list=ls())
require(rgdal)
require(raster)
require(foreach)
require(ENMeval)
require(randomForest)
require(DescTools)
require(plyr)
require(tidyverse)
require(dplyr)
require(plyr)
require(caret)
require(gtools)
require(ggplot2)
require(viridis)
require(cowplot)

#Set your working directory.  Your's will be different on your machine.
wd <- "/Users/levisimons/Desktop/Practicum/LASAN/Code"
#wd <- "/home1/alsimons/LASAN"
setwd(wd)
if(wd=="/home1/alsimons/LASAN"){
  rasterOptions(todisk = FALSE)
  options( java.parameters = "-Xmx100g" )
  rasterOptions(todisk=FALSE,maxmemory=2e+11,chunksize=1e+09)
  require(rJava)
  #To deal with java issues
  .jinit()
  require(dismo)
}
if(wd=="~/Desktop/Practicum/LASAN/Code"){
  options(java.parameters = "-Xmx1g" )
  require(rJava)
  #To deal with java issues
  .jinit()
  require(dismo)
}

#Read in the boundaries for LA.  The CRS is EPSG:2229.
LABoundaries <- readOGR("LACityBoundaries.shp")

#Read in filtered environmental map layers generated here: https://github.com/LASanitation/LASAN/blob/main/LAEcologyIndexGeneratorV2.R
env.filtered.further <- read.table("EnvFilteredFurtherV2.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
env.filtered <- c(paste(wd,"/envLayers/",env.filtered.further$env.filtered.further,".tif",sep=""))
nLAEI <- raster("LANativesEcologicalIndexVersion4.tif")
#Stack environmental layers
env.data <- stack(c(env.filtered,nLAEI))
#Get environmental layer names
env.filenames <- env.filtered.further$env.filtered.further

RFEvaluationTotal <- data.frame()
RFImportanceTotal <- data.frame()
RFPartialPlotsTotal <- data.frame()
set.seed(1)
sampleNum <- 200
LoopNum <- 100
for(i in 1:LoopNum){
  # Initialize data containing environmental layer values at presence locations.
  samplevals <- as.data.frame(raster::extract(env.data,as.data.frame(spsample(LABoundaries,n=2*sampleNum,"random"))))
  # Remove rows with missing data.
  samplevals <- samplevals[complete.cases(samplevals),]
  # Resample down to a uniform number of points per iteration.
  samplevals <- samplevals[sample(nrow(samplevals),sampleNum),]
  #Set factor variables.
  envFactors <- c("cal_fire","ClimateZones","DLC","DominantCanopyCover","EcoRegion","Ecotopes","FloodPlain","LandCover","LandUse","PublicLandStatus","Totrcv","Vegcover","WHR","WildlandUrbanInterface")
  samplevals[,colnames(samplevals) %in% envFactors] <- lapply(samplevals[,colnames(samplevals) %in% envFactors],factor)
  
  #Construct a training and testing set for the streams subset data.
  group <- kfold(samplevals,5)
  LAEITrain <- samplevals[group!=1,]
  LAEITest <- samplevals[group==1,]
  #Run the random forest model using the training data.
  rf1 <- suppressWarnings(tuneRF(LAEITrain[,c(env.filenames)],LAEITrain$LANativesEcologicalIndex,stepFactor=1,plot=FALSE,doBest=TRUE))
  
  #Evaluate the random forest model.
  prediction <- predict(rf1,LAEITest)
  #Evaluate how well the random forest model can predict the testing data's output.
  RFEvaluation <- data.frame(matrix(nrow=1,ncol=2))
  colnames(RFEvaluation) <- c("r","p")
  #Evaluate how closely correlated the predicted and actual values are.
  erf <- cor.test(prediction,LAEITest$LANativesEcologicalIndex)
  RFEvaluation$r <- erf$estimate
  RFEvaluation$p <- erf$p.value
  #Aggregate evaluations of each data subset's random forest model
  RFEvaluationTotal <- rbind(RFEvaluationTotal,RFEvaluation) 
  #Calculate variable importance within the random forest model
  RFImportance <- as.data.frame(importance(rf1))
  RFImportance <- data.frame(names=row.names(RFImportance),RFImportance)
  RFImportance$names <- as.character(RFImportance$names)
  RFImportanceTotal <- rbind(RFImportanceTotal,RFImportance)
  #Create the data frame to store the partial response plot source data.
  RFPartialPlots <- data.frame(matrix(ncol=1,nrow=nrow(LAEITrain)))
  for(modelVar in dplyr::top_n(RFImportance,31,IncNodePurity)$names){
    tmp <-  as.data.frame(partialPlot(rf1,LAEITrain,x.var=c(modelVar),plot=FALSE))
    colnames(tmp) <- c(modelVar,paste("LAEI",modelVar))
    RFPartialPlots <- plyr::rbind.fill(RFPartialPlots,tmp)
    RFPartialPlots <- RFPartialPlots[,colSums(is.na(RFPartialPlots)) < nrow(RFPartialPlots)]
    RFPartialPlots <- RFPartialPlots[rowSums(is.na(RFPartialPlots)) < ncol(RFPartialPlots),]
    RFPartialPlots <- RFPartialPlots[!duplicated(RFPartialPlots),]
  }
  RFPartialPlotsTotal <- plyr::rbind.fill(RFPartialPlotsTotal,RFPartialPlots)
}

#Coerce character columns to factor.
RFPartialPlotsTotal <- RFPartialPlotsTotal %>% mutate(across(where(is.character), as.factor))

#Generate summary statistics for how well our random forest models perform on a per watershed basis.
#Filter out random forest evaluations which are not significant.
RFEvaluationSignificant <- RFEvaluationTotal[RFEvaluationTotal$p <= 0.0001,]
#Save RF model evaluation statistics.
write.table(RFEvaluationSignificant,paste("RFEvaluationNatives",LoopNum,"Iterations",sampleNum,"Samples.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#Calculate mean Pearson correlation coefficients, and their significance scores,
#Using this method: https://www.tandfonline.com/doi/abs/10.1080/00221309809595548
FisherZInv(mean(FisherZ(RFEvaluationSignificant$r)))
FisherZInv(sd(FisherZ(RFEvaluationSignificant$r)))
print(paste("R=",FisherZInv(mean(FisherZ(RFEvaluationSignificant$r))),"+/-",FisherZInv(sd(FisherZ(RFEvaluationSignificant$r)))))
#Calculate the mean, and standard deviation on the relative importance of RF model variables.
RFImportanceTotal <- ddply(RFImportanceTotal, .(names), summarize,  MeanIncNodePurity=mean(IncNodePurity), SDIncNodePurity=sd(IncNodePurity))
colnames(RFImportanceTotal) <- c("Variable","MeanIncNodePurity","SDIncNodePurity")
#Save relative variable importance values.
write.table(RFImportanceTotal,paste("RFImportanceNatives",LoopNum,"Iterations",sampleNum,"Samples.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
#Select the most important RF model variables.
ImportantVars <- dplyr::top_n(RFImportanceTotal,10,MeanIncNodePurity)$Variable
#Generate partial dependence plots for the RF model against the most important variables.
i=1
for(modelVar in ImportantVars){
  #Subset the partial plot dataframe so only the columns with the variable name and its LAEI response are extracted.
  RFTotal <- RFPartialPlotsTotal[,c(grep(paste("^",modelVar,"$",sep=""),names(RFPartialPlotsTotal)),grep(paste("^LAEI ",modelVar,"$",sep=""),names(RFPartialPlotsTotal)))]
  RFTotal <- RFPartialPlotsTotal[,modelVar==names(RFPartialPlotsTotal) | paste("LAEI ",modelVar,sep="")==names(RFPartialPlotsTotal)]
  #Remove empty rows.
  RFTotal <- RFTotal[complete.cases(RFTotal),]
  #Coerce character columns to factor.
  RFTotal <- RFTotal %>% dplyr::mutate(across(where(is.character), as.factor))
  #Create a violinplot for categorical model variables.
  print(paste(modelVar,"Factor:",sum(sapply(RFTotal,is.factor)),"Numeric:",sum(sapply(RFTotal,is.numeric)),"Character:",sum(sapply(RFTotal,is.character))))
  if(is.factor(RFTotal[,1])){
    RFPlot <- ggplot(RFTotal, aes_string(x=colnames(RFTotal)[1],y=as.name(colnames(RFTotal)[2])))+
      geom_violin(na.rm=T)+
      xlab(colnames(RFTotal)[1])+ylab("nLAEI")+
      theme_bw(base_size=25)+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))+scale_x_discrete(guide=guide_axis(n.dodge=3))
    assign(paste("RFPlotNatives",LoopNum,"Iterations",sampleNum,"SamplesPanel",i,sep=""),RFPlot)
  }
  #Create heatmaps for numerical model variables.
  if(is.numeric(RFTotal[,1])){
    RFPlot <- ggplot(RFTotal, aes_string(x=as.name(colnames(RFTotal)[1]), y=as.name(colnames(RFTotal)[2])))+
      xlab(colnames(RFTotal)[1])+ylab("nLAEI")+
      geom_bin2d(bins = 20)+scale_fill_continuous(type = "viridis")+
      stat_smooth(aes_string(y = as.name(colnames(RFTotal)[2]), fill=as.name(colnames(RFTotal)[2])),method="auto",formula=y~x,color="violet",fill="red",n=1*nrow(RFTotal))+theme_bw(base_size=25)+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15))
    assign(paste("RFPlotNatives",LoopNum,"Iterations",sampleNum,"SamplesPanel",i,sep=""),RFPlot)
  }
  i=i+1
}
#Get the list of plots generated by the above script.
RFPlotList <- ls(pattern="RFPlotNatives(.*?)Iterations(.*?)SamplesPanel(.*?)")
#Output the list of plots as pdfs.
RFPlotsTotal <- mget(RFPlotList)
invisible(mapply(ggsave, file=paste0(names(RFPlotsTotal), ".pdf"), plot=RFPlotsTotal, width=10,height=10,units="in"))
###

RFEvaluationSignificant <- read.table(paste("RFEvaluationNatives",LoopNum,"Iterations",sampleNum,"Samples.txt",sep=""), header=T, sep="\t",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
RFImportanceTotal <- read.table(paste("RFImportanceNatives",LoopNum,"Iterations",sampleNum,"Samples.txt",sep=""), header=T, sep="\t",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")

#Calculate mean Pearson correlation coefficients, and their significance scores,
#Using this method: https://www.tandfonline.com/doi/abs/10.1080/00221309809595548
FisherZInv(mean(FisherZ(RFEvaluationSignificant$r)))
FisherZInv(sd(FisherZ(RFEvaluationSignificant$r)))
