#This script using MaxEnt to identify potential indicator species for Los Angeles.
rm(list=ls())
require(rgdal)
require(rJava)
require(raster)
require(parallel)
require(doParallel)
require(foreach)
require(ENMeval)

#To deal with java issues
.jinit()

#Set your working directory.  Your's will be different on your machine.
wd <- "~/Desktop/Practicum/LASAN/Code"
#wd <- "/home1/alsimons/LASAN"
setwd(wd)
if(wd=="/home1/alsimons/LASAN"){
  rasterOptions(todisk = FALSE)
  options( java.parameters = "-Xmx48g" )
  rasterOptions(todisk=FALSE,maxmemory=2e+11,chunksize=1e+09)
  require(dismo)
}
if(wd=="~/Desktop/Practicum/LASAN/Code"){
  options(java.parameters = "-Xmx1g" )
  require(dismo)
}

#Read in raw GBIF species data.
SpeciesInput <- read.table("0175678-200613084148143.csv", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
SpeciesInput <- SpeciesInput[SpeciesInput$taxonRank=="SPECIES",]
names(SpeciesInput)[names(SpeciesInput)=="decimalLatitude"] <- "latitude"
names(SpeciesInput)[names(SpeciesInput)=="decimalLongitude"] <- "longitude"

#Convert species location data into a spatial dataframe.  Assume a CRS of EPSG:4326.
xy <- SpeciesInput[,c("longitude","latitude")]
xy <- sp::coordinates(xy)
SpeciesLocations <- SpatialPointsDataFrame(coords = xy, data=SpeciesInput)
sp::proj4string(SpeciesLocations) <- sp::CRS('+init=epsg:4326')

#Read in the boundaries for LA.  The CRS is EPSG:2229.
LABoundaries <- readOGR("LACityBoundaries.shp")

#Reproject the species observations into the same CRS as the LA boundaries.
SpeciesLocations <- sp::spTransform(SpeciesLocations, CRS(proj4string(LABoundaries)))

#Clip species locations by the LA city boundaries.
SpeciesLocations <- SpeciesLocations[LABoundaries,]

#Convert LA species observations back to a dataframe.
SpeciesLocations <- as.data.frame(SpeciesLocations)

#Filter LA species observations by time and prevalence.
SpeciesFreq <- as.data.frame(table(SpeciesLocations[SpeciesLocations$year<=2020 & SpeciesLocations$year>=2010,"species"]))
Prevalence <- 30
SpeciesFreq <- SpeciesFreq[SpeciesFreq$Freq >= Prevalence,]
colnames(SpeciesFreq) <- c("species","Freq")
SpeciesFreq$species <- as.character(SpeciesFreq$species)

#Create a list of species above a certain number of observation points.
speciesList <- SpeciesFreq$species

#Get the list of species already evaluated.
speciesDone <- list.files(pattern="MaxentEvaluation(.*?).txt")
speciesDone <- gsub("MaxentEvaluation","",gsub(".txt","",speciesDone))
#Remove the species already evaluated from the next iteration of analysis.
speciesList <-speciesList[!(speciesList %in% speciesDone)]

#Filter out species observation points to only include species with sufficient observational data.
all.presvals <- SpeciesLocations[SpeciesLocations$species %in% speciesList,]

#Specify sampling number for MaxEnt input data.
sampleNum <- 25

#Generate random background points for MaxEnt modeling.
bgPoints <- as.data.frame(spsample(LABoundaries,n=100000,"random"))
#Get all species locations.
all.obs.data <- SpeciesLocations[,c("longitude.1","latitude.1")]
#Remove any potential overlap between random background points and actual observations.
bgPoints <- bgPoints[!(bgPoints$x %in% all.obs.data$decimalLongitude.1) & !(bgPoints$y %in% all.obs.data$decimalLatitude.1),]
colnames(bgPoints) <- c("longitude.1","latitude.1")

#List environmental rasters
env.files <- list.files(pattern=".tif$",full.names=TRUE)
#Stack environmental layers
env.data <- stack(c(env.files))
#Get environmental layer names
env.filenames <- gsub("^./","",gsub(".tif","",env.files))

SDM <- function(i) {
  #Randomly select a species from the list of those not already analyzed.
  #Get the list of species already evaluated.
  speciesDone <- list.files(pattern="MaxentEvaluation(.*?).txt")
  speciesDone <- gsub("MaxentEvaluation","",gsub(".txt","",speciesDone))
  #Remove the species already evaluated from the next iteration of analysis.
  speciesList <-speciesList[!(speciesList %in% speciesDone)]
  species <- speciesList[sample(length(speciesList))[1]]
  #Extract observation locations for a given species.
  obs.data <- SpeciesLocations[SpeciesLocations$species==species,c("longitude.1","latitude.1")]
  
  # Initialize data containing environmental layer values at presence locations.
  presvals <- as.data.frame(raster::extract(stack(env.files),obs.data))
  #Convert NA values to 0.
  presvals[is.na(presvals)] <- 0
  #Remove rows with missing data.
  presvals <- presvals[complete.cases(presvals),]
  #Remove lat/lon
  presvals$longitude.1 <- NULL
  presvals$latitude.1 <- NULL
  #Specify a presence/absence column
  presvals$pa <- 1
  
  #Create a background point set.
  abs.data <- bgPoints[sample(nrow(bgPoints),20*nrow(presvals)),]
  absvals <- as.data.frame(raster::extract(stack(env.files),abs.data))
  #Convert NA values to 0.
  absvals[is.na(absvals)] <- 0
  #Remove rows with missing data.
  absvals <- absvals[complete.cases(absvals),]
  #Specify a presence/absence column
  absvals$pa <- 0
  #Remove lat/lon
  absvals$longitude.1 <- NULL
  absvals$latitude.1 <- NULL
  
  #Merge presence/absence data
  #Set factor variables.
  modelData <- rbind(presvals,absvals)
  modelData$LandUse <- factor(modelData$LandUse,levels=unique(modelData$LandUse))
  modelData$Ecotopes <- factor(modelData$Ecotopes,levels=unique(modelData$Ecotopes))
  modelData$FloodPlain <- factor(modelData$FloodPlain,levels=unique(modelData$FloodPlain))
  modelData$ClimateZones <- factor(modelData$ClimateZones, levels=unique(modelData$ClimateZones))
  modelData$FloodPlain <- factor(modelData$FloodPlain, levels=unique(modelData$FloodPlain))
  modelData$PublicLandStatus <- factor(modelData$PublicLandStatus, levels=unique(modelData$PublicLandStatus))
  modelData$WildlandUrbanInterface <- factor(modelData$WildlandUrbanInterface, levels=unique(modelData$WildlandUrbanInterface))
  modelData$LandCover <- factor(modelData$LandCover, levels=unique(modelData$LandCover))
  modelData$DominantCanopyCover <- factor(modelData$DominantCanopyCover, levels=unique(modelData$DominantCanopyCover))
  modelData$PotentialNaturalVegetation <- factor(modelData$PotentialNaturalVegetation, levels=unique(modelData$PotentialNaturalVegetation))
  modelData$DLC <- factor(modelData$DLC, levels=unique(modelData$DLC))
  modelData$Vegcover <- factor(modelData$Vegcover, levels=unique(modelData$Vegcover))
  modelData$Totrcv <- factor(modelData$Totrcv, levels=unique(modelData$Totrcv))
  modelData$WHR <- factor(modelData$WHR, levels=unique(modelData$WHR))
  modelData$cal_fire <- factor(modelData$cal_fire, levels=unique(modelData$cal_fire))
  
  #Create dataframe to store MaxEnt model evaluation metrics.
  XMEvaluations <- data.frame()
  #Run MaxEnt model using subsamples of data
  for(i in 1:100){
    #Randomly subsample the number of observation and pseudoabsence points.
    presenceSubsample <-  modelData[modelData$pa==1,]
    presenceSubsample <- presenceSubsample[sample(nrow(presenceSubsample),sampleNum),]
    absenceSubsample <- modelData[modelData$pa==0,]
    absenceSubsample <- absenceSubsample[sample(nrow(absenceSubsample),20*sampleNum),]
    modelSubsample <- rbind(presenceSubsample,absenceSubsample)
    #maxent()
    xm <- maxent(x=modelSubsample[,c(env.filenames)],p=modelSubsample$pa,factors=c('Ecotopes','FloodPlain','LandUse','ClimateZones','FloodPlain','PublicLandStatus','WildlandUrbanInterface','LandCover','DominantCanopyCover','PotentialNaturalVegetation','DLC','Vegcover','Totrcv','WHR','cal_fire'))
    #f <- list("Ecotopes"=levels(modelData$Ecotopes),"FloodPlain"=levels(modelData$FloodPlain),"LandUse"=levels(modelData$LandUse),"ClimateZones"=levels(modelData$ClimateZones),"FloodPlain"=levels(modelData$FloodPlain),"PublicLandStatus"=levels(modelData$PublicLandStatus),"WildlandUrbanInterface"=levels(modelData$WildlandUrbanInterface),"LandCover"=levels(modelData$LandCover),"DominantCanopyCover"=levels(modelData$DominantCanopyCover),"PotentialNaturalVegetation"=levels(modelData$PotentialNaturalVegetation),"DLC"=levels(modelData$DLC),"Vegcover"=levels(modelData$Vegcover),"Totrcv"=levels(modelData$Totrcv),"WHR"=levels(modelData$WHR),"cal_fire"=levels(modelData$cal_fire))
    #Evaluate maxent model.
    exm <- suppressWarnings(evaluate(p=modelSubsample[modelSubsample$pa==1,c(env.filenames)],a=modelSubsample[modelSubsample$pa==0,c(env.filenames)],model=xm))
    tmp <- data.frame(matrix(nrow=1,ncol=5))
    colnames(tmp) <- c("species","AUC","SEDI","r","p")
    tmp$species <- species
    tmp$AUC <- exm@auc
    #SEDI: https://natureconservation.pensoft.net/article/33918/
    a <- mean(exm@TPR,na.rm=T)
    b <- mean(exm@FPR,na.rm=T)
    c <- mean(exm@FNR,na.rm=T)
    d <- mean(exm@TNR,na.rm=T)
    H <- a/(a+c)
    F <- b/(b+d)
    tmp$SEDI <- (log(F)-log(H)-log(1-F)+log(1-H))/(log(F)+log(H)+log(1-F)+log(1-H))
    tmp$r <- exm@cor
    tmp$p <- exm@pcor
    XMEvaluations <- rbind(XMEvaluations,tmp)
    print(paste(species,"auc:",tmp$AUC,"SEDI:",tmp$SEDI,"cor:",tmp$r,"p:",tmp$p))
  }
  write.table(XMEvaluations,paste("MaxentEvaluation",species,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
}

#Run MaxEnt evaluations on all of the species.
lapply(1:length(speciesList),SDM)

#Get the list of species already evaluated.
speciesDone <- list.files(pattern="MaxentEvaluation(.*?).txt")
##Summarize all evaluations so each row is a unique species with all of its associated maxent model evaluations.
XMEvaluationsTotal <- data.frame()
for(speciesEval in speciesDone){
  speciesEvalInput <- read.table(speciesEval, header=T, sep="\t",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
  tmp <- speciesEvalInput %>% dplyr::summarise_if(is.numeric,mean,na.rm=T)
  tmp$species <- unique(speciesEvalInput$species)
  XMEvaluationsTotal <- rbind(XMEvaluationsTotal,tmp)
}

#Output MaxEnt model evaluation summary to a single file.
write.table(XMEvaluationsTotal,"LAIndicatorTaxa.txt",quote=FALSE,sep="\t",row.names = FALSE)
