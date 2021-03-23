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
#Filter observations by positional uncertainty.
SpeciesInput <- SpeciesInput[SpeciesInput$coordinateUncertaintyInMeters <= 10,]

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

#Get list of Los Angeles county native plants from CalFlora.
LAPlants <- read.table("CalPlants.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
#Get list of Los Angeles county native animals as developed for the LA city biodiversity index.
LAAnimals <- read.table("CalAnimals.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
#Create a merged LA natives list.
LASpecies <- rbind(LAPlants,LAAnimals)

#Filter data to only contain LA natives.
SpeciesLocations <- SpeciesLocations[SpeciesLocations$species %in% LASpecies$species,]

#Filter LA species observations by time and prevalence.
SpeciesFreq <- as.data.frame(table(SpeciesLocations[SpeciesLocations$year<=2020 & SpeciesLocations$year>=2010,"species"]))
Prevalence <- 30
SpeciesFreq <- SpeciesFreq[SpeciesFreq$Freq >= Prevalence,]
colnames(SpeciesFreq) <- c("species","Freq")
SpeciesFreq$species <- as.character(SpeciesFreq$species)

#Create a list of species above a certain number of observation points.
speciesList <- SpeciesFreq$species

#Get the list of species already evaluated.
speciesDone <- list.files(pattern="MaxentEvaluationNatives(.*?).txt")
speciesDone <- gsub("MaxentEvaluationNatives","",gsub(".txt","",speciesDone))
#Remove the species already evaluated from the next iteration of analysis.
speciesList <-speciesList[!(speciesList %in% speciesDone)]

#Filter out species observation points to only include species with sufficient observational data.
all.presvals <- SpeciesLocations[SpeciesLocations$species %in% speciesList,]

#Specify sampling number for MaxEnt input data.
sampleNum <- 25

#Remove environmental layers with a high degree of multicollinearity.
#Only run this once to generate a list of layers with low multicollinearity.
#List environmental rasters
#env.files <- list.files(path=paste(wd,"/envLayers",sep=""),pattern=".tif$",full.names=TRUE)
#Stack environmental layers
#env.data <- stack(c(env.files))
#Get environmental layer names
#env.filenames <- gsub(paste("^",wd,"/envLayers/",sep=""),"",gsub(".tif","",env.files))
#env.filtered <- removeCollinearity(env.data,nb.points = 10000,sample.points = T,select.variables = T,multicollinearity.cutoff = 0.6)
#env.filtered <- as.data.frame(env.filtered)
#write.table(env.filtered,"EnvFilteredV1.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Create a filtered raster stack.
env.filtered <- read.table("EnvFilteredV1.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
env.filtered <- c(paste(wd,"/envLayers/",env.filtered$env.filtered,".tif",sep=""))
env.data <- stack(c(env.filtered))

SDM <- function(i) {
  #Randomly select a species from the list of those not already analyzed.
  #Get the list of species already evaluated.
  speciesDone <- list.files(pattern="MaxentEvaluationNatives(.*?).txt")
  speciesDone <- gsub("MaxentEvaluationNatives","",gsub(".txt","",speciesDone))
  #Remove the species already evaluated from the next iteration of analysis.
  speciesList <-speciesList[!(speciesList %in% speciesDone)]
  species <- speciesList[sample(length(speciesList))[1]]
  #Extract observation locations for a given species.
  obs.data <- SpeciesLocations[SpeciesLocations$species==species,c("longitude.1","latitude.1")]
  
  # Initialize data containing environmental layer values at presence locations.
  presvals <- as.data.frame(raster::extract(env.data,obs.data))
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
  absvals <- as.data.frame(raster::extract(env.data,as.data.frame(spsample(LABoundaries,n=20*nrow(presvals),"random"))))
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
  envFactors <- c("cal_fire","ClimateZones","DLC","DominantCanopyCover","EcoRegion","Ecotopes","FloodPlain","LandCover","LandUse","PublicLandStatus","Totrcv","Vegcover","WHR","WildlandUrbanInterface")
  modelData[,colnames(modelData) %in% envFactors] <- lapply(modelData[,colnames(modelData) %in% envFactors],factor)
  env.filtered <- gsub(paste("^",wd,"/envLayers/",sep=""),"",gsub(".tif","",env.filtered))
  modelData <- modelData[,c(env.filtered,"pa")]
  
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
    xm <- maxent(x=modelSubsample[,env.filtered],p=modelSubsample$pa,factors=names(Filter(is.factor, modelData)))
    #Get relative importance of environmental layers.
    XMImportance <- as.data.frame(t(var.importance(xm)))
    colnames(XMImportance) <- as.character(unlist(XMImportance[1,]))
    XMImportance <- as.data.frame(XMImportance[2,])
    XMImportance[] <- lapply(XMImportance, function(x) as.numeric(as.character(x)))
    #Evaluate maxent model.
    exm <- suppressWarnings(evaluate(p=modelSubsample[modelSubsample$pa==1,c(env.filtered)],a=modelSubsample[modelSubsample$pa==0,c(env.filtered)],model=xm))
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
    tmp <- cbind(tmp,XMImportance)
    XMEvaluations <- rbind(XMEvaluations,tmp)
    print(paste(species,"auc:",tmp$AUC,"SEDI:",tmp$SEDI,"cor:",tmp$r,"p:",tmp$p))
  }
  write.table(XMEvaluations,paste("MaxentEvaluationNatives",species,".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
}

#Run MaxEnt evaluations on all of the species.
lapply(1:length(speciesList),SDM)

##Get the list of species already evaluated.
speciesDone <- list.files(pattern="MaxentEvaluationNatives(.*?).txt")
##Summarize all evaluations so each row is a unique species with all of its associated maxent model evaluations.
XMEvaluationsTotal <- data.frame()
for(speciesEval in speciesDone){
  speciesEvalInput <- read.table(speciesEval, header=T, sep="\t",as.is=T,skip=0,fill=T,quote="\"",check.names=F,encoding = "UTF-8")
  tmp <- speciesEvalInput %>% dplyr::summarise_if(is.numeric,mean,na.rm=T)
  tmp$species <- unique(speciesEvalInput$species)
  XMEvaluationsTotal <- rbind(XMEvaluationsTotal,tmp)
}

#Output MaxEnt model evaluation summary to a single file.
write.table(XMEvaluationsTotal,"LAIndicatorTaxaNatives.txt",quote=FALSE,sep="\t",row.names = FALSE)

XMEvaluationsTotal <- read.table("LAIndicatorTaxaNatives.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
