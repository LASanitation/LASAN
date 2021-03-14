rm(list=ls())
require(rgdal)
require(raster)
require(foreach)
require(ENMeval)

#Set your working directory.  Your's will be different on your machine.
wd <- "~/Desktop/Practicum/LASAN/Code"
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

#Generate random background points for MaxEnt modeling.
bgPoints <- as.data.frame(spsample(LABoundaries,n=100000,"random"))
#Get all species locations.
all.obs.data <- SpeciesLocations[,c("longitude.1","latitude.1")]
#Remove any potential overlap between random background points and actual observations.
bgPoints <- bgPoints[!(bgPoints$x %in% all.obs.data$decimalLongitude.1) & !(bgPoints$y %in% all.obs.data$decimalLatitude.1),]
colnames(bgPoints) <- c("longitude.1","latitude.1")

#List environmental rasters
env.files <- list.files(path=paste(wd,"/envLayers",sep=""),pattern=".tif$",full.names=TRUE)
#Stack environmental layers
env.data <- stack(c(env.files))
#Get environmental layer names
env.filenames <- gsub("^./","",gsub(".tif","",env.files))

#Read in the list of species, and their MaxEnt model evaluations, as generated from here: https://github.com/LASanitation/LASAN/blob/main/LAIndicatorTaxaV2.R
speciesList <- read.table("LAIndicatorTaxaNatives.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
#Filter this list so only the species with the highest SEDI scores are retained as indicators.
speciesList <- dplyr::top_n(speciesList,100,SEDI)$species

SDM <- function(i) {
  #Randomly select a species from the list of those not already analyzed.
  #Get the list of species already evaluated.
  speciesDone <- list.files(pattern="PredictionNatives(.*?).tif$",full.names=T)
  speciesDone <- gsub("^./PredictionNatives","",gsub(".tif","",speciesDone))
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
  
  #Get geographic extent for running the SDM.
  testExtent <- extent(env.data$Aspect)
  #Run MaxEnt model
  xm <- maxent(x=modelData[,c(env.filenames)],p=modelData$pa,factors=c('Ecotopes','FloodPlain','LandUse','ClimateZones','FloodPlain','PublicLandStatus','WildlandUrbanInterface','LandCover','DominantCanopyCover','PotentialNaturalVegetation','DLC','Vegcover','Totrcv','WHR','cal_fire'))
  f <- list("Ecotopes"=levels(modelData$Ecotopes),"FloodPlain"=levels(modelData$FloodPlain),"LandUse"=levels(modelData$LandUse),"ClimateZones"=levels(modelData$ClimateZones),"FloodPlain"=levels(modelData$FloodPlain),"PublicLandStatus"=levels(modelData$PublicLandStatus),"WildlandUrbanInterface"=levels(modelData$WildlandUrbanInterface),"LandCover"=levels(modelData$LandCover),"DominantCanopyCover"=levels(modelData$DominantCanopyCover),"PotentialNaturalVegetation"=levels(modelData$PotentialNaturalVegetation),"DLC"=levels(modelData$DLC),"Vegcover"=levels(modelData$Vegcover),"Totrcv"=levels(modelData$Totrcv),"WHR"=levels(modelData$WHR),"cal_fire"=levels(modelData$cal_fire))
  #Evaluate maxent model.
  exm <- suppressWarnings(evaluate(modelData[modelData$pa==1,c(env.filenames)],modelData[modelData$pa==0,c(env.filenames)],xm))
  #Evaluate probability threshold of species detection
  bc.threshold <- threshold(x = exm, stat = "spec_sens")
  
  r <- raster::predict(env.data, xm,extent=testExtent,na.rm=TRUE,inf.rm=TRUE,factors=f,progress='text')
  #Convert prediction probability raster to a presence/absence prediction.
  rPA <- r > bc.threshold
  writeRaster(rPA,paste("PredictionNatives",species,".tif",sep=""),overwrite=T)
  #Summarize all map layers into a single layer representing the LA Ecological Index.
  files=list.files(pattern="PredictionNatives(.*?).tif$",full.names=T)
  if(length(files)==100){
    rs <- raster::stack(files)
    rs1 <- raster::calc(rs,sum,na.rm=T)
    writeRaster(rs1,"summaryNatives100.tif") 
  }
}

#Run MaxEnt evaluations on all of the species.
lapply(1:length(speciesList),SDM)
