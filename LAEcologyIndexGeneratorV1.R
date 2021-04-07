rm(list=ls())
require(rgdal)
require(raster)
require(foreach)
require(ENMeval)
require(spThin)

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

set.seed(1)
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

#Filter LA species data by prevalence.
SpeciesLocations <- SpeciesLocations[SpeciesLocations$species %in% names(table(SpeciesLocations$species))[table(SpeciesLocations$species) >= 30],]

#Spatially thin species observations.
SpeciesThinned <- data.frame()
for(taxon in unique(SpeciesLocations$species)){
  TaxonLocations <- SpeciesLocations[SpeciesLocations$species==taxon,]
  TaxonThinned <- as.data.frame(thin(TaxonLocations,lat.col="latitude",long.col="longitude",spec.col="species",thin.par=0.5,reps=100,write.files=F,write.log.file=F,locs.thinned.list.return=T,out.dir=wd)[100])
  SpeciesThinnedTmp <- TaxonLocations[TaxonLocations$latitude %in% TaxonThinned$Latitude & TaxonLocations$longitude %in% TaxonThinned$Longitude,]
  if(nrow(SpeciesThinnedTmp) >= 30){
    SpeciesThinned <- rbind(SpeciesThinned,SpeciesThinnedTmp)
  }
}

#Read in the list of species, and their MaxEnt model evaluations, as generated from here: https://github.com/LASanitation/LASAN/blob/main/LAIndicatorTaxaV1.R
XMEvaluationsTotal <- read.table("LAIndicatorTaxa.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")

#Merge in evaluation data with species observations.
SpeciesThinned <- dplyr::left_join(SpeciesThinned,XMEvaluationsTotal)

#Create a list so only the species with the highest SEDI scores are retained as indicators.
tmp <- SpeciesThinned[,c("SEDI","species")]
tmp <- tmp[!duplicated(tmp),]
speciesList <- dplyr::top_n(tmp,100,SEDI)$species
SpeciesThinned <- SpeciesThinned[SpeciesThinned$species %in% speciesList,]

#This filtered environmental list was generated here: https://github.com/LASanitation/LASAN/blob/main/LAIndicatorTaxaV1.R
env.filtered <- read.table("EnvFilteredV1.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")

#Filter environmental list further by their mean relative importance in 
#maxent models for all indicator species.
XMLayersSummary <- as.data.frame(colMeans(SpeciesThinned[,env.filtered$env.filtered]))
colnames(XMLayersSummary) <- c("MeanImportance")
XMLayersSummary$variable <- row.names(XMLayersSummary)
env.filtered.further <- XMLayersSummary[XMLayersSummary$MeanImportance >= 3,c("variable")]
write.table(as.data.frame(env.filtered.further),"EnvFilteredFurtherV1.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Create a filtered raster stack.
env.filtered <- c(paste(wd,"/envLayers/",env.filtered.further,".tif",sep=""))
env.data <- stack(c(env.filtered))

#Read in bias raster, generated here: https://github.com/LASanitation/LASAN/blob/main/LABiasRasterV1.R
dens.ras <- raster("SpeciesBiasV1.tif")

#Remove the species already evaluated from the next iteration of analysis.
speciesDone <- list.files(pattern="Prediction(.*?).tif$",full.names=T)
speciesList <-speciesList[!(speciesList %in% speciesDone)]

SDM <- function(i) {
  #Randomly select a species from the list of those not already analyzed.
  set.seed(as.numeric(Sys.time()))
  #Get the list of species already evaluated.
  speciesDone <- list.files(pattern="Prediction(.*?).tif$",full.names=T)
  speciesDone <- gsub("^./Prediction","",gsub(".tif","",speciesDone))
  #Remove the species already evaluated from the next iteration of analysis.
  speciesList <-speciesList[!(speciesList %in% speciesDone)]
  species <- speciesList[sample(length(speciesList))[1]]
  
  #Extract observation locations for a given species.
  obs.data <- SpeciesThinned[SpeciesThinned$species==species,c("longitude.1","latitude.1")]
  
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
  #Create background points using bias background layer.
  bg <- as.data.frame(xyFromCell(dens.ras, sample(which(!is.na(values(dens.ras))), 20*nrow(presvals), prob=values(dens.ras)[!is.na(values(dens.ras))])))
  colnames(bg) <- c("longitude.1","latitude.1")
  absvals <- as.data.frame(raster::extract(env.data,bg))
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
  #env.filtered <- gsub(paste("^",wd,"/envLayers/",sep=""),"",gsub(".tif","",env.filtered))
  #modelData <- modelData[,c(env.filtered,"pa")]
  
  #Get geographic extent for running the SDM.
  testExtent <- extent(env.data[[1]])
  #Run MaxEnt model
  xm <- maxent(x=modelData[,c(env.filtered.further)],p=modelData$pa,factors=names(Filter(is.factor, modelData)))
  #f <- list("Ecotopes"=levels(modelData$Ecotopes),"FloodPlain"=levels(modelData$FloodPlain),"LandUse"=levels(modelData$LandUse),"ClimateZones"=levels(modelData$ClimateZones),"FloodPlain"=levels(modelData$FloodPlain),"PublicLandStatus"=levels(modelData$PublicLandStatus),"WildlandUrbanInterface"=levels(modelData$WildlandUrbanInterface),"LandCover"=levels(modelData$LandCover),"DominantCanopyCover"=levels(modelData$DominantCanopyCover),"PotentialNaturalVegetation"=levels(modelData$PotentialNaturalVegetation),"DLC"=levels(modelData$DLC),"Vegcover"=levels(modelData$Vegcover),"Totrcv"=levels(modelData$Totrcv),"WHR"=levels(modelData$WHR),"cal_fire"=levels(modelData$cal_fire))
  #Evaluate maxent model.
  exm <- suppressWarnings(evaluate(modelData[modelData$pa==1,c(env.filtered.further)],modelData[modelData$pa==0,c(env.filtered.further)],xm))
  #Evaluate probability threshold of species detection
  bc.threshold <- threshold(x = exm, stat = "spec_sens")
  
  r <- dismo::predict(object=xm,x=env.data,extent=testExtent,progress='text')
  #Convert prediction probability raster to a presence/absence prediction.
  rPA <- r > bc.threshold
  writeRaster(rPA,paste("Prediction",species,".tif",sep=""),overwrite=T)
  #Summarize all map layers into a single layer representing the LA Ecological Index.
  files=list.files(pattern="Prediction(.*?).tif$",full.names=T)
  if(length(files)==100){
    rs <- raster::stack(files)
    rs1 <- raster::calc(rs,sum,na.rm=T)
    writeRaster(rs1,"summary100.tif") 
  }
}

#Run MaxEnt evaluations on all of the species.
lapply(1:length(speciesList),SDM)
