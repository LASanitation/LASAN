rm(list=ls())
require(rgdal)
require(raster)
require(dismo) # interface with MaxEnt
require(raster) # spatial data manipulation
require(MASS) # for 2D kernel density function
require(magrittr) # for piping functionality, i.e., %>%
require(maptools) # reading shapefiles

#Set your working directory.  Your's will be different on your machine.
wd <- "/Users/levisimons/Desktop/Practicum/LASAN/Code"
#wd <- "/home1/alsimons/LASAN"
setwd(wd)

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

#Filter LA species observations by time and prevalence.
SpeciesFreq <- as.data.frame(table(SpeciesLocations[SpeciesLocations$year<=2020 & SpeciesLocations$year>=2010,"species"]))
Prevalence <- 30
SpeciesFreq <- SpeciesFreq[SpeciesFreq$Freq >= Prevalence,]
colnames(SpeciesFreq) <- c("species","Freq")
SpeciesFreq$species <- as.character(SpeciesFreq$species)
SpeciesLocations <- SpeciesLocations[SpeciesLocations$species %in% unique(SpeciesFreq$species),]

#Create a filtered raster stack.
env.filtered <- read.table("EnvFilteredV1.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote="", encoding = "UTF-8")
env.filtered <- c(paste(wd,"/envLayers/",env.filtered$env.filtered,".tif",sep=""))
env.data <- stack(c(env.filtered))

#Create a bias file for illustrating the observation density bias in species observations.
#This is for downstream use in Maxent modeling for generating background points.
occur.ras <- rasterize(SpeciesLocations[,c("longitude.1","latitude.1")],env.data$Aspect,1)
occur.ras.tmp <- raster::crop(occur.ras,extent(LABoundaries))
occur.ras.LA <- raster::mask(occur.ras.tmp,LABoundaries)
dens <- kde2d(SpeciesLocations$longitude.1, SpeciesLocations$latitude.1, n = c(nrow(occur.ras.LA), ncol(occur.ras.LA)))
dens.ras <- raster(dens)
crs(dens.ras) <- CRS(proj4string(LABoundaries))
writeRaster(dens.ras, "SpeciesBiasV1tmp.tif")
#To clip and align this raster run the following command in GDAL:
# gdalwarp -srcnodata -3.39999999999999996e+38 -dstnodata -32768 -of GTiff -co COMPRESS=LZW -t_srs EPSG:2229 -cutline LACityBoundaries.shp -crop_to_cutline -co BIGTIFF=YES -tap -tr 30 30 SpeciesBiasV1tmp.tif SpeciesBiasV1.tif
# Output file: SpeciesBiasV1.tif

#Remove environmental layers with a high degree of multicollinearity.
#List environmental rasters
env.files <- list.files(path=paste(wd,"/envLayers",sep=""),pattern=".tif$",full.names=TRUE)
#Stack environmental layers
env.data <- stack(c(env.files))
#Get environmental layer names
env.filenames <- gsub(paste("^",wd,"/envLayers/",sep=""),"",gsub(".tif","",env.files))
#Using r=0.5 as the cutoff: https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1600-0587.2010.06229.x
env.filtered <- removeCollinearity(env.data,nb.points = 10000,sample.points = T,select.variables = T,multicollinearity.cutoff = 0.5)
env.filtered <- as.data.frame(env.filtered)
#Save filtered list of environmental variables.
write.table(env.filtered,"EnvFilteredV1.txt",quote=FALSE,sep="\t",row.names = FALSE)