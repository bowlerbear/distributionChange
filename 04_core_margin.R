#script to pull out core and marginal areas
library(raster)

### get data ####

#run script 01

### subset data #####

#just to test below code
modelSummaries$PA <- sapply(modelSummaries$mean, function(x) rbinom(1,1,x))

#add on lon and lat
modelSummaries$lon <- mtbsDF$lon[match(modelSummaries$MTB, mtbsDF$MTB)]
modelSummaries$lat <- mtbsDF$lat[match(modelSummaries$MTB, mtbsDF$MTB)]
modelSummaries <- subset(modelSummaries, !is.na(lon) & !is.na(lat))

#just take first and last year
modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))

### lists ####

allspecies <- sort(unique(modelSummaries_Limits$Species))

allyears <- sort(unique(modelSummaries_Limits$Year))

utmProj <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

### example species #####


getCoreRegions <- function(myspecies,myyear){

speciesSummary <- modelSummaries_Limits %>%
                  filter(Species == myspecies) %>%
                  filter(Year == myyear)

speciesPixels <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                        data = data.frame(PA=speciesSummary$PA),
                                        tolerance = 0.6,
                                        proj4string = crs("+proj=longlat +datum=WGS84"))

#make into a raster
speciesRaster <- raster(speciesPixels)

#for each raster cell, get mean class of surroundings
## in a 3-by-3 window surrounding each focal cell 
rmean <- focal(speciesRaster, 
               w = matrix(1, ncol=3, nrow=3), 
               fun=mean, na.rm=TRUE)
#plot(rmean)

#define core as all those fully surrounded by presences
speciesRasterCore <- rmean
speciesRasterCore[speciesRasterCore < 1] <- 0.5 #indicator for marginal sites
speciesRasterCore[is.na(speciesRaster)] <- NA
speciesRasterCore[speciesRaster==0] <- 0 #absent sites
plot(speciesRasterCore)

#convert back into a data frame
coreDF <- as.data.frame(speciesRasterCore, xy=TRUE)
coreDF$Obs <- as.data.frame(speciesRaster)[,1] #original observation
coreDF$cellNu <- cellFromXY(speciesRaster, coreDF[,c("x","y")])
coreDF <- subset(coreDF, !is.na(layer)) #remove sites beyond german border

#clean marginal information
coreDF$Core <- ifelse(coreDF$layer==1,"core",
                      ifelse(coreDF$layer==0.5, "marginal", "absent"))

#add on another data
coreDF$Year <- myyear
coreDF$Species <- myspecies

return(coreDF)

}

### apply function ###

getCoreRegions(allspecies[1],allyears[1])


#other possibilities

#define core and marginal areas
#core area is defined as all cells that are further away from the edge of each patch than a specified edge depth

#https://github.com/r-spatialecology/landscapemetrics/blob/main/R/show_cores.R

show_patches(speciesRaster)
show_cores(speciesRaster) #inside areas beyond an edge

#explore possibilities
show_lsm(speciesRaster, class="1", what = "lsm_p_area")
show_lsm(speciesRaster, class="1",what = "lsm_p_core")
show_lsm(speciesRaster, class="1",what = "lsm_p_cai", edge_depth=10)
show_lsm(speciesRaster, class="1",what = "lsm_p_frac")
show_lsm(speciesRaster, class="1",what = "lsm_p_contig")
show_lsm(speciesRaster, class="1",what = "lsm_p_contig")
show_lsm(speciesRaster, class="1",what = "lsm_p_para")
show_lsm(speciesRaster, class="1",what = "lsm_p_ncore")

#Get data frame
rasterCore <- spatialize_lsm(speciesRaster, what = "lsm_p_core")
rasterCai <- spatialize_lsm(speciesRaster, what = "lsm_p_cai")

plot(rasterCai$layer_1$lsm_p_cai)
plot(rasterCore$layer_1$lsm_p_core)
