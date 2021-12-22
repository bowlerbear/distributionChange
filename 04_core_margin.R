#script to pull out core and marginal areas

library(landscapemetrics)
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

allYears <- sort(unique(modelSummaries_Limits$Year))

utmProj <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

### example species #####

speciesSummary <- subset(modelSummaries_Limits,
                         Species = allspecies[1] & Year == 2010)
speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                        data = data.frame(PA=speciesSummary$PA),
                                        tolerance = 0.6,
                                        proj4string = crs("+proj=longlat +datum=WGS84"))

#make into a raster
speciesRaster <- raster(speciesRaster)

#convert to utm projection
speciesRaster <- projectRaster(speciesRaster, crs=utmProj, method="ngb")
plot(speciesRaster)

#define core and marginal areas
#core area is defined as all cells that are further away from the edge of each patch than a specified edge depth

#https://github.com/r-spatialecology/landscapemetrics/blob/main/R/show_cores.R

show_patches(speciesRaster)
show_cores(speciesRaster) #inside areas beyond an edge

show_lsm(speciesRaster,what = "lsm_p_area")
show_lsm(speciesRaster,what = "lsm_p_core")
show_lsm(speciesRaster,what = "lsm_p_cai", edge_depth=1)
show_lsm(speciesRaster,what = "lsm_p_frac")


lsm_p_core(speciesRaster)
temp <- lsm_p_cai(speciesRaster)

rasterCore <- spatialize_lsm(speciesRaster, what = "lsm_p_core")
rasterCai <- spatialize_lsm(speciesRaster, what = "lsm_p_cai")

plot(rasterCai$layer_1$lsm_p_cai)
plot(rasterCore$layer_1$lsm_p_core)

myRaster <- spatialize_lsm(speciesRaster, what = "lsm_p_area")
plot(myRaster$layer_1$lsm_p_area)

speciesRaster[speciesRaster==0] <- NA
