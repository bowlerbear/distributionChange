#script to pull out core and marginal areas
library(raster)

### get data ####

#run script 01

source("05_core_functions.R")

### subset data #####

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

### apply function #####

temp <- getCoreRegions(allspecies[5])

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(100,1,x))
PA_matrix <- do.call(rbind,PAs)

#add on data
modelSummaries_Limits$PA <- PA_matrix[,1] 
modelSummaries_Limits$Core <- temp$Core[match(modelSummaries_Limits$MTB, temp$MTB)]

#summarize
coreSummary <- modelSummaries_Limits %>%
                  filter(!is.na(Core)) %>%
                  group_by(Core, Year) %>%
                  summarize(occ = sum(PA), total= length(PA)) %>%
                  mutate(prop = occ/total)


### appendix ####
#other possibilities

#define core and marginal areas
#core area is defined as all cells that are further away from the edge of each patch than a specified edge depth

#https://github.com/r-spatialecology/landscapemetrics/blob/main/R/show_cores.R

# show_patches(speciesRaster)
# show_cores(speciesRaster) #inside areas beyond an edge
# 
# #explore possibilities
# show_lsm(speciesRaster, class="1", what = "lsm_p_area")
# show_lsm(speciesRaster, class="1",what = "lsm_p_core")
# show_lsm(speciesRaster, class="1",what = "lsm_p_cai", edge_depth=10)
# show_lsm(speciesRaster, class="1",what = "lsm_p_frac")
# show_lsm(speciesRaster, class="1",what = "lsm_p_contig")
# show_lsm(speciesRaster, class="1",what = "lsm_p_contig")
# show_lsm(speciesRaster, class="1",what = "lsm_p_para")
# show_lsm(speciesRaster, class="1",what = "lsm_p_ncore")
# 
# #Get data frame
# rasterCore <- spatialize_lsm(speciesRaster, what = "lsm_p_core")
# rasterCai <- spatialize_lsm(speciesRaster, what = "lsm_p_cai")
# 
# plot(rasterCai$layer_1$lsm_p_cai)
# plot(rasterCore$layer_1$lsm_p_core)
