#run the modelSummaries_stan

#http://datazone.birdlife.org/species/spcredcrit

#Severely fragmented - Severely fragmented refers to the situation where increased extinction risks to the #species result from the fact that most individuals within a species are found in small and relatively #isolated subpopulations. These small subpopulations may go extinct, with a reduced probability of #recolonisation.

# https://r-spatialecology.github.io/landscapemetrics/
# https://r-spatialecology.github.io/landscapemetrics/articles/articles/comparing_tools.html
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.04617

#install.packages("landscapemetrics")
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

# speciesSummary <- subset(modelSummaries_Limits, 
#                          Species = allspecies[1] & Year == 2010)
# speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
#                                         data = data.frame(PA=speciesSummary$PA),
#                                         tolerance = 0.6,
#                                         proj4string = crs("+proj=longlat +datum=WGS84"))
# 
# #make into a raster
# speciesRaster <- raster(speciesRaster)
# 
# #convert to utm projection
# speciesRaster <- projectRaster(speciesRaster, crs=utmProj, method="ngb")
# plot(speciesRaster)

### function ####

#make raster stack for a species - layers represent each year
#just take first and last year
getFragStats <- function(species, modelSummaries_Limits){
  
    #data for 1990
    speciesSummary <- modelSummaries_Limits %>%
      filter(Species==species) %>%
      filter(Year == 1990)
    
    speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                            data = data.frame(PA=speciesSummary$PA),
                                            tolerance = 0.6,
                                            proj4string = crs("+proj=longlat +datum=WGS84"))
    
    #make into a raster
    r1 <- raster(speciesRaster)
    
    #data for 2016
    speciesSummary <- modelSummaries_Limits %>%
      filter(Species==species) %>%
      filter(Year == 2016)
    
    speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                            data = data.frame(PA=speciesSummary$PA),
                                            tolerance = 0.6,
                                            proj4string = crs("+proj=longlat +datum=WGS84"))
    
    #make into a raster
    r2 <- raster(speciesRaster)
    

  #make stack  
  speciesStack <- stack(list(r1,r2))
  speciesStack <- projectRaster(speciesStack, crs=utmProj, method="ngb")
  
  #calc metrics
  calculate_lsm(speciesStack, 
                what = c("lsm_l_ta",
                         "lsm_c_pland",
                         "lsm_c_clumpy"),
                full_name = TRUE) %>%
    add_column(Year = c(rep(1990,5),rep(2016,5))) %>%
    add_column(Species = species)
}


### apply ####

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(5,1,x))
PA_matrix <- do.call(rbind,PAs)

#apply it
fragDF <- lapply(allspecies,function(x){
                 applyFragStats(x, modelSummaries_Limits, summary = "annual")
})
fragDF <- do.call(rbind,fragDF)

### Appendix: metric play ####

#https://cran.r-project.org/web/packages/landscapemetrics/vignettes/getstarted.html

#a patch being defined as contiguous cells belonging to the same land-cover class

# # general structure
# lsm_"level"_"metric"
# 
# # Patch level
# ## lsm_p_"metric"
# lsm_p_enn()
# 
# # Class level
# ## lsm_c_"metric"
# lsm_c_enn()
# 
# # Landscape level
# ## lsm_p_"metric"
# lsm_l_enn()

#what is available

list_lsm()

list_lsm()$name

list_lsm(level = 'patch')

list_lsm(level = "class")$type

list_lsm(name = "total core area")

list_lsm(level = "landscape")$type

#visualizations

show_patches(speciesRaster, class="all")

show_cores(speciesRaster, class="all") #inside areas beyond an edge

show_lsm(speciesRaster, what = "lsm_p_contig")

#spatialize_lsm(speciesRaster, what = "lsm_p_contig")

#correlations

metrics <- calculate_lsm(speciesRaster, what = c("class"))
correlations <- calculate_correlation(metrics)
show_correlation(data = correlations, method = "pearson")

#apply a function
lsm_c_clumpy(speciesRaster)
lsm_c_clumpy(speciesStack)

#appply multiple functions
calculate_lsm(speciesRaster, 
              what = c("lsm_l_ta",
                       "lsm_c_pland",
                       "lsm_c_clumpy"),
              full_name = TRUE)