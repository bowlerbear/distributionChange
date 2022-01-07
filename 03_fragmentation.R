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

### apply ####

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(5,1,x))
PA_matrix <- do.call(rbind,PAs)

#apply it
fragDF <- lapply(allspecies[-1],function(x){
                 applyFragStats(x, modelSummaries_Limits, summary = "annual")
})
fragDF <- do.call(rbind,fragDF)

### plotting ####

fragDF <- fragDF %>%
  filter(class == 1) %>%
  filter(name == "clumpiness index") %>%
  group_by(Species) %>%
  mutate(change = log(medianMetric[Year==2016]/medianMetric[Year==1990]))
save(fragDF, file="outputs/clumpiChanges.rds")

fragDF%>%
  ggplot(aes(x=Year,y=medianMetric,group=Species,colour=change))+
  scale_colour_gradient2(low="green",high="purple")+ 
  geom_point(alpha=0.2)+
  geom_line(size=2,alpha=0.5)+
  scale_x_continuous(breaks=c(1990,2016),labels=c(1990,2016))+
  ylab("Clumpiness")

ggsave("plots/clumpiChange.png")

### rel with AOCC

fragDF <- lapply(allspecies[-1],function(x){
  applyFragStats(x, modelSummaries_Limits, summary = "change")
})
fragChanges <- do.call(rbind,fragDF)

areaChanges <- readRDS("outputs/areaChanges.rds")

allChanges <- inner_join(fragChanges,areaChanges,
                         by=c("species"),
                         suffix = c("_frag","_area"))

ggplot(data = allChanges,
       aes(x = medianChange_area,y = medianChange_frag)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange_frag,ymax = upperChange_frag)) + 
  geom_errorbarh(aes(xmin = lowerChange_area, xmax = upperChange_area))+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in AOCC") + ylab("Change in clumpiness")+
  theme_few()

ggsave("plots/eoccChange_vs_clumpiChange.png")


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

# #what is available
# 
# list_lsm()
# 
# list_lsm()$name
# 
# list_lsm(level = 'patch')
# 
# list_lsm(level = "class")$type
# 
# list_lsm(name = "total core area")
# 
# list_lsm(level = "landscape")$type
# 
# #visualizations
# 
# show_patches(speciesRaster, class="all")
# 
# show_cores(speciesRaster, class="all") #inside areas beyond an edge
# 
# show_lsm(speciesRaster, what = "lsm_p_contig")
# 
# #spatialize_lsm(speciesRaster, what = "lsm_p_contig")
# 
# #correlations
# 
# metrics <- calculate_lsm(speciesRaster, what = c("class"))
# correlations <- calculate_correlation(metrics)
# show_correlation(data = correlations, method = "pearson")
# 
# #apply a function
# lsm_c_clumpy(speciesRaster)
# lsm_c_clumpy(speciesStack)
# 
# #appply multiple functions
# calculate_lsm(speciesRaster, 
#               what = c("lsm_l_ta",
#                        "lsm_c_pland",
#                        "lsm_c_clumpy"),
#               full_name = TRUE)