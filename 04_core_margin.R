#script to pull out core and marginal areas
library(raster)

### get data ####

#run script 01

source("05_core_functions.R")
theme_set(theme_few())

### subset data #####

#add on lon and lat
modelSummaries$lon <- mtbsDF$lon[match(modelSummaries$MTB, mtbsDF$MTB)]
modelSummaries$lat <- mtbsDF$lat[match(modelSummaries$MTB, mtbsDF$MTB)]
modelSummaries <- subset(modelSummaries, !is.na(lon) & !is.na(lat))

#sort species name
modelSummaries$Species <- as.character(sapply(modelSummaries$Species, function(x) strsplit(x,"_")[[1]][1]))

#just take first and last year
modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))

### lists ####

allspecies <- sort(unique(modelSummaries_Limits$Species))

allyears <- sort(unique(modelSummaries_Limits$Year))

utmProj <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

### apply function #####

#define MTBs as core or marginal
coreDF <- allspecies %>%
        map(.,getCoreRegions) %>%
        reduce(rbind)

#get realizations of occupancy
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(100,1,x))
PA_matrix <- do.call(rbind,PAs)

#apply function
coreChanges <- allspecies %>%
                map(~getCoreCalc(., summary="change")) %>%
                reduce(rbind)

saveRDS(coreChanges, file="outputs/coreChanges.rds")

### plotting ####

coreChanges <- readRDS("outputs/coreChanges.rds")

coreChanges %>%
  mutate(Core = fct_relevel(Core, "core", "marginal", "absent")) %>%
  mutate(Species = fct_reorder(Species, desc(medianChange))) %>%
  ggplot()+
  geom_pointrange(aes(x = Species, y = medianChange, ymin = lowerChange, 
                      max = upperChange,
                      colour = Core))+
  scale_color_viridis_d()+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed") +
  ylab("Relative changes")+
  theme(axis.text.y = element_text(size=rel(0.7)))
  

ggsave("plots/coreChange.png")

### relationships ####

hullChanges <- readRDS("outputs/hullChanges.rds")

#comparison to core
allChanges <- coreChanges %>%
              filter(Core=="core") %>%
              rename(species = Species) %>%
              inner_join(.,hullChanges,
                         by=c("species"),
                         suffix = c("_core","_extent"))

ggplot(data = subset(allChanges, medianChange_extent<5),
       aes(x = medianChange_core, y = medianChange_extent)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange_extent,ymax = upperChange_extent)) + 
  geom_errorbarh(aes(xmin = lowerChange_core, xmax = upperChange_core))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in core AOCC") + ylab("Change in EOCC")+
  theme_few()

ggsave("plots/coreChange_vs_aoccChange.png")


#comparison to marginal
allChanges <- coreChanges %>%
  filter(Core=="marginal") %>%
  rename(species = Species) %>%
  inner_join(.,hullChanges,
             by=c("species"),
             suffix = c("_core","_extent"))

ggplot(data = subset(allChanges, medianChange_extent<5),
       aes(x = medianChange_core, y = medianChange_extent)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange_extent,ymax = upperChange_extent)) + 
  geom_errorbarh(aes(xmin = lowerChange_core, xmax = upperChange_core))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in marginal AOCC") + ylab("Change in EOCC")+
  theme_few()

ggsave("plots/marginalChange_vs_aoccChange.png")


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
