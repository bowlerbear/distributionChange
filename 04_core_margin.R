#script to pull out core and marginal areas
library(raster)
theme_set(theme_few())

### get data ####

source("01_getModels.R")
source("05_core_functions.R")

#get realizations of occupancy
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(100,1,x))
PA_matrix <- do.call(rbind,PAs)

### define core or marginal ####

coreDF <- allspecies %>%
        map(.,getCoreRegions) %>%
        reduce(rbind)

#apply function to get changes
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
  ylab("Relative change in AOO")+
  theme(axis.text.y = element_text(size=rel(0.5)))
  
ggsave("plots/coreChange.png")

### relationships ####

coreChanges <- readRDS("outputs/coreChanges.rds")
areaChanges <- readRDS("outputs/areaChanges.rds")

#comparison to core
allChanges <- coreChanges %>%
              filter(Core=="core") %>%
              rename(species = Species) %>%
              inner_join(.,areaChanges,
                         by=c("species"),
                         suffix = c("_core","_total"))

ggplot(data = allChanges,
       aes(x = medianChange_core, y = medianChange_total)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange_total,ymax = upperChange_total)) + 
  geom_errorbarh(aes(xmin = lowerChange_core, xmax = upperChange_core))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in core AOCC") + ylab("Change in total AOCC")+
  theme_few()

ggsave("plots/coreChange_vs_aoccChange.png")

#comparison to marginal
allChanges <- coreChanges %>%
  filter(Core=="marginal") %>%
  rename(species = Species) %>%
  inner_join(.,arealChanges,
             by=c("species"),
             suffix = c("_core","_total"))

ggplot(data = allChanges,
       aes(x = medianChange_core, y = medianChange_total)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange_total,ymax = upperChange_total)) + 
  geom_errorbarh(aes(xmin = lowerChange_core, xmax = upperChange_core))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in marginal AOCC") + ylab("Change in EOCC")+
  theme_few()

ggsave("plots/marginalChange_vs_aoccChange.png")

### summary statistics ####

nrow(coreChanges_wide)

nrow(subset(coreChanges_wide, absent>0 & core >0 & marginal>0))

nrow(subset(coreChanges_wide, absent<0 & core <0 & marginal<0))

sd(coreChanges_wide$core)
sd(coreChanges_wide$marginal)
sd(coreChanges_wide$absent)

### more plots #####

#also merge with area changes
coreChanges_full <- coreChanges %>%
                      rename(species = Species) %>%
                      inner_join(.,areaChanges,
                      by=c("species"),
                      suffix = c("_core","_total")) %>%
                      mutate(Core = fct_relevel(Core, "core", "marginal", "absent"))

ggplot(data = coreChanges_full,
       aes(y = medianChange_core, x = medianChange_total, colour=Core)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = lowerChange_total,xmax = upperChange_total, colour=Core)) + 
  geom_errorbar(aes(ymin = lowerChange_core, ymax = upperChange_core, colour=Core))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  geom_smooth(method="gam", aes(colour=Core), se=FALSE)+
  scale_color_viridis_d()+
  ylab("Change in regional AOO") + xlab("Change in total AOO")+
  theme_few()+
  theme(legend.position = "top")

ggsave("plots/marginalChange_vs_coreChange.png", width = 4, height = 4)

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
