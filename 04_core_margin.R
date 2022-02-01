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

coreDF <- allspecies[-63] %>%
        map(.,getCoreRegions) %>%
        reduce(rbind)

coreSummary <- coreDF %>%
                group_by(Species) %>%
                summarize(nuMarginal = sum(Core=="absent" | Core=="marginal"),
                        nuCore = sum(Core=="core"),
                        sum = length(Core),
                        ratio = nuMarginal/nuCore)%>%
                rename(species = "Species")

saveRDS(coreSummary, file="outputs/coreSummary.rds")

#annual
coreAnnual <- allspecies[-63] %>%
  map(~getCoreCalc(., summary="annual")) %>%
  reduce(rbind)
saveRDS(coreAnnual, file="outputs/coreAnnual.rds")


#apply function to get changes
coreChanges <- allspecies[-63] %>%
                map(~getCoreCalc(., summary="change")) %>%
                reduce(rbind)

saveRDS(coreChanges, file="outputs/coreChanges.rds")
#problem with "Sympetrum meridionale"

### plotting ####

coreChanges <- readRDS("outputs/coreChanges.rds")

(g1 <- coreChanges %>%
  mutate(Core = fct_relevel(Core, "core", "marginal")) %>%
  mutate(Species = fct_reorder(Species, desc(medianChange))) %>%
  ggplot()+
  geom_pointrange(aes(x = Species, y = medianChange, ymin = lowerChange, 
                      max = upperChange,
                      fill = Core),
                      shape=21)+
  scale_fill_brewer("Region",type="div")+
  #coord_flip()+
  geom_hline(yintercept=0, linetype="dashed") +
  ylab("Relative change in AOO")+
  theme(axis.text.x = element_blank(),
        legend.position="top"))
  
ggsave("plots/coreChange.png")

nrow(coreChanges_wide)

nrow(subset(coreChanges_wide, absent>0 & core >0 & marginal>0))

nrow(subset(coreChanges_wide, absent<0 & core <0 & marginal<0))

sd(coreChanges_wide$core)
sd(coreChanges_wide$marginal)
sd(coreChanges_wide$absent)

### relationships ####

coreChanges <- readRDS("outputs/coreChanges.rds")
areaChanges <- readRDS("outputs/areaChanges.rds")

coreChanges_full <- coreChanges %>%
                      rename(species = Species) %>%
                      inner_join(.,areaChanges,
                      by=c("species"),
                      suffix = c("_core","_total")) %>%
                      mutate(Core = fct_relevel(Core, "core", "marginal"))

(g2 <- ggplot(data = coreChanges_full,
       aes(y = medianChange_core, x = medianChange_total)) + 
  geom_point(aes(fill=Core), shape=21) + 
  geom_errorbarh(aes(xmin = lowerChange_total,xmax = upperChange_total, colour=Core)) + 
  geom_errorbar(aes(ymin = lowerChange_core, ymax = upperChange_core, colour=Core))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  geom_smooth(method="gam", aes(colour=Core), se=FALSE)+
  scale_fill_brewer("Region",type="div")+
  scale_colour_brewer("Region",type="div")+
  ylab("Change in regional AOO") + xlab("Change in total AOO")+
  theme_few()+
  theme(legend.position = "none"))

ggsave("plots/marginalChange_vs_coreChange.png", width = 4, height = 4)

### core summary ####

coreSummary <- readRDS("outputs/coreSummary.rds")
coreSummary$Direction <- areaChanges$Direction[match(coreSummary$species,
                                                      areaChanges$species)]
ggplot(coreSummary)+
  geom_boxplot(aes(x=Direction,y=ratio))+
  scale_y_log10()#not much difference!!
#yeah!


hist(coreSummary$nuMarginal)
hist(coreSummary$nuCore)
summary(coreSummary$nuMarginal)
summary(coreSummary$nuCore)

#summary(coreSummary$nuMarginal/coreSummary$sum)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.02043 0.40442 0.79976 0.66281 0.94268 0.99633 

#summary(coreSummary$nuCore/coreSummary$sum)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003672 0.057323 0.200239 0.337195 0.595579 0.979571 

### col vs ext rate ####

colextDF <- allspecies[-63] %>%
  map(~getColExt(.)) %>%
  reduce(rbind)

saveRDS(colextDF, file="outputs/colextDF.rds")

#boxplots
colextDF <- readRDS("outputs/colextDF.rds")

(g3 <- ggplot(colextDF)+
        geom_boxplot(aes(x=Type, y=medianProp, fill=Core))+
        ylab("Proportion of grids")+xlab("Occupancy change")+
        scale_fill_brewer("Region",type="div"))


colextDF$Direction <- areaChanges$Direction[match(colextDF$species,
                                                  areaChanges$species)]

(g4 <- ggplot(colextDF)+
    geom_boxplot(aes(x=Direction, y=medianProp, fill=Core))+
    ylab("Proportion of grids")+xlab("Region")+
    facet_wrap(~Type)+
    scale_fill_brewer("Region",type="div")+
    theme(legend.position="none"))

plot_grid(g1,g2,g4,nrow=3)

ggsave("plots/Fig.4.png",width=6,height=8)

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
