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

### lists ####

allspecies <- sort(unique(modelSummaries_Limits$Species))

allYears <- sort(unique(modelSummaries_Limits$Year))

utmProj <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(50,1,x))
PA_matrix <- do.call(rbind,PAs)

### fragmentation ####

#change
fragChange <- lapply(allspecies[-c(1,64)],function(x){
  applyFragStats(x, modelSummaries_Limits)}) %>% 
  reduce(rbind) %>%
  filter(class == 1) %>%
  filter(name == "clumpiness index") 

saveRDS(fragChange, file="outputs/clumpiChange.rds")

#annual
fragAnnual <- lapply(allspecies[-c(1,64)],function(x){
  applyFragStats(x, modelSummaries_Limits, summary = "annual")}) %>% 
  reduce (rbind) %>%
  filter(class == 1) %>%
  filter(name == "clumpiness index") %>%
  group_by(Species) %>%
  filter(Species!="Gomphus flavipes") %>%
  mutate(change = (log(medianMetric[Year==2016])- log((medianMetric[Year==1990]))))

saveRDS(fragAnnual, file="outputs/clumpiAnnual.rds")

#plotting
fragAnnual %>%
  filter(medianMetric > (-1)) %>%
  ggplot(aes(x=Year,y=medianMetric,group=Species,colour=change))+
  scale_colour_gradient2(low="green",high="purple")+ 
  geom_point()+
  geom_line(size=2)+
  scale_x_continuous(breaks=c(1990,2016),labels=c(1990,2016))+
  ylab("Clumpiness")

ggsave("plots/clumpiChange.png")


### saturation #####

#need to sort out units

saturationDF <- lapply(allspecies[-c(1,64)],function(x){
  applySaturation(x, modelSummaries_Limits, summary = "change")
}) %>% reduce(rbind)

saveRDS(saturationDF, file="outputs/satuationChange.rds")

### relationship with AOCC ####

areaChanges <- readRDS("outputs/areaChanges.rds")

allChanges <- fragChanges %>%
  rename(species = "Species") %>%
  dplyr::select(species, change) %>%
  filter(!duplicated(species)) %>%
  ungroup() %>%
  inner_join(.,areaChanges,
             by=c("species"),
             suffix = c("_frag","_area"))

g1 <- ggplot(data = allChanges,
             aes(x = medianChange,y = change)) + 
  geom_point() + 
  geom_hline(yintercept=0,linetype="dashed")+
  stat_smooth(method="gam")+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in AOC") + ylab("Change in clumpiness")
g1

ggsave("plots/eoccChange_vs_clumpiChange.png")


allChanges <- saturationDF  %>%
  inner_join(.,areaChanges,
             by=c("species"),
             suffix = c("_sat","_area"))

g2 <- ggplot(data = allChanges,
             aes(x = medianChange_area,y = medianChange_sat)) + 
  geom_point() + 
  geom_hline(yintercept=0,linetype="dashed")+
  stat_smooth(method="gam")+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in AOC") + ylab("Change in saturation")

g2


plot_grid(g1,g2,ncol=2)

ggsave("plots/clumpiness_vs_saturation.png",width=8, height=3)

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