#http://datazone.birdlife.org/species/spcredcrit

#Severely fragmented - Severely fragmented refers to the situation where increased extinction risks to the #species result from the fact that most individuals within a species are found in small and relatively #isolated subpopulations. These small subpopulations may go extinct, with a reduced probability of #recolonisation.

# https://r-spatialecology.github.io/landscapemetrics/
# https://r-spatialecology.github.io/landscapemetrics/articles/articles/comparing_tools.html
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.04617

#install.packages("landscapemetrics")
library(landscapemetrics)
library(raster)
library(ggthemes)
library(cowplot)

theme_set(theme_few())

#as a fragmentation measure, we use the clumpiness index

### get data ####

source("01_getModels.R")
source("05_core_functions.R")

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(50,1,x))
PA_matrix <- do.call(rbind,PAs)

### fragmentation ####

# #of the distributions at each time step, and then the change of it
# 
# #change
# fragChange <- lapply(allspecies[-c(1,64)],function(x){
#   applyFragStats(x, modelSummaries_Limits)}) %>% 
#   reduce(rbind) %>%
#   filter(class == 1) %>%
#   filter(name == "clumpiness index") 
# 
# saveRDS(fragChange, file="outputs/clumpiChange.rds")
# 
# #annual
# fragAnnual <- lapply(allspecies[-c(1,64)],function(x){
#   applyFragStats(x, modelSummaries_Limits, summary = "annual")}) %>% 
#   reduce (rbind) %>%
#   filter(class == 1) %>%
#   filter(name == "clumpiness index") %>%
#   group_by(Species) %>%
#   filter(Species!="Gomphus flavipes") %>%
#   mutate(change = (log(medianMetric[Year==2016])- log((medianMetric[Year==1990]))))
# 
# saveRDS(fragAnnual, file="outputs/clumpiAnnual.rds")
# 
# #plotting
# fragAnnual %>%
#   filter(medianMetric > (-1)) %>%
#   ggplot(aes(x=Year,y=medianMetric,group=Species,colour=change))+
#   scale_colour_gradient2(low="green",high="purple")+ 
#   geom_point()+
#   geom_line(size=2)+
#   scale_x_continuous(breaks=c(1990,2016),labels=c(1990,2016))+
#   ylab("Clumpiness")
# 
# ggsave("plots/clumpiChange.png")

### change fragmentation ####

#first get the change - then get the fragmentation of that

fragChanges <- lapply(allspecies,function(x){
  applyFragChange(x, modelSummaries_Limits)}) %>%
  reduce(rbind) %>%
  rename(species = "Species")
saveRDS(fragChanges, file="outputs/clumpiChange_change.rds")


#boxplots
fragChanges <- readRDS("outputs/clumpiChange_change.rds")
fragChanges$Direction <- areaChanges$Direction[match(fragChanges$species,
                                                    areaChanges$species)]
  
(g1 <- fragChanges %>%
    filter(medianChange>0) %>%
    filter(!is.na(Direction)) %>% #check why we need this
    ggplot(aes(x = Direction, y = medianChange))+
    geom_pirate(aes(colour = Direction),bars=FALSE)+
    xlab("Direction of change") + ylab("Clumpiness of change (log-scale)")+
    scale_y_log10()+
    scale_colour_manual(values=c("deeppink","dodgerblue4")))

saveRDS(g1, "plots/boxplot_fragChange_change.rds")

### relationship with AOCC ####

areaChanges <- readRDS("outputs/areaChanges.rds")

allChanges <- fragChanges %>%
  ungroup() %>%
  inner_join(.,areaChanges,
             by=c("species"),
             suffix = c("_frag","_area"))

(g2 <- allChanges %>%
  filter(medianChange_frag>0) %>%
  filter(class =="change") %>%
  #filter(class %in% c("increase","decrease")) %>%
  ggplot(aes(x = medianChange_area, y = medianChange_frag)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange_frag,ymax = upperChange_frag)) + 
  geom_errorbarh(aes(xmin = lowerChange_area, xmax = upperChange_area)) +
  stat_smooth(method="gam", se=FALSE, colour="black")+
  #geom_hline(yintercept=0,linetype="dashed")+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in AOO") + ylab("Clumpiness of change"))

ggsave("plots/eoccChange_vs_clumpiChange_change.png",
       width = 7.5, height = 3.5)

### fig. 4 ####

plot_grid(g1,g2,nrow=2)
ggsave("plots/Fig.4.png",
       width = 6.5, height = 6)

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
