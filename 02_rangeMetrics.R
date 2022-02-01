#the script asks whether species about:
#range filler or range expander, and 
#change in number of occupied sites versus change in areas of convex hull

#https://help.natureserve.org/biotics/content/record_management/Element_Files/Element_Ranking/ERANK_Definitions_of_Extent_of_Occurrence_and_Area_of_Occupancy.htm
#The range extent will be overestimated by using a minimum convex polygon (also called a convex hull) to calculate the area. In these cases, the α-hull is recommended. The α-hull can be estimated by making a Delauney triangulation of the points in a sample (connect the points with lines, constrained so that no lines intersect), and then deleting lines that are longer than two times the average line length. The range extent is then the sum of enclosed areas. 
#https://www.ala.org.au/spatial-portal-help/aoo/

#Joppa, L. N., Butchart, S. H. M., Hoffmann, M., Bachman, S. P., Akçakaya, H. R., Moat, J. F., Böhm, M., #Holland, R. A., Newton, A., Polidoro, B. and Hughes, A. (2016), Impact of alternative metrics on #estimates of extent of occurrence for extinction risk assessment. Conservation Biology, 30: 362–370. doi:10.1111/cobi.12591

#dealing with outliers
#https://search.r-project.org/CRAN/refmans/CoordinateCleaner/html/cc_outl.html

### libraries ####

library(sf)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(ggpirate)

theme_set(theme_few())
options(scipen=10000)

### get data ####

source("01_getModels.R")
source("05_core_functions.R")

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(100,1,x))
PA_matrix <- do.call(rbind,PAs)

### area ####

#area
areaAnnual <- lapply(allspecies, function(x){
  applyRangeArea(x,modelSummaries_Limits, summary="annual")}) %>%
  reduce(rbind) 
saveRDS(areaAnnual,file="outputs/areaAnnual.rds")

# #check useful transformations
# areaAnnual <- areaAnnual %>%
#                 filter(Year %in% c(1990,2016)) %>%
#                 select(Species,Year,medianArea) %>%
#                 pivot_wider(names_from="Year", values_from="medianArea") %>%
#                 janitor::clean_names()
# 
# temp <- areaAnnual %>%
#           mutate(log_ratio = log(x2016/(x1990+1)),
#                   log_ratio2 = log(x2016/(x1990+as.numeric(meanArea))),
#                   perc_change = (x2016-x1990)/(x1990+as.numeric(meanArea)))


areaChanges <- lapply(allspecies, function(x){
  applyRangeArea(x,modelSummaries_Limits)})%>%
  reduce(rbind)
saveRDS(areaChanges,file="outputs/areaChanges.rds")

#plot
areaChanges %>%
  mutate(species = fct_reorder(species, desc(medianChange))) %>%
  ggplot()+
  geom_pointrange(aes(x = species, y = medianChange, ymin = lowerChange, max = upperChange))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  ylab("Relative change in AOCC")+
  theme(axis.text.y = element_text(size=rel(0.6)))

ggsave("plots/areaChange.png")

#cal some summary statistics
nrow(subset(areaChanges, medianChange>0))
nrow(subset(areaChanges, medianChange<0))

summary(subset(areaChanges, medianChange>0)$medianChange)
summary(subset(areaChanges, medianChange<0)$medianChange)

#boxplots
areaChanges <- readRDS("outputs/areaChanges.rds")
areaChanges$Direction <- ifelse(areaChanges$medianChange>0,
                                "Winners","Losers")
                            
(fig3b <- ggplot(areaChanges, aes(x = Direction, y = abs(medianChange)))+
  geom_pirate(aes(colour = Direction),bars=FALSE)+
  xlab("Direction of change") + ylab("|Change in AOO| (log-scale)")+
  scale_y_log10()+
  scale_colour_manual(values=c("deeppink","dodgerblue4")))

saveRDS(fig2b, "plots/boxplot_areaChange.rds")

### concave hull ####

#hull - annual
hullAnnual <- lapply(allspecies, function(x){
  applyConcaveMan(x,modelSummaries_Limits, summary="annual")})%>%
  reduce(rbind)
saveRDS(hullAnnual,file="outputs/concavehullAnnual.rds")


#check useful transformations
# hullAnnual <- hullAnnual %>%
#                 filter(Year %in% c(1990,2016)) %>%
#                 select(Species,Year,medianArea) %>%
#                 pivot_wider(names_from="Year", values_from="medianArea") %>%
#                 janitor::clean_names()
# 
# temp <- hullAnnual %>%
#           mutate(log_ratio = log(x2016/(x1990+1)),
#                   log_ratio2 = log(x2016/(x1990+as.numeric(meanArea))),
#                   perc_change = (x2016-x1990)/(x1990+as.numeric(meanArea)))


#hull-changes
hullChanges <- lapply(allspecies, function(x){
  applyConcaveMan(x,modelSummaries_Limits)})%>%
  reduce(rbind) %>%
  mutate(species = fct_reorder(species, desc(medianChange))) 
saveRDS(hullChanges,file="outputs/concavehullChanges.rds")

#plot
hullChanges %>%
  #filter(medianChange < 10) %>%
  ggplot()+
  geom_pointrange(aes(x = species, y = medianChange, ymin = lowerChange, 
                      max = upperChange))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed") +
  ylab("Relative change in EOCC")+
  theme(axis.text.y = element_text(size=rel(0.7)))

ggsave("plots/hullChange.png")

#cal some statistics
nrow(subset(hullChanges, medianChange>0))
nrow(subset(hullChanges, medianChange<0))

summary(subset(hullChanges, medianChange>0)$medianChange)
summary(subset(hullChanges, medianChange<0)$medianChange)

#boxplots
hullChanges$Direction <- areaChanges$Direction[match(hullChanges$species,
                                                   areaChanges$species)]

(fig3b.eocc <- ggplot(hullChanges, aes(x = Direction, y = abs(medianChange)))+
    geom_pirate(aes(colour = Direction), bars=FALSE)+
    xlab("Direction of change") + ylab("|Change in EOO| (log-scale)")+
    scale_y_log10()+
    scale_colour_manual(values=c("deeppink","dodgerblue4")))

saveRDS(fig2b.eocc, "plots/boxplot_hullChange.rds")

### alpha hull ####

hullChanges <- lapply(allspecies, function(x){
  applyAlphaHull(x,modelSummaries_Limits)
  print(x)
  }) %>%
  reduce(rbind) %>%
  mutate(species = fct_reorder(species, desc(medianChange))) 

### MCP hull ####

# MCP = minimum convex hull

#change
hullChanges <- lapply(allspecies, function(x){
  applyMCPHull(x,modelSummaries_Limits)
})%>%
  reduce(rbind) %>%
  mutate(species = fct_reorder(species, desc(medianChange))) 
saveRDS(hullChanges,file="outputs/mcphullChanges.rds")

#annual
hullAnnual <- lapply(allspecies, function(x){
  applyMCPHull(x,modelSummaries_Limits, summary="annual")
})%>%
  reduce(rbind)
saveRDS(hullAnnual,file="outputs/mcpahullAnnual.rds")

### compare hulls ####

hullChanges1 <- readRDS("outputs/concavehullChanges.rds")
hullChanges2 <- readRDS("outputs/mcphullChanges.rds")

bothChanges <- inner_join(hullChanges1,hullChanges2,
                         by=c("species"),
                         suffix = c("_1","_2"))
qplot(medianChange_1, medianChange_2, data=bothChanges)
#strongly correlated -2 is smaller!

### saturation #####

saturationChanges <- lapply(allspecies,function(x){
  applySaturation(x, modelSummaries_Limits, summary = "change")
}) %>% reduce(rbind)

saveRDS(saturationChanges, file="outputs/saturationChanges.rds")

#boxplots
saturationChanges$Direction <- areaChanges$Direction[match(saturationChanges$species,
                                                     areaChanges$species)]

(g1 <- ggplot(saturationChanges, aes(x = Direction, y = abs(medianChange)))+
    geom_pirate(aes(colour = Direction),bars=FALSE)+
    xlab("Direction of change") + ylab("|Change in saturation|")+
    scale_colour_manual(values=c("deeppink","dodgerblue4")))

saveRDS(g1, "plots/boxplot_saturationChange.rds")

### relationships ####

#AOO vs EOO
hullChanges <- readRDS("outputs/concavehullChanges.rds")

allChanges <- inner_join(hullChanges,areaChanges,
                         by=c("species"),
                         suffix = c("_extent","_area"))

(fig3c <- ggplot(data = allChanges,
       aes(x = medianChange_area,y = medianChange_extent)) + 
  geom_point() + 
  geom_smooth(method="gam",se=FALSE, colour="grey") +
  geom_errorbar(aes(ymin = lowerChange_extent,ymax = upperChange_extent)) + 
  geom_errorbarh(aes(xmin = lowerChange_area, xmax = upperChange_area))+
  scale_colour_viridis_c()+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in AOO") + ylab("Change in EOO"))

ggsave("plots/eoccChange_vs_aoccChange.png")

#AOO vs Saturation
saturationChanges <- readRDS("outputs/saturationChanges.rds")

allChanges <- saturationChanges  %>%
  inner_join(.,areaChanges,
             by=c("species"),
             suffix = c("_sat","_area"))

(fig3c.eocc <- ggplot(data = subset(allChanges, medianChange_sat>(-6)),
              aes(x = medianChange_area,y = medianChange_sat)) + 
    geom_point() + 
    geom_smooth(method="gam", se=FALSE, colour="black") +
    geom_errorbar(aes(ymin = lowerChange_sat,ymax = upperChange_sat)) + 
    geom_errorbarh(aes(xmin = lowerChange_area, xmax = upperChange_area))+
    geom_hline(linetype="dashed",yintercept=0)+
    geom_vline(linetype="dashed",xintercept=0)+
    xlab("Change in AOO") + ylab("Change in saturation (AOO/EOO)")+
    theme_few())

### regression ####

#Deming regression or orthogonal regression
library(deming)

#get error (approximate for the most)
allChanges$sd_extent <- (allChanges$upperChange_extent - allChanges$medianChange_extent)/2
allChanges$sd_area <- (allChanges$upperChange_area - allChanges$medianChange_area)/2

#fit models to each group
model_winners <- deming(medianChange_extent ~ exp(medianChange_area), 
       xstd = sd_area, ystd = sd_extent, 
       data = subset(allChanges, medianChange_area > 0 & medianChange_extent < 10))

model_losers <- deming(medianChange_extent ~ medianChange_area, 
                        xstd = sd_area, ystd = sd_extent, 
                        data = subset(allChanges, medianChange_area<0 & medianChange_extent<10))

#https://stackoverflow.com/questions/68684079/restrict-range-of-geom-ablineslope-somenumber

g1 +
  # geom_abline(aes(xmax=5),slope = 1,color='red') +
  geom_function(fun=Vectorize(function(x) {
    if(x > 0)
      return(NA)
    else
      return(x * model_losers$coefficients[2] + model_losers$coefficients[1])
  }), color = viridis::viridis(2)[1]) +
  
  # geom_abline(aes(xmin=5),slope = 1,intercept = 3,color='blue') +
  geom_function(fun=Vectorize(function(x) {
    if(x > 0)
      return(exp(x) * model_winners$coefficients[2] + model_winners$coefficients[1])
    else
      return(NA)
  }), color = viridis::viridis(2)[2])


ggsave("plots/eoccChange_vs_aoccChange.png", width = 7, height = 4)


#test difference between slopes
#model_winners
#model_losers
#?
#deming allow only one predictor

### species ts ####

#pick species
myspecies <- "Crocothemis erythraea"

hullSpecies <- applyConcaveMan(myspecies,modelSummaries, summary="annual")
g1 <- ggplot(hullSpecies)+
  geom_pointrange(aes(x = Year, y = medianArea, ymin = lowerArea, ymax = upperArea))+
  ylab("EOCC")

areaSpecies <- applyRangeArea(myspecies,modelSummaries, summary="annual")
mult <- as.numeric(meanArea)
g2 <- ggplot(areaSpecies)+
  geom_pointrange(aes(x = Year, y = medianArea *mult, ymin = lowerArea*mult, ymax = upperArea*mult))+
  ylab("AOCC")

p1_CE <- cowplot::plot_grid(g1,g2,nrow=1)

myspecies <- "Sympetrum danae"

#same code as above
p1_SD <- cowplot::plot_grid(g1,g2,nrow=1)


#combine both together
cowplot::plot_grid(p1_CE,p1_SD,
                   nrow=2,
                   align="hv",
                   axis = "tblr",
        labels=c("a) C. erythraea","b) S. danae"),
        hjust=c(-0.35,-0.5))

ggsave("plots/example_ts.png")


### nationwide-mean ####

#occupancy area
(aocc <- modelSummaries_Limits %>%
  group_by(Species,Year) %>%
  summarise(area= sum(mean)*as.numeric(meanArea))%>%
  mutate(area = area/1000000000) %>%
  group_by(Species) %>%
  mutate(change=(area[Year==2016]-area[Year==1990]))%>%
  ungroup()%>%
  ggplot(aes(x=Year,y=area,group=Species,colour=change))+
  scale_colour_gradient2(low="deeppink",mid="grey99",high="dodgerblue4")+
  geom_point(size=4,alpha=0.9)+
  geom_line(size=2,alpha=0.5)+
  scale_x_continuous(breaks=c(1990,2016),labels=c(1990,2016))+
  theme(legend.position = "top") +
  ylab("Area of occupancy"))

#extent of occupancy

hullAnnual <- readRDS("outputs/concavehullAnnual.rds")

(eocc <- hullAnnual %>%
  group_by(Species) %>%
  mutate(medianArea = medianArea/1000000000,
         change=(medianArea[Year==2016]-medianArea[Year==1990])) %>%
  ungroup()%>%
  ggplot(aes(x=Year,y=medianArea,group=Species,colour=change))+
  scale_colour_gradient2(low="deeppink",mid="grey97",high="dodgerblue4")+
  geom_point(size=4,alpha=0.9)+
  geom_line(size=2,alpha=0.5)+
  scale_x_continuous(breaks=c(1990,2016),labels=c(1990,2016))+
  theme(legend.position = "top") +
  ylab("Extent of occupancy"))


### fig. 3 ####

#plot together
top <- plot_grid(aocc,eocc,ncol=2)
middle <- plot_grid(fig3b,fig3b.eocc,ncol=2)
bottom <- plot_grid(fig3c,fig3c.eocc,ncol=2)

plot_grid(top,middle,bottom,nrow=3)

ggsave("plots/Fig.3.png",width=7.5, height=9)

### end ####