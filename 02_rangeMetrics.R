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

theme_set(theme_few())

source("05_core_functions.R")

### get data ####

#run script 01
modelSummaries$Species <- sapply(modelSummaries$Species, function(x)strsplit(x,"_")[[1]][1])

#all those less than 5% should be removed.

modelSummaries$mean <- ifelse(modelSummaries$mean<0.05,0,modelSummaries$mean)

#also remove rare species
speciesSummary <- modelSummaries %>%
                    dplyr::group_by(Species) %>%
                    dplyr::summarise(medOccu = median(mean)) %>%
                    arrange(medOccu) %>%
                    filter(medOccu>0.05)

modelSummaries <- subset(modelSummaries, Species %in% speciesSummary$Species)
sort(unique(modelSummaries$Species))

### GIS data ####

mtbs <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",layer="MTB_25832")
#Projected CRS: ETRS89 / UTM zone 32N
equalProj <- st_crs(mtbs)
area <- st_area(mtbs)
meanArea <- mean(area)

#get germany outline
germanOutline <- raster::getData(name='GADM', country='DE',level=0)
germanOutline <- st_as_sf(germanOutline)
germanOutline <- st_transform(germanOutline,equalProj)

### subset data #####

#subset to first and last year
modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))

allspecies <- sort(unique(modelSummaries$Species))

allYears <- sort(unique(modelSummaries$Year))

### area ####

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(100,1,x))
PA_matrix <- do.call(rbind,PAs)

#area
areaChanges <- lapply(allspecies, function(x){
  applyRangeArea(x,modelSummaries_Limits)})
areaChanges  <- do.call(rbind,areaChanges)
saveRDS(areaChanges,file="outputs/areaChanges.rds")

#plot
areaChanges %>%
  mutate(species = fct_reorder(species, desc(medianChange))) %>%
  ggplot()+
  geom_pointrange(aes(x = species, y = medianChange, ymin = lowerChange, max = upperChange))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  ylab("Relative area change")+
  theme(axis.text.y = element_text(size=rel(0.7)))

ggsave("plots/areaChange.png")


### lat extent ####

#lat extent
rangeExtents <- lapply(allspecies, function(x){
  applyRangeExtent(x,modelSummaries_Limits)
  })
rangeExtents  <- do.call(rbind,rangeExtents)
saveRDS(rangeExtents,file="outputs/rangeExtents.rds")

#plot
rangeExtents %>%
  pivot_longer(contains("_"),names_to="type", values_to="change") %>%
  separate(type, c("quantile","extent")) %>%
  filter(quantile=="median") %>%
  mutate(extent = ifelse(extent=="Max","upper","lower")) %>%
  ggplot()+
  geom_col(aes(x = species, y = change, fill=extent))+
  scale_fill_brewer("Position",type="qual") +
  ylab("Latitudinal extent change (m)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme(axis.text.y = element_text(size=rel(0.6)))
  
ggsave("plots/latitudeChange.png")

#maximum y changes the most

### hull ####

#hull
hullChanges <- lapply(allspecies, function(x){
  applyConcaveMan(x,modelSummaries_Limits)})
hullChanges  <- do.call(rbind,hullChanges)
save(hullChanges,file="outputs/hullChanges.rds")

#plot
hullChanges %>%
  mutate(species = fct_reorder(species, desc(medianChange))) %>%
  ggplot()+
  geom_pointrange(aes(x = species, y = medianChange, ymin = lowerChange, 
                      max = upperChange))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed") +
  ylab("Relative hull change")+
  theme(axis.text.y = element_text(size=rel(0.7)))

ggsave("plots/hullChange.png")

### relationships ####

allChanges <- inner_join(hullChanges,areaChanges,
                         by=c("species"),
                         suffix = c("_extent","_area"))

ggplot(data = subset(allChanges, medianChange_extent<3),
       aes(x = medianChange_area,y = medianChange_extent)) + 
  geom_point() + 
  #geom_text(data=subset(allChanges,abs(medianChange_extent)>2),
  #          aes(x = medianChange_area,y = medianChange_extent,label=Species))+
  geom_errorbar(aes(ymin = lowerChange_extent,ymax = upperChange_extent)) + 
  geom_errorbarh(aes(xmin = lowerChange_area, xmax = upperChange_area))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in AOCC") + ylab("Change in EOCC")+
  theme_few()

ggsave("plots/eoccChange_vs_aoccChange.png")

### species ts ####

#get realizations
PAs <- lapply(modelSummaries$mean, function(x) rbinom(100,1,x))
PA_matrix <- do.call(rbind,PAs)

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

#same as above
p1_SD <- cowplot::plot_grid(g1,g2,nrow=1)

cowplot::plot_grid(p1_CE,p1_SD,
                   nrow=2,
                   align="hv",
                   axis = "tblr",
        labels=c("a) C. erythraea","b) S. danae"),
        hjust=c(-0.35,-0.5))

ggsave("plots/example_ts.png")

### ecoregion analysis ####

#clean ecoregion
modelSummaries_Limits$Naturraum <- gsub("Naturräume Deutschland/","",
                                         modelSummaries_Limits$Coarsenaturraum)

modelSummaries_Limits$Naturraum[modelSummaries_Limits$Naturraum=="Alpen"]<- "Alpen / vorland"
modelSummaries_Limits$Naturraum[modelSummaries_Limits$Naturraum=="Alpenvorland"]<- "Alpen / vorland"
table(modelSummaries_Limits$Naturraum)

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(100,1,x))
PA_matrix <- do.call(rbind,PAs)

#apply function
ecoChange <- lapply(allspecies,function(x){
  applyEcoregion(x, modelSummaries_Limits,summary="annual")})
ecoChange <- do.call(rbind,ecoChange)
ecoChange$change <- boot::logit(ecoChange$median2016+0.01)-boot::logit(ecoChange$median1990+0.001)
save(ecoChange,file="outputs/ecoChange.rds")


ecoChange %>%
  ungroup() %>%
  select(Species,change,Naturraum,median2016,median1990) %>%
  pivot_longer(starts_with("median"),names_to="Year",values_to="value") %>%
  mutate(Year = parse_number(Year)) %>%
  ggplot(aes(x=Year,y=value,group=Species,colour=change))+
    scale_colour_gradient2(low="green",high="purple")+
    geom_point(alpha=0.2)+
    geom_line(size=2,alpha=0.5)+
    scale_x_continuous(breaks=c(1990,2016),labels=c(1990,2016))+
    facet_wrap(~Naturraum)+
    ylab("Occupancy proportion")

ggsave("plots/ecoChange.png",height=5.5,width=7)


### nationwide-mean ####

modelSummaries_Limits %>%
  group_by(Species,Year) %>%
  summarise(prop = sum(mean)/length(mean))%>%
  group_by(Species) %>%
  mutate(change=boot::logit(prop[Year==2016]+0.01)-boot::logit(prop[Year==1990]+0.001))%>%
  ungroup()%>%
  ggplot(aes(x=Year,y=prop,group=Species,colour=change))+
  scale_colour_gradient2(low="green",high="purple")+
  geom_point(alpha=0.2)+
  geom_line(size=2,alpha=0.5)+
  scale_x_continuous(breaks=c(1990,2016),labels=c(1990,2016))+
  ylab("Occupancy proportion")

ggsave("plots/nationChange.png")

### state  analysis ####

### end ####