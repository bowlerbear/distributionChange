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

source("splines/stan/spline_functions.R")

### get data ####

#run the modelSummaries_stan

#all those less than 5% should be removed.

modelSummaries$mean <- ifelse(modelSummaries$mean<0.3,0,modelSummaries$mean)

#also remove rare species
speciesSummary <- modelSummaries %>%
                    dplyr::group_by(Species) %>%
                    dplyr::summarise(medOccu = median(mean)) %>%
                    arrange(medOccu) %>%
                    filter(medOccu>0.05)

modelSummaries <- subset(modelSummaries, Species %in% speciesSummary$Species)

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

allYears <- 1990:2016

### application ####

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(100,1,x))
PA_matrix <- do.call(rbind,PAs)

#area
areaChanges <- lapply(allspecies, function(x){
  applyRangeArea(x,modelSummaries_Limits)})
areaChanges  <- do.call(rbind,areaChanges)
saveRDS(areaChanges,file="splines/outputs/areaChanges.rds")

#order species by range changes
areaChanges <- arrange(areaChanges, myMedian)
areaChanges$Species <- factor(areaChanges$species, levels=areaChanges$species)

#plot
ggplot(areaChanges)+
  geom_pointrange(aes(x = Species, y = myMedian, ymin = lower, max = upper))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_few()

#lat extent
rangeExtents <- lapply(allspecies, function(x){
  applyRangeExtent(x,modelSummaries_Limits)
  })
rangeExtents  <- do.call(rbind,rangeExtent)
saveRDS(rangeExtents,file="splines/outputs/rangeExtents.rds")

#order species by range changes
rangeExtents <- arrange(rangeExtents, median_Max)
rangeExtents$Species <- factor(rangeExtents$species, levels=rangeExtents$species)

#plot
ggplot(rangeExtents)+
  geom_pointrange(aes(x = Species, y = median_Max, ymin = lower_Max, max = upper_Max),
                  colour="darkblue")+
  geom_pointrange(aes(x = Species, y = median_Min, ymin = lower_Min, max = upper_Min),
                  colour="darkgreen")+
  ylab("Latitudinal extent change (m)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_few()+
  theme(axis.text.y = element_text(size=rel(0.6)))
  
#maximum y changes the most

#hull
hullChanges <- lapply(allspecies, function(x){
  applyConcaveMan(x,modelSummaries_Limits)})
hullChanges  <- do.call(rbind,hullChanges)
save(hullChanges,file="splines/outputs/hullChanges.rds")

#order species by range changes
hullChanges <- arrange(hullChanges, myMedian)
hullChanges$Species <- factor(hullChanges$species, levels=hullChanges$species)

#plot
ggplot(hullChanges)+
  geom_pointrange(aes(x = Species, y = myMedian, ymin = lower, max = upper))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_few()


### relationships ####

allChanges <- inner_join(hullChanges,areaChanges,
                         by=c("Species","species"),
                         suffix = c("_extent","_area"))

ggplot(data = allChanges,
       aes(x = myMedian_area,y = myMedian_extent)) + 
  geom_point() + 
  geom_text(data=subset(allChanges,abs(myMedian_extent)>2),
            aes(x = myMedian_area,y = myMedian_extent,label=Species))+
  geom_errorbar(aes(ymin = lower_extent,ymax = upper_extent)) + 
  geom_errorbarh(aes(xmin = lower_area, xmax = upper_area))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("area") + ylab("extent")+
  theme_few()


### species ts ####

#get realizations
PAs <- lapply(modelSummaries$mean, function(x) rbinom(10,1,x))
PA_matrix <- do.call(rbind,PAs)

#pick species
myspecies <- "Crocothemis erythraea"

hullSpecies <- applyConcaveMan(myspecies,modelSummaries, summary="annual")
g1 <- ggplot(hullSpecies)+
  geom_pointrange(aes(x = Year, y = medianArea, ymin = lowerArea, ymax = upperArea))+
  ylab("EOCC")+  
  theme_few()

areaSpecies <- applyRangeArea(myspecies,modelSummaries, summary="annual")
mult <- as.numeric(meanArea)
g2 <- ggplot(areaSpecies)+
  geom_pointrange(aes(x = Year, y = myMedian *mult, ymin = lower*mult, ymax = upper*mult))+
  ylab("AOCC") + 
  theme_few()

#packing
hullSpecies$packing <- hullSpecies$medianArea/areaSpecies$myMedian
g3 <- ggplot(hullSpecies)+
  geom_point(aes(x = Year, y = packing))+
  ylab("Packing") + 
  theme_few()

p1_CE <- cowplot::plot_grid(g1,g2,g3,nrow=1)

myspecies <- "Sympetrum danae"

#same as above
p1_SD <- cowplot::plot_grid(g1,g2,g3,nrow=1)

cowplot::plot_grid(p1_CE,p1_SD,
                   nrow=2,
                   align="hv",
                   axis = "tblr",
        labels=c("a) C. erythraea","b) S. danae"))


### ecoregion analysis ####

#clean ecoregion
modelSummaries_Limits$Naturraum <- gsub("Naturräume Deutschland/","",
                                         modelSummaries_Limits$Coarsenaturraum)

modelSummaries_Limits$Naturraum[modelSummaries_Limits$Naturraum=="Alpen"]<- "Alpen / vorland"
modelSummaries_Limits$Naturraum[modelSummaries_Limits$Naturraum=="Alpenvorland"]<- "Alpen / vorland"

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(500,1,x))
PA_matrix <- do.call(rbind,PAs)

#apply function
ecoChange <- lapply(allspecies,function(x)applyEcoregion(x,
                                                         modelSummaries_Limits, 
                                                         summary="annual"))
ecoChange <- do.call(rbind,ecoChange)
ecoChange$change <- boot::logit(ecoChange$median2016)-boot::logit(ecoChange$median1990)

g1 <- ggplot(ecoChange) +
  geom_boxplot(aes(x=Naturraum,y=median1990))+
  coord_flip()+
  theme_few()
  
g2 <- ggplot(ecoChange) +
    geom_boxplot(aes(x=Naturraum,y=median2016))+
    coord_flip()+
    theme_few()

cowplot::plot_grid(g1,g2)

ggplot(ecoChange) +
  geom_boxplot(aes(x=Naturraum,y=change))+
  coord_flip()+
  geom_hline(yintercept=0,linetype="dashed")+
  theme_few()

### state  analysis ####