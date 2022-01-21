library(tidyverse)
library(ggthemes)
library(sf)

### sMon folder ####

sMonFolder <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects"

### mtbq info ####

load(paste(sMonFolder,"splines/mtbqsDF.RData",sep="/"))
names(mtbqsDF)[2] <- "MTB"
mtbsDF <- subset(mtbqsDF,!duplicated(MTB))

#### choose model directory ####

modelDirectory <- paste(sMonFolder,"model-outputs/Odonata_stan_spline/v17",sep="/")

### get list of models ####

stanFiles <- list.files(modelDirectory) %>% str_subset("m_fit")

### psi predictions ####

#function to apply to each file
readStanModel <- function(file, get="psi"){
  
    readRDS(paste(modelDirectory,file,sep="/")) %>%
    as.data.frame() %>%
    add_column(Param = row.names(.)) %>%
    dplyr::filter(str_detect(Param, get)) %>% 
    dplyr::filter(str_detect(Param, "beta_psi", negate = TRUE)) %>%
    add_column(File = rep(file,nrow(.))) %>%
    dplyr::mutate(Species = gsub("m_fit_summary_spacetime_","",File)) %>%
    dplyr::mutate(Species = gsub(".rds","",Species))
}

modelSummaries <- stanFiles %>%
  map_df(readStanModel) %>%
  mutate(siteIndex = parse_number(Param))

#complete site and year information 
siteInfo_NAs <- readRDS(paste(sMonFolder,"splines/siteInfo_NAs.rds",sep="/")) %>%
    dplyr::select(!c(Species,SpeciesOrig)) %>%
    dplyr::filter(type!="extension")

#merge
modelSummaries <- inner_join(modelSummaries,siteInfo_NAs, by="siteIndex")
modelSummaries$Species <- sapply(modelSummaries$Species, function(x)strsplit(x,"_")[[1]][1])

#### check time-series ####

#analyse nationwide time series for each species
nuMTBs <- length(unique(siteInfo_NAs$MTB))
myspecies <- sort(unique(modelSummaries$Species))

annualTS <- modelSummaries %>%
            dplyr::group_by(Species,Year) %>%
            dplyr::summarise(total = sum(mean)/nuMTBs)

ggplot(annualTS)+
  geom_line(aes(x=Year, y=total))+
  facet_wrap(~Species)+
  theme_bw()

### processing ####

#remove rare species
speciesSummary <- modelSummaries %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(medOccu = median(mean)) %>%
  arrange(medOccu) %>%
  filter(medOccu>0)

#we lose 6 species
modelSummaries <- subset(modelSummaries, Species %in% speciesSummary$Species)

sort(unique(modelSummaries$Species))

#remove S. flaveolum - migratory and fluctuates alot!
modelSummaries <- subset(modelSummaries, Species!="Sympetrum flaveolum")

#weird maps at the moment - check later
modelSummaries <- subset(modelSummaries, !Species %in% c("Anax ephippiger",
                                                         "Boyeria irene"))

#all those occupancies than 5% should be absent
modelSummaries$mean <- ifelse(modelSummaries$mean<0.05,0,modelSummaries$mean)

### GIS data ####

mtbs <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",layer="MTB_25832")
#Projected CRS: ETRS89 / UTM zone 32N
equalProj <- st_crs(mtbs)
area <- st_area(mtbs)
meanArea <- mean(area)
totalArea <- as.numeric(sum(area))

#get germany outline
germanOutline <- raster::getData(name='GADM', country='DE',level=0) 

germanOutline <-   st_as_sf(germanOutline) %>%
  st_transform(.,equalProj)

#add on lon and lat
modelSummaries$lon <- mtbsDF$lon[match(modelSummaries$MTB, mtbsDF$MTB)]
modelSummaries$lat <- mtbsDF$lat[match(modelSummaries$MTB, mtbsDF$MTB)]
modelSummaries <- subset(modelSummaries, !is.na(lon) & !is.na(lat))

### subset data #####

#subset to first and last year
modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))

allspecies <- sort(unique(modelSummaries$Species))

allYears <- sort(unique(modelSummaries$Year))

utmProj <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


