library(tidyverse)
library(ggthemes)

### sMon folder ####

sMonFolder <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects"

### mtbq info ####

load(paste(sMonFolder,"splines/mtbqsDF.RData",sep="/"))
names(mtbqsDF)[2] <- "MTB"
mtbsDF <- subset(mtbqsDF,!duplicated(MTB))

#### choose model directory ####

#spatial-temporal and k = 7 + expanded grid
modelDirectory <- paste(sMonFolder,"model-outputs/Odonata_stan_spline/v11",sep="/")

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

