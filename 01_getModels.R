library(tidyverse)
library(ggthemes)

#sMon folder
sMonFolder <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects"

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
  
  temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
  
  #get rows for the parameter of interest
  temp$Param <- row.names(temp)
  temp <- subset(temp, grepl(get,temp$Param))
  temp <- subset(temp, !grepl("beta_psi",temp$Param))
                 
  #get species name
  temp$File <- file
  temp$Species <- gsub("m_fit_summary_spacetime_","",temp$File)
  temp$Species <- gsub(".rds","",temp$Species)
  
  return(temp)
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

