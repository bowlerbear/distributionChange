#maps of spatial change

library(sf)
library(tmap)

#tmaptools::palette_explorer()

### mtb grid #####

mtbGrid <- siteInfo_NAs %>%
            filter(!duplicated(MTB))

#mtb shapefile
mtbs <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                layer="MTB_25832") %>%
        rename(MTB = "Value")


### annual richness ####

speciesrichness <- modelSummaries_Limits %>% 
                    group_by(Year,MTB) %>%
                    summarise(richness=sum(mean)) %>%
                    merge(mtbs,.) %>%
                    janitor::clean_names()


sr <- tm_shape(speciesrichness)+
  tm_fill(col="richness", palette="YlOrRd",
          legend.is.portrait = FALSE,
          n = 10)+
  tm_facets(by="year")+
  tm_layout(legend.outside = TRUE, 
            legend.outside.position = "top",
            legend.position = c(0.20, 0.05))
sr

tmap_save(sr, "plots/maps_speciesrichness.png",height=6,width=10)

### increase vs decrease species ####

speciesrichnessChange <- modelSummaries_Limits %>% 
  dplyr::select(mean, Species, MTB, Year) %>%
  dplyr::group_by(MTB, Species) %>%
  pivot_wider(names_from="Year", values_from="mean") %>%
  ungroup() %>%
  janitor::clean_names() %>%
  dplyr::mutate(increase = ifelse(x2016>x1990,1,0),
                decrease= ifelse(x1990>=x2016,1,0)) %>%
  rename(MTB = "mtb")

#where do species increase
speciesrichnessIncrease <- speciesrichnessChange %>% 
  group_by(MTB) %>%
  summarise(increase=sum(x2016[increase==1])-x1990[increase==1]) %>%
  merge(mtbs,.)

sr_I <- tm_shape(speciesrichnessIncrease)+
  tm_fill(col="increase", palette="YlOrRd",
          style="cont")

#where do species decrease
speciesrichnessDecrease <- speciesrichnessChange %>% 
  group_by(MTB) %>%
  summarise(decrease=sum(x2016[decrease==1])-x1990[decrease==1]) %>%
  merge(mtbs,.)

sr_D <- tm_shape(speciesrichnessDecrease)+
  tm_fill(col="decrease", palette="YlOrRd",
          style="cont")

sr_C <- tmap_arrange(sr_I, sr_D)
sr_C

tmap_save(sr_C, "plots/maps_speciesrichnessChange.png",height=6,width=10)

### plot example species ####

exampleData <- modelSummaries_Limits %>% 
                filter(Species %in% c("Anax imperator",
                                      "Crocothemis erythraea",
                                      "Lestes sponsa",
                                      "Sympetrum danae")) %>%
                rename(Occupancy = "mean") 

(g1 <- ggplot(exampleData) +
  geom_point(aes(x = x_MTB, y = y_MTB, colour = Occupancy)) +
  facet_grid(Year ~ Species) +
  theme_map() +
  scale_colour_viridis_c(option = "A", direction = 1))

ggsave("plots/maps_examplespecies.png",height=6,width=10)


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
  applyEcoregion(x, modelSummaries_Limits,summary="annual")})%>%
  reduce(rbind)
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

### lat extent ####

#lat extent
rangeExtents <- lapply(allspecies, function(x){
  applyRangeExtent(x,modelSummaries_Limits)
})%>%
  reduce(rbind)
saveRDS(rangeExtents,file="outputs/rangeExtents.rds")

#plot
rangeExtents %>%
  filter(species %in% hullChanges$species[hullChanges$medianChange<10]) %>%
  pivot_longer(contains("_"),names_to="type", values_to="change") %>%
  separate(type, c("quantile","extent")) %>%
  filter(quantile=="median") %>%
  mutate(extent = ifelse(extent=="Max","upper","lower")) %>%
  mutate(species = factor(species, levels(factor(hullChanges$species[hullChanges$medianChange<10])))) %>%
  ggplot()+
  geom_col(aes(x = species, y = change, fill=extent))+
  scale_fill_brewer("Position",type="qual") +
  ylab("Latitudinal extent change (m)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme(axis.text.y = element_text(size=rel(0.6)))

ggsave("plots/latitudeChange.png")


modelSummaries_Limits %>%
  group_by(Species,Year) %>%
  summarise(medianExtent = weighted.mean(y_MTB, mean)) %>%
  mutate(medianExtent = medianExtent/1000) %>%
  group_by(Species) %>%
  mutate(change=medianExtent[Year==2016]-medianExtent[Year==1990]) %>%
  ungroup()%>%
  ggplot(aes(x=Year,y=medianExtent,group=Species,colour=change))+
  scale_colour_gradient2(low="green",high="purple")+
  geom_point(alpha=0.2)+
  geom_line(size=2,alpha=0.5)+
  scale_x_continuous(breaks=c(1990,2016),labels=c(1990,2016))+
  theme(legend.position = "top") +
  ylab("Extent position")
g2

ggsave("plots/nationextentChange.png")


### case study species ####

#change in range extent - contract, expand
#only change in AOO

areaChanges <- readRDS("outputs/areaChanges.rds")
hullChanges <- readRDS("outputs/concavehullChanges.rds")

allChanges <- inner_join(hullChanges,areaChanges,
                         by=c("species"),
                         suffix = c("_extent","_area"))
ggplot(data = allChanges,
              aes(x = medianChange_area,y = medianChange_extent)) + 
    geom_text(aes(label=species))

#Change in AOO and EOO
#Increase - Aeshna affinis
#Decrease - Nehalennia speciosa

#Change in AOO but not EOO
#Increase - Anax parthenope
#Decrease - Sympetrum pedemontanum

modelSummaries_Limits %>% 
  filter(Species %in% c("Sympetrum pedemontanum")) %>%
  rename(Occupancy = "mean") %>%
  ggplot() +
  geom_point(aes(x = x_MTB, y = y_MTB, colour = Occupancy)) +
  facet_wrap(~ Year) +
  theme_map() +
  scale_colour_viridis_c(option = "A", direction = 1)+
  theme(legend.position="none")


#Change in saturation
#Increase - 
#Decrease - 

#Core change
#Increase in marginal - 
#Decrease in core - 

#clumpiness of change
#Clumpy  - 
#Non-clumpy - 

### end ####