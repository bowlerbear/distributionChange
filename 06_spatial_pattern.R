#map of spatial change

library(sf)
library(tmap)

tmaptools::palette_explorer()

### mtb grid #####

mtbGrid <- siteInfo_NAs %>%
            filter(!duplicated(MTB))

#mtb shapefile
mtbs <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                layer="MTB_25832") %>%
        rename(MTB = "Value")


### change in richness ####

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

### community turnover ####


### plot example species ###

exampleData <- modelSummaries_Limits %>% 
                filter(Species %in% c("Anax imperator",
                                      "Crocothemis erythraea",
                                      "Lestes sponsa",
                                      "Sympetrum danae")) %>%
                rename(Occupancy = "mean") 

g1 <- ggplot(exampleData) +
  geom_point(aes(x = x_MTB, y = y_MTB, colour = Occupancy)) +
  facet_grid(Year ~ Species) +
  theme_map() +
  scale_colour_viridis_c(direction = -1)

ggsave(g1,"plots/maps_examplespecies.png",height=6,width=10)
