# A play at making a simulation...
# Will come back and comment better later after I know what the heck I'm doing.

# packages
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(concaveman)
library(nngeo)

# make a spatial grid for simulations
# corners
box <- data.frame(lat=c(53.3886, 53.45407, 48.2305, 48.3183),
                  lng=c(9.60172, 16.14958, 8.98649, 17.95133))
a <- c(53.3886, 9.60172)
b <- c(53.45407, 16.14958)
c <- c(48.2305, 8.98649)
d <- c(48.3183, 17.95133)

# make bounding box
bb <- box %>%
  st_as_sf(coords=c("lng", "lat"), crs=4326)

plot(bb)


# make grids
# can mess with grid size a bit more down the road
grids <- bb %>%
  st_make_grid(cellsize=0.1, square=FALSE) %>%
  st_as_sf() %>%
  mutate(index=1:nrow(.)) %>%
  mutate(centroid_lat=st_coordinates(st_centroid(.))[,2]) %>%
  mutate(centroid_lng=st_coordinates(st_centroid(.))[,1]) %>%
  mutate(area_m2=as.numeric(st_area(.)))

# quick test
ggplot()+
  geom_sf(data=bb)+
  geom_sf(data=grids, aes(color=centroid_lat))


# now we have a 'landscape' set up that can easily be changed
# we can start thinking about the simulation a little bit
# we want to vary the following I think:
# gain or loss
# percent gain or loss
# marginal core or random

# start with a 'high EOO and high AOO' scenario
scenario_1_function_random <- function(magnitude, prop_change, prop_start=.9){
  
  # top level 'boot function'
  # to get multiple runs of the different starting parameters
  boot_function <- function(draw){
    
    message(paste0("Running draw number ", draw))
    
    # first simulate a starting landscape
    # number of grids to sample
    grid_count <- round(nrow(grids)*prop_start)
    
    # occupied cells by index
    occupied_cells <- grids %>%
      sample_n(grid_count) %>%
      .$index
    
    # create a 'starting distribution'
    distribution_0 <- grids %>%
      mutate(presence=ifelse(index %in% occupied_cells, 1, 0))
    
    # visualize
    ggplot()+
      geom_sf(data=distribution_0, aes(fill=as.factor(presence)))
    
    # quick check
    sum(distribution_0$presence)==length(occupied_cells)
    
    # now simulate a step change
    # distribution_0=time0
    # and need to simulate
    # a change through time
    # put a switch in here if a run is a 'gain'
    # or a 'loss'
    cells_change <- if(magnitude=="decrease"){
      
      distribution_0 %>%
        dplyr::filter(presence==1) %>%
        sample_n(round(distribution_0 %>%
                         dplyr::filter(presence==1) %>%
                         nrow(.)*prop_change))
      
    } else {
      
      distribution_0 %>%
        dplyr::filter(presence==0) %>%
        sample_n(round(distribution_0 %>%
                         dplyr::filter(presence==0) %>%
                         nrow(.)*prop_change))
      
    }
    
    # now create a second distribution dataframe/sf object
    # for the second timepoint
    distribution_1 <- if(magnitude=="decrease"){
      
      distribution_0 %>%
        mutate(presence=ifelse(index %in% cells_change$index, 0, presence))
      
    } else {
      
      distribution_0 %>%
        mutate(presence=ifelse(index %in% cells_change$index, 1, presence))
      
    }
    
    # quick check
    sum(distribution_0$presence)
    sum(distribution_1$presence)
    
    # now calculate stuff on the 'distribution' dataframe
    # I'm just doing my own simple stuff for now! Diana will have to adapt to what
    # she did in her other scripts. I had a quick look but wasn't totally sure
    # and I'm just doing to do extent and area for now
    extent_0 <- distribution_0 %>%
      st_set_geometry(NULL) %>%
      dplyr::filter(presence==1) %>%
      st_as_sf(., coords=c("centroid_lng", "centroid_lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric()
    
    extent_1 <- distribution_1 %>%
      st_set_geometry(NULL) %>%
      dplyr::filter(presence==1) %>%
      st_as_sf(., coords=c("centroid_lng", "centroid_lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric()
    
    area_0 <- distribution_0 %>%
      st_set_geometry(NULL) %>%
      dplyr::filter(presence==1) %>%
      .$area_m2 %>%
      sum()
    
    area_1 <- distribution_1 %>%
      st_set_geometry(NULL) %>%
      dplyr::filter(presence==1) %>%
      .$area_m2 %>%
      sum()
    
    # pop into a dataframe
    summary_df <- data.frame(extent=c(extent_0, extent_1),
                             area=c(area_0, area_1),
                             time=c("time_0", "time_1")) %>%
      mutate(magnitude=magnitude) %>%
      mutate(prop_change=prop_change) %>%
      mutate(run_number=draw)
    
    return(summary_df)
    
  }
  
  results <- bind_rows(lapply(c(1:10), boot_function))

}

# now run this function twice
# once for 'increase'
# once for 'decrease'
# and each time from .1 to .9 prop change
increase_results <- bind_rows(lapply(seq(0.1, 0.9, by=.1), function(x){scenario_1_function_random("increase", x)})) 

decrease_results <- bind_rows(lapply(seq(0.1, 0.9, by=.1), function(x){scenario_1_function_random("decrease", x)})) 

scenario_1_results <- increase_results %>%
  bind_rows(decrease_results)

# do a little calculations
# and summarize over runs for now
temp <- scenario_1_results %>%
  pivot_wider(names_from=time, values_from=c(area, extent)) %>%
  mutate(delta_area=area_time_1-area_time_0) %>%
  mutate(delta_extent=extent_time_1-extent_time_0) %>%
  mutate(delta_area_percent=(delta_area/area_time_0)*100) %>%
  mutate(delta_extent_percent=(delta_extent/extent_time_0)*100) %>%
  group_by(magnitude, prop_change) %>%
  summarize(delta_area_percent=mean(delta_area_percent),
            delta_extent_percent=mean(delta_extent_percent))


ggplot(temp, aes(x=delta_area_percent, y=delta_extent_percent, color=prop_change))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="black", linetype="dashed")+
  geom_vline(xintercept=0, color="black", linetype="dashed")+
  scale_color_viridis_c(option="D")+
  ggtitle("Scenario 1: High Extent and High Area starting")



#############################################################
#############################################################
######### Now do the opposite case
######## Start with low EOO and low AOO
#############################################################
scenario_2_function_random <- function(magnitude, prop_change, prop_start=.1){
  
  # top level 'boot function'
  # to get multiple runs of the different starting parameters
  boot_function <- function(draw){
    
    message(paste0("Running draw number ", draw))
    
    # first simulate a starting landscape
    # number of grids to sample
    grid_count <- round(nrow(grids)*prop_start)
    
    # now we want to randomly sample 1 grid cell
    # and we'll then calculate the nearest neighbors from that
    # single grid cell and then sample ONLY within the nearest neighbors
    # but need to ensure the number of nearest neighbors > the number of grids to be samples
    
    # one random cell
    starting_cell <- grids %>%
      sample_n(1) 
    
    ggplot()+
      geom_sf(data=grids)+
      geom_sf(data=starting_cell, fill="red")
    
    # now get the nearest neighbors
    # first calculate k
    # we'll allow k to be 10% greater than the
    # number of samples that need to be collected
    # this percentage is what can be changed to control the level of 'smallness' we want in a starting distribution
    k <- round(grid_count+grid_count*.1)
    
    nn <- st_nn(starting_cell, grids, k=k)
    
    # occupied cells by index
    occupied_cells <- grids %>%
      dplyr::filter(index %in% nn[[1]]) %>%
      sample_n(grid_count) %>%
      .$index
    
    # create a 'starting distribution'
    distribution_0 <- grids %>%
      mutate(presence=ifelse(index %in% occupied_cells, 1, 0))
    
    # visualize
    ggplot()+
      geom_sf(data=distribution_0, aes(fill=as.factor(presence)))
    
    # quick check
    sum(distribution_0$presence)==length(occupied_cells)
    
    # now simulate a step change
    # distribution_0=time0
    # and need to simulate
    # a change through time
    # put a switch in here if a run is a 'gain'
    # or a 'loss'
    cells_change <- if(magnitude=="decrease"){
      
      distribution_0 %>%
        dplyr::filter(presence==1) %>%
        sample_n(round(distribution_0 %>%
                         dplyr::filter(presence==1) %>%
                         nrow(.)*prop_change))
      
    } else {
      
      distribution_0 %>%
        dplyr::filter(presence==0) %>%
        sample_n(round(distribution_0 %>%
                         dplyr::filter(presence==0) %>%
                         nrow(.)*prop_change))
      
    }
    
    # now create a second distribution dataframe/sf object
    # for the second timepoint
    distribution_1 <- if(magnitude=="decrease"){
      
      distribution_0 %>%
        mutate(presence=ifelse(index %in% cells_change$index, 0, presence))
      
    } else {
      
      distribution_0 %>%
        mutate(presence=ifelse(index %in% cells_change$index, 1, presence))
      
    }
    
    # quick check
    sum(distribution_0$presence)
    sum(distribution_1$presence)
    
    # now calculate stuff on the 'distribution' dataframe
    # I'm just doing my own simple stuff for now! Diana will have to adapt to what
    # she did in her other scripts. I had a quick look but wasn't totally sure
    # and I'm just doing to do extent and area for now
    extent_0 <- distribution_0 %>%
      st_set_geometry(NULL) %>%
      dplyr::filter(presence==1) %>%
      st_as_sf(., coords=c("centroid_lng", "centroid_lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric()
    
    extent_1 <- distribution_1 %>%
      st_set_geometry(NULL) %>%
      dplyr::filter(presence==1) %>%
      st_as_sf(., coords=c("centroid_lng", "centroid_lat"), crs=4326) %>%
      concaveman(., concavity=1) %>%
      st_area() %>%
      as.numeric()
    
    area_0 <- distribution_0 %>%
      st_set_geometry(NULL) %>%
      dplyr::filter(presence==1) %>%
      .$area_m2 %>%
      sum()
    
    area_1 <- distribution_1 %>%
      st_set_geometry(NULL) %>%
      dplyr::filter(presence==1) %>%
      .$area_m2 %>%
      sum()
    
    # pop into a dataframe
    summary_df <- data.frame(extent=c(extent_0, extent_1),
                             area=c(area_0, area_1),
                             time=c("time_0", "time_1")) %>%
      mutate(magnitude=magnitude) %>%
      mutate(prop_change=prop_change) %>%
      mutate(run_number=draw)
    
    return(summary_df)
    
  }
  
  results <- bind_rows(lapply(c(1:10), boot_function))
  
}

# now run this function twice
# once for 'increase'
# once for 'decrease'
# and each time from .1 to .9 prop change
increase_results <- bind_rows(lapply(seq(0.1, 0.9, by=.1), function(x){scenario_2_function_random("increase", x)})) 

decrease_results <- bind_rows(lapply(seq(0.1, 0.9, by=.1), function(x){scenario_2_function_random("decrease", x)})) 

scenario_2_results <- increase_results %>%
  bind_rows(decrease_results)

# do a little calculations
# and summarize over runs for now
temp2 <- scenario_2_results %>%
  pivot_wider(names_from=time, values_from=c(area, extent)) %>%
  mutate(delta_area=area_time_1-area_time_0) %>%
  mutate(delta_extent=extent_time_1-extent_time_0) %>%
  mutate(delta_area_percent=(delta_area/area_time_0)*100) %>%
  mutate(delta_extent_percent=(delta_extent/extent_time_0)*100) %>%
  group_by(magnitude, prop_change) %>%
  summarize(delta_area_percent=mean(delta_area_percent),
            delta_extent_percent=mean(delta_extent_percent))


ggplot(temp2, aes(x=delta_area_percent, y=delta_extent_percent, color=prop_change))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_hline(yintercept=0, color="black", linetype="dashed")+
  geom_vline(xintercept=0, color="black", linetype="dashed")+
  scale_color_viridis_c(option="D")+
  ggtitle("Scenario 2: Low Extent and Low Area starting")




