### range area (AOO) ####

getRangeArea <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  dat <- subset(speciesData, PA == 1)[,c("Species","PA","x_MTB","y_MTB","Year","MTB")]
  
  dat %>%
    dplyr::group_by(Species,Year) %>%
    dplyr::summarise(nuGrids = length(unique(MTB))) %>%
    dplyr::ungroup() %>%
    tidyr::complete(Year = allYears) %>%
    dplyr::mutate(Species = species) %>%
    dplyr::mutate(nuGrids = ifelse(is.na(nuGrids),0,nuGrids)) %>%
    dplyr::mutate(sumGrids = nuGrids * as.numeric(meanArea))
  
}


applyRangeArea <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getRangeArea(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  
  #summarise change
  if(summary == "change"){
  temp %>%
    select(Year,Species,sumGrids,simNu) %>%
    tidyr::pivot_wider(everything(),names_from = Year,values_from=sumGrids) %>%
    janitor::clean_names() %>%  
    tidyr::complete() %>%
    dplyr::mutate(x1990 = ifelse(is.na(x1990),as.numeric(meanArea),x1990)) %>%
    dplyr::mutate(x2016 = ifelse(is.na(x2016),as.numeric(meanArea),x2016)) %>%
    dplyr::mutate(x1990 = ifelse(x1990==0,as.numeric(meanArea),x1990)) %>%
    dplyr::mutate(x2016 = ifelse(x2016==0,as.numeric(meanArea),x2016)) %>%
    dplyr::mutate(change = log(x2016/x1990)) %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(medianChange = median(change), 
                     lowerChange = quantile(change, 0.025),
                     upperChange = quantile(change, 0.975))
    
  }else if (summary == "annual") {
    
  #summarise annual
   temp %>%
     dplyr::group_by(Species,Year) %>%
     dplyr::summarise(medianArea = median(sumGrids), 
                      lowerArea = quantile(sumGrids, 0.025),
                      upperArea = quantile(sumGrids, 0.975))
    
  }
  
}

### latitudinal extents ####

getRangeExtents <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  dat <- subset(speciesData, PA == 1)[,c("Species","PA","x_MTB","y_MTB","Year","MTB")]
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  dat %>%
    dplyr::group_by(Species,Year) %>%
    dplyr::summarise(max_Y = max(y_MTB), 
                     min_Y = min(y_MTB)) %>%
    dplyr::ungroup()
  
}

applyRangeExtent <- function(species, modelSummaries_Limits, summary="change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getRangeExtents(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  if((sum(temp$Year==1990)) == (sum(temp$Year==2016))){
    
  #summarise change
  if(summary =="change"){
  temp %>%
    tidyr::pivot_wider(everything(),names_from = Year,values_from=c(max_Y,min_Y)) %>%
    janitor::clean_names() %>%  
    tidyr::complete() %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(change_maxY = (max_y_2016 - max_y_1990), 
                  change_minY = (min_y_2016 - min_y_1990))  %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(median_Max = median(change_maxY), 
                     lower_Max = quantile(change_maxY, 0.025),
                     upper_Max = quantile(change_maxY, 0.975),
                     median_Min = median(change_minY), 
                     lower_Min = quantile(change_minY, 0.025),
                     upper_Min = quantile(change_minY, 0.975))
    
  }else if(summary =="annual"){
    
  #summarise annual
   temp %>%
     dplyr::group_by(Species,Year) %>%
     dplyr::summarise(median_Max = median(max_Y), 
                      lower_Max = quantile(max_Y, 0.025),
                      upper_Max = quantile(max_Y, 0.975),
                      median_Min = median(min_Y), 
                      lower_Min = quantile(min_Y, 0.025),
                      upper_Min = quantile(min_Y, 0.975))
  }
  }
  
}

### extent (EOCC) ####

getConvexHull <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeHull <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
      ch <- chull(ydat)
      coords <- ydat[c(ch, ch[1]), ] 
      plot(ydat, pch=19, main = species)
      lines(coords, col="red")
      sp_poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), 
                                                       ID=1)),
                                     proj4string=sp::CRS("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
      
      sp_poly_cut <- raster::intersect(sp_poly,germanOutline)#make sense or not??
      rangeSize <- raster::area(sp_poly_cut)
      
      return(rangeSize)
    }else {
      return(0)
    }})
  
  data.frame(Species = species, Year = allYears, rangeHull = rangeHull) 
}


getAlphaHull <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeAlpha <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
      
      dist <- (max(ydat$y_MTB)-min(ydat$y_MTB))/3
      ah <- alphahull::ahull(ydat, alpha = dist)
      #plot(ah, main = species)
      alphahull::areaahull(ah)
      
    }else {
      return(0)
    }})
  
  data.frame(Species = species, Year = allYears, rangeAlpha = rangeAlpha) 
  
}


getMCP <- function(species, modelSummaries_Limits){
  
  require(adehabitatHR)
  require(sp)
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeHull <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>4){
      
      mycoords <- ydat[,c("x_MTB","y_MTB")]
      coordinates(mycoords) <- c("x_MTB","y_MTB")
      proj4string(mycoords) <- CRS("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
      ch <- mcp(mycoords, percent=99, unout="m2")
      rangeSize <- ch@data$area
      return(rangeSize)
    }else {
      return(0)
    }})
  
  data.frame(Species = species, Year = allYears, rangeHull = rangeHull) 
}

getConcaveMan <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeMan <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
      
      ydat <- sf::st_as_sf(ydat, coords =c("x_MTB", "y_MTB"),crs = 25832)
      sp_poly <- concaveman::concaveman(ydat)
      germanOutline_sf <- st_transform(germanOutline, st_crs(ydat))
      sp_poly_cut <- st_intersection(sp_poly, germanOutline_sf)
      #plot(sp_poly_cut)
      return(as.numeric(st_area(sp_poly)))#m2 units
    }else {
      return(0)
    }})
  
  data.frame(Species = species, Year = allYears, rangeMan = rangeMan)
  
}

applyConcaveMan <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getConcaveMan(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
  temp %>%
    tidyr::pivot_wider(everything(),names_from = Year, values_from = rangeMan) %>%
    janitor::clean_names() %>%  
    dplyr::mutate(change = log((x2016+as.numeric(meanArea))/(x1990+as.numeric(meanArea))))  %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(medianChange = median(change), 
                     lowerChange = quantile(change, 0.025),
                     upperChange = quantile(change, 0.975))
    
  }else if(summary == "annual"){
    
    #summarise annual
     temp %>%
       dplyr::group_by(Species,Year) %>%
       dplyr::summarise(medianArea = median(rangeMan), 
                        lowerArea = quantile(rangeMan, 0.025),
                        upperArea = quantile(rangeMan, 0.975))
    
  }
  
}


applyAlphaHull <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getAlphaHull(species, modelSummaries_Limits)
    #print(i)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
    temp %>%
      tidyr::pivot_wider(everything(),names_from = Year, values_from = rangeAlpha) %>%
      janitor::clean_names() %>%  
      dplyr::mutate(change = log((x2016+10)/(x1990+10)))  %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change,na.rm=T), 
                       lowerChange = quantile(change, 0.025,na.rm=T),
                       upperChange = quantile(change, 0.975,na.rm=T))
    
  }else if(summary == "annual"){
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year) %>%
      dplyr::summarise(medianArea = median(rangeAlpha,na.rm=T), 
                       lowerArea = quantile(rangeAlpha, 0.025,na.rm=T),
                       upperArea = quantile(rangeAlpha, 0.975,na.rm=T))
    
  }
  
}


applyMCPHull <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getMCP(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
    temp %>%
      tidyr::pivot_wider(everything(),names_from = Year, values_from = rangeHull) %>%
      janitor::clean_names() %>%  
      dplyr::mutate(change = log10((x2016+100000)/(x1990+100000)))  %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change), 
                       lowerChange = quantile(change, 0.025),
                       upperChange = quantile(change, 0.975))
    
  }else if(summary == "annual"){
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year) %>%
      dplyr::summarise(medianArea = median(rangeHull), 
                       lowerArea = quantile(rangeHull, 0.025),
                       upperArea = quantile(rangeHull, 0.975))
    
  }
  
}

### spatial var ####

getChangeSD <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  speciesData %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(nuStable = sum(PA[Year==1990]==PA[Year==2016]),
                     nu1990 = sum(PA[Year==1990]),
                     nu2016 = sum(PA[Year==2016]),
                     nuStable0 = sum(PA[Year==1990]==0 & PA[Year==2016]==0),
                     nuStable1 = sum(PA[Year==1990]==1 & PA[Year==2016]==1),
                     nuIncrease = sum(PA[Year==1990]==0 & PA[Year==2016]==1),
                     nuDecrease = sum(PA[Year==1990]==1 & PA[Year==2016]==0),
                     total = length(PA[Year==1990]),
                     totalChange = abs(sum(PA[Year==2016])-sum(PA[Year==1990]))) %>%
    #Site turnover is defined as the probability that a 
    #randomly chosen occupied site is newly occupied
    dplyr::mutate(Turnover = abs(log(nuIncrease)-log(nuDecrease)),
                  probNew = nuIncrease/nu2016,
                  Change = (nuIncrease + nuDecrease)/(nuStable0 + nuStable1),
                  Diff = ifelse(nuIncrease>nuDecrease,
                                totalChange/nuIncrease,totalChange/nuDecrease),
                  probStay = nuDecrease/nu1990 + nuIncrease/(1-nu1990))
  
}


applyChangeSD <-function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getChangeSD(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  }) %>%
    reduce(rbind) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(medianChange = quantile(probStay,0.5),
                     lowerChange = quantile(probStay,0.025),
                     upperChange = quantile(probStay,0.975))
  
}


### clumpiness ####

#calculate the local occupancy change and get the clumpiness of the change

getFragChange <- function(species, modelSummaries_Limits){
  
  #data for 1990
  speciesSummary <- modelSummaries_Limits %>%
    filter(Species==species) %>%
    filter(Year == 1990)
  
  speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$PA),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #make into a raster
  r1 <- raster(speciesRaster)
  
  #data for 2016
  speciesSummary <- modelSummaries_Limits %>%
    filter(Species==species) %>%
    filter(Year == 2016)
  
  speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$PA),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #make into a raster
  r2 <- raster(speciesRaster)
  
  
  #make stack  
  speciesStack <- stack(list(r1,r2))
  speciesStack <- projectRaster(speciesStack, crs=utmProj, method="ngb")
  speciesChange <- calc(speciesStack, function(x) x[2]-x[1])
  
  #calc change metrics
  changes <- calculate_lsm(speciesChange, 
                what = c("lsm_c_pland",
                         "lsm_c_clumpy"),
                full_name = TRUE) %>%
    mutate(class = ifelse(class==0, "no change",
                          ifelse(class==1, "increase", "decrease"))) %>%
    add_column(Species = species)
  
  #also absolute change
  speciesAbsChange <- calc(speciesChange, function(x) ifelse(x==0,0,1))
  
  calculate_lsm(speciesAbsChange, 
                           what = c("lsm_c_pland",
                                    "lsm_c_clumpy"),
                           full_name = TRUE) %>%
    filter(class == 1) %>%
    mutate(class = "change") %>%
    add_column(Species = species) %>%
    bind_rows(.,changes)
  
}


applyFragChange <- function(species, modelSummaries_Limits){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getFragChange(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise across sims
  temp %>%
      dplyr::filter(metric == "clumpy" & !is.na(value)) %>%
      dplyr::select(Species,class,value,simNu) %>%
      dplyr::group_by(Species,class) %>%
      dplyr::summarise(medianChange = median(value), 
                       lowerChange = quantile(value, 0.025),
                       upperChange = quantile(value, 0.975))
    
}


#get the clumpiness of the annual distributions, and then the change in the clumpiness
getFragStats <- function(species, modelSummaries_Limits){
  
  #data for 1990
  speciesSummary <- modelSummaries_Limits %>%
    filter(Species==species) %>%
    filter(Year == 1990)
  
  speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$PA),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #make into a raster
  r1 <- raster(speciesRaster)
  
  #data for 2016
  speciesSummary <- modelSummaries_Limits %>%
    filter(Species==species) %>%
    filter(Year == 2016)
  
  speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$PA),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #make into a raster
  r2 <- raster(speciesRaster)
  
  
  #make stack  
  speciesStack <- stack(list(r1,r2))
  speciesStack <- projectRaster(speciesStack, crs=utmProj, method="ngb")
  
  #calc metrics
  calculate_lsm(speciesStack, 
                what = c("lsm_l_ta",
                         "lsm_c_pland",
                         "lsm_c_clumpy"),
                full_name = TRUE) %>%
    add_column(Year = c(rep(1990,5),rep(2016,5))) %>%
    add_column(Species = species)
}

applyFragStats <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getFragStats(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
    temp %>%
      dplyr::filter(class == 1 & metric == "clumpy" & !is.na(value)) %>%
      dplyr::select(Species,Year,value,simNu) %>%
      tidyr::pivot_wider(everything(),names_from = Year, values_from = value) %>%
      janitor::clean_names() %>%  
      dplyr::mutate(change = x2016/x1990)  %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change), 
                       lowerChange = quantile(change, 0.025),
                       upperChange = quantile(change, 0.975))
    
  }else if(summary == "annual"){
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year,class, metric,name,type,function_name) %>%
      dplyr::summarise(medianMetric = median(value,na.rm=T), 
                       lowerMetric = quantile(value, 0.025, na.rm=T),
                       upperMetric = quantile(value, 0.975, na.rm=T))
    
  }
  
}


### core regions ####

getCoreRegions <- function(myspecies){
  
  speciesSummary <- modelSummaries %>%
    filter(Species == myspecies) %>%
    filter(Year %in% 1990:1995) %>%
    dplyr::group_by(Species, MTB, lon, lat) %>%
    dplyr::summarise(meanpsi = median(mean)) %>%
    dplyr::mutate(normalpsi = (meanpsi-min(.$meanpsi))/(max(.$meanpsi)-min(.$meanpsi)))
  
  speciesPixels <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$normalpsi),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #make into a raster
  speciesRaster <- raster(speciesPixels)
  
  #for each raster cell, get mean class of surroundings
  ## in a 3-by-3 window surrounding each focal cell 
  rmean <- focal(speciesRaster, 
                 w = matrix(1, ncol=3, nrow=3), 
                 fun=mean, na.rm=TRUE)
  #plot(rmean)
  
  #define core as all those surrounded by 50% presences
  speciesRasterCore <- rmean
  speciesRasterCore[speciesRasterCore > 0.5] <- 1 #indicator for core sites
  speciesRasterCore[speciesRasterCore < 0.5] <- 0.5 #indicator for marginal sites
  speciesRasterCore[is.na(speciesRaster)] <- NA
  speciesRasterCore[speciesRaster < 0.2] <- 0 #absent sites
  plot(speciesRasterCore)
  
  
  #plot example
  # library(tmap)
  #  t1 <- tm_shape(speciesRaster)+
  #    tm_raster(style="cont",legend.show=FALSE)
  #  t2 <- tm_shape(speciesRasterCore)+
  #    tm_raster(legend.show=FALSE)
  #  tmap_arrange(t1,t2,ncol=2)
  
  
  #convert back into a data frame
  coreDF <- as.data.frame(speciesRasterCore, xy=TRUE)
  coreDF$Obs <- as.data.frame(speciesRaster)[,1] #original observation
  coreDF$cellNu <- cellFromXY(speciesRaster, coreDF[,c("x","y")])
  coreDF <- subset(coreDF, !is.na(layer)) #remove sites beyond german border
  
  #clean marginal information
  coreDF$Core <- ifelse(coreDF$layer==1,"core",
                        ifelse(coreDF$layer==0.5, "marginal", "absent"))
  
  #add on other data
  coreDF$Species <- myspecies
  
  
  #map to mtbqs
  mtbPixels <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                      data = data.frame(MTB = as.numeric(speciesSummary$MTB)),
                                      tolerance = 0.6,
                                      proj4string = crs("+proj=longlat +datum=WGS84"))
  
  mtbRaster <- raster(mtbPixels)
  mtbDF <- as.data.frame(mtbRaster, xy=TRUE)
  mtbDF$cellNu <- cellFromXY(mtbRaster, mtbDF[,c("x","y")])
  coreDF <- left_join(coreDF, mtbDF[,3:4], by="cellNu")
  
  return(coreDF)
  
}

getCoreCalc <- function(myspecies, summary="annual"){
  
  coreDF_species <- subset(coreDF, Species== myspecies)
  modelSummaries_Limits$Core <- coreDF_species$Core[match(modelSummaries_Limits$MTB, 
                                                          coreDF_species$MTB)]
  modelSummaries_Limits$Core[modelSummaries_Limits$Core=="absent"] <- "marginal"
  
  #for each realization do the following
  temp <- lapply(1:dim(PA_matrix)[2], function(i){
    modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Species==myspecies) %>%
      group_by(Core, Year) %>%
      summarize(occ = sum(PA), total= length(PA)) %>%
      mutate(prop = occ/total) %>%
      add_column(simNu = i)
  }) %>% bind_rows() %>% ungroup() 
  
  
  if(summary=="annual"){
    
    temp %>%
      dplyr::group_by(Core,Year) %>%
      dplyr::summarize(medianProp = quantile(prop,0.5), 
                       lowerProp = quantile(prop,0.025),
                       upperProp = quantile(prop,0.975)) %>%
      add_column(Species = myspecies)
    
  }else if(summary=="change"){
    
    temp %>% 
      dplyr::select(Core, Year, prop, simNu) %>%
      dplyr::group_by(Core,Year, simNu) %>%
      tidyr::pivot_wider(names_from="Year", values_from="prop") %>%
      janitor::clean_names(case = "title") %>%
      dplyr::mutate(change = boot::logit(X2016+0.01) - boot::logit(X1990+0.01)) %>%
      dplyr::mutate(change = ifelse(is.infinite(change),0, change)) %>%
      dplyr::mutate(change = ifelse(is.na(change),0, change)) %>%
      dplyr::group_by(Core) %>%
      dplyr::summarize(medianChange = quantile(change,0.5), 
                       lowerChange = quantile(change,0.025),
                       upperChange = quantile(change,0.975)) %>%
      add_column(Species = myspecies)
    
  }
  
}

getColExt <- function(myspecies){
  
  #add on core region info
  coreDF_species <- subset(coreDF, Species== myspecies)
  modelSummaries_Limits$Core <- coreDF_species$Core[match(modelSummaries_Limits$MTB, 
                                                          coreDF_species$MTB)]
  modelSummaries_Limits$Core[modelSummaries_Limits$Core=="absent"] <- "marginal"
  
  #colonizations
  colTemp <- lapply(1:dim(PA_matrix)[2], function(i){
    
    emptySites <- modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Species==myspecies) %>%
      filter(Year==1990 & PA==0) %>%
      pull(MTB)
    
    modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Year==2016) %>%
      filter(Species==myspecies) %>%
      filter(MTB %in% emptySites) %>%
      group_by(Core) %>%
      summarize(occ = sum(PA), total= length(PA)) %>%
      mutate(prop = occ/total) %>%
      add_column(simNu = i)
  }) %>% 
    bind_rows() %>% 
    ungroup() %>%
    add_column(Type = "Colonization")
  
  
  #extinctions
  extTemp <- lapply(1:dim(PA_matrix)[2], function(i){
    
    occuSites <- modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Species==myspecies) %>%
      filter(Year==1990 & PA==1) %>%
      pull(MTB)
    
    modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Year==2016) %>%
      filter(Species==myspecies) %>%
      filter(MTB %in% occuSites) %>%
      group_by(Core) %>%
      summarize(occ = sum(PA), total= length(PA)) %>%
      mutate(prop = 1-occ/total) %>%
      add_column(simNu = i)
  }) %>% 
    bind_rows() %>% 
    ungroup() %>%
    add_column(Type = "Extinction")
  
  bind_rows(colTemp,extTemp) %>%
    group_by(Core,Type) %>%
    summarise(medianProp = median(prop),
              sdProp = sd(prop),
              lowerProp = quantile(prop,0.025),
              upperProp = quantile(prop, 0.975)) %>%
    group_by(Type) %>%
    mutate(ratio = boot::logit(medianProp[Core=="core"])-boot::logit(medianProp[Core=="marginal"])) %>%
    ungroup() %>%
    add_column(species = myspecies)
  
}

### spatial variability ####

applySaturation <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    
    #get hull
    out <- getConcaveMan(species, modelSummaries_Limits)
    
    #add on number of occupied grids
    out2 <- getRangeArea(species, modelSummaries_Limits) %>%
              filter(Year %in% c(1990,2016))
    out$sumGrids <- out2$sumGrids
    
    #deal with zeros - assume one grid occupied
    out$rangeMan[is.na(out$rangeMan)] <- as.numeric(meanArea)
    out$sumGrids[is.na(out$sumGrids)] <- as.numeric(meanArea)
    out$rangeMan[out$rangeMan==0] <- as.numeric(meanArea)
    out$sumGrids[out$sumGrids==0] <- as.numeric(meanArea)
    
    out$saturation <- out$sumGrids/out$rangeMan
    out$saturation[is.infinite(out$saturation)] <- 0
    out$sim <- i
    return(out[,c("Species","Year","saturation","sim")])
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
    temp %>%
      tidyr::pivot_wider(everything(),names_from = Year, values_from = saturation) %>%
      janitor::clean_names() %>%  
      dplyr::mutate(change = (boot::logit(x2016) - boot::logit(x1990-0.01)))  %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change,na.rm=T), 
                       lowerChange = quantile(change, 0.025,na.rm=T),
                       upperChange = quantile(change, 0.975,na.rm=T))
    
  }else if(summary == "annual"){
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year) %>%
      dplyr::summarise(medianArea = median(saturation), 
                       lowerArea = quantile(saturation, 0.025),
                       upperArea = quantile(saturation, 0.975))
    
  }
  
}



### ecoregion analysis ####

applyEcoregion <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i]
    
    modelSummaries_Limits %>%
      filter(Species==species) %>%
      dplyr::group_by(Naturraum,Year) %>%
      dplyr::summarise(nu=sum(PA),total=length(PA)) %>%
      add_column(Species=species, simNu=i)
    
  })
  
  temp <- do.call(rbind,temp)
  
  
  #summarise
  if(summary=="change"){
    temp %>%
      dplyr::group_by(Naturraum,Species,simNu) %>%
      dplyr::summarize(change = (nu[Year==2016]-nu[Year==1990]+1)/(total[Year==1990])) %>%
      dplyr::group_by(Species,Naturraum) %>%
      dplyr::summarise(medianChange = quantile(change,0.5),
                       lowerChange = quantile(change,0.025),
                       upperChange = quantile(change,0.975))
    
    
  }else if(summary=="annual"){
    temp %>%
      dplyr::group_by(Species,Naturraum,simNu) %>%
      dplyr::summarize(prop1990 = nu[Year==1990]/total[Year==1990],
                       prop2016 = nu[Year==2016]/total[Year==2016]) %>%
      dplyr::group_by(Species,Naturraum) %>%
      dplyr::summarise(median1990 = quantile(prop1990,0.5),
                       lower1990 = quantile(prop1990,0.025),
                       upper1990 = quantile(prop1990,0.975),
                       median2016 = quantile(prop2016,0.5),
                       lower2016 = quantile(prop2016,0.025),
                       upper2016 = quantile(prop2016,0.975))
    
  }
  
}

