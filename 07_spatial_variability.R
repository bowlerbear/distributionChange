### get data ####

source("01_getModels.R")
source("05_core_functions.R")

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(50,1,x))
PA_matrix <- do.call(rbind,PAs)

### turnover ####

#Jaccard turnover index

#turnover - probability of a random selected patch occupied at time 1, 
#is still occupied at time 2
#of patches occupied at time 1, how many are occupied at time 2

#null is that total change is driven all by increase or decrease 

### get SDs ####

sdChanges <- lapply(allspecies,function(x){
  applyChangeSD(x, modelSummaries_Limits)
}) %>% 
  reduce(rbind)%>%
  rename(species="Species")
saveRDS(sdChanges, file="outputs/sdChanges.rds")

#boxplots

sdChanges$Direction <- areaChanges$Direction[match(sdChanges$species,
                                                           areaChanges$species)]

(g1 <- ggplot(sdChanges, aes(x = Direction, y = medianChange)) +
    geom_pirate(aes(colour = Direction, fill = Direction),bars=FALSE) +
    xlab("Direction of change") + ylab("Turnover") +
    scale_colour_brewer(type="qual", direction=-1) +
    scale_fill_brewer(type="qual", direction=-1))

### relationships ####

allChanges <- inner_join(sdChanges,areaChanges,
                         by=c("species"),
                         suffix = c("_sd","_area"))

(g1 <- ggplot(data = allChanges,
              aes(x = medianChange_area,y = medianChange_sd)) + 
    geom_point() + 
    #geom_smooth(method="gam") +
    geom_errorbar(aes(ymin = lowerChange_sd,ymax = upperChange_sd)) + 
    geom_errorbarh(aes(xmin = lowerChange_area, xmax = upperChange_area))+
    scale_colour_viridis_c()+
    geom_hline(linetype="dashed",yintercept=0)+
    geom_vline(linetype="dashed",xintercept=0)+
    xlab("Change in AOO") + ylab("Turnover"))


### boxplots ####

p1 <- readRDS("plots/boxplot_saturationChange.rds")
p2 <- readRDS("plots/boxplot_areaChange.rds")
p3 <- readRDS("plots/boxplot_hullChange.rds")
p4 <- readRDS("plots/boxplot_fragChange_change.rds")
plot_grid(p2,p3,p1,p4,nrow=2)
