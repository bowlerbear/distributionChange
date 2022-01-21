### get data ####

source("01_getModels.R")
source("05_core_functions.R")

#get realizations
PAs <- lapply(modelSummaries_Limits$mean, function(x) rbinom(50,1,x))
PA_matrix <- do.call(rbind,PAs)

### saturation #####

#need to sort out units

saturationDF <- lapply(allspecies[-c(1,64)],function(x){
  applySaturation(x, modelSummaries_Limits, summary = "change")
}) %>% reduce(rbind)

saveRDS(saturationDF, file="outputs/satuationChange.rds")

### plotting ####

allChanges <- saturationDF  %>%
  inner_join(.,areaChanges,
             by=c("species"),
             suffix = c("_sat","_area"))

g2 <- ggplot(data = allChanges,
             aes(x = medianChange_area,y = medianChange_sat)) + 
  geom_point() + 
  geom_hline(yintercept=0,linetype="dashed")+
  stat_smooth(method="gam")+
  geom_vline(linetype="dashed",xintercept=0)+
  xlab("Change in AOC") + ylab("Change in saturation")

g2