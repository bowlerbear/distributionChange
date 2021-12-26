library(AHMbook)
library(raster)
library(RandomFields)

#using random fiels

#original function
?simOccSpatial
dat <- simOccSpatial(nsurveys = 3, mean.psi = 0.6, beta = c(0, 0),
                    mean.p = 1, alpha = c(0, 0), 
                    sample.size = 500, 
                    variance.RF = 1, theta.RF = 10,
                    seeds = c(10, 100), 
                    show.plots = TRUE, verbose = TRUE)

#modified function
my_simOccSpatial <- function(mean.psi = 0.6, nsites = 2500,latEffect = 0,
                              variance.RF = 1, theta.RF = 10, 
                              show.plots = TRUE){
  
  s <- simExpCorrRF(variance = variance.RF, theta = theta.RF, 
                    size = sqrt(nsites),
                    show.plots = show.plots)
  
  #parameters to predict occupancy
  beta0 <- qlogis(mean.psi)
  lpsi <- beta0 + c(s$field) + (s$grid[,2]-25) * latEffect
  psi <- plogis(lpsi)
  z <- rbinom(n = nsites, 1, psi)

  if (show.plots) {
    oldpar <- par(mfrow = c(1, 2), mar = c(5, 8, 5, 2), cex.lab = 1.5, 
                  "cex.main")
    
    oldAsk <- devAskNewPage(ask = dev.interactive(orNone = TRUE))
    on.exit({
      par(oldpar)
      devAskNewPage(oldAsk)
    })
    
    tryPlot <- try({
      
      r <- raster::rasterFromXYZ(data.frame(x = s$grid[,1], y = s$grid[, 2], 
                                            z = c(s$field)))
      raster::plot(r, col = topo.colors(20), axes = FALSE, 
                   box = FALSE, main = "Spatial effect (neg. exp. corr.)")
      
      r <- raster::rasterFromXYZ(data.frame(x = s$grid[,1], y = s$grid[,2], 
                                            z = z))
      raster::plot(r, col = c("white", "black"), axes = FALSE, 
                   box = FALSE, main = "Presence/absence (z)")
      

    }, silent = TRUE)
    
    if (inherits(tryPlot, "try-error")) 
      tryPlotError(tryPlot)
    
  }
  
}

#test it
my_simOccSpatial(mean.psi = 0.2, nsites = 2500,latEffect = -0.1,
                 variance.RF = 1, theta.RF = 20)

#see https://cran.r-project.org/web/packages/RandomFields/RandomFields.pdf


s <- simExpCorrRF(variance = variance.RF, theta = theta.RF, 
                  size = sqrt(nsites),
                  show.plots = show.plots)

simExpCorrRF <- function (variance = 1, theta = 1, size = 50, seed = NA, show.plots = TRUE){
  
  step <- 1
  x <- seq(1, size, step)
  y <- seq(1, size, step)
  grid <- cbind(x = rep(x, each = size), y = y)
  
  RFoptions(seed = seed)
  
  RFmodel = RMexp(var = variance, scale = theta) 
  
  field <- RFsimulate(RFmodel, x = x, y = y, grid = TRUE)@data$variable1
  
  RFoptions(seed = NA)
  
  if (show.plots) {
    oldpar <- par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), "cex.main")
    on.exit(par(oldpar))
    tryPlot <- try({
      dis <- seq(0.01, 20, by = 0.01)
      corr <- exp(-dis/theta)
      plot(dis, corr, type = "l", xlab = "Distance", ylab = "Correlation", 
           ylim = c(0, 1), col = "blue", lwd = 2)
      text(0.8 * max(dis), 0.8, labels = paste("theta:", 
                                               theta))
      par(mar = c(3, 2, 5, 1))
      raster::plot(rasterFromXYZ(cbind(grid, field)), col = topo.colors(20), 
                   main = paste("Gaussian random field with \n negative exponential correlation (theta =", 
                                theta, ")"), cex.main = 1, legend = FALSE, 
                   box = FALSE)
      box()
    }, silent = TRUE)
    if (inherits(tryPlot, "try-error")) 
      tryPlotError(tryPlot)
  }
  return(list(variance = variance, theta = theta, size = size, 
              seed = seed, field = field, grid = grid))
}


#using gstat

#http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/
library(gstat)
xy <- expand.grid(1:50, 1:50)
names(xy) <- c('x','y')

#if z is linearly dependent on x and y use the formula z~x+y
g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, 
                 beta=1, 
                 model=vgm(psill=5, range=50, model='Exp'), nmax=50)

#large the range, the more coarse
#the larger the nmax, the more coarse

yy <- predict(g.dummy, newdata=xy, nsim=4)
gridded(yy) = ~x+y
spplot(yy)

#a model with a trend in one dimension
g.dummy <- gstat(formula=z~1+y, locations=~x+y, dummy=T, beta=c(0,0.1), 
                 model=vgm(psill=1, range=50, model='Exp'), nmax=50)
yy <- predict(g.dummy, newdata=xy, nsim=4)
summary(yy)
gridded(yy) = ~x+y
spplot(yy)




