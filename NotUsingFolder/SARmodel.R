rm(list=ls())

## -- Required Packages --:
source("~/Desktop/Tentative Results-Proj2/R scripts/requiredPackages.R")
source("~/Desktop/Tentative Results-Proj2/R scripts/riskFunc.R")
source("~/Desktop/Tentative Results-Proj2/R scripts/funcPlotImagePlot.R")
source("~/Desktop/Tentative Results-Proj2/R scripts/Convert-2Dto3D.R")
source("~/Desktop/Tentative Results-Proj2/R scripts/marginalTransformation.R")
source("~/Desktop/Tentative Results-Proj2/R scripts/16x16gridBoxes.R")
source("~/Desktop/Tentative Results-Proj2/R scripts/funcGen_r_Pareto.R")
source("~/Desktop/Tentative Results-Proj2/R scripts/ParameterConfiguration.R") 
source("~/Desktop/Tentative Results-Proj2/R scripts/VarFuncDev(G_G_bar).R")
source("~/Desktop/Tentative Results-Proj2/R scripts/r_Pareto_SimFunc.R")
source("~/calcExtremogram.R")

load('~/Downloads/redseatemperature_extra.rdata')
load("~/Downloads/redseatemperature.rdata")                 

## -------------------------- Script -------------------------------------------

## Spatial domain of the SST:
print(dim(loc)) # dim: 16703 2

## Extract the summer obs (July, August, September):
summer_months <- month %in% c(7,8,9)
data  <- data[summer_months, ]
month <- month[summer_months]
year  <- year[summer_months]
time  <- time[summer_months]

### ---- Code: for removing the trend ----- ###
# Separate the date part and the time zone offset
library('lubridate')
date_part <- substr(time, 1, 10)
timezone_offset <- substr(time, 12, 14)

# Parse the date part
parsed_date <- ymd(date_part)

# Convert the timezone offset to hours and create a duration object
offset_hours <- as.numeric(timezone_offset)
duration_offset <- hours(offset_hours)

# Adjust the parsed date by adding the duration offset
parsed_datetime <- parsed_date + duration_offset # we worked witj parsed datetime for trend fit
### ---------------------------------------------- ###


## Extract observations from the southern part of interest:
idx  <- loc[ ,"lat"] < 20.2 & loc[,"lat"] > 15.75
data <- data[ , idx]
loc  <- loc[idx, ] # 6239    2 

## Transform locations:
loc[, "lon"] <- loc[, "lon"] * 1.04
loc[, "lat"] <- loc[, "lat"] * 1.11

## Rename the full data so that we can access it in each of the following sections:
# 2852 is the daily summer observation for 31 years 
# 6239 is the spatial locations
full_data <- data  
print(dim(full_data)) # 2852 6239


full_loc  <- loc
print(dim(full_loc)) # 6239 2


par(mfrow=c(1,1))
plot(full_loc[,1], 
     full_loc[,2],
     col= alpha('black', 0.2),
     pch = 21,
     xlab = 'lon',
     ylab = 'lat',
     main = 'Focusing on Southern Region')
cat("Number of locations in each field in the (regular) Red Sea data set:",
    length(idx))

## Function to fit a linear trend and compute residuals
## Function to fit a linear trend and compute residuals
fit_linear_trend <- function(dat_val, time_vals=parsed_date) {
  # Fit the linear model
  linear_model <- lm(dat_val ~ time_vals)
  
  # Calculate the residuals (errors)
  residuals <- dat_val - linear_model$fitted.values
  
  # Return the residuals
  return(residuals)
}

## Detrend the temporal component
tic()
Detrended_data <- suppressWarnings(apply(data, 2, fit_linear_trend))
toc()

## Function to fit marginal GPD and transform the marginals to standard Pareto:
# ( to be noted, only the standard Pareto margins will be applicable with GM risk function):
XStarMargins <- function(univariateData, thrQ = 0.95)
{
  size <- length(univariateData)
  empiricalMargins <- rank(univariateData)/(size+1)
  
  empiricalMargins[univariateData > quantile(univariateData, thrQ, na.rm=TRUE)] <- NA
  
  invisible(capture.output(fit <- gpd.fit(univariateData,
                                          threshold=quantile(univariateData,thrQ, na.rm=TRUE))$mle))
  scale <- fit[1]
  shape <- fit[2]
  
  # computing the probability
  empiricalMargins[is.na(empiricalMargins)] <- 1 - pgpd(univariateData[univariateData>quantile(univariateData, thrQ)],
                                                        quantile(univariateData, thrQ),
                                                        scale,
                                                        shape,
                                                        lower.tail=FALSE)*(1-thrQ) 
  
  
  # Calculate Standard Pareto  margins
  XPMargins <- 1/(1 - empiricalMargins)
  
  ## convert the standard Pareto to exponential margins
  #expMargins <- log(XPMargins)
  
  return(XPMargins) # returning standard Pareto marginals
}

## -- Step 0: Transform the marginals standard Pareto marginals --:
# fit standard Pareto to marginal of the spatial process
tic()
X_star <- suppressWarnings(apply(Detrended_data, 2, XStarMargins))
toc()

dim(X_star) # 2852 6239
# 11.426 sec elapsed

## -- Step 1: Defining the Red Sea Surface Temperature on grid in 3D array format --:
# Extracting the unique values of longitude and latitude values
lon <- full_loc[, 1]
lat <- full_loc[, 2]
unique_lon <- sort(unique(lon))
unique_lat <- sort(unique(lat))

grid_matrix <- array(NA, c(length(unique_lon), length(unique_lat), nrow(data)))
dim(grid_matrix)
for (i in 1:nrow(full_loc)) {
  if(i %% 125 == 0) {
    print(i)
  }
  lon_index <- which(unique_lon == full_loc[i, 1])
  lat_index <- which(unique_lat == full_loc[i, 2])
  grid_matrix[lon_index, lat_index, ] <- X_star[,i]  # or some other value or identifier
}  
dim(grid_matrix) # 113 89 2852

## -- Step 1: Divide the spatial domain to 16x16 grids --:
# Define the sequence for latitude; it's fixed for all grids
grid_size <- 16
obs16x16Boxes <- list()
gridPoints16x16Boxes <- list()
countSpatialPoints16x16Grids <- list() # count of spatial points falling in 16x16 grids

lat_indices <- seq(1, length(unique_lat), by = grid_size-1)
print(lat_indices)

# Loop over each latitude start point
for(lat_start in lat_indices)
{
  lat_end <- min(lat_start + grid_size - 1, length(unique_lat))
  lat_indices_subset <- lat_start:lat_end
  for(lon_start in seq(1, length(unique_lon), by = grid_size-1))
  {
    lon_end <- min(lon_start + grid_size - 1, length(unique_lon))
    lon_indices <- lon_start:lon_end
    
    # Define the lon/lat range for the current grid box
    lon_range <- unique_lon[lon_indices]
    lat_range <- unique_lat[lat_indices_subset]
    
    # Check if there are any points within the current grid box
    points_in_box <- full_loc[,1] >= min(lon_range) & full_loc[,1] <= max(lon_range) & 
      full_loc[,2] >= min(lat_range) & full_loc[,2] <= max(lat_range)
    
    # Count the points in the box
    points_count <- sum(points_in_box)
    
    if (any(points_in_box))
    { # Add the grid only if there are points inside it
      obs16x16Boxes[[length(obs16x16Boxes)+1]] <- grid_matrix[lon_indices, lat_indices_subset, ]
      gridPoints16x16Boxes[[length(gridPoints16x16Boxes)+1]] <- list('lon' = unique_lon[lon_indices],
                                                                     'lat' = unique_lat[lat_indices_subset])
      countSpatialPoints16x16Grids[[length(countSpatialPoints16x16Grids)+1]] <- points_count
      cat('Print point count in','box', i, ' : ', points_count,'\n')
    }
  }
}

# Plotting the all possible 16x16 grid boxes across the Red Sea spatial domain:
# png("~/Desktop/16x16GridBox.png",
#     units = "in",
#     width = 12,
#     height = 10,
#     res = 200)
# par(mfrow = c(1,1),
#     mar = c(4, 4, 2, 2),
#     oma = c(2, 2.5, 2, 2),
#     cex.main = 1.5)

plot(full_loc[,1],
     full_loc[,2], 
     col = alpha('black', 0.1),
     pch = 19, 
     cex = 2.5,
     yaxt='n',
     xaxt='n',
     cex.axis = 2, 
     cex.lab = 2.5, 
     cex.main = 3,
     xlab = '',
     ylab = '',
     main = '')
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext('Lat',
      side=2,
      las=2,
      cex=2.5,
      line=3.5)
mtext('Lon',
      side=1,
      cex=2.5,
      line=3.5)
mtext('Spatial Domain of the Red Sea',
      side=3,
      cex=3,
      line=1.5)

# Overlay the 16x16 boxes
for (box in gridPoints16x16Boxes){
  # Assuming box is a list with 'lon' and 'lat' as vectors
  lon_range <- range(box$lon)
  lat_range <- range(box$lat)
  
  # Draw the rectangle for the 16x16 box
  rect(lon_range[1],
       lat_range[1],
       lon_range[2],
       lat_range[2],
       border = 'magenta',
       lwd = 4)
}
# dev.off()
print(length(gridPoints16x16Boxes))

# Scatter plots for checking the grid count across the 16x16 grid boxes:
# png("~/Desktop/CountsFalling16x16Grids.png",
#     units = "in",
#     width = 10,
#     height = 10,
#     res = 200)
plot(x = 1:length(countSpatialPoints16x16Grids), 
     y = countSpatialPoints16x16Grids, 
     pch =  19,
     cex = 1.5,
     cex.axis = 1.5, 
     cex.lab = 1.5, 
     col = alpha('grey', 0.8),
     xlab = '16x16 grid by index',
     ylab = 'Grid box point count')
abline(h=256, lwd=3, col= 'magenta')
# dev.off()

## -- Step 2: Defining r-Pareto process on the 16x16 tiles --:
threshold_across_boxes <- rep(NA)
R_vals <- list()
rPareto_across_boxes <- list()
log_W_process_across_boxes <- list()
extreme_days <- list()
obs_at_boxes_2d <- list()
for(i in 1:length(gridPoints16x16Boxes))
{
  obs_at_boxes <- obs16x16Boxes[[i]]
  dim_ar <- dim(obs_at_boxes)
  cat('Print dim of the 3D array',i,':', dim_ar)
  # 16   16 2852 
  
  # Defining 3 dimensional matrix to 2D array:
  X_star_obs <- matrix(obs_at_boxes, 
                       nrow=dim_ar[1]*dim_ar[2],
                       ncol=dim_ar[3])
  cat('  Print dim of the 3D array:', dim(X_star_obs))
  
  risk_averages <- apply(X_star_obs,
                         2, 
                         geometricMean)
  print(length(risk_averages))
  # Note: each column corresponds to daily observation across the space
  
  # Note: we can vary the threshold based on box-size
  threshold <- quantile(risk_averages, 0.97)
  threshold_across_boxes[i] <-  threshold
  
  extreme_ds <- which(risk_averages > threshold) # finding the extreme days
  extreme_days[[length(extreme_days)+1]] <- extreme_ds
  
  print(length(extreme_ds))
  
  X_star_extremes <- (X_star_obs[, extreme_ds])
  obs_at_boxes_2d[[length(rPareto_across_boxes)+1]]  <- X_star_obs
  
  R <- risk_averages[extreme_ds]
  Y_r <- X_star_extremes/threshold
  # using the spectral decomposition:
  log_W <- log(Y_r) - log(R)
  
  n <- length(rPareto_across_boxes)
  rPareto_across_boxes[[n+1]] <- Y_r # standard Pareto margins
  R_vals[[n+1]] <- R
  log_W_process_across_boxes[[n+1]] <- log_W
}
## ---------- End: Data processing and defining 16x16 grid tiles ---------- ##

## -- EDA of the transformed r-Pareto process with exponential marginals --:
# Check for validation:
# 1) fit Matern covariance across the tiles
# 2) check the fitted versus empricial variogram
# 3) Transform the spatial process with White Noise and check for model fit

## First, we randomly select tiles: let do 3 for EDA
set.seed(222)
rand_tiles <- round(runif(n=4, min=1, max=35))
cat('Rand tiles indexes:', rand_tiles, '\n')
# Rand tiles indexes: 33 3 18 1 

## -- Variogram of the log_W process for randomly selected tiles --:
## Tile 1: 33
coor1 <- expand.grid(gridPoints16x16Boxes[[rand_tiles[1]]]$lon, 
                     gridPoints16x16Boxes[[rand_tiles[1]]]$lat)
dim(coor1) # 256x2

log_W_1 <- log_W_process_across_boxes[[rand_tiles[1]]]
dim(log_W_1) # 256x86

# Fit Gaussian process for all the realizations all together:
## Fit LK model to get an idea about the
obj1 <- LatticeKrig(coor1, log_W_1, findAwght=TRUE)
print(obj1)

# fit_GP_tile33 <- spatialProcess(coor1,
#                                 log_W_1,
#                                 na.rm=TRUE,
#                                 collapseFixedEffect = FALSE,
#                                 smoothness=1,
#                                 mKrig.args=list(m=1)) # we fit linear mean model, with intercept and slope terms
# 
# 
# cov_parameter_tile1 <- as.numeric(fit_GP_tile33$MLESummary[6:9])
# fixed_constant_parameter_tile1  <- as.matrix(fit_GP_tile33$d)

## Fit Gaussian process to log_W process for Tile 11:
cov_parameter_tile1_dM <- matrix(NA,
                                 nrow=ncol(log_W_1),
                                 ncol=4)
fixed_constant_parameter_tile1_dM <- rep(NA)
# fixed_constant_parameter_tile1_dM  <- matrix(NA, 
#                                              nrow=ncol(log_W_1),
#                                              ncol=1)

# Covariance matrix varies across the space and time:
for(i in 1:ncol(log_W_1))
{
  y_process <- log_W_1[,i]
  print(i)
  
  tryCatch({
    fit_GP <- spatialProcess(coor1,
                             y_process,
                             na.rm=TRUE,
                             smoothness=1, # Changed from 0.5 to 1
                             mKrig.args=list(m=1))
    cov_parameter_tile1_dM[i,] <- as.numeric(fit_GP$MLESummary[6:9])
    fixed_constant_parameter_tile1_dM[i] <- as.numeric(fit_GP$d)
  }, error = function(e) {
    cat("Error occurred for i=", i, ", j=", j, ": ", conditionMessage(e), "\n")}
  )
}

## Variogram function for fitted model:
vTrue <- function(d, sigma2, range_val, nu)
{
  # vario_vals <- sigma2*(1-exp(-d/range_val))
  vario_vals <- sigma2*(1-Matern(d=d, range=range_val, smoothness = nu))
  return(vario_vals)
}

## Randomly selected day:
set.seed(4)
rand_days <- round(runif(n=4, 1, 86))
print(rand_days)

## Log W process for randomly selected days for Tile 1:
imgPlt_temp11 <- imagePlot(coor1,
                           log_W_1[, rand_days[1]]) # day 1
imgPlt_temp12 <- imagePlot(coor1,
                           log_W_1[, rand_days[2]]) # day 2
imgPlt_temp13 <- imagePlot(coor1,
                           log_W_1[, rand_days[3]]) # day 3
imgPlt_temp14 <- imagePlot(coor1,
                           log_W_1[, rand_days[4]]) # day 4

## Day 1
rdist1 <- rdist(coor1) # dim: 224 224
sigma2_1 <- cov_parameter_tile1_dM[rand_days[1], 2]
tau1 <- cov_parameter_tile1_dM[rand_days[1], 1]
rho1 <- cov_parameter_tile1_dM[rand_days[1], 3]

sigma_mat11 <- (tau1^2) + sigma2_1* Matern(d=rdist1, 
                                           range=rho1,
                                           smoothness=1) # tau^2 is added for nugget effect
L1 <- t(chol(sigma_mat11)) # check for lower triangular
L_inv1 <- t(L1)%*%solve(sigma_mat11)

mean_wN11 <- as.numeric(fixed_constant_parameter_tile1_dM[rand_days[1]])
wN11 <- L_inv1%*% (log_W_1[,rand_days[1]])
imgPlt_wn11 <- imagePlot(coor1, wN11)
set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp11[[1]],
           imgPlt_temp11[[2]],
           imgPlt_temp11[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn11[[1]],
           imgPlt_wn11[[2]],
           imgPlt_wn11[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))

## Further validation:
# QQplot:
stats::qqnorm((wN11-mean(wN11))/sd(wN11))
abline(0,1, col='magenta', lty=2)

# Variogram:
bins <- seq( 0, 0.8,length.out=30)
dGrid<- seq(0, 0.8, length.out=200)
lookBins11 <- vgram(coor1, log_W_1[,rand_days[1]], N=30)
vMean11 <- lookBins11$stats["mean",]
binCenter11 <- lookBins11$centers

bplot.xy(lookBins11$d,
         lookBins11$vgram,
         breaks=lookBins11$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[1])),
         xlab='lag',
         #xaxt='n',
         #yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,8))
points(binCenter11, vMean11, col="red")
lines(dGrid,
      (cov_parameter_tile1[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile1[2],
                                        range_val = cov_parameter_tile1[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile1_dM[rand_days[1], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile1_dM[rand_days[1], 2],
                                                         range_val=cov_parameter_tile1_dM[rand_days[1], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-1-I.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp11[[1]],
           imgPlt_temp11[[2]],
           imgPlt_temp11[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn11[[1]],
           imgPlt_wn11[[2]],
           imgPlt_wn11[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[1]),
                          ", ", 
                          "Ext Day" == .(rand_days[1])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()

png("~/Desktop/QQ_variogram_Tile-1-I.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins11$d,
         lookBins11$vgram,
         breaks = lookBins11$breaks,
         outline = FALSE, 
         axes = FALSE,
         xlab='',
         ylab= '', 
         ylim = c(0,8))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter11,
       vMean11,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile1[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile1[2],
                                        range_val = cov_parameter_tile1[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile1_dM[rand_days[1], 1])^2+ vTrue(dGrid,
                                                         sigma2 = cov_parameter_tile1_dM[rand_days[1], 2],
                                                         range_val = cov_parameter_tile1_dM[rand_days[1], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN11-mean(wN11))/sd(wN11),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[1]),
                          ", ", 
                          "Ext Day" == .(rand_days[1])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()


## Day 2
sigma2_2 <- cov_parameter_tile1_dM[rand_days[2], 2]
tau2 <- cov_parameter_tile1_dM[rand_days[2], 1]
rho2 <- cov_parameter_tile1_dM[rand_days[2], 3]

sigma_mat12 <- (tau2^2) + sigma2_2* Matern(d=rdist1, 
                                           range=rho2,
                                           smoothness=1) # tau^2 is added for nugget effect
L2 <- t(chol(sigma_mat12)) # check for lower triangular
L_inv2 <- t(L2)%*%solve(sigma_mat12)

mean_wN12 <- as.numeric(fixed_constant_parameter_tile1_dM[rand_days[2]])
wN12 <- L_inv2%*%(log_W_1[,rand_days[2]])
imgPlt_wn12 <- imagePlot(coor1, wN12)
set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp12[[1]],
           imgPlt_temp12[[2]],
           imgPlt_temp12[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn12[[1]],
           imgPlt_wn12[[2]],
           imgPlt_wn12[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN12-mean(wN12))/sd(wN12))
abline(0,1, col='magenta', lty=2)

# Variogram:
bins <- seq( 0, 0.8,length.out=30)
dGrid<- seq(0, 0.8, length.out=200)
lookBins12 <- vgram(coor1, log_W_1[,rand_days[2]], N=30)
vMean12 <- lookBins12$stats["mean",]
binCenter12 <- lookBins12$centers

bplot.xy(lookBins12$d,
         lookBins12$vgram,
         breaks=lookBins12$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[2])),
         xlab='lag',
         #xaxt='n',
         #yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2)
points(binCenter12, vMean12, col="red")
lines(dGrid,
      (cov_parameter_tile1[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile1[2],
                                        range_val = cov_parameter_tile1[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile1_dM[rand_days[2], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile1_dM[rand_days[2], 2],
                                                         range_val=cov_parameter_tile1_dM[rand_days[2], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-1-II.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp12[[1]],
           imgPlt_temp12[[2]],
           imgPlt_temp12[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn12[[1]],
           imgPlt_wn12[[2]],
           imgPlt_wn12[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[1]),
                          ", ", 
                          "Ext Day" == .(rand_days[2])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()


png("~/Desktop/QQ_variogram_Tile-1-II.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins12$d,
         lookBins12$vgram,
         breaks = lookBins12$breaks,
         outline = FALSE, 
         axes = FALSE,
         xlab='',
         ylab= '')
#ylim = c(0,8))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter12,
       vMean12,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile1[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile1[2],
                                        range_val = cov_parameter_tile1[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile1_dM[rand_days[2], 1])^2+ vTrue(dGrid,
                                                         sigma2 = cov_parameter_tile1_dM[rand_days[2], 2],
                                                         range_val = cov_parameter_tile1_dM[rand_days[2], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN12-mean(wN12))/sd(wN12),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[1]),
                          ", ", 
                          "Ext Day" == .(rand_days[2])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()

## Day 3
sigma2_3 <- cov_parameter_tile1_dM[rand_days[3], 2]
tau3 <- cov_parameter_tile1_dM[rand_days[3], 1]
rho3 <- cov_parameter_tile1_dM[rand_days[3], 3]

sigma_mat13 <- (tau3^2) + sigma2_3* Matern(d=rdist1, 
                                           range=rho3,
                                           smoothness=1) # tau^2 is added for nugget effect
L3 <- t(chol(sigma_mat13)) # check for lower triangular
image(L3) # L1 is the cholesky factor

L_inv3 <- t(L3)%*%solve(sigma_mat13)
image(L_inv3)

# fixed_constant_3 <- as.numeric(fixed_constant_parameter_tile1_dM[rand_days[3]])
# mean_wN13 <- fixed_constant_3[1] 
wN13 <- L_inv1%*% (log_W_1[,rand_days[3]])
imgPlt_wn13 <- imagePlot(coor1, wN13)
set.panel(1,2)
par(mar=c(4,4,3,2))
image.plot(imgPlt_temp13[[1]],
           imgPlt_temp13[[2]],
           imgPlt_temp13[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn13[[1]],
           imgPlt_wn13[[2]],
           imgPlt_wn13[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN13-mean(wN13))/sd(wN13))
abline(0,1, col='magenta', lty=2)

# Variogram:
lookBins13 <- vgram(coor1, log_W_1[,rand_days[3]], N=30)
vMean13 <- lookBins13$stats["mean",]
binCenter13 <- lookBins13$centers

bplot.xy(lookBins13$d,
         lookBins13$vgram,
         breaks=lookBins13$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[2])),
         xlab='lag',
         xaxt='n',
         yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,4))
points(binCenter13, vMean13, col="red")
lines(dGrid,
      (cov_parameter_tile1[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile1[2],
                                        range_val = cov_parameter_tile1[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile1_dM[rand_days[3], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile1_dM[rand_days[3], 2],
                                                         range_val=cov_parameter_tile1_dM[rand_days[3], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)

## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-1-III.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp13[[1]],
           imgPlt_temp13[[2]],
           imgPlt_temp13[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn13[[1]],
           imgPlt_wn13[[2]],
           imgPlt_wn13[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[1]),
                          ", ", 
                          "Ext Day" == .(rand_days[3])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()


png("~/Desktop/QQ_variogram_Tile-1-III.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins13$d,
         lookBins13$vgram,
         breaks=lookBins13$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab= '',
         ylim=c(0,8))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter13,
       vMean13,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile1[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile1[2],
                                        range_val = cov_parameter_tile1[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile1_dM[rand_days[3], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile1_dM[rand_days[3], 2],
                                                         range_val=cov_parameter_tile1_dM[rand_days[3], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN13-mean(wN13))/sd(wN13),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[1]),
                          ", ", 
                          "Ext Day" == .(rand_days[3])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()


## Day 4
rdist1 <- rdist(coor1) # dim: 224 224
sigma2_4 <- cov_parameter_tile1_dM[rand_days[4], 2]
tau4 <- cov_parameter_tile1_dM[rand_days[4], 1]
rho4 <- cov_parameter_tile1_dM[rand_days[4], 3]

sigma_mat14 <- (tau4^2) + sigma2_4* Matern(d=rdist1, 
                                           range=rho4,
                                           smoothness=1) # tau^2 is added for nugget effect
L4 <- t(chol(sigma_mat14)) # check for lower triangular
L_inv4 <- t(L4)%*%solve(sigma_mat14)

mean_wN14 <- as.numeric(fixed_constant_parameter_tile1_dM[rand_days[4]])
wN14 <- L_inv1%*% (log_W_1[,rand_days[4]])
imgPlt_wn14 <- imagePlot(coor1, wN14)
set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp14[[1]],
           imgPlt_temp14[[2]],
           imgPlt_temp14[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn14[[1]],
           imgPlt_wn14[[2]],
           imgPlt_wn14[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN14-mean(wN14))/sd(wN14))
abline(0,1, col='magenta', lty=2)

# Variogram:
bins <- seq( 0, 0.8,length.out=30)
dGrid<- seq(0, 0.8, length.out=200)
lookBins14 <- vgram(coor1, log_W_1[,rand_days[4]], N=30)
vMean14 <- lookBins14$stats["mean",]
binCenter14 <- lookBins14$centers

bplot.xy(lookBins14$d,
         lookBins14$vgram,
         breaks=lookBins14$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[4])),
         xlab='lag',
         #xaxt='n',
         #yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,3))
points(binCenter14, vMean14, col="red")
lines(dGrid,
      (cov_parameter_tile1[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile1[2],
                                        range_val = cov_parameter_tile1[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile1_dM[rand_days[4], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile1_dM[rand_days[4], 2],
                                                         range_val=cov_parameter_tile1_dM[rand_days[4], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-1-IV.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp14[[1]],
           imgPlt_temp14[[2]],
           imgPlt_temp14[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn14[[1]],
           imgPlt_wn14[[2]],
           imgPlt_wn14[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[1]),
                          ", ", 
                          "Ext Day" == .(rand_days[4])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()

png("~/Desktop/QQ_variogram_Tile-1-IV.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins14$d,
         lookBins14$vgram,
         breaks = lookBins14$breaks,
         outline = FALSE, 
         axes = FALSE,
         xlab='',
         ylab= '', 
         ylim = c(0,3))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter14,
       vMean14,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile1[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile1[2],
                                        range_val = cov_parameter_tile1[3],
                                        nu=0.5),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile1_dM[rand_days[4], 1])^2+ vTrue(dGrid,
                                                         sigma2 = cov_parameter_tile1_dM[rand_days[4], 2],
                                                         range_val = cov_parameter_tile1_dM[rand_days[4], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN14-mean(wN14))/sd(wN14),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[1]),
                          ", ", 
                          "Ext Day" == .(rand_days[1])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()


## -- Tile2: 3 --:
coor2 <- expand.grid(gridPoints16x16Boxes[[rand_tiles[2]]]$lon, 
                     gridPoints16x16Boxes[[rand_tiles[2]]]$lat)
dim(coor2) # 256x2

log_W_2 <- log_W_process_across_boxes[[rand_tiles[2]]]
dim(log_W_2) # 256x86

# Fit Gaussian process for all the realizations all together:
fit_GP_tile2 <- spatialProcess(coor2,
                               log_W_2,
                               na.rm=TRUE,
                               collapseFixedEffect = FALSE,
                               smoothness=1,
                               mKrig.args=list(m=1)) # just intercept term for the fixed part of the model


cov_parameter_tile2 <- as.numeric(fit_GP_tile2$MLESummary[6:9])
fixed_constant_parameter_tile2 <- (fit_GP_tile2$d)

## Fit Gaussian process to log_W process for Tile 11:
cov_parameter_tile2_dM <- matrix(NA,
                                 nrow=ncol(log_W_2),
                                 ncol=4)
fixed_constant_parameter_tile2_dM <- rep(NA)

# Covariance matrix varies across the space and time:
for(i in 1:ncol(log_W_2))
{
  y_process <- log_W_2[,i]
  print(i)
  
  tryCatch({
    fit_GP <- spatialProcess(coor2,
                             y_process,
                             na.rm=TRUE,
                             smoothness=1, # Changed from 0.5 to 1
                             mKrig.args=list(m=1))
    cov_parameter_tile2_dM[i,] <- as.numeric(fit_GP$MLESummary[6:9])
    fixed_constant_parameter_tile2_dM[i] <- as.numeric(fit_GP$d)
  }, error = function(e) {
    cat("Error occurred for i=", i, ": ", conditionMessage(e), "\n")}
  )
}

## Variogram function for fitted model:
vTrue <- function(d, sigma2, range_val, nu)
{
  # vario_vals <- sigma2*(1-exp(-d/range_val))
  vario_vals <- sigma2*(1-Matern(d=d, range=range_val, smoothness = nu))
  return(vario_vals)
}

## Randomly selected day:
set.seed(222)
rand_days <- round(runif(n=4, 1, 86))
print(rand_days)

## Log W process for randomly selected days for Tile 1:
imgPlt_temp21 <- imagePlot(coor2,
                           log_W_2[, rand_days[1]]) # day 1
imgPlt_temp22 <- imagePlot(coor2,
                           log_W_2[, rand_days[2]]) # day 2
imgPlt_temp23 <- imagePlot(coor2,
                           log_W_2[, rand_days[3]]) # day 3
imgPlt_temp24 <- imagePlot(coor2,
                           log_W_2[, rand_days[4]]) # day 4

## Day 1
rdist2 <- rdist(coor2) # dim: 224 224
sigma2 <- cov_parameter_tile2_dM[rand_days[1], 2]
tau <- cov_parameter_tile2_dM[rand_days[1], 1]
rho <- cov_parameter_tile2_dM[rand_days[1], 3]

sigma_mat21 <- (tau^2) + sigma2* Matern(d=rdist2, 
                                        range=rho,
                                        smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat21)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat21)

mean_wN21 <- as.numeric(fixed_constant_parameter_tile2_dM[rand_days[1]])
wN21 <- L_inv%*%(log_W_2[ ,rand_days[1]])
imgPlt_wn21 <- imagePlot(coor2, wN21)

set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp21[[1]],
           imgPlt_temp21[[2]],
           imgPlt_temp21[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn21[[1]],
           imgPlt_wn21[[2]],
           imgPlt_wn21[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN21-mean(wN21))/sd(wN21))
abline(0,1, col='magenta', lty=2)

# Variogram:
lookBins21 <- vgram(coor2, log_W_2[,rand_days[1]], N=30)
vMean21 <- lookBins21$stats["mean",]
binCenter21 <- lookBins21$centers

bplot.xy(lookBins21$d,
         lookBins21$vgram,
         breaks=lookBins21$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[1])),
         xlab='lag',
         #xaxt='n',
         #yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,3))
points(binCenter21, vMean21, col="red")
lines(dGrid,
      (cov_parameter_tile2[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile2[2],
                                        range_val=cov_parameter_tile2[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile2_dM[rand_days[1], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile2_dM[rand_days[1], 2],
                                                         range_val=cov_parameter_tile2_dM[rand_days[1], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-2-I.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp21[[1]],
           imgPlt_temp21[[2]],
           imgPlt_temp21[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn21[[1]],
           imgPlt_wn21[[2]],
           imgPlt_wn21[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[2]),
                          ", ", 
                          "Ext Day" == .(rand_days[1])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()

png("~/Desktop/QQ_variogram_Tile-2-I.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins21$d,
         lookBins21$vgram,
         breaks = lookBins21$breaks,
         outline = FALSE, 
         axes = FALSE,
         xlab='',
         ylab= '', 
         ylim = c(0,3))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter21,
       vMean21,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile2[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile2[2],
                                        range_val = cov_parameter_tile2[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile2_dM[rand_days[1], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile2_dM[rand_days[1], 2],
                                                         range_val=cov_parameter_tile2_dM[rand_days[1], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN21-mean(wN21))/sd(wN21),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main=bquote(paste("Tile" == .(rand_tiles[2]),
                        ", ", 
                        "Ext Day" == .(rand_days[1])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()


## Day 2
sigma2 <- cov_parameter_tile2_dM[rand_days[2], 2]
tau <- cov_parameter_tile2_dM[rand_days[2], 1]
rho <- cov_parameter_tile2_dM[rand_days[2], 3]

sigma_mat22 <- (tau^2) + sigma2* Matern(d=rdist2, 
                                        range=rho,
                                        smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat22)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat22)

mean_wN22 <- as.numeric(fixed_constant_parameter_tile2_dM[rand_days[2]])
wN22 <- L_inv%*% (log_W_2[,rand_days[2]])
imgPlt_wn22 <- imagePlot(coor2, wN22)

set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp22[[1]],
           imgPlt_temp22[[2]],
           imgPlt_temp22[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn22[[1]],
           imgPlt_wn22[[2]],
           imgPlt_wn22[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN22-mean(wN22))/sd(wN22))
abline(0,1, col='magenta', lty=2)

# Variogram:
lookBins22 <- vgram(coor2, log_W_2[,rand_days[2]], N=30)
vMean22 <- lookBins22$stats["mean",]
binCenter22 <- lookBins22$centers

bplot.xy(lookBins22$d,
         lookBins22$vgram,
         breaks=lookBins22$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[2])),
         xlab='lag',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2)
points(binCenter22, vMean22, col="red")
lines(dGrid,
      (cov_parameter_tile2[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile2[2],
                                        range_val=cov_parameter_tile2[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile2_dM[rand_days[2], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile2_dM[rand_days[2], 2],
                                                         range_val=cov_parameter_tile2_dM[rand_days[2], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-2-II.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp22[[1]],
           imgPlt_temp22[[2]],
           imgPlt_temp22[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn22[[1]],
           imgPlt_wn22[[2]],
           imgPlt_wn22[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[2]),
                          ", ", 
                          "Ext Day" == .(rand_days[2])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()

png("~/Desktop/QQ_variogram_Tile-2-II.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins22$d,
         lookBins22$vgram,
         breaks=lookBins22$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab= '',
         ylim=c(0,2))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter22,
       vMean22,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile2[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile2[2],
                                        range_val=cov_parameter_tile2[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile2_dM[rand_days[2], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile2_dM[rand_days[2], 2],
                                                         range_val=cov_parameter_tile2_dM[rand_days[2], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN22-mean(wN22))/sd(wN22),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[2]),
                          ", ", 
                          "Ext Day" == .(rand_days[2])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()

## Day 3
sigma2 <- cov_parameter_tile2_dM[rand_days[3], 2]
tau <- cov_parameter_tile2_dM[rand_days[3], 1]
rho <- cov_parameter_tile2_dM[rand_days[3], 3]

sigma_mat23 <- (tau^2) + sigma2* Matern(d=rdist2, 
                                        range=rho,
                                        smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat23)) # check for lower triangular
image(L) # L1 is the cholesky factor

L_inv <- t(L)%*%solve(sigma_mat23)
image(L_inv)

# fixed_constant_3 <- as.numeric(fixed_constant_parameter_tile1_dM[rand_days[3]])
# mean_wN13 <- fixed_constant_3[1] 
wN23 <- L_inv%*% (log_W_2[,rand_days[3]])
imgPlt_wn23 <- imagePlot(coor2, wN23)
set.panel(1,2)
par(mar=c(4,4,3,2))
image.plot(imgPlt_temp23[[1]],
           imgPlt_temp23[[2]],
           imgPlt_temp23[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn23[[1]],
           imgPlt_wn23[[2]],
           imgPlt_wn23[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN23-mean(wN23))/sd(wN23))
abline(0,1, col='magenta', lty=2)

# Variogram:
lookBins23 <- vgram(coor2, log_W_2[,rand_days[3]], N=30)
vMean23 <- lookBins23$stats["mean",]
binCenter23 <- lookBins23$centers

bplot.xy(lookBins23$d,
         lookBins23$vgram,
         breaks=lookBins23$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[3])),
         xlab='lag',
         xaxt='n',
         yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,4))
points(binCenter23, vMean23, col="red")
lines(dGrid,
      (cov_parameter_tile2[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile2[2],
                                        range_val = cov_parameter_tile2[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile2_dM[rand_days[3], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile2_dM[rand_days[3], 2],
                                                         range_val=cov_parameter_tile2_dM[rand_days[3], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)

## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-2-III.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp23[[1]],
           imgPlt_temp23[[2]],
           imgPlt_temp23[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn23[[1]],
           imgPlt_wn23[[2]],
           imgPlt_wn23[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[2]),
                          ", ", 
                          "Ext Day" == .(rand_days[3])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()

png("~/Desktop/QQ_variogram_Tile-2-III.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins23$d,
         lookBins23$vgram,
         breaks=lookBins23$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab= '',
         ylim=c(0,8))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter23,
       vMean23,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile2[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile2[2],
                                        range_val=cov_parameter_tile2[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile2_dM[rand_days[3], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile2_dM[rand_days[3], 2],
                                                         range_val=cov_parameter_tile2_dM[rand_days[3], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN23-mean(wN23))/sd(wN23),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[2]),
                          ", ", 
                          "Ext Day"== .(rand_days[3])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()


## Day 4
rdist2 <- rdist(coor2) # dim: 224 224
sigma2 <- cov_parameter_tile2_dM[rand_days[4], 2]
tau <- cov_parameter_tile2_dM[rand_days[4], 1]
rho <- cov_parameter_tile2_dM[rand_days[4], 3]

sigma_mat24 <- (tau^2) + sigma2* Matern(d=rdist2, 
                                        range=rho,
                                        smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat24)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat24)

mean_wN24 <- as.numeric(fixed_constant_parameter_tile2_dM[rand_days[4]])
wN24 <- L_inv%*% (log_W_2[,rand_days[4]])
imgPlt_wn24 <- imagePlot(coor2, wN24)

set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp24[[1]],
           imgPlt_temp24[[2]],
           imgPlt_temp24[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn24[[1]],
           imgPlt_wn24[[2]],
           imgPlt_wn24[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN24-mean(wN24))/sd(wN24))
abline(0,1, col='magenta', lty=2)

# Variogram:
bins <- seq( 0, 0.8,length.out=30)
dGrid<- seq(0, 0.8, length.out=200)
lookBins24 <- vgram(coor1, log_W_2[,rand_days[4]], N=30)
vMean24 <- lookBins24$stats["mean",]
binCenter24 <- lookBins24$centers

bplot.xy(lookBins24$d,
         lookBins24$vgram,
         breaks=lookBins24$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[4])),
         xlab='lag',
         #xaxt='n',
         #yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,3))
points(binCenter24, vMean24, col="red")
lines(dGrid,
      (cov_parameter_tile2[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile2[2],
                                        range_val=cov_parameter_tile2[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile2_dM[rand_days[4], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile2_dM[rand_days[4], 2],
                                                         range_val=cov_parameter_tile2_dM[rand_days[4], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-2-IV.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp24[[1]],
           imgPlt_temp24[[2]],
           imgPlt_temp24[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn24[[1]],
           imgPlt_wn24[[2]],
           imgPlt_wn24[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main=bquote(paste("Tile" == .(rand_tiles[2]),
                        ", ", 
                        "Ext Day" == .(rand_days[4])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()

png("~/Desktop/QQ_variogram_Tile-2-IV.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins24$d,
         lookBins24$vgram,
         breaks=lookBins24$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab='', 
         ylim=c(0,3))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter24,
       vMean24,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile2[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile2[2],
                                        range_val=cov_parameter_tile2[3],
                                        nu=0.5),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile2_dM[rand_days[4], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile2_dM[rand_days[4], 2],
                                                         range_val=cov_parameter_tile2_dM[rand_days[4], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN24-mean(wN24))/sd(wN24),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[2]),
                          ", ", 
                          "Ext Day" == .(rand_days[4])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()



## -- Tile3: 15 --
coor3 <- expand.grid(gridPoints16x16Boxes[[rand_tiles[3]]]$lon, 
                     gridPoints16x16Boxes[[rand_tiles[3]]]$lat)
dim(coor3) # 256x2

log_W_3 <- log_W_process_across_boxes[[rand_tiles[3]]]
dim(log_W_3) # 256x86

## Clean the NA values:
combined_data <- cbind(coor3, log_W_3)

# Remove rows with any NA values
cleaned_data <- na.omit(combined_data)

# Extract the cleaned coordinates and log_W_3 matrix
coor3_cleaned <- cleaned_data[, 1:2]
log_W_3_cleaned <- cleaned_data[, -c(1:2)]

# Check dimensions
dim(coor3_cleaned)
dim(log_W_3_cleaned)

# Fit Gaussian process for all the realizations all together:
fit_GP_tile3 <- spatialProcess(coor3_cleaned,
                               log_W_3_cleaned,
                               na.rm=TRUE,
                               collapseFixedEffect = FALSE,
                               smoothness=1,
                               mKrig.args=list(m=1)) # we fit linear mean model, with intercept and slope terms

cov_parameter_tile3 <- as.numeric(fit_GP_tile3$MLESummary[6:9])
fixed_constant_parameter_tile3  <- as.matrix(fit_GP_tile3$d)

## Fit Gaussian process to log_W process for Tile 11:
cov_parameter_tile3_dM <- matrix(NA,
                                 nrow=ncol(log_W_3_cleaned),
                                 ncol=4)
fixed_constant_parameter_tile3_dM <- rep(NA)
# Covariance matrix varies across the space and time:
for(i in 1:ncol(log_W_3_cleaned))
{
  y_process <- log_W_3_cleaned[,i]
  print(i)
  
  tryCatch({
    fit_GP <- spatialProcess(coor3_cleaned,
                             y_process,
                             na.rm=TRUE,
                             smoothness=1, # Changed from 0.5 to 1
                             mKrig.args=list(m=1))
    cov_parameter_tile3_dM[i,] <- as.numeric(fit_GP$MLESummary[6:9])
    fixed_constant_parameter_tile3_dM[i] <- as.numeric(fit_GP$d)
  }, error = function(e) {
    cat("Error occurred for i=", i, ": ", conditionMessage(e), "\n")}
  )
}

## Variogram function for fitted model:
vTrue <- function(d, sigma2, range_val, nu)
{
  # vario_vals <- sigma2*(1-exp(-d/range_val))
  vario_vals <- sigma2*(1-Matern(d=d, range=range_val, smoothness = nu))
  return(vario_vals)
}

## Randomly selected day:
set.seed(4)
rand_days <- round(runif(n=4, 1, 86))
print(rand_days)

## Log W process for randomly selected days for Tile 1:
imgPlt_temp31 <- imagePlot(coor3_cleaned,
                           log_W_3_cleaned[, rand_days[1]]) # day 1
imgPlt_temp32 <- imagePlot(coor3_cleaned,
                           log_W_3_cleaned[, rand_days[2]]) # day 2
imgPlt_temp33 <- imagePlot(coor3_cleaned,
                           log_W_3_cleaned[, rand_days[3]]) # day 3
imgPlt_temp34 <- imagePlot(coor3_cleaned,
                           log_W_3_cleaned[, rand_days[4]]) # day 4

## Day 1
rdist3 <- rdist(coor3_cleaned) # dim: 117 117
sigma2 <- cov_parameter_tile3_dM[rand_days[1], 2]
tau <- cov_parameter_tile3_dM[rand_days[1], 1]
rho <- cov_parameter_tile3_dM[rand_days[1], 3]

sigma_mat31 <- (tau^2) + sigma2* Matern(d=rdist3,
                                        range=rho,
                                        smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat31)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat31)

mean_wN31 <- as.numeric(fixed_constant_parameter_tile3_dM[rand_days[1]])
wN31 <- L_inv%*% (log_W_3_cleaned[,rand_days[1]])
imgPlt_wn31 <- imagePlot(coor3_cleaned, wN31)
set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp31[[1]],
           imgPlt_temp31[[2]],
           imgPlt_temp31[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn31[[1]],
           imgPlt_wn31[[2]],
           imgPlt_wn31[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN31-mean(wN31))/sd(wN31))
abline(0,1, col='magenta', lty=2)

# Variogram:
lookBins31 <- vgram(coor3_cleaned, log_W_3_cleaned[,rand_days[1]], N=30)
vMean31 <- lookBins31$stats["mean",]
binCenter31 <- lookBins31$centers

bplot.xy(lookBins31$d,
         lookBins31$vgram,
         breaks=lookBins31$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[1])),
         xlab='lag',
         #xaxt='n',
         #yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,3))
points(binCenter31, vMean31, col="red")
lines(dGrid,
      (cov_parameter_tile3[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile3[2],
                                        range_val=cov_parameter_tile3[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile3_dM[rand_days[1], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile3_dM[rand_days[1], 2],
                                                         range_val=cov_parameter_tile3_dM[rand_days[1], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-3-I.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp31[[1]],
           imgPlt_temp31[[2]],
           imgPlt_temp31[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn31[[1]],
           imgPlt_wn31[[2]],
           imgPlt_wn31[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[1])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()

png("~/Desktop/QQ_variogram_Tile-3-I.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins31$d,
         lookBins31$vgram,
         breaks=lookBins31$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab= '', 
         ylim = c(0,3))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter31,
       vMean31,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile3[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile3[2],
                                        range_val=cov_parameter_tile3[3],
                                        nu=0.5),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile3_dM[rand_days[1], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile3_dM[rand_days[1], 2],
                                                         range_val=cov_parameter_tile3_dM[rand_days[1], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN31-mean(wN31))/sd(wN31),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[1])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()


## Day 2
sigma2 <- cov_parameter_tile3_dM[rand_days[2], 2]
tau <- cov_parameter_tile3_dM[rand_days[2], 1]
rho <- cov_parameter_tile3_dM[rand_days[2], 3]

sigma_mat32 <- (tau^2) + sigma2* Matern(d=rdist3, # 117x117
                                        range=rho,
                                        smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat32)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat32)

mean_wN32 <- as.numeric(fixed_constant_parameter_tile3_dM[rand_days[2]])
wN32 <- L_inv%*% (log_W_3_cleaned[,rand_days[2]])
imgPlt_wn32 <- imagePlot(coor3_cleaned, wN32)

set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp32[[1]],
           imgPlt_temp32[[2]],
           imgPlt_temp32[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn32[[1]],
           imgPlt_wn32[[2]],
           imgPlt_wn32[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN32-mean(wN32))/sd(wN32))
abline(0,1, col='magenta', lty=2)

# Variogram:
lookBins32 <- vgram(coor3_cleaned, log_W_3_cleaned[,rand_days[2]], N=30)
vMean32 <- lookBins32$stats["mean",]
binCenter32 <- lookBins32$centers

bplot.xy(lookBins32$d,
         lookBins32$vgram,
         breaks=lookBins32$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[2])),
         xlab='lag',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2)
points(binCenter32, vMean32, col="red")
lines(dGrid,
      (cov_parameter_tile3[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile3[2],
                                        range_val=cov_parameter_tile3[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile3_dM[rand_days[2], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile3_dM[rand_days[2], 2],
                                                         range_val=cov_parameter_tile3_dM[rand_days[2], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-3-II.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp32[[1]],
           imgPlt_temp32[[2]],
           imgPlt_temp32[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn32[[1]],
           imgPlt_wn32[[2]],
           imgPlt_wn32[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[2])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()


png("~/Desktop/QQ_variogram_Tile-3-II.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins32$d,
         lookBins32$vgram,
         breaks=lookBins32$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab= '',
         ylim=c(0,6))
#ylim = c(0,8))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter32,
       vMean32,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile3[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile3[2],
                                        range_val=cov_parameter_tile3[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile3_dM[rand_days[2], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile3_dM[rand_days[2], 2],
                                                         range_val=cov_parameter_tile3_dM[rand_days[2], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN32-mean(wN32))/sd(wN32),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[2])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()

## Day 3
sigma2 <- cov_parameter_tile3_dM[rand_days[3], 2]
tau <- cov_parameter_tile3_dM[rand_days[3], 1]
rho <- cov_parameter_tile3_dM[rand_days[3], 3]

sigma_mat33 <- (tau^2) + sigma2* Matern(d=rdist3,
                                        range=rho,
                                        smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat33)) # check for lower triangular
image(L) # L1 is the cholesky factor

L_inv <- t(L)%*%solve(sigma_mat33)
image(L_inv)

# fixed_constant_3 <- as.numeric(fixed_constant_parameter_tile1_dM[rand_days[3]])
# mean_wN13 <- fixed_constant_3[1] 
wN33 <- L_inv%*% (log_W_3_cleaned[,rand_days[3]])
imgPlt_wn33 <- imagePlot(coor3_cleaned, wN33)
set.panel(1,2)
par(mar=c(4,4,3,2))
image.plot(imgPlt_temp33[[1]],
           imgPlt_temp33[[2]],
           imgPlt_temp33[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn33[[1]],
           imgPlt_wn33[[2]],
           imgPlt_wn33[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN33-mean(wN33))/sd(wN33))
abline(0,1, col='magenta', lty=2)

# Variogram:
lookBins33 <- vgram(coor3_cleaned, log_W_3_cleaned[,rand_days[3]], N=30)
vMean33 <- lookBins33$stats["mean",]
binCenter33 <- lookBins33$centers

bplot.xy(lookBins33$d,
         lookBins33$vgram,
         breaks=lookBins33$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[3])),
         xlab='lag',
         xaxt='n',
         yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,4))
points(binCenter33, vMean33, col="red")
lines(dGrid,
      (cov_parameter_tile3[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile3[2],
                                        range_val = cov_parameter_tile3[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile3_dM[rand_days[3], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile3_dM[rand_days[3], 2],
                                                         range_val=cov_parameter_tile3_dM[rand_days[3], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)

## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-3-III.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp33[[1]],
           imgPlt_temp33[[2]],
           imgPlt_temp33[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn33[[1]],
           imgPlt_wn33[[2]],
           imgPlt_wn33[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[3])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()


png("~/Desktop/QQ_variogram_Tile-3-III.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins33$d,
         lookBins33$vgram,
         breaks=lookBins33$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab= '',
         ylim=c(0,3))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter33,
       vMean33,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile3[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile3[2],
                                        range_val=cov_parameter_tile3[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile3_dM[rand_days[3], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile3_dM[rand_days[3], 2],
                                                         range_val=cov_parameter_tile3_dM[rand_days[3], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN33-mean(wN33))/sd(wN33),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[3])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()


## Day 4
rdist3 <- rdist(coor3_cleaned) # dim: 224 224
sigma2 <- cov_parameter_tile3_dM[rand_days[4], 2]
tau <- cov_parameter_tile3_dM[rand_days[4], 1]
rho <- cov_parameter_tile3_dM[rand_days[4], 3]

sigma_mat34 <- (tau^2) + sigma2* Matern(d=rdist3, 
                                        range=rho,
                                        smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat34)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat34)


mean_wN34 <- as.numeric(fixed_constant_parameter_tile3_dM[rand_days[4]])
wN34 <- L_inv%*% (log_W_3_cleaned[,rand_days[4]])
imgPlt_wn34 <- imagePlot(coor3_cleaned, wN34)
set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp34[[1]],
           imgPlt_temp34[[2]],
           imgPlt_temp34[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn34[[1]],
           imgPlt_wn34[[2]],
           imgPlt_wn34[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
## Further validation:
# QQplot:
stats::qqnorm((wN34-mean(wN34))/sd(wN34))
abline(0,1, col='magenta', lty=2)

# Variogram:
lookBins34 <- vgram(coor3_cleaned, log_W_3_cleaned[,rand_days[4]], N=30)
vMean34 <- lookBins34$stats["mean",]
binCenter34 <- lookBins34$centers

bplot.xy(lookBins34$d,
         lookBins34$vgram,
         breaks=lookBins34$breaks,
         outline=FALSE, 
         main=bquote('Ext day' ~ .(rand_days[4])),
         xlab='lag',
         #xaxt='n',
         #yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,3))
points(binCenter34, vMean34, col="red")
lines(dGrid,
      (cov_parameter_tile3[1])^2+ vTrue(dGrid,
                                        sigma2 = cov_parameter_tile3[2],
                                        range_val = cov_parameter_tile3[3],
                                        nu=1),
      col='blue',
      lwd=2,
      lty=2)
lines(dGrid,
      (cov_parameter_tile3_dM[rand_days[4], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile3_dM[rand_days[4], 2],
                                                         range_val=cov_parameter_tile3_dM[rand_days[4], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN-3-IV.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp34[[1]],
           imgPlt_temp34[[2]],
           imgPlt_temp34[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn34[[1]],
           imgPlt_wn34[[2]],
           imgPlt_wn34[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[4])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-2.5)
dev.off()

png("~/Desktop/QQ_variogram_Tile-3-IV.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins34$d,
         lookBins34$vgram,
         breaks=lookBins34$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab='', 
         ylim = c(0,3))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter34,
       vMean34,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter_tile3[1])^2+ vTrue(dGrid,
                                        sigma2=cov_parameter_tile3[2],
                                        range_val=cov_parameter_tile3[3],
                                        nu=1),
      col='blue',
      lwd=3,
      lty=2)
lines(dGrid,
      (cov_parameter_tile3_dM[rand_days[4], 1])^2+ vTrue(dGrid,
                                                         sigma2=cov_parameter_tile3_dM[rand_days[4], 2],
                                                         range_val=cov_parameter_tile3_dM[rand_days[4], 3], 
                                                         nu=1),
      col='magenta', 
      lwd=3, 
      lty=2)
## QQ plot
stats::qqnorm((wN34-mean(wN34))/sd(wN34),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[4])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()
##########################################################################

## -- Fitting Gaussian Process on the entire spatial domain --:
risk_averages <- apply(X_star,
                       1, 
                       geometricMean)
print(length(risk_averages)) # length 2852
threshold <- quantile(risk_averages,
                      0.97)

extreme_ds <- which(risk_averages > threshold) # finding the extreme days
print(length(extreme_ds))

X_star_extremes <- (X_star[extreme_ds,])
R <- risk_averages[extreme_ds]
Y_r <- X_star_extremes/threshold

# Using the spectral decomposition:
log_W <- log(Y_r) - log(R)
dim(log_W) #dimension: 2852   86


## Randomly selected day:
set.seed(333)
r_days <- round(runif(n=3, min=1, max=nrow(log_W)))
print(r_days)

## Getting image plot:
set.panel(1,3)
tempG1 <- imagePlot(loc[1:2000, ], log_W[r_days[1],1:2000])
image.plot(tempG1[[1]],
           tempG1[[2]],
           tempG1[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
tempG2 <- imagePlot(loc[1:2000, ], log_W[r_days[2],1:2000])
image.plot(tempG2[[1]],
           tempG2[[2]],
           tempG2[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
tempG3 <- imagePlot(loc[1:2000, ], log_W[r_days[3],1:2000])
image.plot(tempG3[[1]],
           tempG3[[2]],
           tempG3[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))

# Picking random days: 
# Day 1
fit_GP_all1 <- spatialProcess(loc[1:2000, ],
                              log_W[r_days[1],1:2000],
                              na.rm=TRUE,
                              smoothness=1,
                              mKrig.args=list(m=1)) # we

cov_parameter1 <- as.numeric(fit_GP_all1$MLESummary[6:9])
fixed_constant1  <- as.matrix(fit_GP_all1$d)


# Day 2
fit_GP_all2 <- spatialProcess(loc[1:2000, ],
                              log_W[r_days[2], 1:2000],
                              na.rm=TRUE,
                              smoothness=1,
                              mKrig.args=list(m=1)) # we

cov_parameter2 <- as.numeric(fit_GP_all2$MLESummary[6:9])
fixed_constant2  <- as.matrix(fit_GP_all2$d)

# Day 3
fit_GP_all3 <- spatialProcess(loc[1:2000, ],
                              log_W[r_days[3], 1:2000],
                              na.rm=TRUE,
                              smoothness=1,
                              mKrig.args=list(m=1)) # we

cov_parameter3 <- as.numeric(fit_GP_all3$MLESummary[6:9])
fixed_constant3  <- as.matrix(fit_GP_all3$d)

# common distance matrix
rdist_common <- rdist(loc[1:2000, ]) # dim: 224 224

# Day 1
sigma2 <- cov_parameter1[2]
tau <- cov_parameter1[1]
rho <- cov_parameter1[3]

sigma_mat1 <- (tau^2) + sigma2* Matern(d=rdist_common, 
                                       range=rho,
                                       smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat1)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat1)


mean_wN1 <- as.numeric(fixed_constant1)
wN1 <- L_inv%*% (log_W[r_days[1],1:2000])

imgPlt_wn1 <- imagePlot(loc[1:2000, ], wN1)
imgPlt_temp1 <- imagePlot(loc[1:2000, ], log_W[r_days[1],1:2000])

set.panel(1,2)
par(mar=c(2,2,2,1.5))
image.plot(imgPlt_temp1[[1]],
           imgPlt_temp1[[2]],
           imgPlt_temp1[[3]],
           main='Spatial Process',
           xlab='lon',
           ylab='lat',
           col = hcl.colors(12, "YlOrRd", rev = TRUE))
image.plot(imgPlt_wn1[[1]],
           imgPlt_wn1[[2]],
           imgPlt_wn1[[3]],
           main='Approx White Noise',
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev = TRUE))


## Further validation:
# QQplot:
stats::qqnorm(wN1-mean(wN1)/sd(wN1))
abline(0,1, col='magenta', lty=2)

# Variogram:
bins <- seq( 0, 4, length.out=30)
dGrid<- seq(0, 4, length.out=200)
lookBins1 <- vgram(loc[1:2000, ],  log_W[r_days[1],1:2000], N=30)
vMean1 <- lookBins1$stats["mean",]
binCenter1 <- lookBins1$centers

bplot.xy(lookBins1$d,
         lookBins1$vgram,
         breaks=lookBins1$breaks,
         outline=FALSE, 
         #main=bquote('Ext day' ~ .(rand_days[4])),
         xlab='lag',
         #xaxt='n',
         #yaxt='n',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,3))
points(binCenter1, vMean1, col="red")
lines(dGrid,
      (cov_parameter1[1])^2+ vTrue(dGrid,
                                   sigma2=cov_parameter1[2],
                                   range_val=cov_parameter1[3], 
                                   nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN_biggerSpatialDomainDay1.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp1[[1]],
           imgPlt_temp1[[2]],
           imgPlt_temp1[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn1[[1]],
           imgPlt_wn1[[2]],
           imgPlt_wn1[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Ext Day" == .(r_days[1])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-0.5)
dev.off()

png("~/Desktop/QQ_variogram_biggerSpatialDomainDay1.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins1$d,
         lookBins1$vgram,
         breaks=lookBins1$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab='', 
         ylim = c(0,3))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter1,
       vMean1,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter1[1])^2+ vTrue(dGrid,
                                   sigma2=cov_parameter1[2],
                                   range_val=cov_parameter1[3],
                                   nu=1),
      col='magenta',
      lwd=3,
      lty=2)
## QQ plot
stats::qqnorm((wN1-mean(wN1))/sd(wN1),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[4])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()

# Day 2
sigma2 <- cov_parameter2[2]
tau <- cov_parameter2[1]
rho <- cov_parameter2[3]

sigma_mat2 <- (tau^2) + sigma2* Matern(d=rdist_common, 
                                       range=rho,
                                       smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat2)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat2)

mean_wN2 <- as.numeric(fixed_constant2)
wN2 <- L_inv%*% (log_W[r_days[2],1:2000])

imgPlt_wn2 <- imagePlot(loc[1:2000, ], wN2)
imgPlt_temp2 <- imagePlot(loc[1:2000, ], log_W[r_days[2], 1:2000])
lookBins2 <- vgram(loc[1:2000, ],  log_W[r_days[2], 1:2000], N=30)

vMean2 <- lookBins2$stats["mean",]
binCenter2 <- lookBins2$centers

bplot.xy(lookBins2$d,
         lookBins2$vgram,
         breaks=lookBins2$breaks,
         outline=FALSE, 
         xlab='lag',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,0.8))
points(binCenter2, vMean2, col="red")
lines(dGrid,
      (cov_parameter2[1])^2+ vTrue(dGrid,
                                   sigma2=cov_parameter2[2],
                                   range_val=cov_parameter2[3], 
                                   nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN_biggerSpatialDomainDay2.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp2[[1]],
           imgPlt_temp2[[2]],
           imgPlt_temp2[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn2[[1]],
           imgPlt_wn2[[2]],
           imgPlt_wn2[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Ext Day" == .(r_days[2])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-0.5)
dev.off()

png("~/Desktop/QQ_variogram_biggerSpatialDomainDay2.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins2$d,
         lookBins2$vgram,
         breaks=lookBins2$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab='', 
         ylim = c(0,0.8))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter2,
       vMean2,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter2[1])^2+ vTrue(dGrid,
                                   sigma2=cov_parameter2[2],
                                   range_val=cov_parameter2[3],
                                   nu=1),
      col='magenta',
      lwd=3,
      lty=2)
## QQ plot
stats::qqnorm((wN2-mean(wN2))/sd(wN2),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[4])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()


# Day 3
sigma3 <- cov_parameter3[2]
tau <- cov_parameter3[1]
rho <- cov_parameter3[3]

sigma_mat3 <- (tau^2) + sigma3* Matern(d=rdist_common, 
                                       range=rho,
                                       smoothness=1) # tau^2 is added for nugget effect
L <- t(chol(sigma_mat3)) # check for lower triangular
L_inv <- t(L)%*%solve(sigma_mat3)

mean_wN3 <- as.numeric(fixed_constant3)
wN3 <- L_inv%*% (log_W[r_days[3],1:2000])

imgPlt_wn3 <- imagePlot(loc[1:2000, ], wN3)
imgPlt_temp3 <- imagePlot(loc[1:2000, ], log_W[r_days[3], 1:2000])
lookBins3 <- vgram(loc[1:2000, ],  log_W[r_days[3], 1:2000], N=30)

vMean3 <- lookBins3$stats["mean",]
binCenter3 <- lookBins3$centers

bplot.xy(lookBins3$d,
         lookBins3$vgram,
         breaks=lookBins3$breaks,
         outline=FALSE, 
         xlab='lag',
         ylab=expression(gamma), 
         cex.main=1.5,  # Increase main title font size
         cex.lab=1.4,   # Increase x and y axis labels font size
         cex.axis=1.2,
         ylim=c(0,2))
points(binCenter3, vMean3, col="red")
lines(dGrid,
      (cov_parameter3[1])^2+ vTrue(dGrid,
                                   sigma2=cov_parameter3[2],
                                   range_val=cov_parameter3[3], 
                                   nu=1),
      col='magenta', 
      lwd=2, 
      lty=2)


## Plots to generate:
png("~/Desktop/Comparing_log_W_with_WN_biggerSpatialDomainDay3.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# title_expression <- bquote('Ext Day' == .(rand_days[1]))
image.plot(imgPlt_temp3[[1]],
           imgPlt_temp3[[2]],
           imgPlt_temp3[[3]],
           cex.main=3,
           main='log W process',
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
# title_expression <- title_expression <- bquote('Ext Day' == .(rand_days[2]))
image.plot(imgPlt_wn3[[1]],
           imgPlt_wn3[[2]],
           imgPlt_wn3[[3]],
           main='White Noise',
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.main=3,
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2)
title_main = bquote(paste("Ext Day" == .(r_days[3])))
mtext(title_main, side=3, line= 5, cex = 3.5, adj =-0.5)
dev.off()

png("~/Desktop/QQ_variogram_biggerSpatialDomainDay3.png",
    units = "in",
    width = 20,
    height = 9,
    res = 200)
par(mfrow = c(1,2),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6))
bplot.xy(lookBins3$d,
         lookBins3$vgram,
         breaks=lookBins3$breaks,
         outline=FALSE, 
         axes=FALSE,
         xlab='',
         ylab='', 
         ylim = c(0,1.6))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('lag',
      side=1,
      cex=2.8,
      line=4)
mtext(expression(gamma),
      side=2,
      cex=2.8,
      line=4)
mtext('Variogram',
      side=3,
      cex=2.8,
      line=2.5)
points(binCenter3,
       vMean3,
       col="red",
       cex =2, 
       yaxt='n',
       xaxt='n')
lines(dGrid,
      (cov_parameter3[1])^2+ vTrue(dGrid,
                                   sigma2=cov_parameter3[2],
                                   range_val=cov_parameter3[3],
                                   nu=1),
      col='magenta',
      lwd=3,
      lty=2)
## QQ plot
stats::qqnorm((wN3-mean(wN3))/sd(wN3),
              yaxt='n',
              xaxt='n',
              main='',
              xlab='',
              ylab= '')
abline(0,1, lty=2, lwd=2, col='magenta')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
title_main = bquote(paste("Tile" == .(rand_tiles[3]),
                          ", ", 
                          "Ext Day" == .(rand_days[4])))
mtext('QQ plot',
      side=3,
      cex=2.8,
      line=2.5)
dev.off()

## -- First, we take the simple r-Pareto process in logarithmic scale --:
set.seed(123)
rand_days1 <- round(runif(n=3, min=1, max= 86))
cat('Rand day indexes: ', rand_days1, '\n')

set.seed(111)
rand_days2 <- round(runif(n=3, min=1, max= 86))
cat('Rand day indexes: ', rand_days2,'\n')

set.seed(222)
rand_days3 <- round(runif(n=3, min=1, max= 86))
cat('Rand day indexes: ', rand_days3,'\n')

set.seed(333)
rand_days4 <- round(runif(n=3, min=1, max= 86))
cat('Rand day indexes: ', rand_days4,'\n')

# Tile 1
# day 1
temp11 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[1]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[1]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[1]]]) [ ,rand_days1[1]]) 

image.plot(temp11[[1]],
           temp11[[2]],
           temp11[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[1]]])[ ,rand_days1[1]])
vario_temp1 <- vgram(expand.grid(gridPoints16x16Boxes[[rand_tiles[1]]]$lon, 
                                 gridPoints16x16Boxes[[rand_tiles[1]]]$lat),
                     (log_W_process_across_boxes[[rand_tiles[1]]]) [ ,rand_days1[1]],
                     N=10)

plot()

# day 2
temp12 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[1]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[1]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[1]]]) [ ,rand_days1[2]]) 

image.plot(temp12[[1]],
           temp12[[2]],
           temp12[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[1]]])[ ,rand_days1[2]])

# day 3
temp13 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[1]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[1]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[1]]]) [ ,rand_days1[3]]) 

image.plot(temp13[[1]],
           temp13[[2]],
           temp13[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[1]]])[ ,rand_days1[3]])

# Tile 2
# day 1
temp21 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[2]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[2]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[2] ]]) [ ,rand_days2[1]]) 

image.plot(temp21[[1]],
           temp21[[2]],
           temp21[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[2]]])[ ,rand_days2[1]])

# day 2
temp22 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[2]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[2]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[2]]]) [ ,rand_days2[2]]) 

image.plot(temp22[[1]],
           temp22[[2]],
           temp22[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[2]]])[ ,rand_days2[2]])

# day 3
temp23 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[2]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[2]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[2]]]) [ ,rand_days2[3]]) 

image.plot(temp23[[1]],
           temp23[[2]],
           temp23[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[2]]])[ ,rand_days2[3]])

# Tile 3
# day 1
temp31 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[3]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[3]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[3]]])[ ,rand_days3[1]]) 

image.plot(temp31[[1]],
           temp31[[2]],
           temp31[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[3]]])[ ,rand_days3[1]])

# day 2
temp32 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[3]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[3]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[3]]]) [ ,rand_days3[2]]) 

image.plot(temp32[[1]],
           temp32[[2]],
           temp32[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[3]]])[ ,rand_days3[2]])

# day 3
temp33 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[3]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[3]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[3]]]) [ ,rand_days3[3]]) 

image.plot(temp33[[1]],
           temp33[[2]],
           temp33[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[3]]])[ ,rand_days3[3]])


# Tile 4
# day 1
temp41 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[4]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[4]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[4]]]) [ ,rand_days4[1]]) 

image.plot(temp41[[1]],
           temp41[[2]],
           temp41[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[4]]])[ ,rand_days4[1]])

# day 2
temp42 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[4]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[4]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[4]]]) [ ,rand_days4[2]]) 

image.plot(temp42[[1]],
           temp42[[2]],
           temp42[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[4]]])[ ,rand_days4[2]])

# day 3
temp43 <- imagePlot(expand.grid(gridPoints16x16Boxes[[rand_tiles[4]]]$lon, 
                                gridPoints16x16Boxes[[rand_tiles[4]]]$lat),
                    (log_W_process_across_boxes[[rand_tiles[4]]]) [ ,rand_days4[3]]) 

image.plot(temp43[[1]],
           temp43[[2]],
           temp43[[3]],
           xlab='lon',
           ylab='lat',
           col=hcl.colors(12, "YlOrRd", rev=TRUE))
summary((log_W_process_across_boxes[[rand_tiles[4]]])[ ,rand_days4[3]])

## -- Creating Plots--:
png("~/Desktop/Fig-Tile1-GP-from-RSST.png",
    units="in",
    width=27,
    height=8,
    res=200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste('Ext day' == .(rand_days1[1])))
image.plot(temp11[[1]],
           temp11[[2]],
           temp11[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2.5)
title_side <- bquote(paste('Tile' == .(rand_tiles[1])))
mtext(title_side, side= 2,  outer=TRUE, cex=2.5)

title_expression <- bquote(paste('Ext day' == .(rand_days1[2])))
image.plot(temp12[[1]],
           temp12[[2]],
           temp12[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis =2.5)
title_expression <- bquote(paste('Ext day' == .(rand_days1[3])))
image.plot(temp13[[1]],
           temp13[[2]],
           temp13[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis =2.5)
dev.off()

png("~/Desktop/Fig-Tile2-GP-from-RSST.png",
    units="in",
    width=27,
    height=8,
    res=200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste('Ext day' == .(rand_days2[1])))
image.plot(temp21[[1]],
           temp21[[2]],
           temp21[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis =2.5)
title_side <- bquote(paste('Tile' == .(rand_tiles[2])))
mtext(title_side, side= 2,  outer=TRUE, cex=2.5)

title_expression <- bquote(paste('Ext day' == .(rand_days2[2])))
image.plot(temp22[[1]],
           temp22[[2]],
           temp22[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis =2.5)
title_expression <- bquote(paste('Ext day' == .(rand_days2[3])))
image.plot(temp23[[1]],
           temp23[[2]],
           temp23[[3]],
           cex.main=3,
           main=title_expression,
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2.5)
dev.off()

png("~/Desktop/Fig-Tile3-GP-from-RSST.png",
    units="in",
    width=27,
    height=8,
    res=200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste('Ext day' == .(rand_days3[1])))
image.plot(temp31[[1]],
           temp31[[2]],
           temp31[[3]],
           cex.main=3,
           main=title_expression,
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2.5)
title_side <- bquote(paste('Tile' == .(rand_tiles[3])))
mtext(title_side, side= 2,  outer=TRUE, cex=2.5)

title_expression <- bquote(paste('Ext day' == .(rand_days3[2])))
image.plot(temp32[[1]],
           temp32[[2]],
           temp32[[3]],
           cex.main=3,
           main=title_expression,
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2.5)
title_expression <- bquote(paste('Ext day' == .(rand_days3[3])))
image.plot(temp33[[1]],
           temp33[[2]],
           temp33[[3]],
           cex.main=3,
           main=title_expression,
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2.5)
dev.off()

png("~/Desktop/Fig-Tile4-GP-from-RSST.png",
    units="in",
    width=27,
    height=8,
    res=200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste('Ext day' == .(rand_days4[1])))
image.plot(temp41[[1]],
           temp41[[2]],
           temp41[[3]],
           cex.main=3,
           main=title_expression,
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2.5)
title_side <- bquote(paste('Tile' == .(rand_tiles[4])))
mtext(title_side, side= 2,  outer=TRUE, cex=2.5)

title_expression <- bquote(paste('Ext day' == .(rand_days4[2])))
image.plot(temp42[[1]],
           temp42[[2]],
           temp42[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis =2.5)
title_expression <- bquote(paste('Ext day' == .(rand_days4[3])))
image.plot(temp43[[1]],
           temp43[[2]],
           temp43[[3]],
           cex.main=3,
           main=title_expression,
           xlab='',
           ylab='',
           col=hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args=list(cex.axis=2.5),
           cex.axis=2.5)
dev.off()

## Variogram plot for randomly selected tiles:
tile1 <- rand_tiles[1]
loc_across_tile1 <- expand.grid(gridPoints16x16Boxes[[tile1]]$lon, 
                                gridPoints16x16Boxes[[tile1]]$lat)
dim(loc_across_tile1)


geoR::variog(loc_across_tile1, log_W_process_across_boxes[[tile1]], option=('bin'))



# # Calculate Empirical Variogram:
# library('mev')
# extremo(dat = rPareto_across_boxes[[rand_tile]][ ,rand_days[1]], margp = 0.96, 
#         coord = as.matrix(expand.grid(gridPoints16x16Boxes[[rand_tile]]$lon, 
#                             gridPoints16x16Boxes[[rand_tile]]$lat)))
# 
# variogram_tile21 <- vgram(expand.grid(gridPoints16x16Boxes[[rand_tile]]$lon, 
#                                       gridPoints16x16Boxes[[rand_tile]]$lat),
#                           log((rPareto_across_boxes[[rand_tile]])[ ,rand_days[1]]),
#                           N=15, lon.lat=TRUE)
# 
# boxplot(variogram_tile21$d, variogram_tile21$d)

# loc_temp <- gridPoints16x16Boxes[[2]]
# loc <- expand.grid(loc_temp$lon, loc_temp$lat)
# dim(loc)
# y <- as.matrix(obs_at_boxes_2d[[2]])
# dim(y)
# 
# # temp <- extremo(dat=y, margp=0.95, coord=loc)
# # extremo <- function(dat, margp, coord){
# #  #stopifnot(isTRUE(all(margp>=0, margp< 1, length(margp)==1, nrow(coord)==ncol(dat))))
# #  # Local quantile - threshold data
#   
# dat <- as.matrix(dat)
#   margthresh <- apply(dat, 2, quantile, margp, na.rm = TRUE)
#   dat <- t(t(dat)-margthresh)
#   
#   # Keep only instances where there is at least one exceedance
#   dimat <- rdist(coord)
#   excind <- apply(dat, 1, function(x){isTRUE(max(x, na.rm = TRUE) > 0)}) #avoid NA
#   
#   dat <- dat[excind,]
#   res <- matrix(0, ncol=4, nrow=choose(ncol(dat), 2))
#   b <- 0L
#   for(i in 1:(ncol(dat)-1)){
#     for(j in (i+1):ncol(dat)){
#       b <- b + 1L
#       subdat <- na.omit(dat[,c(i,j)])
#       if(length(subdat) > 0){
#         res[b, ] <- c(i,j, nrow(subdat), mean(I(subdat[,2] > 0) * I(subdat[,1] > 0))/(0.5*mean(I(subdat[,1] > 0)) + 0.5*mean(I(subdat[,2] > 0))) )
#       } else{
#         res[b, ] <- NA
#       }
#     }
#   }
#   res <- na.omit(res)
#   ellips <- list(...)
#   ellips$y <- res[,4]
#   ellips$x <- apply(res, 1, function(x){dimat[x[1], x[2]]})
#   return(invisible(cbind(dist = ellips$x, prob = ellips$y)))
# }

# ## Calculating extremogram for the r-Pareto process 
# extremogram <- function(loc, y, thres=0.96, id=NULL, d=NULL, lon.lat=FALSE, 
#                         dmax=NULL, N=NULL, breaks=NULL, prettyBins=FALSE){
#   
#   margthresh <- apply(y, 2, quantile, probs=thres, na.rm=TRUE)
#   dat <- t(t(y)- margthresh)
#   
#   loc_temp <- gridPoints16x16Boxes[[1]]
#   loc <- expand.grid(loc_temp$lon, loc_temp$lat)
#   dim(loc)
#   
#   # If nearest neighbor indices are missing, create all possible pairs.
#   if(is.null(id)){
#     n <- nrow(loc)
#     is = rep(1:n, n)
#     js = rep(1:n, rep(n, n))
#     ind <- is > js
#     id <- cbind(is, js)[ind, ]
#   }
#   
#   # If distances are missing, calculate them
#   if(is.null(d)){
#     loc <- as.matrix(loc)
#     if(lon.lat){
#       d <- rdist.earth.vec(loc[id[,1],], loc[id[,2],]) # Result in miles, not meters
#     }
#     else{
#       d <- rdist.vec(loc[id[,1],], loc[id[,2],])
#     }
#   }
#   
#   # Center the columns by their mean and get row means if y is a matrix
#   colMeans <- apply(y, 2, mean, na.rm=TRUE)
#   yCntr = sweep(y, 2, colMeans) 
#   y1Cntr = yCntr[id[,1],]
#   y2Cntr = yCntr[id[,2],]
#   
#   # Compute the extremogram
#   ext <- rowMeans(cbind(y1Cntr * y2Cntr), na.rm = TRUE)
#   
#   # Information for the returned object
#   call <- match.call()
#   if (is.null(dmax)) {
#     dmax <- max(d)
#   }
#   od <- order(d)
#   d <- d[od]
#   ext <- ext[od]
#   ind <- d <= dmax & !is.na(ext)
#   
#   # Add a binned extremogram if breaks are supplied
#   out <- list(d = d[ind], extremogram = ext[ind], call = call)
#   if (!is.null(breaks) | !is.null(N)) {
#     out <- c(out, stats.bin(d[ind], ext[ind], N = N, breaks = breaks,
#                             prettyBins = prettyBins))
#   }
#   class(out) <- c("extremogram", class(out))
#   out
# }

## --- Fitting Gaussian process to the tiles -- :
## Initialize the array to store results
set.seed(123)
t <- dim(rPareto_across_boxes[[1]])[2]
print(t)

cov_parameter_est <- array(NA, dim=c(35, t, 4)) # tau, sigma2, rho, EDF
dim(cov_parameter_est)

location_parameter_est <- matrix(NA, nrow=35, ncol=t)
print(dim(location_parameter_est))

# Loop through each i
tic()
for (i in 1:35) {
  cat("Box index:", i, "\n")
  loc <- expand.grid(gridPoints16x16Boxes[[i]]$lon,
                     gridPoints16x16Boxes[[i]]$lat)
  
  # T_star <- rPareto_across_boxes[[i]]
  T_star <- log_W_process_across_boxes[[i]]
  # Fitting Gaussian Processes for individual time points
  for (j in 1:t) {
    print(j)
    
    tryCatch({
      # Fitting Gaussian processes to the tiles
      fit_Gauss <- spatialProcess(loc,
                                  T_star[,j],
                                  na.rm=TRUE,
                                  smoothness = 0.5,
                                  mKrig.args=list(m=1))  # Constant part of the model
      
      cov_parameter_est[i, j, ] <- as.numeric(fit_Gauss$MLESummary[6:9])
      location_parameter_est[i, j] <- as.numeric(fit_Gauss$d)
      
    }, error = function(e) {
      cat("Error occurred for i=", i, ", j=", j, ": ", conditionMessage(e), "\n")
    })
  }
}
toc()

## Plot the fixed slope:
# randomly selected tiles
set.seed(123)
rand_tiles <- round(runif(n=13, min =0, max=length(rPareto_across_boxes)))
rand_tiles <- unique(rand_tiles)
print(rand_tiles)

set.panel(3,3)
par(mar=c(2, 1.5, 1.5, 3))
hist(location_parameter_est[rand_tiles[1],],
     probability=TRUE,
     #xlab= 'values',
     ylab='probability',
     main='Tile 11')
hist(location_parameter_est[rand_tiles[2],],
     probability=TRUE,
     #ylab='probability',
     main='Tile 28')
hist(location_parameter_est[rand_tiles[3],],
     probability=TRUE,
     #ylab='probability',
     main='Tile 15')
mtext('Historgram of the location parameter in the model',
      adj=1.5,
      line=3)
hist(location_parameter_est[rand_tiles[4],],
     probability=TRUE,
     ylab='probability',
     main='Tile 31')
hist(location_parameter_est[rand_tiles[5],],
     probability=TRUE, 
     main='Tile 33')
hist(location_parameter_est[rand_tiles[6],],
     probability=TRUE,
     main='Tile 3')
hist(location_parameter_est[rand_tiles[7],],
     probability=TRUE, 
     ylab='probability',
     xlab='values',
     main='Tile 19')
hist(location_parameter_est[rand_tiles[8],],
     probability=TRUE,
     xlab='values',
     main ='Tile 20')
hist(location_parameter_est[rand_tiles[9],],
     probability=TRUE,
     xlab='values',
     main='Tile 17')

## -- Plotting the parameter estimates --
# Randomly selected Tilew
set.panel(3,3)
par(mar=c(2, 1.5, 1.5, 3))
hist(cov_parameter_est[rand_tiles[1], ,1],
     probability = T,
     xlab='',
     ylab='',
     main='Tile 11')
hist(cov_parameter_est[rand_tiles[2], ,1],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 28')
mtext(expression(tau),
      line=2,
      cex=1.8)
hist(cov_parameter_est[rand_tiles[3], ,1],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 15')

hist(cov_parameter_est[rand_tiles[4], ,1],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 31')
hist(cov_parameter_est[rand_tiles[5], ,1],
     probability=TRUE, 
     xlab='',
     ylab='',
     main='Tile 33')
hist(cov_parameter_est[rand_tiles[6], , 1],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 3')
hist(cov_parameter_est[rand_tiles[7], , 1],
     probability=TRUE, 
     xlab='',
     ylab='',
     main='Tile 19')
hist(cov_parameter_est[rand_tiles[8], ,1],
     probability=TRUE,
     xlab='',
     ylab='',
     main ='Tile 20')
hist(cov_parameter_est[rand_tiles[9], ,1],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 17')


# parameter estimate: sigma2
set.panel(3,3)
par(mar=c(2, 1.5, 1.5, 3))
hist(cov_parameter_est[rand_tiles[1], ,2],
     probability = T,
     xlab='',
     ylab='',
     main='Tile 11')
hist(cov_parameter_est[rand_tiles[2], ,2],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 28')
mtext(expression(sigma^2),
      line=2,
      cex=1.8)
hist(cov_parameter_est[rand_tiles[3], ,2],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 15')

hist(cov_parameter_est[rand_tiles[4], ,2],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 31')
hist(cov_parameter_est[rand_tiles[5], ,2],
     probability=TRUE, 
     xlab='',
     ylab='',
     main='Tile 33')
hist(cov_parameter_est[rand_tiles[6], , 2],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 3')
hist(cov_parameter_est[rand_tiles[7], , 2],
     probability=TRUE, 
     xlab='',
     ylab='',
     main='Tile 19')
hist(cov_parameter_est[rand_tiles[8], ,2],
     probability=TRUE,
     xlab='',
     ylab='',
     main ='Tile 20')
hist(cov_parameter_est[rand_tiles[9], ,2],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 17')

# parameter estimate: Range parameter
set.panel(3,3)
par(mar=c(2, 1.5, 1.5, 3))
hist(cov_parameter_est[rand_tiles[1], ,3],
     probability = T,
     xlab='',
     ylab='',
     main='Tile 11')
hist(cov_parameter_est[rand_tiles[2], ,3],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 28')
mtext(expression(rho),
      line=2,
      cex=1.8)
hist(cov_parameter_est[rand_tiles[3], ,3],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 15')
hist(cov_parameter_est[rand_tiles[4], ,3],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 31')
hist(cov_parameter_est[rand_tiles[5], ,3],
     probability=TRUE, 
     xlab='',
     ylab='',
     main='Tile 33')
hist(cov_parameter_est[rand_tiles[6], , 3],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 3')
hist(cov_parameter_est[rand_tiles[7], , 3],
     probability=TRUE, 
     xlab='',
     ylab='',
     main='Tile 19')
hist(cov_parameter_est[rand_tiles[8], ,3],
     probability=TRUE,
     xlab='',
     ylab='',
     main ='Tile 20')
hist(cov_parameter_est[rand_tiles[9], ,3],
     probability=TRUE,
     xlab='',
     ylab='',
     main='Tile 17')

# Loop through each j
# for (j in 1:143) {
#   if (j %% 50 == 0) {
#     cat("j:", j, "\n")
#   }
#   tryCatch({
#     fit_Gauss <- spatialProcess(loc, log_W[, j])
#     fit_Gaussian[i, j, ] <- as.numeric(fit_Gauss$MLESummary)[6:8]
#   }, error = function(e) {
#     cat("Error occurred for i=", i, ", j=", j, ": ", conditionMessage(e), "\n")
#     next
#   })
#}


## Tiles 1:
r_Pareto_obs_1 <- rPareto_across_boxes[[rand_tiles[1]]][,1] # 256x86
loc_tile_1 <- expand.grid(gridPoints16x16Boxes[[ rand_tiles[1] ]]$lon, 
                          gridPoints16x16Boxes[[ rand_tiles[1] ]]$lat)

variogram_1 <- vgram(loc_tile_1, r_Pareto_obs_1)
set.panel(1,1)
plot(variogram_1$d, variogram_1$vgram, col = alpha('gray', 0.1))
bplot.xy(variogram_1$d, variogram_1$vgram)


chi96[i,j]<-chi96[j,i]<-mean(U[,i]>0.96 & U[,j]>0.96) / (1-0.96)







fit_Gaussian_part2 <- fit_Gaussian[14:35, ,]
fit_Gaussian_part1 <- fit_Gaussian[1:13, ,]
dim(fit_Gaussian_part1)

boxplot(fit_Gaussian[,,1], main = expression(tau^2))
boxplot(fit_Gaussian[,,2], main = expression(sigma^2))
boxplot(fit_Gaussian[,,3], main = expression(rho))
plot(threshold_across_boxes,
     main = 'Threshold across the 16x16 grids',
     ylab= 'u')

## Check for NA values in each matrix
no_na_indices <- sapply(rPareto_across_boxes, function(x) !any(is.na(x)))
no_na_rPareto <- rPareto_across_boxes[no_na_indices]
no_na_16x16GridBox <- gridPoints16x16Boxes[no_na_indices]

## Save the list to an RDS file
#saveRDS(no_na_rPareto, "~/Desktop/rParetoRRST.rds")


## Step 5: Simulate extreme process:  log(u) + G_i - G_bar - Var(G_i- G_bar)
# Parameter configuration:
# # n = n_train + n_valid
# n_train <- 8000
# n_valid <- 2000
# # n_test <- 10000
# n <- n_train + n_valid 
# print(n)
# 
range_vals_train <- 10^seq(log10(0.2), log10(4), length.out = 6)
smoothness_vals_train <- seq(0.7, 2, length.out = 6)
u_vals_train <- 10^seq(log10(2.8), log10(5), length.out = 6)
sigma2_vals_train <- 10^seq(log10(1), log10(8), length.out = 6)

Theta_vals <-expand.grid(range_vals_train, smoothness_vals_train, u_vals_train, sigma2_vals_train)
dim(Theta_vals)

# colnames(Theta_vals) <- c('rho', 'nu', 'sigma2', 'u')
# dim(Theta_vals) # 32768 3

# # Fixed training parameter configuration for the moment:
# Theta_train <- Theta_vals[1:n_train,]
# dim(Theta_train)
# save(Theta_train, file="~/Desktop/theta_train.RData")
# 
# Theta_valid <- Theta_vals[(n_train+1):n,]
# dim(Theta_valid)
# save(Theta_valid, file="~/Desktop/theta_valid.RData")
# 
# range_vals_test <- 10^runif(n_test, log10(1.1), log10(10))
# smoothness_vals_test <- runif(n_test, 0.5, 3)
# sigma2_vals_test <- 10^runif(n_test, log10(1), log10(10))
# threshold_vals_test <- 10^runif(n_test, log10(2.2), log10(4.2))
# Theta_vals_test <- cbind(range_vals_test, smoothness_vals_test, sigma2_vals_test, threshold_vals_test)
# dim(Theta_vals_test)
# save(Theta_vals_test, file="~/Desktop/theta_test.RData")
# 

# Validation set:
# Theta_valid <- parameters_grid(K=20)
# dim(Theta_valid) # 5832 3
#save(Theta_valid, file="~/Desktop/theta_valid_sim.RData")

# # Test set:
# Theta_test <- parameters_grid(K=10)
# dim(Theta_test) # 1000  3
#save(Theta_test, file="~/Desktop/theta_test_sim.RData")

# ## Save the 3D array observations in RData format:
# save(Theta_train, file="~/Desktop/Tentative Results-Proj2/theta_train_sim.RData")

## Function defined to calculate the Var(G-G_bar):
Var_G_G_bar <- function(sigma_mat){
  d <- nrow(sigma_mat)
  # Covariance matrix is symmetric, 
  # therefore, rowSumsMatrix = columnSumsMatrix
  rowSumsSigma <- rowSums(sigma_mat)
  rowSumsMatrix <- do.call(cbind, replicate(d, rowSumsSigma, simplify=FALSE))
  
  # Compute overall sum once
  totalSumSigma <- sum(sigma_mat)
  
  # Calculate the variance directly without creating additional matrices
  a <- diag(sigma_mat) - (2*diag(rowSumsMatrix/d)) + (totalSumSigma/(d^2))
  
  return(a)
}

simulate_extremes <- function(spatial_grid, Theta, m, seed)
{
  # set.seed(seed)
  ## Spatial coordinates defining the grid
  coord <- spatial_grid
  D_matrix <- rdist(coord) # 1024 1024
  
  # Initialize a list to store Y_r matrices
  Y_r_matrix <- array(NA, 
                      dim = c(nrow(Theta), m, ncol(D_matrix)))
  
  for(i in 1:nrow(Theta)){
    # if(i%%25==0){
    #   print(i)}
    print(i)
    
    range <- Theta[i,1]
    smooth <- Theta[i,2]
    sigma2 <- Theta[i,3]
    u <- Theta[i,4] # threshold of the process
    
    nugget <- 0.0
    cov_mat <- fields::Matern(D_matrix, range = range, nu = smooth) + diag(nugget, nrow = nrow(D_matrix))
    L <- chol(cov_mat)
    
    # variance of G-G_bar, in vector format:
    var_G_dev_mean <- Var_G_G_bar(cov_mat)
    
    for(j in 1:m){
      set.seed(seed + j)
      # Gaussian Process
      log_R <- rexp(1) + u # to account for the extremal intensity
      z <- t(L) %*% rnorm(nrow(coord))
      G <- (z - c(mean(z)) - 0.5*c(var_G_dev_mean))  # spatail variability
      log_Y_r <- log_R + G  
      Y_r_matrix[i, j, ] <- log_Y_r # converting the exponential marginals to standard Pareto
    }
  }
  return(Y_r_matrix)
}

## -- Simulating the r-Pareto process --:
# for grid 1000: a = -0.9, b = 1.05
# library('qmap')
# Y_r_process_RedSeaSST_afterQM <- array(NA, dim = c(9, 143, 6239))
# for(i in 1:6239)
# {
#   print(i)
#   for(j in 1:9)
#   {
#     print(j)
#     p_m <- Y_r_process_RedSeaSST[j, ,i]
#     p_o <- as.numeric(r_pareto_process_acrs_RS[,i])
#     
#     fitQM <- fitQmapPTF(p_o, p_m, transfun = 'linear')
#     qm3 <- doQmapPTF(p_m , fitQM)
#     Y_r_process_RedSeaSST_afterQM[j, , i] <-   qm3
#   }
# }
# rspace <- round(runif(n=10, min = 1, max =  6239))

## Generate r-Pareto process in log scale
# gridbox_16by16 <- expand.grid(1:16, 1:16)
range_vals_train <- 10^seq(log10(0.1), log10(5), length.out = 10)
smoothness_vals_train <- seq(0.5, 2, length.out = 10)
sigma2_vals_train <- 10^seq(log10(0.1), log10(5), length.out = 10)
# u_vals_train <- 10^seq(log10(2.8), log10(5), length.out = 6)
u_vals_train <- threshold_across_boxes[[2]]

Theta_vals <-expand.grid(range_vals_train, smoothness_vals_train, u_vals_train, sigma2_vals_train)
dim(Theta_vals)

tic()
Y_r_process_RedSeaSST <- simulate_extremes(expand.grid(gridPoints16x16Boxes[[2]]$lon, 
                                                       gridPoints16x16Boxes[[2]]$lat),
                                           Theta_vals,
                                           m=143,
                                           seed=123)
toc()


temp <- imagePlot(expand.grid(gridPoints16x16Boxes[[2]]$lon, 
                              gridPoints16x16Boxes[[2]]$lat),
                  (Y_r_process_RedSeaSST)[20,3,])
image.plot(temp[[1]],
           temp[[2]],
           temp[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
summary(c(Y_r_process_RedSeaSST))


library('qmap')
Y_r_process_RedSeaSST_afterQM <- array(NA, dim = c(1000, 143, 256))
est_intercept <- matrix(NA, nrow=216, ncol=256)
est_slope  <- matrix(NA, nrow=216, ncol=256)
for(i in 1:216)
{
  print(i)
  for(j in 1:256)
  {
    print(j)
    p_m <- Y_r_process_RedSeaSST[i, ,j]
    p_o <- (rPareto_across_boxes[[2]])[j,]
    
    # Check for NA values in p_o
    if (any(is.na(p_o))) {
      next  # Skip this iteration
    }
    
    fitQM <- fitQmapPTF(p_o, p_m, transfun = 'linear')
    qm3 <- doQmapPTF(p_m , fitQM)
    Y_r_process_RedSeaSST_afterQM[i, , j] <- qm3
    
    # fit lm to get the a and b:
    regression <- lm(qm3 ~ p_m)
    intercept <- coef(regression)[1]
    slope <- coef(regression)[2]
    
    est_slope[i, j] <- slope
    est_intercept[i, j] <- intercept
    
  }
}
hist(c(est_slope))
hist(c(est_intercept))



## General simulation 
gridbox_16by16 <- expand.grid(1:16, 1:16)
tic()
Y_r_process_RedSeaSST <- simulate_extremes(gridbox_16by16,
                                           Theta_vals,
                                           m=143,
                                           seed=123)
toc()
dim(Y_r_process_RedSeaSST)
summary(c(Y_r_process_RedSeaSST))
summary(c(rPareto_across_boxes[[4]], rPareto_across_boxes[[6]], rPareto_across_boxes[[8]]))

# library('qmap')
Y_r_process_RedSeaSST_afterQM <- array(NA, dim = c(1296, 143, 256))
est_slope <- array(NA, dim=c(1296, 14, 256)) 
est_intercept <- array(NA, dim=c(1296,14,256)) 

for(i in 1:1296) {
  print(i)
  for(j in 1:14) {
    print(j)
    for(k in 1:256) {
      p_m <- Y_r_process_RedSeaSST[i, ,k]
      p_o <- as.matrix(no_na_rPareto[[j]])[k,]
      
      # Check if p_o contains only NA values
      if (all(is.na(p_o))) {
        next  # Skip to the next iteration of k
      }
      
      fitQM <- parametric_quantile_mapping(p_o, p_m)
      Y_r_process_RedSeaSST_afterQM[i, ,k] <- fitQM$corrected_model
      est_slope[i, j, k] <- fitQM$slope
      est_intercept[i, j ,k] <- fitQM$intercept
    }
  }
}


hist(c(est_slope))
hist(c(est_intercept))

dim(Y_r_process_RedSeaSST_afterQM)
temp <- imagePlot(expand.grid(1:16,  1:16),
                  (Y_r_process_RedSeaSST[8,4,]))
image.plot(temp[[1]],
           temp[[2]],
           temp[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))

temp <- imagePlot(expand.grid(gridPoints16x16Boxes[[4]]$lon, 
                              gridPoints16x16Boxes[[4]]$lat),
                  (rPareto_across_boxes[[4]])[,2])
image.plot(temp[[1]],
           temp[[2]],
           temp[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))

parametric_quantile_mapping <- function(observed, model) {
  # Calculate quantiles
  obs_quantiles <- quantile(observed, probs = seq(0, 1, length.out = 100), na.rm=TRUE)
  model_quantiles <- quantile(model, probs = seq(0, 1, length.out = 100), na.rm=TRUE)
  
  # Fit linear regression
  regression <- lm(obs_quantiles ~ model_quantiles)
  intercept <- coef(regression)[1]
  slope <- coef(regression)[2]
  
  # Apply linear transformation
  corrected_model <- intercept + slope * model
  
  return(list(corrected_model = corrected_model, intercept = intercept, slope = slope))
}

temp <- imagePlot(gridbox_16by16,
                  Y_r_process_RedSeaSST[12^2,9,])
image.plot(temp[[1]], 
           temp[[2]],
           temp[[3]],
           col = hcl.colors(12, "YlOrRd", rev = TRUE))

# 1207.649
save(Y_r_process_RedSeaSST, file="~/Desktop/Train_fields_rep10.RData")
# dim(Y_r_process_RedSeaSST) # 27   10 6239


tic()
Y_r_process_RedSeaSST_valid <- simulate_extremes(gridbox_16by16,
                                                 Theta_valid,
                                                 m=10,
                                                 seed=123)
toc()
dim(Y_r_process_RedSeaSST_valid)

# 1207.649
save(Y_r_process_RedSeaSST_valid, file="~/Desktop/Valid_fields_rep10.RData")

tic()
Y_r_process_RedSeaSST_test <- simulate_extremes(gridbox_16by16,
                                                Theta_vals_test,
                                                m=10,
                                                seed=123)
toc()
dim(Y_r_process_RedSeaSST_test)

# 1207.649
save(Y_r_process_RedSeaSST_test, file="~/Desktop/Test_fields_rep10.RData")


# Latest-April 13, 2024
###  --- Not required ------------------
p01 <- quantile(Y_r_process_RedSeaSST[2,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_RedSeaSST[2,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_RedSeaSST[2,1,], p99))
imgPlt_temp1 <- imagePlot(spatial_loc, clipped_data)

summary(c(Y_r_process_RedSeaSST)) 
## Red Sea data exploration after selecting the extremes days: 
imgPlt_temp1 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,1,])
image.plot(imgPlt_temp1[[1]],
           imgPlt_temp1[[2]],
           imgPlt_temp1[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp2 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,2, ])
image.plot(imgPlt_temp2[[1]],
           imgPlt_temp2[[2]],
           imgPlt_temp2[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp3 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,3,])
image.plot(imgPlt_temp3[[1]],
           imgPlt_temp3[[2]],
           imgPlt_temp3[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp4 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,3,])
image.plot(imgPlt_temp4[[1]],
           imgPlt_temp4[[2]],
           imgPlt_temp4[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp5 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,5,])
image.plot(imgPlt_temp5[[1]],
           imgPlt_temp5[[2]],
           imgPlt_temp5[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp6 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,6])
image.plot(imgPlt_temp6[[1]],
           imgPlt_temp6[[2]],
           imgPlt_temp6[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp7 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,7,])
image.plot(imgPlt_temp7[[1]],
           imgPlt_temp7[[2]],
           imgPlt_temp7[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp8 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,8,])
image.plot(imgPlt_temp8[[1]],
           imgPlt_temp8[[2]],
           imgPlt_temp8[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp9 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[1,9,])
image.plot(imgPlt_temp9[[1]],
           imgPlt_temp9[[2]],
           imgPlt_temp9[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE)) 
imgPlt_temp9 <- imagePlot(spatial_loc, Y_r_process_RedSeaSST[9,10,])
image.plot(imgPlt_temp9[[1]],
           imgPlt_temp9[[2]],
           imgPlt_temp9[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE)) 
# Plots for storage:
png("~/Desktop/Tentative Results-Proj2/Fig-sim-RedSeaSST-Panel1-nu-1.5.png",
    units = "in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[1,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[1,3])))
image.plot(imgPlt_temp1[[1]],
           imgPlt_temp1[[2]],
           imgPlt_temp1[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[2,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[2,3])))
image.plot(imgPlt_temp2[[1]],
           imgPlt_temp2[[2]],
           imgPlt_temp2[[3]],
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[3,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[3,3])))
image.plot(imgPlt_temp3[[1]],
           imgPlt_temp3[[2]],
           imgPlt_temp3[[3]],
           main = title_expression, 
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_main <- bquote(paste(nu == .(1.5)))
mtext(title_main, line = -0.5, outer=TRUE, cex=2.5)
dev.off()
png("~/Desktop/Tentative Results-Proj2/Fig-sim-RedSeaSST-Panel2-nu-1.5.png",
    units = "in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[4,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[4,3])))
image.plot(imgPlt_temp4[[1]],
           imgPlt_temp4[[2]],
           imgPlt_temp4[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)

title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[5,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[5,3])))
image.plot(imgPlt_temp5[[1]],
           imgPlt_temp5[[2]],
           imgPlt_temp5[[3]],
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[6,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[6,3])))
image.plot(imgPlt_temp6[[1]],
           imgPlt_temp6[[2]],
           imgPlt_temp6[[3]],
           main = title_expression, 
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
dev.off()
png("~/Desktop/Tentative Results-Proj2/Fig-sim-RedSeaSST-Panel3-nu-1.5.png",
    units = "in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[7,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[7,3])))
image.plot(imgPlt_temp7[[1]],
           imgPlt_temp7[[2]],
           imgPlt_temp7[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)

title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[8,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[8,3])))
image.plot(imgPlt_temp8[[1]],
           imgPlt_temp8[[2]],
           imgPlt_temp8[[3]],
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1.5[9,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1.5[9,3])))
image.plot(imgPlt_temp9[[1]],
           imgPlt_temp9[[2]],
           imgPlt_temp9[[3]],
           main = title_expression, 
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
dev.off()

#  -- Extracting parameter for smoothness 1 --:
index_nu_1 <- which(parameter_vals[,2]==smoothness_vals[1])

parameter_vals_nu_1 <- parameter_vals[index_nu_1,]
dim(parameter_vals_nu_1)

r_Pareto_process_nu_1 <- Y_r_process_RedSeaSST[parameter_index_nu_1, ,]
dim(r_Pareto_process_nu_1)

## Red Sea data exploration after selecting the extremes days: 
imgPlt_temp1 <- imagePlot(loc, r_Pareto_process_nu_1[1,8,])
image.plot(imgPlt_temp1[[1]],
           imgPlt_temp1[[2]],
           imgPlt_temp1[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp2 <- imagePlot(loc, r_Pareto_process_nu_1[2,7,])
image.plot(imgPlt_temp2[[1]],
           imgPlt_temp2[[2]],
           imgPlt_temp2[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp3 <- imagePlot(loc, r_Pareto_process_nu_1[3,1,])
image.plot(imgPlt_temp3[[1]],
           imgPlt_temp3[[2]],
           imgPlt_temp3[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp4 <- imagePlot(loc, r_Pareto_process_nu_1[4,5,])
image.plot(imgPlt_temp4[[1]],
           imgPlt_temp4[[2]],
           imgPlt_temp4[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp5 <- imagePlot(loc, r_Pareto_process_nu_1[5,6,])
image.plot(imgPlt_temp5[[1]],
           imgPlt_temp5[[2]],
           imgPlt_temp5[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp6 <- imagePlot(loc, r_Pareto_process_nu_1[6,1,])
image.plot(imgPlt_temp6[[1]],
           imgPlt_temp6[[2]],
           imgPlt_temp6[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp7 <- imagePlot(loc, r_Pareto_process_nu_1[7,7,])
image.plot(imgPlt_temp7[[1]],
           imgPlt_temp7[[2]],
           imgPlt_temp7[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp8 <- imagePlot(loc, r_Pareto_process_nu_1[8,8,])
image.plot(imgPlt_temp8[[1]],
           imgPlt_temp8[[2]],
           imgPlt_temp8[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp9 <- imagePlot(loc, r_Pareto_process_nu_1[9,9,])
image.plot(imgPlt_temp9[[1]],
           imgPlt_temp9[[2]],
           imgPlt_temp9[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE)) 


# Plots for storage:
png("~/Desktop/Tentative Results-Proj2/Fig-sim-RedSeaSST-Panel1-nu-1.png",
    units = "in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[1,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[1,3])))
image.plot(imgPlt_temp1[[1]],
           imgPlt_temp1[[2]],
           imgPlt_temp1[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)

title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[2,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[2,3])))
image.plot(imgPlt_temp2[[1]],
           imgPlt_temp2[[2]],
           imgPlt_temp2[[3]],
           # main = as.expression(title_expression),
           # main = paste('Day', randLoc[2]),
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[3,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[3,3])))
image.plot(imgPlt_temp3[[1]],
           imgPlt_temp3[[2]],
           imgPlt_temp3[[3]],
           # main = paste('Day', randLoc[3]),
           main = title_expression, 
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_main <- bquote(paste(nu == .(1)))
mtext(title_main, line = -0.5, outer=TRUE, cex=2.5)
dev.off()

# Plots for storage:
png("~/Desktop/Tentative Results-Proj2/Fig-sim-RedSeaSST-Panel2-nu-1.png",
    units = "in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[4,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[4,3])))
image.plot(imgPlt_temp4[[1]],
           imgPlt_temp4[[2]],
           imgPlt_temp4[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)

title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[5,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[5,3])))
image.plot(imgPlt_temp5[[1]],
           imgPlt_temp5[[2]],
           imgPlt_temp5[[3]],
           # main = as.expression(title_expression),
           # main = paste('Day', randLoc[2]),
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[6,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[6,3])))
image.plot(imgPlt_temp6[[1]],
           imgPlt_temp6[[2]],
           imgPlt_temp6[[3]],
           # main = paste('Day', randLoc[3]),
           main = title_expression, 
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
# title_main <- bquote(paste(nu == .(1.5)))
# mtext(title_main, line = -0.5, outer=TRUE, cex=2.5)
dev.off()


# Plots for storage:
png("~/Desktop/Tentative Results-Proj2/Fig-sim-RedSeaSST-Panel3-nu-1.png",
    units = "in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[7,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[7,3])))
image.plot(imgPlt_temp7[[1]],
           imgPlt_temp7[[2]],
           imgPlt_temp7[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)

title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[8,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[8,3])))
image.plot(imgPlt_temp8[[1]],
           imgPlt_temp8[[2]],
           imgPlt_temp8[[3]],
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(parameter_vals_nu_1[9,1]),
                                 ", ", 
                                 tau^2 == .(parameter_vals_nu_1[9,3])))
image.plot(imgPlt_temp9[[1]],
           imgPlt_temp9[[2]],
           imgPlt_temp9[[3]],
           main = title_expression, 
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
dev.off()

## -- SİMULATE EXTREME SAMPLES FOR THE TRAINING-VALID-TEST SET --: 
# Spatial extreme process: X_star - threshold = model
# u threshold - computed on the standard Pareto scale
# function simulate_extremes() returns the simulated extreme process
# with threshold u in standard Pareto margins
coord <- expand.grid(1:32, 1:32)
index_tau2_0.1 <- which(Theta_train[,3]==0.1) # 1024
index_tau2_0.5 <- which(Theta_train[,3]==unique(Theta_train[,3])[15]) # 1024
index_tau2_0.7 <- which(Theta_train[,3]==unique(Theta_train[,3])[21]) # 1024
working_index <- c(1, 14337, 20481, 2, 14338, 20482, 3, 14339, 20483)
Theta_train_work <- Theta_train[working_index,]

tic()
Y_r_process_train <- simulate_extremes(coord,
                                       # Theta_train,
                                       Theta_train_work,
                                       m=10,
                                       seed=123,
                                       u=log(threshold))
toc()

## Red Sea data exploration after selecting the extremes days: 
# imgPlt_temp1 <- imagePlot(loc, exp(X_star_extremes[randLoc[1],]))
# Calculate the 0.01 and 0.99 percentiles
p01 <- quantile(Y_r_process_train[1,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[1,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[1,1,], p99))
imgPlt_temp1 <- imagePlot(coord, clipped_data)

image.plot(imgPlt_temp1[[1]],
           imgPlt_temp1[[2]],
           imgPlt_temp1[[3]],
           main = paste('Day', randLoc[1]),
           col=hcl.colors(12, "YlOrRd", rev = TRUE))

p01 <- quantile(Y_r_process_train[2,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[2,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[2,1,], p99))
# imgPlt_temp2 <- imagePlot(loc, exp(X_star_extremes[randLoc[2],]))
imgPlt_temp2 <- imagePlot(coord, clipped_data)
image.plot(imgPlt_temp2[[1]],
           imgPlt_temp2[[2]],
           imgPlt_temp2[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))


p01 <- quantile(Y_r_process_train[3,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[3,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[3,1,], p99))
# imgPlt_temp3 <- imagePlot(loc, exp(X_star_extremes[randLoc[3],]))
imgPlt_temp3 <- imagePlot(coord, clipped_data)
image.plot(imgPlt_temp3[[1]],
           imgPlt_temp3[[2]],
           imgPlt_temp3[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))


p01 <- quantile(Y_r_process_train[4,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[4,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[4,1,], p99))
imgPlt_temp4 <- imagePlot(coord, clipped_data)
# imgPlt_temp4 <- imagePlot(loc, exp(X_star_extremes[randLoc[4],]))
image.plot(imgPlt_temp4[[1]],
           imgPlt_temp4[[2]],
           imgPlt_temp4[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))


p01 <- quantile(Y_r_process_train[5,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[5,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[4,1,], p99))
imgPlt_temp5 <- imagePlot(coord, clipped_data)
# imgPlt_temp5 <- imagePlot(loc, exp(X_star_extremes[randLoc[5],]))
image.plot(imgPlt_temp5[[1]],
           imgPlt_temp5[[2]],
           imgPlt_temp5[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))


p01 <- quantile(Y_r_process_train[6,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[6,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[6,1,], p99))
# imgPlt_temp6 <- imagePlot(loc, exp(X_star_extremes[randLoc[6],]))
imgPlt_temp6 <- imagePlot(coord, clipped_data)
image.plot(imgPlt_temp6[[1]],
           imgPlt_temp6[[2]],
           imgPlt_temp6[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))


p01 <- quantile(Y_r_process_train[7,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[7,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[7,1,], p99))
imgPlt_temp7 <- imagePlot(coord, clipped_data)
# imgPlt_temp7 <- imagePlot(loc, exp(X_star_extremes[randLoc[7],]))
image.plot(imgPlt_temp7[[1]],
           imgPlt_temp7[[2]],
           imgPlt_temp7[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))


p01 <- quantile(Y_r_process_train[8,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[8,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[8,1,], p99))
imgPlt_temp8 <- imagePlot(coord, clipped_data)
#imgPlt_temp8 <- imagePlot(loc, exp(X_star_extremes[randLoc[8],]))
image.plot(imgPlt_temp8[[1]],
           imgPlt_temp8[[2]],
           imgPlt_temp8[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))


p01 <- quantile(Y_r_process_train[9,1,],
                probs=0.01,
                na.rm=TRUE)
p99 <- quantile(Y_r_process_train[9,1,],
                probs=0.99,
                na.rm=TRUE)
clipped_data <- pmax(p01, pmin(Y_r_process_train[9,1,], p99))
imgPlt_temp9 <- imagePlot(coord, clipped_data)
# imgPlt_temp9 <- imagePlot(loc, exp(X_star_extremes[randLoc[9],]))
image.plot(imgPlt_temp9[[1]],
           imgPlt_temp9[[2]],
           imgPlt_temp9[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE)) 

# -- Plots to plot the r-Pareto process extracted from the Red Sea SST --:
png("~/Desktop/Tentative Results-Proj2/Fig-sim32x32-Panel1-clip.png",
    units = "in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste(theta == .(Theta_train_work[1,1] ),
                                 ", ", 
                                 nu == .(Theta_train_work[1,2] )))
image.plot(imgPlt_temp1[[1]],
           imgPlt_temp1[[2]],
           imgPlt_temp1[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(Theta_train_work[2,1] ),
                                 ", ", 
                                 nu == .(Theta_train_work[2,2] )))
image.plot(imgPlt_temp2[[1]],
           imgPlt_temp2[[2]],
           imgPlt_temp2[[3]],
           # main = as.expression(title_expression),
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(Theta_train_work[3,1] ),
                                 ", ", 
                                 nu == .(Theta_train_work[3,2] )))
image.plot(imgPlt_temp3[[1]],
           imgPlt_temp3[[2]],
           imgPlt_temp3[[3]],
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(tau^2 == .(0.1)))
mtext(title_expression, line = -0.9, outer=TRUE, cex=2.5)
dev.off()

# -- Plots to plot the r-Pareto process extracted from the Red Sea SST --:
png("~/Desktop/Tentative Results-Proj2/Fig-sim32x32-Panel2-clip.png",
    units ="in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
title_expression <- bquote(paste(theta == .(Theta_train_work[4,1] ),
                                 ", ", 
                                 nu == .(Theta_train_work[4,2] )))
image.plot(imgPlt_temp4[[1]],
           imgPlt_temp4[[2]],
           imgPlt_temp4[[3]],
           cex.main = 3,
           main = title_expression,
           xlab = '',
           ylab = '',
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(Theta_train_work[5,1] ),
                                 ", ", 
                                 nu == .(Theta_train_work[5,2] )))
image.plot(imgPlt_temp5[[1]],
           imgPlt_temp5[[2]],
           imgPlt_temp5[[3]],
           # main = as.expression(title_expression),
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(theta == .(Theta_train_work[6,1] ),
                                 ", ", 
                                 nu == .(Theta_train_work[6,2] )))
image.plot(imgPlt_temp6[[1]],
           imgPlt_temp6[[2]],
           imgPlt_temp6[[3]],
           main = title_expression,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
title_expression <- bquote(paste(tau^2 == .(0.5)))
mtext(title_expression, line = -0.9, outer=TRUE, cex=2.5)
dev.off()

png("~/Desktop/Tentative Results-Proj2/Fig-ExtremeObsSST-expMargins-Panel3-clip.png",
    units = "in",
    width = 27,
    height = 8,
    res = 200)
par(mfrow = c(1,3),
    mar = c(4, 4, 6, 5),
    oma = c(2, 4, 2.5, 6),
    cex.main = 1.5)
# par(mfrow = c(1,3),
#     mar = c(4, 4, 5, 5),
#     oma = c(2, 4, 0.4, 6),
#     cex.main = 1.5)
image.plot(imgPlt_temp7[[1]],
           imgPlt_temp7[[2]],
           imgPlt_temp7[[3]],
           main = paste('Day', randLoc[7]),
           cex.main = 3,
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           axis.args = list(cex.axis=2.5),
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           cex.axis = 2)
image.plot(imgPlt_temp8[[1]],
           imgPlt_temp8[[2]],
           imgPlt_temp8[[3]],
           main = paste('Day', randLoc[8]),
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           axis.args = list(cex.axis=2.5),
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           cex.axis = 2)
image.plot(imgPlt_temp9[[1]],
           imgPlt_temp9[[2]],
           imgPlt_temp9[[3]],
           main = paste('Day', randLoc[9]),
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           cex.main = 3,
           col = hcl.colors(12, "YlOrRd", rev = TRUE),
           axis.args = list(cex.axis=2.5),
           cex.axis = 2)
dev.off()

# 34.833 sec elapsed
# dim of Y_r process:  216   10 1024
## Save the 3D array observations in RData format:
save(Y_r_process,
     file="~/Desktop/Tentative Results-Proj2/simulateTrain.RData")

# Load the simulated Red SST:
load("~/Desktop/Tentative Results-Proj2/simulateTrain.RData")
# Step 2: List the objects to see what was loaded
loaded_objects <- ls()
print(loaded_objects)



# To generate the test set:
tic()
Y_r_process_test <- simulate_extremes(coord,
                                      Theta_test,
                                      m=10,
                                      seed=123,
                                      u=log(threshold))
toc()
# timing:
save(Y_r_process_test, file="~/Desktop/Tentative Results-Proj2/simulateTest.RData")

tic()
Y_r_process_valid <- simulate_extremes(coord,
                                       Theta_valid,
                                       m=10,
                                       seed=123,
                                       u=log(threshold))
toc()
# timing: 1705.072 secs
# save(Y_r_process_valid, file="~/Desktop/Tentative Results-Proj2/simulateValid.RData")

imgPlt_temp <- imagePlot(coord, (Y_r_process[1,1,]))
image.plot(imgPlt_temp[[1]],
           imgPlt_temp[[2]],
           imgPlt_temp[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp <- imagePlot(coord, (Y_r_process[1,2,]))
image.plot(imgPlt_temp[[1]],
           imgPlt_temp[[2]],
           imgPlt_temp[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp <- imagePlot(coord, (Y_r_process[1,4,]))
image.plot(imgPlt_temp[[1]],
           imgPlt_temp[[2]],
           imgPlt_temp[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))
imgPlt_temp <- imagePlot(coord, (Y_r_process[1,5,]))
image.plot(imgPlt_temp[[1]],
           imgPlt_temp[[2]],
           imgPlt_temp[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))

imgPlt_temp <- imagePlot(coord, (Y_r_process[1,10,]))
image.plot(imgPlt_temp[[1]],
           imgPlt_temp[[2]],
           imgPlt_temp[[3]],
           col=hcl.colors(12, "YlOrRd", rev = TRUE))

# Step 4 & 5: Model fittingl.  and simulation (assuming X_star is a matrix)
## Model: 
# log(X_star) = threshold  +  G_i- average(G) - Var(G_i-average(G) 
# where G_i ~ Gaussian distribution with stationary increaments with mean zero and
# and variance V(G_i- average(G))
coord <- expand.grid(1:32, 1:32)

# Step 6: Backtransform (assuming a function `reverse_transform` exists)
# we need to transform the observation from simulation, such that marginals are 
# coming from the standard Pareto:
simulated_original_scale <- lapply(simulated_extremes, reverse_transform)


## Declustering: Not required at the moment
# ## Finding the u-threshold :
# ## at least 5 observation every year is captured
# ## check for spatial de-clustering and threshold selection
# t=0
# store_quantiles <- rep(NA)
# for(i in 1: 31)
# {
#   print(i)
#   
#   store_quantiles[i] <- quantile(risk_t[(1:92)+t], 0.90)
#   print(quantile(risk_t[(1:92)+t], 0.90))
#   
#   u_temp <- store_quantiles[i]
#   
#   print(which(risk_t[(1:92)+t]>u_temp))
#   t <- t+ 92
# }
# ## Check the threshold value:
# u_threshold <- quantile(store_quantiles, 0.90)
# print(u_threshold)

## -- Creating 3D - array  of the spatial process with spatial locations:
array_r_Y <- array(NA, dim = c(num_lon, num_lat, num_time_points))
dim(array_r_Y)

lon_vals <- unique(loc[,1])
lat_vals <- unique(loc[,2])

## Loop through each location
for (d in 1:num_locations) {
  if(d %% 625 == 0){
    print(d)
  }
  # Find the indices of x and y for the current location
  ix <- match(loc[d, 1], x)
  iy <- match(loc[d, 2], y)
  
  # Assign the values from log_Y_r to the corresponding position in array_r_Y
  array_r_Y[ix, iy, ] <- log_Y_r[ , d]
}
dim(array_r_Y) # 113  89  97




# x <- unique(window[,1])
# y <- unique(window[,2])
# 
# ix <- match(window[,1], x)
# iy <- match(window[,2], y)
# 
# n1 <- max(ix)
# print(n1)
# 
# n2 <- max(iy)
# print(n2)
# 
# mat <- matrix(NA, nrow = n1, ncol = n2)
# dim(mat)
# 
# mat[cbind(ix, iy)] <- spatial_obs[1,]
# image.plot(mat,col = hcl.colors(12, "YlOrRd", rev = TRUE))

## Generate r-Pareto process: 
simulate_r_Pareto <- function(theta, alpha=5, u_threshold){
  sig2 <- 1 # marginal variance
  nu <- 1 # smoothness parameter
  cov_matrix <- sig2 * (fields::Matern(D_matrix, range=theta, nu=nu))
  # vTrue <- 1- sort(unique(cov_matrix))
  L <- chol(cov_matrix)
  
  log_R <- rexp(1)
  condition_met <- FALSE
  Y_r <- 0
  k<- 1
  while (!condition_met)
  {
    print(k)
    ## Generate the GP realization
    z <- as.numeric(L%*% rnorm(16))
    log_W <- (z - mean(z))
    log_Y_r <- log_W + log_R
    Y_r <- exp(log_Y_r)
    Y_r <- (Y_r-min(Y_r))/(max(Y_r)-min(Y_r))
    
    # Check the condition
    condition_met <- mean(Y_r) > 1
    R <- geometricMean(Y_r)
    log_R <- log(R)
    k <- k+1
  }
}

## Defining parameter configuration:
## for training and validation set compatible for Red Sea data
n <- 10
# n_train <- 2000
# n_valid <- 2000
range_vals <- seq(2, 10, length.out= n)
# tau_vals <- exp(runif(n, log(0.0001), log(1))) # nugget effect
#tau_vals <-  exp(seq(2, 10, length.out= n)


