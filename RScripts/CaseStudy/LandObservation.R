rm(list=ls())

library(tictoc)
library(fields)
library(sp)
library(sf)
library(maps)
library(maptools)
library(terra)
library(rnaturalearth)
library(ggplot2)
library(ncdf4)

# remotes::install_version("maptools", version = "1.1-7", repos = "http://cran.us.r-project.org")
# remotes::install_version("rgdal", version = "1.1-1", repos = "http://cran.us.r-project.org")

source("~/Desktop/GEV-SAR/RScripts/CaseStudy/RCM-dataLoad.R")

# Extract Lambert Conformal Projection parameters
proj <- ncatt_get(nc_file, "Lambert_Conformal")

# Print extracted projection attributes
print(proj)

# Extract specific attributes
standard_parallels <- proj$standard_parallel
central_longitude <- proj$longitude_of_central_meridian
latitude_origin <- proj$latitude_of_projection_origin

## Information:
# precip unit: mm/day
# resolution: 22 degree grids/ 25 km projected grid

# ## -- Annual-Maximum -- :
# # Extract the year from the time variable
# library(lubridate)
# years <- year(time_dates)
# 
# # Initialize an array to store annual max precipitation (lat, lon, years)
# annual_maxima <- array(NA,
#                        dim=c(dim(precip_data)[1],
#                              dim(precip_data)[2],
#                              length(unique(years))))
# 
# # Loop through each year and compute annual max at each grid point
# for (i in seq_along(unique(years))) {
#   print(i)
#   year_mask <- (years == unique(years)[i])  # Mask for the current year
#   annual_maxima[,,i] <- apply(precip_data[, ,year_mask],
#                               c(1,2),
#                               max,
#                               na.rm = TRUE)
# }
# # print dimensions of the computed annual maxima
# print(dim(annual_maxima))  # Should be (lat, lon, years)
# save(annual_maxima, file='AnnualMaxima-acrossNA.RData')

## Save the data
setwd("~/Desktop/GEV-SAR/Data")
load('AnnualMaxima-acrossNA.RData')
dim(annual_maxima) # 297 281  31

############ -- Cutting-off the ocean part -- ##################################
dim(lon_data) # 297 281 # orginal domain

# New spatial domain:
lon_data_new <- lon_data  # Copy the original matrix
lat_data_new <- lat_data  # Copy the original matrix

# Create logical masks for longitude and latitude conditions
lon_mask <- (lon_data >= -128 & lon_data <= -64)
lat_mask <- (lat_data >= 16.72 & lat_data <= 72)

# Apply the filter while keeping the matrix shape
lon_data_new[!(lon_mask & lat_mask)] <- NA  # Only keep values satisfying both conditions
lat_data_new[!(lon_mask & lat_mask)] <- NA  # Keep the matrix shape, mask out values

annual_maxima_new <- annual_maxima
annual_maxima_new[!(lon_mask & lat_mask)] <- NA
################################################################################

#################  -- Spatial Domain -- ########################################
set.panel(1,1)
par(oma=c(0,0,0,1))
par(mar=c(3,3,2,2))
plot(lon_data,
     lat_data,
     xlab='Lon',
     ylab='Lat',
     xlim=c(min(lon_data), max(lon_data)),
     ylim=c(min(lat_data), max(lat_data)),
     main='NA-Cordex',
     col='grey')
world <- ne_countries(returnclass = "sv")
terra::plot(world, add=T)
# Add the reprojected world map
#plot(st_geometry(world_proj), col = 'magenta', add = TRUE)
world(add=TRUE, col='magenta')

set.panel(1,1)
plot(lon_data_new,
     lat_data_new,
     xlab='Lon',
     ylab='Lat',
     xlim=c(min(lon_data), max(lon_data)),
     ylim=c(min(lat_data), max(lat_data)),
     main='Focused range',
     col='grey')
world(add=TRUE, col='magenta')

######### Color Palette ########################################################
# jet.colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
## Generate the desired number of colors from this palette
#library(RColorBrewer)
# nbcol <- 50
# color <- jet.colors(nbcol)
#ramp <- colorRamp(c(alpha("blue", 0.4), "white", "red"))

color_palette <- colorRampPalette(c("darkblue", "lightblue", "white", "orange", "red", "darkred"))
################################################################################



#################### -- Data Exploration -- ####################################
png("RCM-annualMax.png",
    units="in", 
    width=24,
    height=8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(5,4,5,8) + 0.3,
    oma=c(1.4,1,2,9))
set.panel(1,3)
bubblePlot(lon_data_new,
           lat_data_new  , 
           (annual_maxima_new[,,1]-median(annual_maxima_new[,,1], na.rm=TRUE))/sd(annual_maxima_new[,,1], na.rm=TRUE),
           col=color_palette,
           zlim=c(-1.2, 12),
           size=1,
           #na.rm=TRUE,
           yaxt='n',
           xaxt='n',
           xlab='',
           ylab= '',
           legend.width=4,
           cex.main=7,
           axis.args=list(cex.axis=2.5,lwd=2))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('Lon',
      side=1,
      cex=1.8,
      line=4.5)
mtext('Lat',
      side=2,
      cex=1.8,
      line=3.5)
world(add=TRUE, col='grey5', lwd=5)

bubblePlot(lon_data_new,
           lat_data_new  , 
           (annual_maxima_new[,,10]-median(annual_maxima_new[,,10], na.rm=TRUE))/sd(annual_maxima_new[,,10], na.rm=TRUE),
           col=color_palette,
           zlim=c(-1.2, 12),
           size=1,
           yaxt='n',
           xaxt='n',
           xlab='',
           ylab= '',
           legend.width=4,
           cex.main=7,
           axis.args=list(cex.axis=2.5,lwd=2))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('Lon',
      side=1,
      cex=1.8,
      line=4.5)
world(add=TRUE, col='grey5', lwd=5)

bubblePlot(lon_data_new,
           lat_data_new  , 
           (annual_maxima_new[,,28]-median(annual_maxima_new[,,28], na.rm=TRUE))/sd(annual_maxima_new[,,28], na.rm=TRUE),
           col=color_palette,
           zlim=c(-1.2, 12),
           size=1,
           yaxt='n',
           xaxt='n',
           xlab='',
           ylab= '',
           legend.width=4,
           cex.main=7,
           axis.args=list(cex.axis=2.5,lwd=2))
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)


mtext('Lon',
      side=1,
      cex=1.8,
      line=4.5)
world(add=TRUE, col='grey5', lwd=5)
dev.off()


## -- NO OVERLAPPING - DEFINING GRID 16X16 GRID ACROSS THE SPATIAL DOMAIN --
annualMaxima16x16Boxes <- list()
gridPoints16x16Boxes <- list()
point_count_list <- list()

grid_size <- 16
step <- 16

n_lon <- dim(lon_data_new)[1]  # 297 longitude points
n_lat <- dim(lat_data_new)[2]  # 281 latitude points

# Loop with overlapping steps
for (i in seq(1, n_lon - grid_size + 1, by = step)) {
  for (j in seq(1, n_lat - grid_size + 1, by = step)) {
    
    # Extract the 16x16 grid
    lon_grid <- lon_data_new[i:(i + grid_size - 1), j:(j + grid_size - 1)]
    lat_grid <- lat_data_new[i:(i + grid_size - 1), j:(j + grid_size - 1)]
    annual_maxima_grid <- annual_maxima_new[i:(i + grid_size - 1), j:(j + grid_size - 1), ]
    
    # Skip grids with NA values
    if (sum(is.na(lon_grid)) > 0 || sum(is.na(lat_grid)) > 0) {
      next
    }
    
    # Calculate grid center coordinates
    center_lon <- mean(lon_grid, na.rm = TRUE)
    center_lat <- mean(lat_grid, na.rm = TRUE)
    
    # Store results (including center)
    point_count_list[[length(point_count_list) + 1]] <- sum(!is.na(lon_grid))
    annualMaxima16x16Boxes[[length(annualMaxima16x16Boxes) + 1]] <- annual_maxima_grid
    gridPoints16x16Boxes[[length(gridPoints16x16Boxes) + 1]] <- list(
      lon = lon_grid,
      lat = lat_grid,
      center = c(center_lon, center_lat)  # Track center of the grid
    )
    
    # Optional: Plot the grid and center
    plot(lon_grid, lat_grid, xlab = 'Lon', ylab = 'Lat', main = paste('Grid:', i, j))
    points(center_lon, center_lat, col = "red", pch = 19)  # Mark center
  }
}
# 657 overlapping tiles

## convert the list into array:
RCMannualMax16x16 <- array(NA, 
                           dim=c(length(annualMaxima16x16Boxes), 
                                 16, 16, 30))
for(i in 1:length(annualMaxima16x16Boxes)){
  
  RCMannualMax16x16[i, , , ] <- (annualMaxima16x16Boxes[[i]])[, ,1:30]
  
}

centerLocs <- matrix(NA,
                     nrow=length(annualMaxima16x16Boxes),
                     ncol=2)
for(i in 1:length(annualMaxima16x16Boxes))
{
  temp <- gridPoints16x16Boxes[[i]] 
  centerLocs[i, ] <- c(temp$center)
}

# ## Save the array as a binary file
# dim(RCMannualMax16x16)
# setwd('~/Desktop')
# save(RCMannualMax16x16, file='RCMannualMax16x16.RData')


## -- Generate spatial plot with 16x16 grid across the focused spatial domain --:
# Define grid size and step size for overlap
grid_size <- 16
step <- 16  # 50% overlap (step = grid_size / 2)


# Get spatial dimensions
n_lon <- dim(lon_data_new)[1]
n_lat <- dim(lat_data_new)[2]

png("Spatial-Domain.png",
    units="in",
    width=24,
    height=8,
    res=200)
par(mfrow=c(1,3),
    mar=c(5,4,5,6) + 0.3,
    oma=c(1.4,1.5,2,4))
plot(lon_data,
     lat_data,
     xlab='',
     ylab='',
     yaxt='n',
     xaxt='n',
     cex=2,
     col='grey')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('NA-Cordex',
      side=3,
      cex=3,
      line=3)
mtext('Lon',
      side=1,
      cex=2,
      line=5)
mtext('Lat',
      side=2,
      cex=2,
      line=4)
#world(add=TRUE, col='magenta', lwd=3)
world(add=TRUE, col='grey5', lwd=5)

plot(lon_data_new, 
     lat_data_new,
     xlim=c(min(c(lon_data)), max(c(lon_data))),
     ylim=c(min(c(lat_data)), max(c(lat_data))),
     yaxt='n',
     xaxt='n',
     xlab='',
     ylab= '',
     cex=2,
     col='grey')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('Focused Region',
      side=3,
      cex=3,
      line=3)
mtext('Lon',
      side=1,
      cex=2,
      line=5)
mtext('Lat',
      side=2,
      cex=2,
      line=4)
world(add=TRUE, col='grey5', lwd=5)
#world(add=TRUE, col='magenta', lwd=3)

# Set up the plot
plot(lon_data_new, 
     lat_data_new,
     #xlim=c(min(c(lon_data)), max(c(lon_data))),
     ylim=c(min(c(lat_data)), max(c(lat_data))),
     yaxt='n',
     xaxt='n',
     xlab='',
     ylab= '',
     cex=2,
     col='grey')
axis(2,
     las=1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.5,
     lwd=2,
     padj=0.9)
mtext('With 16x16 tiles',
      side=3,
      cex=3,
      line=3)
mtext('Lon',
      side=1,
      cex=2,
      line=5)
mtext('Lat',
      side=2,
      cex=2,
      line=4)
world(add=TRUE, col='grey5', lwd=5)
#world(add=TRUE, col='magenta', lwd=3)

# Loop over the domain to plot the grid boxes
for (i in seq(1, n_lon - grid_size + 1, by = step)) {
  for (j in seq(1, n_lat - grid_size + 1, by = step)) {
    # Extract the 16×16 grid for lon/lat
    lon_grid <- lon_data_new[i:(i + grid_size - 1), j:(j + grid_size - 1)]
    lat_grid <- lat_data_new[i:(i + grid_size - 1), j:(j + grid_size - 1)]
    
    # Skip NA grids
    if (sum(is.na(lon_grid)) > 0 || sum(is.na(lat_grid)) > 0) {
      next
    }
    
    # Determine grid box corners
    lon_min <- min(lon_grid, na.rm=TRUE)
    lon_max <- max(lon_grid, na.rm=TRUE)
    lat_min <- min(lat_grid, na.rm=TRUE)
    lat_max <- max(lat_grid, na.rm=TRUE)
    
    # Draw the grid box using segments (4 boundary lines)
    segments(lon_min, lat_min, lon_max, lat_min, col='blue', lwd=2.5) # Bottom line
    segments(lon_min, lat_max, lon_max, lat_max, col='blue', lwd=2.5) # Top line
    segments(lon_min, lat_min, lon_min, lat_max, col='blue', lwd=2.5) # Left line
    segments(lon_max, lat_min, lon_max, lat_max, col='blue', lwd=2.5) # Right line
    
    # Calculate and plot the grid center
    lon_center <- mean(c(lon_min, lon_max))
    lat_center <- mean(c(lat_min, lat_max))
    points(lon_center, lat_center, col="red", pch=16, cex=1.5)  # Grid center
  }
}
#world(add=TRUE, col='magenta', lwd=3)  # Add world boundaries
# world(add=TRUE, col='grey5', lwd=5)
dev.off()
################################################################################


########################### --- BIAS-CORRECTION --- ############################
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Proj2-Data2Evaluate.R")

## Loading the training parameter and test true parameters --:
shape_true_train <- trainParameters[ , 1]
kappa2_true_train <- (trainParameters[ , 2]-4)
tau2_true_train <- (trainParameters[ , 3]) # 50000

# Load necessary libraries
library(splines)   # For smoothing splines

## Bias-correction:
## Log Kappa:
log_kappa_true <- 0.5 * log(kappa2_true_train)
log_kappa_pred <- 0.5 * log(kappa2_pred_train)

## Fit a smoothing spline to residuals
fit_bC_kappa <- smooth.spline(x=log_kappa_pred,
                              y=log_kappa_true,
                              spar = 0.8)  # Adjust 'spar' for smoothness

## Extract fitted residuals from the spline model
bclogKappa <- predict(fit_bC_kappa,
                      x=log_kappa_pred)$y  # Smoothed residuals


## Log Tau:
log_tau_true <- 0.5 * log(tau2_true_train)
log_tau_pred <- 0.5 * log(tau2_pred_train)

## Fit a smoothing spline to residuals
fit_bC_tau <- smooth.spline(x=log_tau_pred,
                            y=log_tau_true,
                            spar = 0.8)  # Adjust 'spar' for smoothness

## Extract fitted residuals from the spline model
bclogTau <- predict(fit_bC_tau,
                    x=log_tau_pred)$y  # Smoothed residuals

## Xi
xi_true <- shape_true_train
xi_pred <- shape_pred_train

# # Fit a smoothing spline to residuals
fit_bC_xi <- smooth.spline(x=xi_pred,
                           y=xi_true,
                           spar = 0.8)  # Adjust 'spar' for smoothness

## Extract fitted residuals from the spline model
bcXi <- predict(fit_bC_xi,
                x=xi_pred)$y  # Smoothed residuals



## Plot residuals after smoothing
bplot.xy(xi_true,
         bcXi,
         pch = 19,
         col = alpha('grey', 0.4),
         xlab = expression('True log('~tau~')'),
         ylab = "Estimated",
         main = expression('Before bias correciton log('~tau~')')
)
abline(0, 1, col = 'magenta', lwd = 2)  # Reference line at 0

# # Plot residuals after smoothing
bplot.xy(log_tau_true,
         bclogTau,
         pch = 19,
         col = alpha('grey', 0.4),
         xlab = expression('True log('~kappa~')'),
         ylab = "Estimated",
         main = expression('After bias correction log('~kappa~')')
)
abline(0, 1, col = 'magenta', lwd = 2)  # Reference line at 0


# # Plot residuals after smoothing
bplot.xy(log_kappa_true,
         bclogKappa,
         pch = 19,
         col = alpha('grey', 0.4),
         xlab = expression('True log('~kappa~')'),
         ylab = "Estimated",
         main = expression('After bias correction log('~kappa~')')
)
abline(0, 1, col = 'magenta', lwd = 2)  # Reference line at 0
################################################################################


## -- CNN ESTIMATE -- ##########################################################
## Loading the estimated parameter values across the RCM  focused region:
setwd('~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy')

## Calling the training estimates:
file_path_rcm <- path.expand("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/CaseStudy/est_parameter_test_RCM-NoOverlaps.npy") # NoOverlaps.npy

## -- Import numpy and load the file -- ##
np <- import("numpy")
matRCM <- np$load(file_path_rcm)

## -- Load the NumPy file -- ##
cnnEstRCM <- matRCM
dim(cnnEstRCM)  # 100x100x3


## Bias:
bclogKappaCNN <- predict(fit_bC_kappa,
                         x=0.5*cnnEstRCM[,2])$y 
bclogTauCNN <- predict(fit_bC_tau,
                      x=0.5*cnnEstRCM[,3])$y 
bcXiCNN <- predict(fit_bC_xi,
                   x=cnnEstRCM[,1])$y 

hist(exp(bclogKappaCNN)+4)
hist(exp(bclogTauCNN))
hist(bcXiCNN)

## Make the surface more smooth:
colorTable<- designer.colors(11, 
                             c( "blue","white", "red") )
fitTPSXi <- Tps(centerLocs,
                bcXiCNN, 
                df=20)
fitTPSlogKappa <- Tps(centerLocs,
                      bclogKappaCNN,
                      df=20)
fitTPSbclogTauCNN <- Tps(centerLocs,
                         bclogTauCNN,
                         df=20)

################################################################################
library(RColorBrewer)
jet.colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

# Generate the desired number of colors from this palette
nbcol <- 50
color <- jet.colors(nbcol)
ramp <- colorRamp(c("blue", "white", "red"))

setwd('~/Desktop')
png("CS-Estimates.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,10.5,10) + 0.3,
    oma=c(0.4,2.5,3,8))
surface(fitTPSXi,
        col=color,
        #viridis(100),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
world(add=TRUE, col='grey5', lwd=5)
#world(add=TRUE, col='magenta3', lwd=5)
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=3,
      cex=4,
      line=4)
mtext('Lon',
      side=1,
      cex=3,
      line=7)
mtext('Lat',
      side=2,
      cex=3,
      line=5)


surface(fitTPSlogKappa,
        col=color,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
world(add=TRUE, col='grey5', lwd=5)
#world(add=TRUE, col='magenta3', lwd=5)
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression('log('~kappa~')'),
      side=3,
      cex=3.8,
      line=4)
mtext('Lon',
      side=1,
      cex=3,
      line=7)
mtext('Lat',
      side=2,
      cex=3,
      line=5)

surface(fitTPSbclogTauCNN,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        col=color,
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression('log('~tau~')'),
      side=3,
      cex=3.8,
      line=4)
mtext('Lon',
      side=1,
      cex=3,
      line=7)
mtext('Lat',
      side=2,
      cex=3,
      line=5)
world(add=TRUE, col='grey5', lwd=5)
#world(add=TRUE,  col='magenta3',  lwd=5)
dev.off()


### Using bubble plots instead: 
setwd('~/Desktop')
png("CS-Estimates-I.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,10.5,10) + 0.3,
    oma=c(0.4,2.5,3,8))
bubblePlot(fitTPSXi$x,
           fitTPSXi$y,
           fitTPSXi$z,
        col=color,
        #viridis(100),
        yaxt='n',
        xaxt='n',
        xlab='',
        size=2,
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
world(add=TRUE, col='grey5', lwd=5)
#world(add=TRUE, col='magenta3', lwd=5)
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=3,
      cex=4,
      line=4)
mtext('Lon',
      side=1,
      cex=3,
      line=7)
mtext('Lat',
      side=2,
      cex=3,
      line=5)


surface(fitTPSlogKappa,
        col=color,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
world(add=TRUE, col='grey5', lwd=5)
#world(add=TRUE, col='magenta3', lwd=5)
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression('log('~kappa~')'),
      side=3,
      cex=3.8,
      line=4)
mtext('Lon',
      side=1,
      cex=3,
      line=7)
mtext('Lat',
      side=2,
      cex=3,
      line=5)

surface(fitTPSbclogTauCNN,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        col=color,
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression('log('~tau~')'),
      side=3,
      cex=3.8,
      line=4)
mtext('Lon',
      side=1,
      cex=3,
      line=7)
mtext('Lat',
      side=2,
      cex=3,
      line=5)
world(add=TRUE, col='grey5', lwd=5)
#world(add=TRUE,  col='magenta3',  lwd=5)
dev.off()
################################################################################
# 
# ## Simulate across the grid --:
# source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/simulateSpatialExtFields.R")
# 
# ## -- Step 1: Extract Center Estimates and Parameters -- ##
# # Assuming `parameters` is a matrix 
# parameters <- matrix(nrow = length(gridPoints16x16Boxes), ncol = 3)
# for (i in seq_along(gridPoints16x16Boxes)) {
#   # Extract GEV parameters for the i-th tile
#   parameters[i, ] <- c(bcXiCNN[i],  (exp(bclogKappaCNN[i])+4), (exp(bclogTauCNN[i])) )  # Replace with actual parameter extraction
# }
# dim(parameters)
# 
# tau_vec <- parameters[ ,3]
# tau_vec[tau_vec>0.1] <- 0.1
# parameters[,3] <- tau_vec
# summary(parameters)
# 
# source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/simulateSpatialExtFields.R")
# m <- 100  # Number of ensemble members
# extremalFieldsList <- list()
# tic()
# for (i in seq_along(gridPoints16x16Boxes)) {
#   cat('Loop: ', i, '\n')
#   center <- gridPoints16x16Boxes[[i]]$center
#   s_tile <- cbind(as.vector(gridPoints16x16Boxes[[i]]$lon),
#                   as.vector(gridPoints16x16Boxes[[i]]$lat))
#   extremalFieldsList[[i]] <- generateSpatialExtFields(s_tile, 
#                                                       m, 
#                                                       parameters[i, , drop = FALSE])
# }
# toc()
# 
# # Step 2: Combine the Simulated Fields
# s_tiles <- array(NA, 
#                  dim=c(length(gridPoints16x16Boxes), 256, 2))
# extremalFieldsCombined <- array(NA, 
#                                 dim=c(length(gridPoints16x16Boxes), 256, m))
# for(i in 1:length(gridPoints16x16Boxes)){
#   s_tiles[i, , ] <- cbind(as.vector(gridPoints16x16Boxes[[i]]$lon),
#                           as.vector(gridPoints16x16Boxes[[i]]$lat))
#   
#   extremalFieldsCombined[i, , ] <- extremalFieldsList[[i]]
# }
# 
# simFieldsExtremes <- array(extremalFieldsCombined, 
#                            dim=c(dim(extremalFieldsCombined)[1]*dim(extremalFieldsCombined)[2], 
#                                  dim(extremalFieldsCombined)[3]))
# 
# simLocs <- array(s_tiles, 
#                  dim=c(dim(s_tiles)[1]*dim(s_tiles)[2], 
#                        dim(s_tiles)[3]))
# 
# set.panel(1,2)
# bubblePlot(simLocs,
#            (simFieldsExtremes[,20]-median(simFieldsExtremes[,20]))/IQR(simFieldsExtremes[,20]),
#            size=0.6,
#            zlim=c(-1, 13),
#            main='Simulated',
#            col=color_palette)
# world(add=TRUE, col='grey3', lwd=2)
# 
# 
# bubblePlot(lon_data_new,
#            lat_data_new  , 
#            (annual_maxima_new[,,20]-median(annual_maxima_new[,,20], na.rm=TRUE))/IQR(annual_maxima_new[,,20], na.rm=TRUE),
#            col=color_palette,
#            size=0.8,
#            main='Orginal',
#            xlab='Lon',
#            ylab='Lat',
#            na.rm=TRUE)
# world(add=TRUE, col='grey3', lwd=2,)
# 
# 
# ### Trying to make the surface smooth:
# library(fields)   # For image.smooth() function
# library(akima)    # For interpolation
# 
# # Determine the spatial dimensions
# nrow <- dim(extremalFieldsCombined)[1]
# ncol <- dim(extremalFieldsCombined)[2]
# 
# # Reshape the simulated field into a matrix for smoothing
# simFieldMatrix <- matrix(simFieldsExtremes[,10],
#                          nrow=nrow,
#                          ncol=ncol)
# 
# # Apply Gaussian smoothing
# smoothedSimFields <- image.smooth(simFieldMatrix)
# 
# set.panel(1,2)
# 
# # Convert back to a flattened vector for bubblePlot
# smoothedSimFieldsVec <- as.vector(smoothedSimFields$z)
# 
# # Plot the smoothed simulated field
# bubblePlot(simLocs,
#            (smoothedSimFieldsVec - median(smoothedSimFieldsVec, na.rm=TRUE)) / 
#              IQR(smoothedSimFieldsVec, na.rm=TRUE),
#            size=0.6, 
#            zlim=c(-0.8, 15), 
#            main='Simulated (Smoothed)', 
#            col=color_palette)
# world(add=TRUE, col='grey3', lwd=2)
# 
# # Original field remains unchanged
# bubblePlot(lon_data_new,
#            lat_data_new,  
#            (annual_maxima_new[,,20] - median(annual_maxima_new[,,20], na.rm=TRUE)) / 
#              IQR(annual_maxima_new[,,20], na.rm=TRUE),
#            col=color_palette,
#            size=0.8,
#            main='Original',
#            xlab='Lon',
#            ylab='Lat',
#            na.rm=TRUE)
# world(add=TRUE, col='grey3', lwd=2)
# 
# 


# Set up the plot



png("Spatial-Domain-NA-Cordex.png",
    units="in", 
    width=7,
    height=8,
    res=200)
plot(lon_data_new, 
     lat_data_new,
     ylim=c(min(c(lat_data)), max(c(lat_data))),
     yaxt='n',
     xaxt='n',
     xlab='',
     ylab= '',
     cex=2,
     col='grey')

axis(2, las=1, cex.axis=1.2, lwd=2, padj=0.9)
axis(1, cex.axis=1.2, lwd=2, padj=0.9)

mtext('With 16x16 tiles', side=3, cex=1.2, line=1)
mtext('Lon', side=1, cex=1.2, line=3)
mtext('Lat', side=2, cex=1.2, line=3)
# , line=4)

world(add=TRUE, col='grey5', lwd=2)

# Initialize tile counter
tile_number <- 1  

selected_Tiles <- c(16, 18, 19, 80, 95, 110, 125,138, 141)

# Loop over the domain to plot the grid boxes
for (i in seq(1, n_lon - grid_size + 1, by = step)) {
  for (j in seq(1, n_lat - grid_size + 1, by = step)) {
    
    # Extract the 16×16 grid for lon/lat
    lon_grid <- lon_data_new[i:(i + grid_size - 1), j:(j + grid_size - 1)]
    lat_grid <- lat_data_new[i:(i + grid_size - 1), j:(j + grid_size - 1)]
    
    # Skip NA grids
    if (sum(is.na(lon_grid)) > 0 || sum(is.na(lat_grid)) > 0) {
      next
    }
    
    # Determine grid box corners
    lon_min <- min(lon_grid, na.rm=TRUE)
    lon_max <- max(lon_grid, na.rm=TRUE)
    lat_min <- min(lat_grid, na.rm=TRUE)
    lat_max <- max(lat_grid, na.rm=TRUE)
    
    # Draw the grid box using segments (4 boundary lines)
    segments(lon_min, lat_min, lon_max, lat_min, col='blue', lwd=2.5) # Bottom line
    segments(lon_min, lat_max, lon_max, lat_max, col='blue', lwd=2.5) # Top line
    segments(lon_min, lat_min, lon_min, lat_max, col='blue', lwd=2.5) # Left line
    segments(lon_max, lat_min, lon_max, lat_max, col='blue', lwd=2.5) # Right line
    
    # Calculate grid center
    lon_center <- mean(c(lon_min, lon_max))
    lat_center <- mean(c(lat_min, lat_max))
    
    # Plot grid center
    points(lon_center, lat_center, col="red", pch=16, cex=1)
    # Label tile number
    #text(lon_center, lat_center, labels=tile_number, col="black", cex=0.6, font=2)
    # Label tile number
    if (tile_number %in% selected_Tiles) {
      text(lon_center, lat_center,
           labels = tile_number,
           col    = "black",
           cex    = 1,
           font   = 2)
    }


    # Increment tile number
    tile_number <- tile_number + 1  
    
  }
}
# world(add=TRUE, col='grey5', lwd=2)
dev.off()



