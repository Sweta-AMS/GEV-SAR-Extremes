rm(list=ls())

library(tictoc)
library(fields)
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/CaseStudy/RCM-dataLoad.R")

## Save the data
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/CaseStudy/ColorPalette.R")
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


## -- NO OVERLAPPING - DEFINING 16X16 GRID ACROSS THE SPATIAL DOMAIN --
annualMaxima16x16Boxes <- list()
gridPoints16x16Boxes <- list()
point_count_list <- list()
tile_info_list <- list()  # Store tile number and center

grid_size <- 16
step <- 16
tile_number <- 1  # Initialize tile counter

n_lon <- dim(lon_data_new)[1]  # 297 longitude points
n_lat <- dim(lat_data_new)[2]  # 281 latitude points

# Loop over domain to create tiles
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
    
    # Store tile information
    tile_info_list[[length(tile_info_list) + 1]] <- list(
      tile_id = tile_number,
      center_lon = center_lon,
      center_lat = center_lat
    )
    
    # Store results
    point_count_list[[length(point_count_list) + 1]] <- sum(!is.na(lon_grid))
    annualMaxima16x16Boxes[[length(annualMaxima16x16Boxes) + 1]] <- annual_maxima_grid
    gridPoints16x16Boxes[[length(gridPoints16x16Boxes) + 1]] <- list(
      tile_id = tile_number,  # Store tile number
      lon = lon_grid,
      lat = lat_grid,
      center = c(center_lon, center_lat)  # Track center of the grid
    )
    
    # Increment tile number
    tile_number <- tile_number + 1
  }
}
# 166 non-overlapping tiles

centerLocs <- matrix(NA,
                     nrow=length(annualMaxima16x16Boxes),
                     ncol=2)
for(i in 1:length(annualMaxima16x16Boxes))
{
  temp <- gridPoints16x16Boxes[[i]] 
  centerLocs[i, ] <- c(temp$center)
}


## -- Loading the bias-corrected parameter of the CNNs -- :
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
load('bias-correct-CNN-CS-Parameters.RData')
dim(parameters) # 166 x 3


### -- Simulating the data across tiles --:
m <- 30
observationList <- array(NA, 
                         dim = c(length(gridPoints16x16Boxes),
                                 256,
                                 m))

for (i in 1:length(gridPoints16x16Boxes)) {  # Ensure safe indexing
  
  # Assign the matrix to the correct index
  observationList[i, , ] <- matrix(
    annualMaxima16x16Boxes[[i]][,,1:30], 
    nrow = 256, ncol = 30
  )
}

## Breaking the observation across the I, II, and III division:
observationListI <- observationList[1:66, ,]
observationListII <- observationList[67:116, ,]
observationListIII <- observationList[117:166, ,]

load('full-tileList-CS.RData')
dim(tilesList)

tilesListI <- tilesList[, 1:66, ,]
tilesListII <- tilesList[ ,67:116, ,]
tilesListIII <- tilesList[ ,117:166, ,]

## Simulating the process:
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Pro2-simulateSpatialExtFields.R")
R <- 400
m <- 30  # Number of ensemble members

## Save the data:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
load('full-extremalFieldsList-CS.RData')
dim(extremalFieldsList)

load('full-centerList-CS.RData')
dim(centerList)

# ### -- Quantile Mapping for bias-correction -- :
# n_yrs <- dim(observationList)[3]
# m <- dim(observationList)[3]
# p_vec <- seq(0.01, 0.99, length.out=100)
# 
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
load('biasCorExtremalFields-CS-I.RData')
load('biasCorExtremalFields-CS-II.RData')
load('biasCorExtremalFields-CS-III.RData')

dim(biasCorExtremalFieldsI) # 1:66
dim(biasCorExtremalFieldsII) # 67:(66+50)
dim(biasCorExtremalFieldsIII) # (66+50+1):(66+50+50)

 
### -- Computing qq plot for true versus simulated --:
# Number of quantile points
n_points <- 256 * 30
p_vec <- (1:n_points) / (n_points + 1)
R <- 400

##  QQ data plot
qq_listIII <- list()
tic()
for (i in 1:50) {
  cat("Processing Tile:", i, "\n")

  # Standardize Observations (tile i)
  obs <- observationListIII[i, ,]
  obs <- (obs - median(obs)) /IQR(obs)
  qObs <- quantile(obs, p=p_vec)

  # Store simulated quantiles across all R replications
  qSim_mat <- matrix(NA,
                     nrow = n_points,
                     ncol = R)

  for (r in 1:R) {
    #cat("Processing Realization:", r, "\n")

    # Standardize Simulated Data (tile i, realization r)
    sim <- biasCorExtremalFieldsIII[r, i, ,]
    sim <- (sim - median(sim)) /IQR(sim) # changed from IQR to SD
    qSim_mat[, r] <- quantile(sim, p=p_vec)
  }

  ## Compute Mean and SD-Based Confidence Intervals
  qSim_mean  <- apply(qSim_mat, 1, mean, na.rm = TRUE)  # Compute mean
  qSim_sd    <- apply(qSim_mat, 1, sd, na.rm = TRUE)    # Compute standard deviation

  ## Compute 95% confidence intervals using normal approximation
  qSim_lower <- qSim_mean - 1.96 * qSim_sd
  qSim_upper <- qSim_mean + 1.96 * qSim_sd

  # qSim_lower <- apply(qSim_mat, 1, quantile, 0.025, na.rm = TRUE)  # Compute mean
  # qSim_upper <- apply(qSim_mat, 1, quantile, 0.975, na.rm = TRUE)  # Compute mean
  #
  # Store in data frame
  qq_df <- data.frame(
    Observed = qObs,
    Simulated_Mean = qSim_mean,
    Simulated_Lower  = qSim_lower,
    Simulated_Upper  = qSim_upper,
    Tile = paste("Tile", gridPoints16x16Boxes[[i+66+50]]$tile_id) #+66+50
  )

  # Append to list
  qq_listIII[[i]] <- qq_df
}
toc()
# # # time:  2.78555 mins

for (i in 1:50) {  # Loop over tiles
  cat('i: ', i, '\n')
  plot(qq_listIII[[i]]$Observed ,
       qq_listIII[[i]]$Simulated_Mean,
       pch=19,
       cex=0.6,
       ylim=c(-2, 12),
       xlim=c(-2, 12),
       xlab='Observed',
       ylab='Simulated',
       main=paste0(qq_listIII[[i]]$Tile[1]),
       col=alpha('blue', 0.8))
  abline(0, 1, col="pink3", lwd=3, lty=2)
  lines(qq_listIII[[i]]$Observed,
        qq_listIII[[i]]$Simulated_Lower,
        col='grey',
        lwd=2,
        lty=2)
  lines(qq_listIII[[i]]$Observed,
        qq_listIII[[i]]$Simulated_Upper,
        col='grey',
        lwd=2,
        lty=2)
}


# ## ----- COMPUTE and COMPARE THE MADOGRAM -------------------------------------:
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Madogram-similar2fields.R")
common_breaks <- seq(0, 10, length.out=50)  # 50 distance bins

# # Compute madogram for simulated data
# madogramSimIII <- array(NA,
#                       dim=c(R,50, 49))
# tic()
# for (r in 1:R) {
#   
#   cat('r :',  r, '\n')
#   
#   for (i in 1:50) {
#     
#     cat('i :', i, '\n')
#     value <- (biasCorExtremalFieldsIII[r, i, , ] - median(biasCorExtremalFieldsIII[r, i, , ]))/sd(biasCorExtremalFieldsIII[r, i, , ])
#     
    # fitMado <- madogram(tilesListIII[1,i, ,],
    #                     value,
    #                     N=50,
    #                     breaks=common_breaks)
# 
#     madogramSimIII[r, i, ] <- fitMado$stats['mean',]
#   }
# }
# toc()
# centerPoint <- fitMado$centers
# save(centerPoint, file='centerPoint.RData')
## Close to an hour

load('centerPoint.RData')

# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
# # load('madogramSimI-Median-SD.RData')
# # load('madogramSimII-Median-SD.RData')
# # load('madogramSimIII-Median-SD.RData')

load('madogramSimI-Median-IQR.RData')
load('madogramSimII-Median-IQR.RData')
load('madogramSimIII-Median-IQR.RData')

## -- Observe Madogram --:
madogramObsI <- array(NA,
                      dim=c(66, 49))
for(i in 1:66)
{
  cat('i: ', i, '\n')
  valueObs <- (observationListI[i, , ] - median(observationListI[i, , ]))/IQR(observationListI[i, , ])
  fitObsMado <-  madogram(tilesListI[1, i, ,],
                          valueObs,
                          N=50,
                          breaks=common_breaks)
  madogramObsI[i, ] <- fitObsMado$stats['mean',]
}

madogramObsII <- array(NA,
                      dim=c(50, 49))
for(i in 1:50)
{
  cat('i: ', i, '\n')
  valueObs <- (observationListII[i, , ] - median(observationListII[i, , ]))/IQR(observationListII[i, , ])
  fitObsMado <-  madogram(tilesListII[1, i, ,],
                          valueObs,
                          N=50,
                          breaks=common_breaks)
  madogramObsII[i, ] <- fitObsMado$stats['mean',]
}


madogramObsIII <- array(NA,
                       dim=c(50, 49))
for(i in 1:50)
{
  cat('i: ', i, '\n')
  valueObs <- (observationListIII[i, , ] - median(observationListIII[i, , ]))/IQR(observationListIII[i, , ])
  fitObsMado <-  madogram(tilesListIII[1, i, ,],
                          valueObs,
                          N=50,
                          breaks=common_breaks)
  madogramObsIII[i, ] <- fitObsMado$stats['mean',]
}

# bplot.xy(centerPoint,
#          madogramSimII[2,1, ])
# points(centerPoint,
#        madogramObsII[1,],
#        pch=19,
#        col='pink3')
# bplot.xy(centerPoint,
#          madogramSimIII[2,1, ])
# points(centerPoint,
#        madogramObsIII[1,],
#        pch=19,
#        col='pink3')

## -- Plots generation --:
# # Compute mean and confidence intervals for simulated madogram
madogramSim_meanIII <- apply(madogramSimIII, c(2,3), mean, na.rm=TRUE)
madogramSim_SDIII <- apply(madogramSimIII, c(2,3), sd, na.rm=TRUE)

# # #madogramSim_mean - 1.96 * madogramSim_SD#madogmedianramSim_mean - 1.96 * madogramSim_SD
madogramSim_lowerIII <- madogramSim_meanIII - 1.96 * madogramSim_SDIII
madogramSim_upperIII <- madogramSim_meanIII + 1.96 * madogramSim_SDIII


# Optional: Store plots in a list
for (i in 1:50) {
  cat('Loop: ', i, '\n')

  df <- data.frame(
    Distance = centerPoint,
    MadogramMean = as.vector(madogramSim_meanIII[i, ]),
    MadogramUpper = as.vector(madogramSim_upperIII[i, ]),
    MadogramLower = as.vector(madogramSim_lowerIII[i, ]),
    MadogramObs = as.vector(madogramObsIII[i, ])
  )

  # Create a new plot
  plot(df$Distance,
       df$MadogramMean,
       type = "l",
       col = "blue",
       lwd = 2,
       xlim = c(0, 5.5),
       ylim = c(0.1, 0.8),
       xlab = "Distance",
       ylab = "Madogram",
       main = paste("Tile", gridPoints16x16Boxes[[i+66+50]]$tile_id))
         #paste("Tile:", rand_tiles[i]))

  lines(df$Distance,
        df$MadogramLower,
        col = "grey",
        lty=2,
        line=0.8,
        lwd = 2)
  lines(df$Distance,
        df$MadogramUpper,
        col = "grey",
        lty=2,
        line=0.8,
        lwd = 2)

  # Replot the mean line to appear on top of the ribbon
  lines(df$Distance, df$MadogramMean, col = "blue", lwd=2)

  # Add observed madogram points
  points(df$Distance, df$MadogramObs, col = "pink3", pch = 19, cex = 1)
}

## -- Ratio of Observe Madogram versus Mean Simulated Madogram --: 
madogramSim_meanI <- apply(madogramSimI,
                           c(2,3),
                           mean,
                           na.rm=TRUE)
madogramSim_meanII <- apply(madogramSimII,
                            c(2,3), 
                            mean,
                            na.rm=TRUE)
madogramSim_meanIII <- apply(madogramSimIII,
                             c(2,3),
                             mean,
                             na.rm=TRUE)

madogramSimMean <- matrix(NA,
                          nrow=166, 
                          ncol=49)

madogramSimMean[1:66, ] <- madogramSim_meanI 
madogramSimMean[67:116, ] <- madogramSim_meanII 
madogramSimMean[117:166, ] <- madogramSim_meanIII

madogramObs <- matrix(NA, 
                      nrow=166,
                      ncol=49)
madogramObs[1:66, ] <- madogramObsI
madogramObs[67:116, ] <- madogramObsII
madogramObs[117:166, ] <- madogramObsIII

madogramRatio <- matrix(NA, nrow= 166, ncol=49)
for(i in 1:166) {
  madogramRatio[i, ] <- abs((madogramObs[i,] - madogramSimMean[i,])/madogramObs[i,])
  
}
length(centerPoint) # 49 

bplot.xy(centerPoint,
         madogramRatio,
         ylim=c(-0.001, 0.5),
         cex=.5,
         xlab='Relative Error',
         ylab='Distance',
         main='Madogram')
abline(h=0, col='pink3', lwd=2, lty=2)


# Calculate the 25th and 75th percentiles for each center point (ignoring NA/NaN values)
q5_curve <- apply(madogramRatio, 2, function(x) quantile(x, 0.05, na.rm = TRUE))
q95_curve <- apply(madogramRatio, 2, function(x) quantile(x, 0.95, na.rm = TRUE))

## Calculate the IQR for each center point
# iqr_curve <- q75_curve - q25_curve
#q25_curve - 1.5 * iqr_curve
#q95_curve + 1.5 * iqr_curve

# Identify outliers based on the IQR method
lower_bound <- q5_curve
upper_bound <- q95_curve 


# Remove outliers (points outside the bounds)
madogramRatio_no_outliers <- apply(madogramRatio, 2, function(x) {
  x[x < lower_bound | x > upper_bound] <- NA  # Replace outliers with NA
  return(x)
})

# Plot non-outlier points using bplot.xy
bplot.xy(centerPoint,
         log(madogramRatio_no_outliers),
         #ylim = c(-0.001, 0.48),
         cex = 0.5,
         xlab = 'log(ARE)',
         ylab = 'Distance',
         main = 'Relative Error in Madogram')

# Add a horizontal line at y = 0 (for reference)
abline(h = 0, col = 'pink3', lwd = 3, lty=2)


## -- QQ mapping --:
cat('i: ', i, '\n')
plot(qq_listI[[i]]$Observed ,
     qq_listI[[i]]$Simulated_Mean,
     pch=19,
     cex=0.6,
     ylim=c(-2, 12),
     xlim=c(-2, 12),
     xlab='Observed',
     ylab='Simulated',
     main=paste0(qq_listI[[i]]$Tile[1]),
     col=alpha('blue', 0.8))
abline(0, 1, col="pink3", lwd=3, lty=2)
lines(qq_listIII[[i]]$Observed,
      qq_listIII[[i]]$Simulated_Lower,
      col='grey',
      lwd=2,
      lty=2)
lines(qq_listI[[i]]$Observed,
      qq_listI[[i]]$Simulated_Upper,
      col='grey',
      lwd=2,
      lty=2)

qqSimMean <- matrix(NA, 
                    nrow=166,
                    ncol= 7680)

qqObs <- matrix(NA,
                nrow=166,
                ncol= 7680)
for(i in 1:66) {
  qqObs[i,] <- qq_listI[[i]]$Observed
  qqSimMean[i, ] <- qq_listI[[i]]$Simulated_Mean
}

for(i in 1:50){
  qqObs[i+66,] <- qq_listII[[i]]$Observed
  qqSimMean[i+66, ] <- qq_listII[[i]]$Simulated_Mean
}

for(i in 1:50){
  qqObs[i+66+50,] <- qq_listIII[[i]]$Observed
  qqSimMean[i+66+50, ] <- qq_listIII[[i]]$Simulated_Mean
}


relativeErrorQQ <- matrix(NA, 
                          nrow=166,
                          ncol=length(p_vec))
for(i in 1:166) {
  relativeErrorQQ[i, ] <- abs((qqObs[i,] - qqSimMean[i,])/qqObs[i,])
}

bplot.xy(p_vec, relativeErrorQQ, ylim=c(-0.01, 0.4))
## Remove Inf and NA values before computing summary statistics:

# Filter values that fall within Q1 and Q3
relativeErrorQQFilter <- matrix(NA,
                                nrow=166,
                                ncol=length(p_vec))
for(i in 1:166) {
  value <- relativeErrorQQ[i, ]
  Q5 <- quantile(value, p=0.05)
  Q95 <- quantile(value, p=0.95)
  value[ value < Q1] <- Q5
  value[value > Q3] <- Q95
  relativeErrorQQFilter[i, ] <- value
}

bplot.xy(p_vec, 
         relativeErrorQQFilter,
         cex=0.5,
         ylim=c(-0.4,0.4))
abline(h=0,
       col='pink3', 
       lwd=2,
       lty=2)

## Calculate the 25th and 75th percentiles for each center point (ignoring NA/NaN values)
q5QQ <- apply(relativeErrorQQ, 
               2,
               function(x) quantile(x, 0.05, na.rm = TRUE))
q95QQ <- apply(relativeErrorQQ,
               2, function(x) quantile(x, 0.95, na.rm = TRUE))

# # Calculate the IQR for each center point
# iqrQQ <- q75QQ - q25QQ
# # q25QQ - 1.5 * iqrQQ
# q75QQ + 1.5 * iqrQQ

# Identify outliers based on the IQR method
lower_bound <- q5QQ
upper_bound <- q95QQ

# Remove outliers (points outside the bounds)
qq_no_outliers <- apply(relativeErrorQQ, 2, function(x) {
  x[x < lower_bound | x > upper_bound] <- NA  # Replace outliers with NA
  return(x)
})

# Plot non-outlier points using bplot.xy
bplot.xy(p_vec,
         log(qq_no_outliers),
         #ylim = c(-0.001, 0.48),
         cex=0.5,
         # outlier=FALSE,
         xlab = 'Relative Error',
         ylab = 'Distance',
         main = 'Relative Error in Quantile values')
# Add a horizontal line at y = 0 (for reference)
abline(h = 0, col = 'pink3', lwd = 3, lty=2)

bplot.xy(p_vec, 
         relativeErrorQQFilter,
         cex=0.5,
         ylab = 'Relative Error',
         xlab = 'Probability',
         main='Quantiles',
         ylim=c(-0.4,0.4))
abline(h=0, col='pink3', lwd=2, lty=2)

