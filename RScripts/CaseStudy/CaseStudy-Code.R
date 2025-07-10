rm(list=ls())
library(tictoc)
library(fields)
source("~/Desktop/GEV-SAR/RScripts/CaseStudy/RCM-dataLoad.R")

## Save the data
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
source("ColorPalette.R")

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



## --- NO OVERLAPPING - DEFINING GRID 16X16 GRID ACROSS THE SPATIAL DOMAIN ---:
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
# 166 non-overlapping tiles

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

## Fit a smoothing spline to residuals
fit_bC_xi <- smooth.spline(x=xi_pred,
                           y=xi_true,
                           spar = 0.8)  # Adjust 'spar' for smoothness

## Extract fitted residuals from the spline model
bcXi <- predict(fit_bC_xi,
                x=xi_pred)$y  # Smoothed residuals
################################################################################


## -- CNN ESTIMATE -- ##########################################################
## Loading the estimated parameter values across the RCM  focused region:
setwd('~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy')

## Calling the training estimates:
file_path_rcm <- path.expand("~/Desktop/Research/Proj2-SAR/GEV-SAR/Rscripts/CaseStudy/est_parameter_test_RCM-NoOverlaps.npy") # NoOverlaps.npy

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
################################################################################



############################### -- Validation -- ###############################
# madogram computation
# return period computation: 10-yr, 30-yr 

## -- Step 1: Extract Center Estimates and Parameters -- ##
# Assuming `parameters` is a matrix 
parameters <- matrix(NA, 
                     nrow = length(gridPoints16x16Boxes),
                     ncol = 3)
for (i in seq_along(gridPoints16x16Boxes)) {
  # Extract GEV parameters for the i-th tile
  parameters[i, ] <- c(bcXiCNN[i],  (exp(bclogKappaCNN[i])+4), (exp(bclogTauCNN[i])) )  # Replace with actual parameter extraction
}
dim(parameters)

tau_vec <- parameters[ ,3]
tau_vec[tau_vec>0.1] <- 0.1
parameters[,3] <- tau_vec
summary(parameters)


###  -- Randomly selected tiles --: 
## lets select just 3 right know:
set.seed(123)
rand_tiles <- sample(seq_along(annualMaxima16x16Boxes),
                     size=27,
                     replace=FALSE)
print(rand_tiles)

# set.seed(111)
# rand_tiles <- sample(seq_along(annualMaxima16x16Boxes),
#                      size=66,
#                      replace=FALSE)
# print(rand_tiles)
# # rand_tiles <- c(108) # focused gulf of mexico
# 
# source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/simulateSpatialExtFields.R")
# R <- 1000
# m <- 30  # Number of ensemble members
# 
# extremalFieldsList <- array(NA,
#                             dim=c(R , length(rand_tiles),  256, m))
# tilesList <- array(NA,
#                    dim=c(R , length(rand_tiles), 256, 2))
# 
# centerList <- array(NA,
#                     dim=c(R , length(rand_tiles), 2, 2))
# tic()
# for (r in 1:R) {
#   cat('r: ', r, '\n')
# 
#   for (idx in seq_along(rand_tiles)) {  # Use seq_along for safe indexing
#     i <- rand_tiles[idx]  # Get the actual tile index
#     cat('i: ', i, '\n')
# 
#     if (i > length(gridPoints16x16Boxes)) {
#       stop("Index out of bounds: ", i, " is larger than ",
#            length(gridPoints16x16Boxes))
#     }
# 
#     center <- gridPoints16x16Boxes[[i]]$center
# 
#     s_tile <- cbind(as.vector(gridPoints16x16Boxes[[i]]$lon),
#                     as.vector(gridPoints16x16Boxes[[i]]$lat))
# 
#     tilesList[r, idx, , ] <- s_tile  # Use `idx` instead of `i`
#     centerList[r, idx, , ] <- c(gridPoints16x16Boxes[[i]]$center)
# 
#     extremalFieldsList[r, idx, , ] <- generateSpatialExtFields(s_tile,
#                                                                m,
#                                                                parameters[i, , drop = FALSE])
#   }
# }
# toc()
# # 10.6 hrs
# dim(extremalFieldsList)
# 

## Save the data:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
# save(extremalFieldsList, 
#      file='extremalFieldsList-CS.RData')

load('extremalFieldsList-CS.RData')
dim(extremalFieldsList)

# save(tilesList, 
#      file='tileList-CS.RData')

load('tileList-CS.RData')
dim(tilesList)

# save(centerList, 
#      file='centerList-CS.RData')

load('centerList-CS.RData')
dim(centerList)

### -- Calling for the obsevations --:
m <- 30
observationList <- array(NA, 
                         dim = c(length(rand_tiles),
                                 256,
                                 m))

for (idx in seq_along(rand_tiles)) {  # Ensure safe indexing
  i <- rand_tiles[idx]  # Get the actual tile index
  
  # Assign the matrix to the correct index
  observationList[idx, , ] <- matrix(
    annualMaxima16x16Boxes[[i]][,,1:30], 
    nrow = 256, ncol = 30
  )
}


### -- Quantile Mapping for bias-correction -- :
n_yrs <- dim(observationList)[3]
m <- n_yrs  # Ensure m aligns with n_yrs
p_vec <- seq(0.01, 0.99, length.out=100)

# # Initialize bias-corrected array
# biasCorExtremalFields <- array(NA, 
#                                dim=c(R, 
#                                      length(rand_tiles), 
#                                      256, 
#                                      m))
# 
# tic()
# for (i in seq_along(rand_tiles)) {  # Avoid potential empty index issues
#   cat('Processing tile:', i, '\n')
#   
#   for (t in seq_len(n_yrs)) {  # Use seq_len() for safe iteration
#     cat('Time:', t, '\n')
#     
#     obsProcess <- observationList[i, , t]
#     qObs <- as.numeric(quantile(obsProcess, 
#                                 p=p_vec, 
#                                 na.rm=TRUE))
#     
#     for (r in seq_len(R)) {  # Use seq_len() instead of 1:R
#       cat('Repetition:', r, '\n')
#       
#       simProcess <- extremalFieldsList[r, i, , t]
#       qSim <- as.numeric(quantile(simProcess, p=p_vec, na.rm=TRUE))
#       
#       # Check for degenerate cases (e.g., constant values in qSim)
#       fitQM <- lm(qObs ~ qSim)
#       biasCorExtremalFields[r, i, , t] <- predict(fitQM, 
#                                                   newdata = data.frame(qSim = simProcess))
#     }
#   }
# }
# toc()
# dim(biasCorExtremalFields)
# # elapsed time: 526.125 sec 

setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
# save(biasCorExtremalFields,
#      file='biasCorExtremalFields-CS.RData')
load('biasCorExtremalFields-CS.RData')


set.panel(1,1)
bubblePlot(tilesList[1,1, , ],
           (observationList[rand_tiles[16], , 1] - median(observationList[rand_tiles[16],  ,1]))/IQR(observationList[rand_tiles[16],  ,1]),
           col=color_palette,
           size=0.8,
           zlim=c(-1.8, 3.5),
           main='Orginal',
           xlab='Lon',
           ylab='Lat')
          # na.rm=TRUE)
bubblePlot(tilesList[1,1, , ],
           (biasCorExtremalFields[rand_tiles[16], 1, , 1]-median(biasCorExtremalFields[rand_tiles[16], 1, , 1]))/IQR(biasCorExtremalFields[rand_tiles[16], 1, , 1]),
           size=0.8,
           zlim=c(-1.8, 3.5),
           main='Simulated',
           xlab='Lon',
           ylab='Lat',
           col=color_palette)


# ## ----- COMPUTE and COMPARE THE MADOGRAM -------------------------------------:
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Madogram-similar2fields.R")
common_breaks <- seq(0, 10, length.out=50)  # 50 distance bins

# Compute madogram for simulated data
R <- 1000
madogramSim <- array(NA,
                     dim=c(R, length(rand_tiles), 49))
for (r in 1:R) {

  cat('r :',  r, '\n')

  for (i in seq_along(rand_tiles)) {

    cat('i :', i, '\n')
    value <- (biasCorExtremalFields[r, i, , ] - median(biasCorExtremalFields[r, i, , ]))/IQR(biasCorExtremalFields[r, i, , ])

    fitMado <- madogram(tilesList[r,i, ,],
                        value,
                        N=50,
                        breaks=common_breaks)

    madogramSim[r, i, ] <- fitMado$stats['mean',]
  }
}
centerPoint <- fitMado$centers

# # elapsed time: 
#save(madogramSim, file='madogramSim-List.RData')

madogramObs <- array(NA,
                     dim=c(length(rand_tiles), 49))
for(i in 1:length(rand_tiles))
{
  cat('i: ', i, '\n')
  valueObs <- (observationList[i, , ] - median(observationList[i, , ]))/IQR(observationList[i, , ])
  fitObsMado <-  madogram(tilesList[1,i, ,],
                          valueObs,
                          N=50,
                          breaks=common_breaks)
  madogramObs[i, ] <- fitObsMado$stats['mean',]
}

bplot.xy(centerPoint,
         madogramSim[2,1, ])
points(centerPoint,
       madogramObs[1,],
       pch=19,
       col='pink3')


 ## -- Plots generation --:
## GPT suggested:
library(ggplot2)

# # Compute mean and confidence intervals for simulated madogram
madogramSim_mean <- apply(madogramSim, c(2,3), mean, na.rm=TRUE)
madogramSim_SD <- apply(madogramSim, c(2,3), sd, na.rm=TRUE)

# ## Define upper and lower bounds (mean ± 1.96 * SD for approx 95% CI)
madogramSim_median <- apply(madogramSim, c(2,3), median, na.rm=TRUE)
madogramSim_upper <- apply(madogramSim, c(2,3), function(x) quantile(x, probs = 0.975, na.rm = TRUE))
madogramSim_lower <- apply(madogramSim, c(2,3), function(x) quantile(x, probs = 0.025, na.rm = TRUE))

# # #madogramSim_mean - 1.96 * madogramSim_SD#madogmedianramSim_mean - 1.96 * madogramSim_SD
# madogramSim_lower <- madogramSim_mean - 1.96 * madogramSim_SD
# madogramSim_upper <- madogramSim_mean + 1.96 * madogramSim_SD

# Convert data to long format for ggplot
plot_list <- list()  # Optional: Store plots in a list

for(i in 1:length(rand_tiles)) {
  cat('Loop: ', i, '\n')

  df <- data.frame(
    Distance = centerPoint,
    MadogramMean = as.vector(madogramSim_median[i, ]),
    MadogramUpper = as.vector(madogramSim_upper[i, ]),
    MadogramLower = as.vector(madogramSim_lower[i, ]),
    MadogramObs = as.vector(madogramObs[i, ])
  )

  # Create the plot
  p <- ggplot(df, aes(x = Distance)) +
    geom_ribbon(aes(ymin = MadogramLower, ymax = MadogramUpper),
                fill = "lightblue", alpha = 0.4) +
    geom_line(aes(y = MadogramMean), color = "blue",
              size = 2, linetype = "solid") +
    geom_point(aes(y = MadogramObs),
               color = "pink3", size = 3, shape = 19) +
    labs(title = paste("Tile: ", rand_tiles[i]),
         x = "Distance",
         y = "Madogram") +
    theme_minimal() +
    xlim(0, 6) +
    ylim(0.1, 0.8)+
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15)
    )

  # Print the plot (necessary inside loops)
  print(p)

  # Optional: Store plot in a list
  plot_list[[i]] <- p

  # Optional: Save each plot as an image file
  #ggsave(filename = paste0("madogram_plot_", i, ".png"), plot = p, width = 6, height = 4, dpi = 300)
}

# Optional: Store plots in a list
plot_list <- list()
for (i in 1:length(rand_tiles)) {
  cat('Loop: ', i, '\n')
  
  df <- data.frame(
    Distance = centerPoint,
    MadogramMean = as.vector(madogramSim_median[i, ]),
    MadogramUpper = as.vector(madogramSim_upper[i, ]),
    MadogramLower = as.vector(madogramSim_lower[i, ]),
    MadogramObs = as.vector(madogramObs[i, ])
  )
  
  # Create a new plot
  plot(df$Distance,
       df$MadogramMean,
       type = "l",
       col = "blue",
       lwd = 2,
       xlim = c(0, 6),
       ylim = c(0.1, 0.8),
       xlab = "Distance",
       ylab = "Madogram",
       main = paste("Tile:", rand_tiles[i]))
  
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
  
  # Store plot in list (if needed)
  plot_list[[i]] <- recordPlot()
}



### -- Computing qq plot for true versus simulated --:
library(ggplot2)

# Number of quantile points
n_points <- 256 * 30
p_vec <- (1:n_points) / (n_points + 1)

# Initialize list to store QQ data
qq_list <- list()


tic()
for (i in seq_along(rand_tiles)) {  
  cat("Processing Tile:", i, "\n")
  
  # Standardize Observations (tile i)
  obs <- observationList[i, ,]
  obs <- (obs - median(obs)) / IQR(obs)
  qObs <- quantile(obs, p = p_vec)
  
  # Store simulated quantiles across all R replications
  qSim_mat <- matrix(NA,
                     nrow = n_points,
                     ncol = 1000)
  
  for (r in 1:1000) {  
    cat("Processing Realization:", r, "\n")
    
    # Standardize Simulated Data (tile i, realization r)
    sim <- biasCorExtremalFields[r, i, ,]
    sim <- (sim - median(sim)) / IQR(sim)
    qSim_mat[, r] <- quantile(sim, p = p_vec)
  }
  
  # Compute CI bounds (2.5% and 97.5% percentiles)
  qSim_median <- apply(qSim_mat, 1, median, na.rm = TRUE)
  qSim_lower  <- apply(qSim_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  qSim_upper  <- apply(qSim_mat, 1, quantile, probs = 0.975, na.rm = TRUE)

  # # # Compute Mean and SD-Based Confidence Intervals
  # qSim_mean  <- apply(qSim_mat, 1, mean, na.rm = TRUE)  # Compute mean
  # qSim_sd    <- apply(qSim_mat, 1, sd, na.rm = TRUE)    # Compute standard deviation
  # 
  # # Compute 95% confidence intervals using normal approximation
  # qSim_lower <- qSim_mean - 1.96 * qSim_sd
  # qSim_upper <- qSim_mean + 1.96 * qSim_sd

  # Store in data frame
  qq_df <- data.frame(
    Observed = qObs,
    Simulated_Median = qSim_median,
    Simulated_Lower  = qSim_lower,
    Simulated_Upper  = qSim_upper,
    Tile = paste("Tile", rand_tiles[i])
  )
  
  # Append to list
  qq_list[[i]] <- qq_df
}
toc()
# time:  2.78555 mins



# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
# save( qq_list, file='qq_list_ObsVsSim.RData')

for (i in seq_along(rand_tiles)) {  # Loop over tiles
  cat('i: ', i, '\n')
  plot(qq_list[[i]]$Observed ,
       qq_list[[i]]$Simulated_Median,
       pch=19,
       cex=0.6,
       #lwd=2,
       ylim=c(-4, 14),
       #type='l',
       xlab='Observed',
       ylab='Simulated',
       main=paste0('Tile', rand_tiles[i]),
       col=alpha('blue', 0.8))
  abline(0, 1, col="pink3", lwd=3, lty=2)
  lines(qq_list[[i]]$Observed,
       qq_list[[i]]$Simulated_Lower,
       col='grey',
       lwd=2,
       lty=2)
  lines(qq_list[[i]]$Observed,
        qq_list[[i]]$Simulated_Upper,
        col='grey',
        lwd=2,
        lty=2)
  
  # # Create QQ-Plot with Confidence Intervals
  # p <- ggplot(qq_list[[i]], aes(x = Observed, y = Simulated_Median )) +
  #   geom_ribbon(aes(ymin = Simulated_Lower, ymax = Simulated_Upper),
  #               fill = "lightblue", linewidth=2, alpha = 0.4) +
  #   geom_point(color = "blue", size = 2.5, alpha = 0.8) +
  #  # Confidence Interval
  #   geom_abline(slope = 1, intercept = 0, linetype = "dashed",
  #               color = "pink3", linewidth =2) +
  #   labs(title = paste("QQ-Plot for Tile", rand_tiles[i]),
  #        x = "Observed Quantiles" ,
  #        y = "Simulated Quantiles (Median & CI)") +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  #     axis.title = element_text(size = 15),
  #     axis.text = element_text(size = 15)
  #   )

  # Display the plot
 #  print(p)
}



for (i in seq_along(rand_tiles)) {  # Loop over tiles
  cat('i: ', i, '\n')
  
  # Create QQ-Plot with Confidence Intervals
  p <- ggplot(qq_list[[i]], aes(x = Observed, y = Simulated_Median)) +
    # Confidence Interval Ribbon
    geom_ribbon(aes(ymin = Simulated_Lower, ymax = Simulated_Upper),
                fill = "lightblue", alpha = 0.4) +
    
    # Scatter points (Observed vs Simulated)
    geom_point(color = "blue", size = 2, alpha = 0.8) +
    
    # Confidence Interval Lines (Grey for Lower and Upper Bounds)
    geom_line(aes(y = Simulated_Lower), color = "grey", linetype = "dashed", linewidth = 0.5) +
    geom_line(aes(y = Simulated_Upper), color = "grey", linetype = "dashed", linewidth = 0.5) +
    
    # 1:1 Reference Line (Magenta)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "pink3", linewidth = 2) +
    
    # Labels and Themes
    labs(title = paste("QQ-Plot for Tile", rand_tiles[i]),
         x = "Observed Quantiles",
         y = "Simulated Quantiles (Median & CI)") +
    ylim(-4, 14)+
    #xlim(-4, 14)+
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 21, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18)
    )
  
  # Display the plot
  print(p)
}


