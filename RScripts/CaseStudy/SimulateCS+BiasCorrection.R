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
# 166 non-overlapping tiles

# ## convert the list into array:
# RCMannualMax16x16 <- array(NA, 
#                            dim=c(length(annualMaxima16x16Boxes), 
#                                  16, 16, 30))
# for(i in 1:length(annualMaxima16x16Boxes)){
#   
#   RCMannualMax16x16[i, , , ] <- (annualMaxima16x16Boxes[[i]])[, ,1:30]
#   
# }

centerLocs <- matrix(NA,
                     nrow=length(annualMaxima16x16Boxes),
                     ncol=2)
for(i in 1:length(annualMaxima16x16Boxes))
{
  temp <- gridPoints16x16Boxes[[i]] 
  centerLocs[i, ] <- c(temp$center)
}


# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
# load('RCMannualMax16x16.RData')
# dim(RCMannualMax16x16)


## Loading the bias-corrected parameter of the CNNs
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
load('bias-correct-CNN-CS-Parameters.RData')
dim(parameters)

# save(extremalFieldsList,
#      file='full-extremalFieldsList-CS.RData')

### Simulating the data across tiles:
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

## Simulating the process:
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/simulateSpatialExtFields.R")
R <- 400
m <- 30  # Number of ensemble members

# extremalFieldsList <- array(NA,
#                             dim=c(R , length(annualMaxima16x16Boxes),  256, m))
# tilesList <- array(NA,
#                    dim=c(R , length(annualMaxima16x16Boxes), 256, 2))
# 
# centerList <- array(NA,
#                     dim=c(R , length(annualMaxima16x16Boxes), 2, 2))
# tic()
# for (r in 1:R) {
#   tic()
#   cat('r: ', r, '\n')
#   for (i in 1:length(annualMaxima16x16Boxes)) {  # Use seq_along for safe indexing
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
#     tilesList[r, i, , ] <- s_tile  
#     centerList[r, i, , ] <- c(gridPoints16x16Boxes[[i]]$center)
# 
#     extremalFieldsList[r, i, , ] <- generateSpatialExtFields(s_tile,
#                                                              m,
#                                                              parameters[i, , drop = FALSE])
#   }
#   toc()
# }
# toc()
# # 8.3 hrs
# dim(extremalFieldsList)


## Save the data:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
load('full-extremalFieldsList-CS.RData')
dim(extremalFieldsList)

# save(extremalFieldsList,
#      file='full-extremalFieldsList-CS.RData')
# save(extremalFieldsList,
#      file='extremalFieldsList-CS.RData')
# load('extremalFieldsList-CS.RData')


# save(tilesList,
#      file='full-tileList-CS.RData')
load('full-tileList-CS.RData')
dim(tilesList)

# save(centerList,
#      file='full-centerList-CS.RData')
load('full-centerList-CS.RData')
dim(centerList)

### -- Quantile Mapping for bias-correction -- :
n_yrs <- dim(observationList)[3]
m <- dim(observationList)[3]
p_vec <- seq(0.01, 0.99, length.out=100)

# # Initialize bias-corrected array
# biasCorExtremalFieldsI <- array(NA,
#                                dim=c(R,
#                                      66,
#                                      256,
#                                      m))
# 
# tic()
# for (i in 1:66) {  # Avoid potential empty index issues
#   cat('Processing tile:', i, '\n')
# 
#   for (t in 1:m) {  # Use seq_len() for safe iteration
#     cat('Time:', t, '\n')
#     qObs <- as.numeric(quantile(observationList[i, , t],
#                                 p=p_vec,
#                                 na.rm=TRUE))
# 
#     for (r in seq_len(R)) {  # Use seq_len() instead of 1:R
#       cat('Repetition:', r, '\n')
#       qSim <- as.numeric(quantile(extremalFieldsList[r, i, , t],
#                                   p=p_vec,
#                                   na.rm=TRUE))
# 
#       # Check for degenerate cases (e.g., constant values in qSim)
#       fitQM <- lm(qObs ~ qSim)
#       biasCorExtremalFieldsI[r, i, , t] <- predict(fitQM,
#                                                    newdata = data.frame(qSim = extremalFieldsList[r, i, , t]))
#     }
#   }
# }
# toc()
# dim(biasCorExtremalFieldsI)
# elapsed time: 9 mins

# # Initialize bias-corrected array
# biasCorExtremalFieldsII <- array(NA,
#                                 dim=c(R,
#                                       50,
#                                       256,
#                                       m))
# 
# tic()
# for (i in 1:50) {  # Avoid potential empty index issues
#   cat('Processing tile:', i+66, '\n')
#   
#   for (t in 1:m) {  # Use seq_len() for safe iteration
#     cat('Time:', t, '\n')
#     qObs <- as.numeric(quantile(observationList[i+66, , t],
#                                 p=p_vec,
#                                 na.rm=TRUE))
#     
#     for (r in seq_len(R)) {  # Use seq_len() instead of 1:R
#       cat('Repetition:', r, '\n')
#       qSim <- as.numeric(quantile(extremalFieldsList[r, i+66, , t],
#                                   p=p_vec,
#                                   na.rm=TRUE))
#       
#       # Check for degenerate cases (e.g., constant values in qSim)
#       fitQM <- lm(qObs ~ qSim)
#       biasCorExtremalFieldsII[r, i, , t] <- predict(fitQM,
#                                                     newdata = data.frame(qSim = extremalFieldsList[r, i+66, , t]))
#     }
#   }
# }
# toc()
# dim(biasCorExtremalFieldsII)
# # 7.42605 mins

# Initialize bias-corrected array
biasCorExtremalFieldsIII <- array(NA,
                                  dim=c(R,
                                        50,
                                        256,
                                        m))

tic()
for (i in 1:50) {  # Avoid potential empty index issues
  cat('Processing tile:', i+(66+50), '\n')
  
  for (t in 1:m) {  # Use seq_len() for safe iteration
    cat('Time:', t, '\n')
    qObs <- as.numeric(quantile(observationList[i+(66+50), , t],
                                p=p_vec,
                                na.rm=TRUE))
    
    for (r in seq_len(R)) {  # Use seq_len() instead of 1:R
      cat('Repetition:', r, '\n')
      qSim <- as.numeric(quantile(extremalFieldsList[r, i+(66+50), , t],
                                  p=p_vec,
                                  na.rm=TRUE))
      
      # Check for degenerate cases (e.g., constant values in qSim)
      fitQM <- lm(qObs ~ qSim)
      biasCorExtremalFieldsIII[r, i, , t] <- predict(fitQM, 
                                                     newdata = data.frame(qSim = extremalFieldsList[r, i+(66+50), , t]))
    }
  }
}
toc()
dim(biasCorExtremalFieldsIII)
## 5.64085 mins

setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
# load('biasCorExtremalFields-CS-I.RData')
# dim(biasCorExtremalFieldsI)
# 
# load('biasCorExtremalFields-CS-II.RData')
# dim(biasCorExtremalFieldsII)

save(biasCorExtremalFieldsIII,
     file='biasCorExtremalFields-CS-III.RData')

load('biasCorExtremalFields-CS-I.RData')
load('biasCorExtremalFields-CS-II.RData')
load('biasCorExtremalFields-CS-III.RData')