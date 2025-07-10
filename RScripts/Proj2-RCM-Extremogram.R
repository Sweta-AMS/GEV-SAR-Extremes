rm(list=ls())
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Read-RCM-Data.R")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Madogram-similar2fields.R")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/CaseStudy/RCM-dataLoad.R")

### -- MADOGRAM COMPUTATION -- ###
## Calling out the annual_maxima's:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
load('AnnualMaxima-acrossNA.RData')
dim(annual_maxima) # 297 281  31

annual_precip_max <- array(annual_maxima, 
                           dim=c(dim(annual_maxima)[1] * dim(annual_maxima)[2] ,
                                 dim(annual_maxima)[3]))
dim(annual_precip_max) # 83457x31

# Compute pairwise distances:
lat_flat <- as.vector(lat_data)
print(length(lat_flat))

lon_flat <- as.vector(lon_data)
print(length(lon_flat))

## -- SUBSETTING -- ##
subset_indices <- sample(1:length(lat_flat),
                         1000)  

lat_subset <- lat_flat[subset_indices]
lon_subset <- lon_flat[subset_indices]

annual_precip_subset <- annual_precip_max[subset_indices, ]
dim(annual_precip_subset) # 10000    31

coords_subset <- cbind(lon_subset, 
                       lat_subset)
dim(coords_subset) # 10000     2

dist_matrix <- as.matrix(dist(coords_subset))  # Compute spatial distances
dim(dist_matrix)


### Calling the madogram functions: 
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Madogram-similar2fields.R")
mgram <- madogram(coords_subset,
                  annual_precip_subset[,1])

################### Compute the madogram ######################################:
# Define bin centers (lag distances)
bin_centers <- seq(0,
                   max(dist_matrix),
                   length.out=20)  # 15 distance bins


epsilon <- diff(bin_centers)[1] / 2  # Bin width (half of difference between bins)

# Initialize storage for variogram values
madogram_values <- numeric(length(bin_centers))
length(madogram_values)

# Loop over bins to compute semivariance
madogramMAT <- matrix(NA, nrow=dim(annual_precip_subset)[2],ncol= length(bin_centers))
for(i in 1:31){
  cat('Outer Loop: ', i, '\n' )
  for (b in seq_along(bin_centers)) {
    cat('Inner Loop: ', b, '\n' )
    d <- bin_centers[b]  # Current bin center
  
    # Identify pairs where distance falls within bin range
    bin_mask <- (dist_matrix >= (d - epsilon)) & (dist_matrix < (d + epsilon))
  
    # Get indices of pairs in this bin
    bin_indices <- which(bin_mask,
                         arr.ind=TRUE)
  
    if (nrow(bin_indices) > 0)
    {
      # Compute squared differences of data values for selected pairs
      data <- (annual_precip_subset[,i] - median(annual_precip_subset[,i]))/sd(annual_precip_subset[,i])
      gamma_values <- abs(data[bin_indices[,1]] - data[bin_indices[,2]])
    
      # Compute variogram estimate
      madogram_values [b] <- mean(gamma_values)
    }
  }
  
  madogramMAT[i, ] <- madogram_values
}
bplot.xy(bin_centers, median(madogramMAT))
bplot.xy(bin_centers, median(madogramMAT))
################################################################################


### -- EXTREMOGRAM COMPUTATION -- ###

library(fields)  # For distance computation

# Ensure coordinates are in the correct format
coords_S <- as.matrix(coords_subset)  # (lon, lat)
dim(coords_S)

# Compute pairwise distances
dist_matrix <- rdist.earth(coords_S, miles=FALSE)

# Define fixed bin structure
bin_centers <- seq(0, max(dist_matrix), length.out=15)  # 15 distance bins
epsilon <- diff(bin_centers)[1] / 2  # Bin width (half of difference between bins)

# Initialize storage for extremogram
extremogramMAT <- matrix(NA, 
                         nrow=dim(annual_precip_subset)[2],
                         ncol=length(bin_centers))

# Define fixed threshold (e.g., 95th percentile)
threshold <- quantile(annual_precip_subset,
                      0.95,
                      na.rm=TRUE)  # Global threshold

# Loop over years
for (i in 1:31) {
  cat("Computing Extremogram for Year:", i, "\n")
  
  # Standardize precipitation data for year i
  data <- (annual_precip_subset[, i] - median(annual_precip_subset[, i], na.rm=TRUE)) /
    sd(annual_precip_subset[, i], na.rm=TRUE)
  
  threshold <- quantile(data,
                        0.95,
                        na.rm=TRUE)  # local threshold
  
  # Identify extreme events (1 if extreme, 0 otherwise)
  extreme_events <- data > threshold
  
  # Loop over bins (lag distances)
  for (b in seq_along(bin_centers)) {
    cat('Print')
    d <- bin_centers[b]  # Current bin center
    
    # Identify pairs where distance falls within bin range
    bin_mask <- (dist_matrix >= (d - epsilon)) & (dist_matrix < (d + epsilon))
    
    # Get indices of pairs in this bin
    bin_indices <- which(bin_mask, arr.ind=TRUE)
    
    if (nrow(bin_indices) > 0) {
      # Compute co-occurrence of extreme events for selected pairs
      extreme_pairs <- extreme_events[bin_indices[,1]] & extreme_events[bin_indices[,2]]
      
      # Compute extremogram value
      extremogramMAT[i, b] <- mean(extreme_pairs, na.rm=TRUE)
    }
  }
}

# Compute mean extremogram across years
extremogram_mean <- colMeans(extremogramMAT, na.rm=TRUE)

# Create a dataframe for plotting
extremogram_df <- data.frame(Distance = bin_centers, Extremogram = extremogram_mean)

# Plot the extremogram
library(ggplot2)

ggplot(extremogram_df, aes(x=Distance, y=Extremogram)) +
  geom_point(color="blue") +
  geom_line(color="blue") +
  labs(title="Spatial Extremogram of Extreme Precipitation",
       x="Lag Distance (km)",
       y="Extremogram Value") +
  theme_minimal()







# --- Extract annual maxima for a selected year (e.g., last year in dataset) ---
library(geoR)

# Ensure coordinates are in the correct format
coords_S <- as.matrix(coords_subset)  # Should be a matrix with (lon, lat)
dim(coords_S)

# Define fixed bin structure for all years
bin_lims <- seq(0, 
                max(dist(coords_S)),
                length.out = 25)  # Define 20 bins up to max distance

# Initialize list to store madograms
madogramLst <- list()
madogramMat <- list()

# Compute madogram for each year
for (i in 1:dim(annual_precip_subset)[2]) {
  cat("Computing Madogram for Year:", i,"\n")
  
  # Extract annual max precipitation & standardize
  values <- as.vector(annual_precip_subset[,i])
  values <- (values - median(values, na.rm=TRUE)) / sd(values, na.rm=TRUE) # length 10000
  
  # Compute empirical madogram:
  madogramLst[[i]] <- variog(coords=coords_S,
                             data=values,
                             uvec=bin_lims,
                             max.dist = 100,
                             estimator.type = "modulus")  # Madogram estimator
  madogramMat[[i]] <- madogramLst[[i]]$v
}
dist <-  madogramLst[[1]]$u

### Madogram Matrix
madogramMatrix <- matrix(NA, nrow=31, ncol=17)
for(i in 1:31)
{
  madogramMatrix [i, ] <- madogramMat[[i]]
}
dim(madogramMatrix)

## Print the first madogram to check results:
plot(dist,
     madogramMat[[30]],
     type="b", 
     pch=19,
     col="blue",
     main="Empirical Madogram",
     xlab="Lag Distance",
     ylab="Madogram")

## Madogram: binning across the boxplot
bplot.xy(dist,
         (madogramMatrix))
lines(smooth.spline(dist, apply(madogramMatrix, 2, median)),
      col="magenta",
      lwd=2)

abline(h=0.8, col='blue', lwd=2, lty=2)

################################################################################
#### -- COMPUTING THE EXTREMOGRAM --
extremogramMat <- matrix(NA,
                         nrow=31,
                         ncol=length(bin_lims)-1)
dim(extremogramMat)

# Loop over 31 years of data
for (i in 1:31)
{
  cat("Computing Extremogram for Year:", i, "\n")
  
  # Extract annual precipitation for year i
  values <- as.vector((annual_precip_subset[, i]- median(annual_precip_subset[, i]))/sd(annual_precip_subset[, i]))
  
  # Define extreme threshold (e.g., 95th percentile)
  threshold <- quantile(values,
                        0.85, 
                        na.rm=TRUE)
  
  # Identify extreme events (binary: 1 = extreme, 0 = non-extreme)
  extreme_events<- (values > threshold)
  
  # Compute extremogram for each distance bin
  for (j in 1:(length(bin_lims) - 1))
  {
    print(j)
    # Get indices of pairs within the current distance bin
    bin_mask <- (dist_matrix >= bin_lims[j] & dist_matrix < bin_lims[j+1])
    
    # Extract extreme event pairs for this bin
    extreme_pairs <- extreme_events[row(bin_mask)] & extreme_events[col(bin_mask)]
    
    # Compute extremogram value as P(X > u | Y > u)
    extremogramMat[i, j] <- mean(extreme_pairs, na.rm=TRUE)
  }
}

## Compute the average extremogram across all years
extremogram_mean <- colMeans(extremogramMat, na.rm=TRUE)

# Convert bins to numeric for plotting
bin_midpoints <- (bin_lims[-1] + bin_lims[-length(bin_lims)]) / 2

# Create a dataframe for visualization
extremogram_df <- data.frame(Distance = bin_midpoints,
                             Extremogram = extremogram_mean)

bplot.xy(bin_midpoints,
         extremogramMat)


## Generate samples:
m <- 30
extremalFields <- array(NA,
                        dim=c(n_sims, n, m))
dim(extremalFields)

tic()  # Start timing the entire process
# Initialize result_list before starting the loop
set.seed(123)
for(i in 1:nrow(parameterComb))
{
  # Defining GEV parameters:
  shape_val <- parameterComb[i, 1]
  cat('Xi: ', shape_val, '\n')
  
  a_wght <- parameterComb[i,2]
  cat('Kappa^2 : ', a_wght, '\n')
  
  lambda_val <- parameterComb[i,3]
  cat('lambda^2 : ', lambda_val, '\n')
  
  # Loop across the parameter configuration:
  cat('GEV parameter loop no', i, '\n')
  
  
  # Setup LKrig
  LKinfo <- LKrigSetup(s,
                       a.wght=a_wght,
                       nlevel=1,
                       nu=1,
                       NC=16,  # Changed from 8
                       NC.buffer=4,
                       fixedFunction=NULL)
  
  # Add measurement error
  scaleParameter <- log(1 + (lambda_val))  # Variance of the log-normal is tau2
  locParameter <- (-scaleParameter / 2)
  
  coefSARList <- LKrigSAREvd(LKinfo,
                               loc=1,
                               scale=shape_val,
                               shape=shape_val,
                               M=m,
                               asList=FALSE)
  coefSAR <- coefSARList$coefSAR
    
    # Simulate field
    PHI <- LKrig.basis(s, LKinfo)
    ySim <- (PHI %*% coefSAR)
    
    # Log-normal matrix generation, unique for each iteration
    LogNormalMat <- exp(matrix(rnorm(n*m,
                                     mean=locParameter,
                                     sd=sqrt(scaleParameter)),
                               n, m))  # n x m matrix
    
    # Spatial Extremal Fields with Nugget effect
    Z <- ySim*LogNormalMat
    extremalFields[i, j,  , ] <- Z
  }
  toc()
}
toc()  # End timing for the entire process
dim(extremalFields)
# 10 hrs 
t <- 2
image.plot(matrix(extremalFields[t, 10, ,1], 16, 16))
parameterComb[t,]


