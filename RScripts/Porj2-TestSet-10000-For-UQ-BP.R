## Simulate the LK fields and compare with matern gaussian field: 
## Working on adding the addition channel for the standardized Gaussian marginal for 
rm(list=ls())

## Required script
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")

# source("~/Downloads/LKrigSAREvd.R")

## Location
M <- 16 # changed from 15
sGrid <- list(x=1:M, y=1:M)
s <- make.surface.grid(sGrid)
plot(s)

## Generating the training set    
## training grid size 50 and validation is 20, let's keep testing 10
n <- nrow(s)
print(paste0('Field size: ', n))

## Latest used
n_full <- 1000 # changed 
print(paste0('full size: ', n_full))


## -- Parameter configuration --:
# Test parameter configuration
set.seed(111)
AWghtGrid <- 4 + exp(runif(n=n_full, log(0.005), log(1.5))) # (0.001, 2)
summary(AWghtGrid)
hist(AWghtGrid)

lambdaGrid <- exp(runif(n=n_full, log(0.00001), log(0.05))) # (0.00001, 0.099)
summary(lambdaGrid)
hist(lambdaGrid)

shapeGrid <- runif(n=n_full, (0.2), (0.8)) # (0.01, 0.9)
summary(shapeGrid)
hist(shapeGrid)

parameterComb <- cbind(shapeGrid, AWghtGrid, lambdaGrid)
dim(parameterComb)

plot(parameterComb[,1], parameterComb[,2])
plot(parameterComb[,1], parameterComb[,3])
plot(parameterComb[,2], parameterComb[,3])

### -- GENERATE TRAINING SET --:

# For replication size 100:
m <- 30
R <- 100
extremalFields <- array(NA,
                        dim=c(n_full, R, n, m))
dim(extremalFields)

tic()  # Start timing the entire process
# Initialize result_list before starting the loop
set.seed(123)
for(i in 1:n_full)
{
  cat('GEV parameter loop no', i, '\n')
  
  # Defining GEV parameters:
  shape_val <- shapeGrid[i]
  cat('Xi: ', shape_val, '\n')
  
  a_wght <- AWghtGrid[i] 
  cat('Kappa^2 : ', a_wght, '\n')
  
  lambda_val <- lambdaGrid[i]
  cat('lambda^2 : ', lambda_val, '\n')
  
  # Set a seed for reproducibility, unique for each iteration
  print(paste("Loop", i))
  
  # Setup LKrig
  LKinfo <- LKrigSetup(s,
                       a.wght=a_wght,
                       nlevel=1,
                       nu=1,
                       NC=16,  # Changed from 8
                       NC.buffer=4,
                       fixedFunction=NULL)
  
  # Simulate SAR coefficients
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
  cat('Dim of ySim: ', dim(ySim), '\n')
  
  # Add measurement error
  scaleParameter <- log(1 + (lambda_val))  # Variance of the log-normal is tau2
  locParameter <- (-scaleParameter / 2)
  
  # Set another seed for the log-normal matrix generation, unique for each iteration
  LogNormalMat <- exp(matrix(rnorm(n*m,
                                   mean=locParameter,
                                   sd=sqrt(scaleParameter)),
                             n, m))  # n x m matrix
  
  dim(LogNormalMat)
  
  # Scale each column (if required by your context)
  Z <- ySim*LogNormalMat
  Z_scaled <- scale_matrix(Z)
  
  # Store the scaled matrix in the extremalFields array
  extremalFields[i, , ] <- Z # removed Z_scaled
  extremalFieldsScl[i, , ] <- Z_scaled
}
toc()  # End timing for the entire process

# Shuffle index in R
set.seed(555)
index_shuffle <- sample(n_full)
index_shuffle[1:10]

storeZRep <- extremalFields
dim(storeZRep)

# Parameter
storeParameterRep <- parameterComb
dim(storeParameterRep)

storeZRep <- storeZRep[index_shuffle, , ]
storeParameterRep <- storeParameterRep[index_shuffle, ]
summary(storeParameterRep)
summary(storeZRep)

## Shuffling the dataset:
dim(storeZRep)
dim(storeParameterRep)

plot(storeParameterRep[,1], storeParameterRep[,2])
plot(storeParameterRep[,2], storeParameterRep[,3])
plot(storeParameterRep[,1], storeParameterRep[,3])

## Defining test set: save(storeParameterRepTrain, file="ParameterGrid-Frechet-Train.RData")
dim(storeParameterRep)
save(storeParameterRep, 
     file="ParameterConfiguration-Jan14-Test-Rep30-Latest.RData")


dim(storeZRep)
save(storeZRep, file="storeZRep-Jan14-Test-Rep30-Latest.RData")

t <- 1600
image.plot(matrix(storeZRep[t, ,10], 16, 16))
storeParameterRep[t,]

# Compute the empirical variogram
vgram_result <- vgram(s, storeZRep[1, ,], dmax=10)

# Plot the empirical variogram
plot(vgram_result)



# ## addition sample:
# ## -- Parameter configuration --:
# n_adj <- 1000
# 
# set.seed(123)
# AWghtGridAdj <- 4 + exp(runif(n=n_adj, log(0.001), log(2)))
# summary(AWghtGridAdj)
# hist(AWghtGridAdj)
# 
# set.seed(123)
# lambdaGridAdj <- exp(seq(log(0.0001), log(0.01), length.out=n_adj))
# summary(lambdaGridAdj)
# hist(lambdaGridAdj)
# 
# set.seed(222)
# shapeGridAdj <- runif(n=n_adj, 0.1, 1)
# summary(shapeGridAdj)
# hist(shapeGridAdj)
# 
# parameterCombAdj <- cbind(shapeGridAdj, AWghtGridAdj, lambdaGridAdj)
# dim(parameterCombAdj)
# 
# 
# # For replication size 10:
# m <- 30
# extremalFieldsAdj <- array(NA, dim=c(n_adj, n, m))
# dim(extremalFieldsAdj)
# 
# tic()
# # Initialize result_list before starting the loop
# for(i in 1:length(shapeGridAdj))
# {
#   tic()
#   cat('GEV parameter loop no', i, '\n')
#   
#   # Defining GEV parameters:
#   shape_val <- shapeGridAdj[i]
#   a_wght <- AWghtGridAdj[i] 
#   lambda_val <- lambdaGridAdj[i]
#   
#   set.seed(123 + i ) 
#   print(paste("Loop", i))
#   
#   # Setup LKrig
#   LKinfo <- LKrigSetup(s,
#                        a.wght=a_wght,
#                        nlevel=1,
#                        nu=1,
#                        NC=10,
#                        fixedFunction=NULL)
#   
#   gridCenter <- LKrigLatticeCenters(LKinfo, 
#                                     Level=1)
#   
#   set.panel(1,1)
#   plot(make.surface.grid( gridCenter))
#   points(s, col=alpha('red',0.3), pch=19)
#   
#   # Simulate SAR coefficients
#   coefSARList <- LKrigSAREvd(LKinfo,
#                              loc=0,
#                              scale=1,
#                              shape=shape_val,
#                              M=m,
#                              asList=FALSE)
#   
#   coefSAR <- coefSARList$coefSAR
#   
#   # Simulate field
#   PHI <- LKrig.basis(s, LKinfo)
#   ySim <- (PHI %*% coefSAR)
#   cat('Dim of ySim: ',  dim(ySim), '\n')
#   
#   set.seed(111 + i)
#   UMat <- matrix(runif(n*m), n, m) # 256*200
#   cat('Dim of UMat : ', dim(UMat), '\n')
#   
#   # Transform to Gumbel: Mu = 0 and sigma = 1
#   standardGumbelMat <- -log(-log(UMat))
#   lambda_values <- rep(lambda_val, times = m)
#   
#   # Add measurement error
#   calScale <- sqrt(6*lambda_values/pi^2)
#   euler_constant <- 0.57721
#   calLoc <- -(calScale)*euler_constant
#   
#   # Scaling each column
#   WN <- sweep(standardGumbelMat, 2, calScale, "*")  
#   WN <- sweep(WN, 2, calLoc, "+")  
#   Z <- ySim + WN
#   Z_scaled <- scale_matrix(Z)
#   toc()
#   extremalFieldsAdj[i, ,] <- Z_scaled
# }
# toc()
# 
# # Spatial Field
# storeZRepAdj <- extremalFieldsAdj
# dim(storeZRepAdj)
# 
# # Parameter
# storeParameterRepAdj <- parameterCombAdj
# dim(storeParameterRepAdj)
# 
# # Shuffle index in R
# set.seed(111)
# # n_full <- 122422
# index_shuffle_adj <- sample(n_adj)
# index_shuffle_adj[1:100]
# 
# storeZRepAdj <- storeZRepAdj[index_shuffle_adj, ,]
# storeParameterRepAdj <- storeParameterRepAdj[index_shuffle_adj, ]
# image.plot(matrix(storeZRepAdj[10, ,15], 32,32))
# storeParameterRepAdj[10,]
# save(storeParameterRepAdj, 
#      file="ParameterConfiguration-Adj-PositiveShape-RandomlySample-LatestOct16.RData")
# dim(storeParameterRepAdj)
# 
# save(storeZRepAdj, file="storeZRep30-Adj-PositiveShape-RandomlySample-LatestOct16.RData")
# dim(storeZRepAdj)
# 
# # ## Now combining the adj with the training data:
# # library(abind)
# # 
# # storeZRepTrain <- abind(storeZRep, storeZRepAdj, along = 1) 
# # dim(storeZRepTrain)
# 
# ## 
# 
# # Load necessary packages
# # install.packages("stabledist") # Uncomment if stabledist package is not installed
# # install.packages("evd")        # Uncomment if evd package is not installed
# library(stabledist)
# library(evd)
# 
# # Set parameters for both distributions
# alpha <- 1.5
# c <- 1
# mu <- 0
# sigma <- 1
# 
# # Define x values for plotting
# x <- seq(0.1, 10, length.out = 100)
# 
# # Calculate density values
# positive_stable_density <- dstable(x, alpha, beta = 1, gamma=c, delta=0, pm=0)
# frechet_density <- dgev(x, loc = mu, scale = sigma, shape = 1 / alpha)
# 
# # Plot the densities
# plot(x, positive_stable_density, type = "l", col = "blue", lwd = 2,
#      main = "Positive Stable vs. Fréchet Density",
#      xlab = "x", ylab = "Density")
# lines(x, frechet_density, col = "red", lwd = 2)
# legend("topright", legend = c("Positive Stable", "Fréchet"), col = c("blue", "red"), lwd = 2)
# 
# 
# 
# 
# load("storeZRep30-PositiveShape-RandomlySample-LatestOct29-Train.RData")
# load("ParameterConfiguration-PositiveShape-RandomlySample-LatestOct29-Train.RData")

