## Simulate the LK fields and compare with matern gaussian field: 
## Working on adding the addition channel for the standardized Gaussian marginal for 
rm(list=ls())

## Required script
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")

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
n_sims <- 10000 # changed 
print(paste0('full size: ', n_sims))


## -- Parameter configuration --:
set.seed(222)
AWghtGrid <- 4 + exp(runif(n=n_sims, log(0.01), log(1.8))) #Training range: 0.001-2
summary(AWghtGrid)
hist(AWghtGrid)

lambdaGrid <- 10^(runif(n=n_sims, log10(0.0001), log10(0.05)))  #Training range: 0.00001-0.1
summary(lambdaGrid)
hist(lambdaGrid)

shapeGrid <- runif(n=n_sims, (0.2), (0.8))  #Training range: 0.01-0.9
summary(shapeGrid)
hist(shapeGrid)

parameterComb <- cbind(shapeGrid, AWghtGrid, lambdaGrid)
dim(parameterComb)

## Sanity Check:
plot(parameterComb[,1], parameterComb[,2])
plot(parameterComb[,1], parameterComb[,3])
plot(parameterComb[,2], parameterComb[,3])

# ## -- Parameter configuration --:
# # Load the training parameter configuration
# trainParameters <-load("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/Data/ParameterConfiguration-Jan17-Train-Rep30.RData")
# trainParameters <- storeParameterRep
# dim(trainParameters) # 50000 3
# 
# # Load the test parameter configuration
# testParameters <-load("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/Data/ParameterConfiguration-Jan14-Test-Rep30-Latest.RData")
# testParameters <- storeParameterRep
# dim(testParameters) # 10000 3

# plot(testParameters[,1], 
#      testParameters[,2],
#      pch=19,
#      col=alpha(alpha=0.4, 'black'))
# plot(testParameters[,1],
#      testParameters[,3],
#      pch=19,
#      col=alpha(alpha=0.4, 'black'))
# plot(testParameters[,2],
#      testParameters[,3],
#      pch=19,
#      col=alpha(alpha=0.42, 'blak'))

### -- GENERATE TRAINING SET --:
# For replication size 30:
source("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/RScripts/LKrigSAREvd.R")
m <- 30
extremalFields <- array(NA,
                        dim=c(n_sims, n, m))
dim(extremalFields)

tic()  # Start timing the entire process
# Initialize result_list before starting the loop
set.seed(123)
for (i in 1:n_sims)
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
  
  # Store the scaled matrix in the extremalFields array
  extremalFields[i, , ] <- Z 
}
toc()  # End timing for the entire process

plot(parameterComb[,1], parameterComb[,2])
plot(parameterComb[,2], parameterComb[,3])
plot(parameterComb[,1], parameterComb[,3])

## Defining test set: save(storeParameterRepTrain, file="ParameterGrid-Frechet-Train.RData")
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/Data")
dim(parameterComb)
save(parameterComb, 
     file="ParameterConfiguration-Jan22-Test-Rep30-Latest.RData")


dim(extremalFields)
save(extremalFields, file="storeZRep-Jan22-Test-Rep30-Latest.RData")

t <- 10000
image.plot(matrix(extremalFields[t, ,20], 16, 16))
parameterComb[t,]

# Compute the empirical variogram
vgram_result <- vgram(s, extremalFields[1, ,], dmax=10)

# Plot the empirical variogram
plot(vgram_result)


