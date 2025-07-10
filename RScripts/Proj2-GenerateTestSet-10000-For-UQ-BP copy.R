## Simulate the LK fields and compare with matern gaussian field: 
## Working on adding the addition channel for the standardized Gaussian marginal for 
rm(list=ls())

## Required script
source("~/Desktop/Tentative Results-Proj2/R scripts/requiredPackages.R")
source("~/Downloads/LKrigSimSAR.R")
source("~/LKrigSARPareto.R")
source("~/Desktop/Tentative Results-Proj2/R script/LKrigSAREvd.R")

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
n_full <- 10000 # changed 
print(paste0('full size: ', n_full))


## -- Parameter configuration --:
# Test parameter configuration
set.seed(111)
AWghtGrid <- 4 + exp(runif(n=n_full, log(0.001), log(2)))
summary(AWghtGrid)
hist(AWghtGrid)

lambdaGrid <- exp(runif(n=n_full, log(0.00001), log(0.099)))
summary(lambdaGrid)
hist(lambdaGrid)

shapeGrid <- runif(n=n_full, (0.1), (0.9))
summary(shapeGrid)
hist(shapeGrid)

parameterComb <- cbind(shapeGrid, AWghtGrid, lambdaGrid)
dim(parameterComb)

plot(parameterComb[,1], parameterComb[,2])
plot(parameterComb[,1], parameterComb[,3])
plot(parameterComb[,2], parameterComb[,3])

# n_train <-  50000
# n_valid <-  10000

# n_train <- n_full*0.80
# n_valid <- n_full*0.20

# print(paste0('n_train: ', n_train))
# print(paste0('n_valid: ', n_valid))


### -- GENERATE TRAINING SET --:
# Function to use a min-max scaling
scale_matrix <- function(matrix_temp)
{
  # Step 1: Compute column-wise minimum and maximum
  column_min <- apply(matrix_temp, 2, min)
  column_max <- apply(matrix_temp, 2, max)
  
  # Step 2: Compute the difference between max and min for each column
  column_max_min <- column_max - column_min
  
  # Step 3: Subtract the column-wise minimum from each element in the matrix
  matrix_temp_scale <- sweep(matrix_temp, 2, column_min, "-")
  
  # Step 4: Divide each column by its respective max-min range
  matrix_temp_scaled <- matrix_temp_scale / matrix(rep(column_max_min, each = nrow(matrix_temp)), nrow = nrow(matrix_temp))
  
  # Return the final scaled matrix
  return(matrix_temp_scaled)
}

# For replication size 10:
source("~/Desktop/Tentative Results-Proj2/R script/LKrigSAREvd.R")
m <- 30
extremalFields <- array(NA,
                        dim=c(n_full, n, m))
extremalFieldsScl <- array(NA,
                           dim=c(n_full, n, m))
dim(extremalFields)
dim(extremalFieldsScl)

tic()  # Start timing the entire process
# Initialize result_list before starting the loop
set.seed(123)
for (i in 1:length(shapeGrid))
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

