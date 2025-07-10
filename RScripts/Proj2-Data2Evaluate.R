## this script is to load the training and test parameter configuration-
## both estimated and true for comparison study:
# rm(list=ls())
## Required script
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")

## -- Parameter configuration --:
# Load the training parameter configuration
trainParameters <- load("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Data/ParameterConfiguration-Jan17-Train-Rep30.RData")
trainParameters <- storeParameterRep
dim(trainParameters) # 50000 3

# Load the test parameter configuration
testParameters <-load("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Data/ParameterConfiguration-Jan22-Test-Rep30-Latest.RData")
testParameters <- parameterComb
dim(testParameters) # 10000 3


## Calling the training estimates:
file_path_train <- path.expand("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/trainEstimates-n-30-Latest-Jan17-2025.npy") # Estimates over the train set

## -- Import numpy and load the file -- ##
np <- import("numpy")
matTrain <- np$load(file_path_train)

## -- Load the NumPy file -- ##
cnnEstTrain <- matTrain
dim(cnnEstTrain)  # 100x100x3

## -- Bias correction for shape -- ##
shape_pred_train <- cnnEstTrain[, 1]
kappa2_pred_train <- exp(cnnEstTrain[, 2])
tau2_pred_train <- exp(cnnEstTrain[, 3])


## Calling the test estimates
file_path_test <- path.expand("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/test_estimates_NGP-n-30-Latest-Jan22-2025.npy") 


## -- Import numpy and load the file -- ##
np <- import("numpy")
matTest <- np$load(file_path_test)

## -- Load the NumPy file -- ##
cnnEstTest <- matTest
dim(cnnEstTest)  # 100x100x3

## -- Bias correction for shape -- ##
shape_pred_test <- cnnEstTest[, 1]
kappa2_pred_test <- exp(cnnEstTest[, 2])
tau2_pred_test <- exp(cnnEstTest[, 3])

