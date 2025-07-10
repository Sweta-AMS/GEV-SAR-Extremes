# rm(list=ls())

## -- Required Package library -- :
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Proj2-Data2Evaluate.R")

## -- Loading the true parameters for both training and testing data set --:
shape_true_train <- trainParameters[ , 1]
kappa2_true_train <- (trainParameters[ , 2]-4)
tau2_true_train <- (trainParameters[ , 3]) # 50,000

shape_true_test <- testParameters[ , 1]
kappa2_true_test <- (testParameters[ , 2]-4)
tau2_true_test <- (testParameters[ , 3]) # 10,000

###  -- Simulation study to evaluate results and compare behaviour across the replication size --: 
# Expand the file path
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs")
file_path_30 <- path.expand("test_estimates_NGP-n-30-Latest-Jan22-2025.npy")
file_path_10 <- path.expand("test_estimates_NGP-n-10-Latest-Jan22-2025.npy")
file_path_1 <- path.expand("test_estimates_NGP-n-1-Latest-Jan22-2025.npy")

# Import numpy and load the file
np <- import("numpy")
mat30 <- np$load(file_path_30)
mat10 <- np$load(file_path_10)
mat1 <- np$load(file_path_1)

# Load the NumPy file
pred_parameterRep30 <- as.matrix(mat30)
dim(pred_parameterRep30)

pred_parameterRep10 <- as.matrix(mat10)
dim(pred_parameterRep10)

pred_parameterRep1 <- as.matrix(mat1)
dim(pred_parameterRep1)

##  --- Bias-Correction ---:
# Shape parameter: 
biasXi30 <- (pred_parameterRep30[,1] - shape_true_test)
biasXi10 <- (pred_parameterRep10[,1] - shape_true_test)
biasXi1 <- (pred_parameterRep1[,1] - shape_true_test)

# log Kappa parameter:
biaslnKappa30 <- (0.5*pred_parameterRep30[,2] - 0.5*log(kappa2_true_test))
biaslnKappa10 <- (0.5*pred_parameterRep10[,2] - 0.5*log(kappa2_true_test))
biaslnKappa1 <- (0.5*pred_parameterRep1[,2] - 0.5*log(kappa2_true_test))

# log Kappa parameter:
biaslnTau30 <- (0.5*pred_parameterRep30[,3] - 0.5*log(tau2_true_test))
biaslnTau10 <- (0.5*pred_parameterRep10[,3] - 0.5*log(tau2_true_test))
biaslnTau1 <- (0.5*pred_parameterRep1[,3] - 0.5*log(tau2_true_test))

# Fit smoothing splines with different smoothing parameters:
## Xi
Xi_biasCorrect30 <- smooth.spline(x=shape_true_test,
                                  y=biasXi30,
                                  spar=0.9)

Xi_biasCorrect10 <- smooth.spline(x=shape_true_test,
                                  y=biasXi10,
                                  spar=0.9)

Xi_biasCorrect1 <- smooth.spline(x=shape_true_test,
                                 y=biasXi1,
                                 spar=0.9)
# Predict bias corrections
biasXI30 <- predict(Xi_biasCorrect30,
                    x=shape_true_test)$y
biasXI10 <- predict(Xi_biasCorrect10,
                    x=shape_true_test)$y
biasXI1 <- predict(Xi_biasCorrect1,
                   x=shape_true_test)$y

## log Kappa
lnKappa_biasCorrect30 <- smooth.spline(x=0.5*log(kappa2_true_test),
                                       y=biaslnKappa30,
                                       spar=0.9)

lnKappa_biasCorrect10 <- smooth.spline(x=0.5*log(kappa2_true_test),
                                       y=biaslnKappa10,
                                       spar=0.9)

lnKappa_biasCorrect1 <- smooth.spline(x=0.5*log(kappa2_true_test),
                                      y=biaslnKappa1,
                                      spar=0.9)
# Predict bias corrections
biaslnKappa30 <- predict(lnKappa_biasCorrect30,
                         x=0.5*log(kappa2_true_test))$y

biaslnKappa10 <- predict(lnKappa_biasCorrect10,
                         x=0.5*log(kappa2_true_test))$y
biaslnKappa1 <- predict(lnKappa_biasCorrect1,
                        x=0.5*log(kappa2_true_test))$y

## log Tau
lnTau_biasCorrect30 <- smooth.spline(x=0.5*log(tau2_true_test),,
                                     y=biaslnTau30,
                                     spar=0.9)

lnTau_biasCorrect10 <- smooth.spline(x=0.5*log(tau2_true_test),,
                                     y=biaslnTau10,
                                     spar=0.9)

lnTau_biasCorrect1 <- smooth.spline(x=0.5*log(tau2_true_test),,
                                    y=biaslnTau1,
                                    spar=0.9)
# Predict bias corrections
biaslnTau30 <- predict(lnTau_biasCorrect30,
                       x=0.5*log(tau2_true_test))$y

biaslnTau10 <- predict(lnTau_biasCorrect10,
                       x=0.5*log(tau2_true_test))$y
biaslnTau1 <- predict(lnTau_biasCorrect1,
                      x=0.5*log(tau2_true_test))$y

