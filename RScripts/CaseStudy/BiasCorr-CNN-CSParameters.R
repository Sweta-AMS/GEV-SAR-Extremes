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
                     nrow = 166,
                     ncol = 3)
for (i in 1:166) {
  # Extract GEV parameters for the i-th tile
  parameters[i, ] <- c(bcXiCNN[i],  (exp(bclogKappaCNN[i])+4), (exp(bclogTauCNN[i])) )  # Replace with actual parameter extraction
}
dim(parameters)

tau_vec <- parameters[ ,3]
tau_vec[tau_vec>0.1] <- 0.1
parameters[,3] <- tau_vec
summary(parameters)

#