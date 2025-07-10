## Confidence bound: Xi
## Required script
rm(list=ls())
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Proj2-Data2Evaluate.R")


shape_true_train <- trainParameters[ , 1]
kappa2_true_train <- (trainParameters[ , 2]-4)
tau2_true_train <- (trainParameters[ , 3])

shape_true_test <- testParameters[ , 1]
kappa2_true_test <- (testParameters[ , 2]-4)
tau2_true_test <- (testParameters[ , 3])


### --- Shape Parameter ---:
## Bias-Correction for Shape parameter: 
biasXi <- (shape_pred_train -shape_true_train)

# Fit smoothing splines with different smoothing parameters:
spline_Xi_biasCorrect <- smooth.spline(x=shape_true_train,
                                       y=biasXi,
                                       spar=0.9)

# Predict bias corrections
pred_biasXI <- predict(spline_Xi_biasCorrect,
                       x=shape_true_test)$y

# Plot residuals after smoothing
bplot.xy(shape_true_test,
         (shape_pred_test -shape_true_test) - pred_biasXI,
         pch = 19,
         col = alpha('grey', 0.4),
         xlab = expression('True '~xi), 
         ylab = "Bias",
         main = expression('After bias correction '~xi)
)
abline(h = 0, 
       col = 'magenta',
       lwd = 2) 

### --- log Kappa Parameter ---:
## Bias-Correction for Shape parameter: 
biaslnKappa <- (0.5*log(kappa2_pred_train) - 0.5*log(kappa2_true_train))

# Fit smoothing splines with different smoothing parameters:
spline_lnKappa_biasCorrect <- smooth.spline(x=0.5*log(kappa2_true_train),
                                            y=biaslnKappa,
                                            spar=0.9)

# Predict bias corrections
pred_biaslnKappa  <- predict(spline_lnKappa_biasCorrect,
                             x=0.5*log(kappa2_true_test))$y

# Plot residuals after smoothing
bplot.xy(0.5*log(kappa2_true_test),
         (0.5*log(kappa2_pred_test) - 0.5*log(kappa2_true_test) - pred_biaslnKappa),
         pch = 19,
         col = alpha('grey', 0.4),
         xlab = expression('True ln('~kappa~')'), 
         ylab = "Bias",
         main = expression('After bias correction ln('~kappa~')')
)
abline(h = 0, 
       col = 'magenta',
       lwd = 2) 
# lines(spline_lnKappa_biasCorrect, col = "blue", lty=2, lwd = 2)  # Add smoothed spline line



### --- log Tau Parameter ---:
## Bias-Correction for Shape parameter: 
biaslnTau<- (0.5*log(tau2_pred_train) - 0.5*log(tau2_true_train))

# Fit smoothing splines with different smoothing parameters:
spline_lnTau_biasCorrect <- smooth.spline(x=0.5*log(tau2_true_train),
                                            y=biaslnTau,
                                            spar=0.9)

# Predict bias corrections
pred_biaslnTau  <- predict(spline_lnTau_biasCorrect,
                           x=0.5*log(tau2_true_test))$y

# Plot residuals after smoothing
bplot.xy(0.5*log(tau2_true_test),
         (0.5*log(tau2_pred_test) - 0.5*log(tau2_true_test) - pred_biaslnTau),
         pch = 19,
         col = alpha('grey', 0.4),
         xlab = expression('True ln('~tau~')'), 
         ylab = "Bias",
         main = expression('After bias correction ln('~tau~')')
)
abline(h = 0, 
       col = 'magenta',
       lwd = 2) 

# lines(spline_lnKappa_biasCorrect, col = "blue", lty=2, lwd = 2)  # Add smoothed spline line

