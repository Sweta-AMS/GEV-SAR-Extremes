rm(list=ls())

## -- Required Package library -- :
source("~/Desktop/Tentative Results-Proj2/R scripts/requiredPackages.R")
library(reticulate)
library(MASS)
library(survival)
library(fitdistrplus)
library(RcppCNPy)

## Load true paramater values:
# load("~/ParameterConfiguration-Dec12-Test.RData")
load("~/ParameterConfiguration-Test-n-30.RData")
load("~/storeZRep-Test-n-30.RData")

dim(storeParameterRep) # 5000x3
dim(storeZRep) # 5000x256x30

# Expand the file path
# setwd("~/Desktop/Tentative Results-Proj2/Latest-Results-Jan08-2025/ResultsofTest-5000PCs") 
file_path <- path.expand("~/Desktop/Tentative Results-Proj2/Latest-Results-Jan08-2025/ResultsofTest-5000PCs/test_estimates_NGP-n-30-Latest-Jan08-2025.npy")
# "~/Downloads/testEst-n30-100Rep-FixedLambda.npy"
# R <- 100
# ns <- 100
# Import numpy and load the file
np <- import("numpy")
mat <- np$load(file_path)

# Load the NumPy file
cnnEst <- mat
dim(cnnEst)  # 100x100x3

# Bias correction for shape
shape_pred <- cnnEst[, 1]
logkappa2_pred <- cnnEst[, 2]
logtau2_pred <- cnnEst[, 3]

shape_T <- storeParameterRep[ ,1]
shape_true <- 
  
  #matrix(rep(shape_T, each = R), nrow = 100, byrow = TRUE)
  dim(shape_true)

kappa2_T <- storeParameterRep[ ,2] - 4
kappa2_true <- matrix(rep(kappa2_T, each = R), nrow = 100, byrow = TRUE)
dim(kappa2_true)

plot(c(shape_pred), c(shape_true))
abline(0,1, col='magenta')

plot(c(kappa2_pred), c(log(kappa2_true)))
abline(0,1, col='magenta')

### ----- Bias correction ----- ###
## Shape Parmeter
x <- cbind(c(shape_pred), c(0.5*(kappa2_pred)))
y <- c(shape_true)

tic()
bCShapeModel <- lm(y~x)
toc()

BCshape_pred <- bCShapeModel$fitted.values

plot(c(BCshape_pred), y, pch=19, col=alpha('black', 0.1))
abline(0,1, col='magenta')

## Kappa2 Parmeter
x <- cbind(c(shape_pred), c(0.5*(kappa2_pred)))
y <- c(0.5*log(kappa2_true))

tic()
bClnKappaModel <- lm(y~x)
toc()

BClnKappa_pred <- bClnKappaModel$fitted.values

##  Fitting Quantile Regression for CI of the prediction:
BCshape_pred <- matrix(BCshape_pred , nrow=100)
MBCshape_pred <- apply(BCshape_pred, 1, mean)
length(MBCshape_pred)

BClnKappa_pred <- matrix(BClnKappa_pred , nrow=100)
MBClnKappa_pred <- apply(BClnKappa_pred, 1, mean)
length(MBClnKappa_pred)

## Fit Quantile Regression: 
library(quantreg)
dataShape <- data.frame(true = storeParameterRep[,1],
                        estimated = MBCshape_pred)

datalnKappa <- data.frame(true = 0.5*log(storeParameterRep[,2]-4),
                          estimated = MBClnKappa_pred)

# Fit quantile regression
quantiles <- seq(0.2, 0.9, length.out=30)
# c(0.2, 0.25, 0.5, 0.75, 0.9)
modelsShape <- lapply(quantiles, function(tau) {
  rq(true ~ estimated, data = dataShape, tau = tau)
})

modelslnKappa <- lapply(quantiles, function(tau) {
  rq(true ~ estimated, data = datalnKappa, tau = tau)
})

# Summarize the models
fittedShape <- lapply(modelsShape, fitted)
QRShapeEst <- matrix(NA, nrow=30, ncol=100)
for(i in 1:30){
  QRShapeEst[i, ] <- fittedShape[[i]]
}
dim(QRShapeEst)

# Summarize the models
fittedlnKappa <- lapply(modelslnKappa, fitted)
QRlnKappaEst <- matrix(NA, nrow=30, ncol=100)
for(i in 1:30){
  QRlnKappaEst[i, ] <- fittedlnKappa[[i]]
}
dim(QRlnKappaEst)


## Compute lower and upper bound:
XiCI_lower_bound <- rep(NA, 100)
XiCI_upper_bound <- rep(NA, 100)
for(i in 1:100){
  lower_bd <- quantile(QRShapeEst[,i], 0.025)
  XiCI_lower_bound[i] <- lower_bd
  
  upper_bd <- quantile(QRShapeEst[,i], 0.975)
  XiCI_upper_bound[i] <- upper_bd
}

## Compute lower and upper bound:
lnKappaCI_lower_bound <- rep(NA, 100)
lnKappaCI_upper_bound <- rep(NA, 100)
for(i in 1:100){
  lower_bd <- quantile(QRlnKappaEst[,i], 0.025)
  lnKappaCI_lower_bound[i] <- lower_bd
  
  upper_bd <- quantile(QRlnKappaEst[,i], 0.975)
  lnKappaCI_upper_bound[i] <- upper_bd
}

# Combine data into a data frame
envelope_data <- data.frame(
  MBCshape_pred = MBCshape_pred,
  Response = storeParameterRep[, 1],
  CI_lower =  XiCI_lower_bound,
  CI_upper =  XiCI_upper_bound
)

# ggplot2 envelope plot with mathematical expression
library(ggplot2)

ggplot(envelope_data, aes(x = MBCshape_pred, y = Response)) +
  geom_point(color = "blue", size = 2) +  # Observed responses
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), 
              fill = "gray", alpha = 0.3) +  # Confidence envelope
  geom_line(aes(y = CI_lower), color = "red", linetype = "dashed") +  # Lower CI
  geom_line(aes(y = CI_upper), color = "red", linetype = "dashed") +  # Upper CI
  theme_minimal() +
  labs(
    title = expression(xi),    # Add xi as title
    x = "Estimate",            # x-axis label
    y = "True"                 # y-axis label
  ) +
  theme(
    plot.title = element_text(hjust = 0.45),  # Center the title
    plot.margin = unit(c(2, 2, 2, 2), "cm")  # Top, right, bottom, left margins
  )
# Save the plot as a PNG file
ggsave(
  filename = "CI-shape.png",   # Name of the file
  plot = envelope_plot,             # Plot object to save
  width = 8,                        # Width in inches
  height = 6,                       # Height in inches
  dpi = 300                         # Resolution in dots per inch
)


## Kappa 

# Combine data into a data frame
envelope_data <- data.frame(
  MBCshape_pred = MBClnKappa_pred,
  Response = 0.5*log(storeParameterRep[, 2]-4),
  CI_lower =  lnKappaCI_lower_bound,
  CI_upper =  lnKappaCI_upper_bound
)

# ggplot2 envelope plot with mathematical expression
library(ggplot2)

ggplot(envelope_data, aes(x =  MBClnKappa_pred, y = Response)) +
  geom_point(color = "blue", size = 2) +  # Observed responses
  geom_ribbon(aes(ymin = lnKappaCI_lower_bound, ymax = lnKappaCI_upper_bound), 
              fill = "gray", alpha = 0.3) +  # Confidence envelope
  geom_line(aes(y = CI_lower), color = "red", linetype = "dashed") +  # Lower CI
  geom_line(aes(y = CI_upper), color = "red", linetype = "dashed") +  # Upper CI
  theme_minimal() +
  labs(
    title = expression('ln'~kappa),    # Add xi as title
    x = "Estimate",            # x-axis label
    y = "True"                 # y-axis label
  ) +
  theme(
    plot.title = element_text(hjust = 0.45),  # Center the title
    plot.margin = unit(c(2, 2, 2, 2), "cm")  # Top, right, bottom, left margins
  )

