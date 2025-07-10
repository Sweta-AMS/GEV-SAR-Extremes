rm(list=ls())

## -- Required Package library -- :
source("~/Desktop/Tentative Results-Proj2/R scripts/requiredPackages.R")
library(reticulate)
library(MASS)
library(survival)
library(fitdistrplus)
library(RcppCNPy)

## -- Loading the Trained CNN model in R, creating python environment --:
# Ensure that the required Python packages are installed
# py_install("tensorflow")

## Load true paramater values:
load("~/ParameterConfiguration-Dec12-Test.RData")
dim(storeParameterRep)

## Loading the CNN Model for evaluation: not need at this time
# # Import the necessary Python modules using reticulate
# tf <- import("tensorflow")
# keras <- import("keras")
# 
# ## Check if the file exists
# # model_path <- "~/Downloads/FitSpatialGEVLogNormalCNN30.keras"
# model_path <- "~/Downloads/FitSpatialGEVLogNormalCNNLatest (1).keras"
# if (file.exists(model_path)) {
#   cat("File exists. Proceeding to load the model...\n")
# } else {
#   stop("File not found. Please check the file path.")
# }
# 
# # Use an absolute path if needed
# model_path_abs <- normalizePath(model_path)
# 
# # Load the TensorFlow model
# model <- tryCatch({
#   keras$models$load_model(model_path_abs)
# }, error = function(e) {
#   cat("Error loading the model: ", e$message, "\n")
#   reticulate::py_last_error()
# })
# 
# # Print model summary to verify it loaded correctly
# if (!is.null(model)) {
#   summary <- model$summary()
#   cat(summary)
# } else {
#   cat("Failed to load the model. Check the error message above for details.\n")
# }
################################################################################


# ## Save files: 
# save(extremalFieldsScaled, 
#      file="SpatialFieldFixedNuggetEffect.RData")
# save(extremalFields, 
#      file="OrginalScaleSpatialFieldFixedNuggetEffect.RData")
# save(testPars, file="TestParameters.RData")


##  Initialize cnnEst with correct dimensions
# cnnEst <- array(NA, dim=c(nc, R, 3))
# dim(cnnEst)
# # Iterate through the ns and R dimensions to reshape data and predict using the model
# for(i in 1:nc) {
#   data <- extremalFieldsScaled[i, , ,]
#   data_array <- reshape_data(data)  # Ensure data is reshaped to (R, 16, 16, 30)
#   dim(data_array)
#   
#   cnnEst[i, , ] <- model$predict(data_array)
# }
# Verify the dimensions of cnnEst
# Load the .npy file

# Expand the file path
file_path <- path.expand("~/Downloads/test_estimates_NGP.npy")
# "~/Downloads/testEst-n30-100Rep-FixedLambda.npy"
R <- 100
ns <- 100
# Import numpy and load the file
np <- import("numpy")
mat <- np$load(file_path)

# Load the NumPy file
cnnEst <- mat
dim(cnnEst)  # 100x100x3

# Bias correction for shape
shape_pred <- cnnEst[, ,1]
kappa2_pred <- cnnEst[, ,2]

shape_T <- storeParameterRep[ ,1]
shape_true <- matrix(rep(shape_T, each = R), nrow = 100, byrow = TRUE)
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

##  Fitting TPS for CI of the prediction:
BCshape_pred <- matrix(BCshape_pred , nrow=100)
MBCshape_pred <- apply(BCshape_pred, 1, mean)
length(MBCshape_pred)

BClnKappa_pred <- matrix(BClnKappa_pred , nrow=100)
MBClnKappa_pred <- apply(BClnKappa_pred, 1, mean)
length(MBClnKappa_pred)

ucShape <- Tps(MBCshape_pred, storeParameterRep[,1], df=10)
SEshape <- predictSE(ucShape, MBCshape_pred)


## Shape Pred and CI
CI_lower <- MBCshape_pred - 1.96*SEshape
CI_upper <- MBCshape_pred + 1.96*SEshape
true_shape <- storeParameterRep[,1]

# Combine data into a data frame
envelope_data <- data.frame(
  MBCshape_pred = MBCshape_pred,
  Response = storeParameterRep[, 1],
  CI_lower = CI_lower,
  CI_upper = CI_upper
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
BClnKappa_pred <- matrix(BClnKappa_pred , nrow=100)
MBClnKappa_pred <- apply(BClnKappa_pred, 1, mean)
length(MBClnKappa_pred)

uclnKappa <- Tps(MBClnKappa_pred, 0.5*log(storeParameterRep[,2]-4), df=10)
SElnKappa <- predictSE(uclnKappa, MBClnKappa_pred)

CI_lower <- MBClnKappa_pred - 1.96*SElnKappa
CI_upper <- MBClnKappa_pred + 1.96*SElnKappa
true_lnkappa <- 0.5*log(storeParameterRep[,2]-4)

# Combine data into a data frame
envelope_data <- data.frame(
  MBClnKappa_pred = MBClnKappa_pred,
  Response = 0.5*log(storeParameterRep[,2] - 4),
  CI_lower = CI_lower,
  CI_upper = CI_upper
)

# ggplot2 envelope plot with mathematical expression
library(ggplot2)

ggplot(envelope_data, aes(x=MBClnKappa_pred, y = Response)) +
  geom_point(color = "blue", size = 2) +  # Observed responses
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), 
              fill = "gray", alpha = 0.3) +  # Confidence envelope
  geom_line(aes(y = CI_lower), color = "red", linetype = "dashed") +  # Lower CI
  geom_line(aes(y = CI_upper), color = "red", linetype = "dashed") +  # Upper CI
  theme_minimal() +
  labs(
    title = expression('log'~kappa),    # Add xi as title
    x = "Estimate",            # x-axis label
    y = "True"                 # y-axis label
  ) +
  theme(
    plot.title = element_text(hjust = 0.45),  # Center the title
    plot.margin = unit(c(2, 2, 2, 2), "cm")  # Top, right, bottom, left margins
  )

# Save the plot as a PNG file
ggsave(
  filename = "CI-lnKappa.png",   # Name of the file
  plot = envelope_plot,             # Plot object to save
  width = 8,                        # Width in inches
  height = 6,                       # Height in inches
  dpi = 300                         # Resolution in dots per inch
)

# SElogkappa <- matrix(SElogkappa , nrow=10)
# SElogkappa_mean <- apply(SElogkappa, 1, mean)
# 
# fit_tps_SE_shape <- Tps(parameter_vals,
#                         SEshape_mean,
#                         df=10)
# fit_tps_SE_logkappa <- Tps(parameter_vals,
#                            SElogkappa_mean,
#                            df=10)
# 
# 
# SElnKappa <- matrix(SElnKappaPred, nrow=100)
# SElnKappaU <- apply(SElnKappa, 1, mean)
# 
# ## Mean estimate around the replicate
# bcShapePredMat <- matrix(bcShape_pred, nrow=100)
# shapeEst <- apply(bcShapePredMat, 1, mean)
# length(shapeEst)
# 
# tempSEshape <- matrix(SEShapePred, nrow=100)
# SEshape <- apply(tempSEshape, 1, mean)
# length(SEshape)
# 
# shapeEst
# LB_shape <- shapeEst-1.96*SEshape
# UB_shape <-shapeEst + 1.96*SEshape
# 
# # Combine shapeEst, LB_shape, and UB_shape into a data frame
# plot_data <- data.frame(
#   x = 1:length(storeParameterRep[,1]),  # Index for x-axis
#   shapeEst = storeParameterRep[,1],     # Mean estimate
#   LB_shape = LB_shape,     # Lower bound
#   UB_shape = UB_shape      # Upper bound
# )
# 
# 
# # Load ggplot2 library
# library(ggplot2)
# 
# # Create the CI plot
# p2 <- ggplot(plot_data, aes(x = x, y = shapeEst)) +
#   geom_point(color = "blue", size = 2) +  # Mean estimate points
#   geom_line(color = "red") +             # Line connecting estimates
#   geom_ribbon(aes(ymin = LB_shape, ymax = UB_shape), 
#               alpha = 0.3, fill = "gray") +  # Confidence interval ribbon
#   theme_minimal() +  # Clean theme
#   labs(
#     title = "Mean Estimates with Confidence Intervals",
#     x = "Index",
#     y = "Mean Estimate"
#   ) +
#   theme(
#     plot.title = element_text(hjust = 0.5)  # Center title
#   )
# 
# # Print the plot
# print(p2)
# 
# 
# bClnKappa_pred <- bClnKappaModel$fitted.values
# 
# plot(c(bClnKappa_pred), y, pch=19, col=alpha('black', 0.1))
# abline(0,1, col='magenta')
# 

## BF 
# ## CI:
# plot(c(bClnKappa_pred), y, pch=19, col=alpha('black', 0.1))
# lines( c(bClnKappa_pred), c(bClnKappa_pred) - 1.96*SEbClnKappaPred, col="red1")
# lines( c(bClnKappa_pred), c(bClnKappa_pred) + 1.96*SEbClnKappaPred, col="red1")
# points(c(shape_pred), bcShape_pred + 1.96* SEshape , col=alpha(alpha=0.1, "red"), lty=2)
# points(c(shape_pred), bcShape_pred - 1.96*SEshape , col=alpha(alpha=0.1, "red"), lty=2)
# 
# U <-  bcShape_pred + 1.96* SEshape
# L <- bcShape_pred - 1.96*SEshape
# 
# df <- data.frame(Pred=c(shape_pred),
#                  True=y,
#                  L=L,
#                  U=U)
# 
# 
# png("CI.png",
#     units="in", 
#     width=21,
#     height=10,
#     res=200)
# par(mfrow=c(1,2), 
#     mai=c(18.5,5,12,16),
#     mar=c(8,8,10.5,8) + 0.3,
#     oma=c(0.4,2.5,3,8))
# 
# ggplot(df, aes(x =Pred, y = True)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L))
# 
# 

# ## --- Kappa2 Parmeter ---:
# y <- c(0.5*log(kappa2_true))
# length(y)
# 
# biasCorrectionlogKappa <- Tps(x, y)
# bclogKappa_pred <- predict(biasCorrectionlogKappa)
# plot(c(bclogKappa_pred), y, pch=19)
# abline(0,1, col='magenta')
# 
# SElogkappa <- predictSE(biasCorrectionlogKappa, x)
# 
# points(c(0.5*log(kappa2_pred)), bclogKappa_pred  + 1.96* SElogkappa, col=alpha(alpha=0.1, "red"), lty=2)
# points(c(0.5*log(kappa2_pred)), bclogKappa_pred  - 1.96*SElogkappa, col=alpha(alpha=0.1, "red"), lty=2)
# 
# U <-  bclogKappa_pred + 1.96* SElogkappa
# L <- bclogKappa_pred - 1.96*SElogkappa
# 
# dfkappa <- data.frame(Pred=c(0.5*log(kappa2_pred)),
#                  True=y,
#                  L=L,
#                  U=U)
# 
# 
# ggplot(dfkappa, aes(x=Pred, y=True)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymax = U, ymin = L))
# mtext(expression('95 % CI log'~kappa),
#       side=3,
#       cex=2,
#       line=3)
# 
# dev.off()

## -- RMSE of the estimates -- :
# -- Shape parameter --:
# shapeEst_bc
bcXi <- matrix(BCshape_pred, nrow=100)
dim(bcXi)

square_error_shape <- (bcXi-shape_true)^2
dim(square_error_shape) # 125x50

rmse_shape <- sqrt(apply(square_error_shape, 1, mean))
length((rmse_shape))

# bplot.xy(shape_true, rmse_shape)

# kappa2Est_bc
bclogKappa <- matrix(BClnKappa_pred, nrow=100)
dim(bclogKappa)

square_error_kappa <- (0.5*bclogKappa - 0.5*log(kappa2_true))^2
dim(square_error_kappa) # 125x50

rmse_log_kappa <- sqrt(apply(square_error_kappa, 1, mean))
length((rmse_log_kappa))

bplot.xy(0.5*log(storeParameterRep[,2]-4), rmse_log_kappa)


## RMSE: 
parameter_vals <- cbind(storeParameterRep[,1], 0.5*log(storeParameterRep[,2]-4))
dim(parameter_vals)

fit_tps_rmse_shape <- Tps(parameter_vals,
                          rmse_shape,
                          df=10)
fit_tps_rmse_logkappa <- Tps(parameter_vals,
                             rmse_log_kappa,
                             df=10)

#################### ----- ########################
setwd("~/Desktop")
png("RMSE-Latest.png",
    units="in", 
    width=21,
    height=10,
    res=200)
par(mfrow=c(1,2), 
    mai=c(18.5,5,12,16),
    mar=c(8,8,10.5,8) + 0.3,
    oma=c(0.4,2.5,3,8))
surface(fit_tps_rmse_shape,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=1.5,
        cex.main=2,
        axis.args=list(cex.axis=1.5,lwd=2))
axis(2,
     las=1,
     cex.axis=2,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=1,
      cex=2.5,
      line=5,
      las=1)
mtext(expression('log '~kappa),
      side=2,
      cex=2.5,
      las=3,
      line=4)
mtext(expression('RMSE('~xi~')'),
      side=3,
      cex=2,
      line=3)
surface(fit_tps_rmse_logkappa,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=2,
        cex.main=2,
        axis.args=list(cex.axis=1.5,lwd=2))
axis(2,
     las=1,
     cex.axis=2,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=1,
      cex=2.5,
      line=5,
      las=1)
mtext(expression('log'~kappa),
      side=2,
      cex=2.5,
      las=3,
      line=4)
mtext(expression('RMSE(log '~kappa~')'),
      side=3,
      cex=2,
      line=3)
dev.off()


