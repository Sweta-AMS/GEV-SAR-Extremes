# rm(list=ls())

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
load("ParameterConfiguration-Dec12-Test.RData")
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
file_path <- path.expand("~/Downloads/testEst-n30-100Rep-FixedLambda.npy")
R <- 100
ns <- 20
# Import numpy and load the file
np <- import("numpy")
mat <- np$load(file_path)

# Load the NumPy file
cnnEst <- mat
dim(cnnEst)  # 400x100x3

# Bias correction for shape
shape_pred <- cnnEst[, ,1]
kappa2_pred <- cnnEst[, ,2]

shape_true <- storeParameterRep[ ,1]
shape_true <- matrix(rep(shape_true, each = R), nrow = 400, byrow = TRUE)
dim(shape_true)

kappa2_true <- storeParameterRep[ ,2] - 4
kappa2_true <- matrix(rep(kappa2_true, each = R), nrow = 400, byrow = TRUE)
dim(kappa2_true)

plot(c(shape_pred), c(shape_true))
abline(0,1)

plot(c(kappa2_pred), c(kappa2_true))
abline(0,1)

x <- cbind(c(shape_pred), c(0.5*log(kappa2_pred)))
dim(x)

y <- c(shape_true)
length(y)

## Bias correction:
## --- Shape Parmeter ---:
biasCorrectionShape <- Tps(x, y, df=10)
bcShape_pred <- predict(biasCorrectionShape)
plot(c(bcShape_pred), y, pch=19)
abline(0,1, col='magenta')

SEshape <- predictSE(biasCorrectionShape, x)

points(c(shape_pred), bcShape_pred + 1.96* SEshape , col=alpha(alpha=0.1, "red"), lty=2)
points(c(shape_pred), bcShape_pred - 1.96*SEshape , col=alpha(alpha=0.1, "red"), lty=2)

U <-  bcShape_pred + 1.96* SEshape
L <- bcShape_pred - 1.96*SEshape

df <- data.frame(Pred=c(shape_pred),
                 True=y,
                 L=L,
                 U=U)


png("CI.png",
    units="in", 
    width=21,
    height=10,
    res=200)
par(mfrow=c(1,2), 
    mai=c(18.5,5,12,16),
    mar=c(8,8,10.5,8) + 0.3,
    oma=c(0.4,2.5,3,8))

ggplot(df, aes(x =Pred, y = True)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L))



## --- Kappa2 Parmeter ---:
y <- c(0.5*log(kappa2_true))
length(y)

biasCorrectionlogKappa <- Tps(x, y, df=10)
bclogKappa_pred <- predict(biasCorrectionlogKappa)
plot(c(bclogKappa_pred), y, pch=19)
abline(0,1, col='magenta')

SElogkappa <- predictSE(biasCorrectionlogKappa, x)

points(c(0.5*log(kappa2_pred)), bclogKappa_pred  + 1.96* SElogkappa, col=alpha(alpha=0.1, "red"), lty=2)
points(c(0.5*log(kappa2_pred)), bclogKappa_pred  - 1.96*SElogkappa, col=alpha(alpha=0.1, "red"), lty=2)

U <-  bclogKappa_pred + 1.96* SElogkappa
L <- bclogKappa_pred - 1.96*SElogkappa

dfkappa <- data.frame(Pred=c(0.5*log(kappa2_pred)),
                 True=y,
                 L=L,
                 U=U)


ggplot(dfkappa, aes(x=Pred, y=True)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L))
mtext(expression('95 % CI log'~kappa),
      side=3,
      cex=2,
      line=3)

dev.off()

## -- RMSE of the estimates -- :
# -- Shape parameter --:
# shapeEst_bc
bcXi <- matrix(bcShape_pred, nrow=10)
dim(bcXi)

square_error_shape <- (bcXi - shape_true)^2
dim(square_error_shape) # 125x50

rmse_shape <- sqrt(apply(square_error_shape, 1, mean))
length((rmse_shape))

bplot.xy(testPars[,1], rmse_shape)

# kappa2Est_bc
bclogKappa <- matrix(bclogKappa_pred, nrow=10)
dim(bclogKappa)

square_error_kappa <- (0.5*bclogKappa - 0.5*log(kappa2_true))^2
dim(square_error_kappa) # 125x50

rmse_log_kappa <- sqrt(apply(square_error_kappa, 1, mean))
length((rmse_log_kappa))

bplot.xy(testPars[,2]-4, rmse_log_kappa)


## RMSE: 

parameter_vals <- cbind(shapeGrid, 0.5*log(AWghtGrid-4))
dim(parameter_vals)
fit_tps_rmse_shape <- Tps(parameter_vals,
                          rmse_shape,
                          df=10)
fit_tps_rmse_logkappa <- Tps(parameter_vals,
                             rmse_log_kappa,
                             df=10)

#################### ----- ########################
setwd("~/Desktop")
png("RMSE.png",
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

## SE
SEshape <- matrix(SEshape , nrow=10)
SEshape_mean <- apply(SEshape, 1, mean)

SElogkappa <- matrix(SElogkappa , nrow=10)
SElogkappa_mean <- apply(SElogkappa, 1, mean)

fit_tps_SE_shape <- Tps(parameter_vals,
                        SEshape_mean,
                        df=10)
fit_tps_SE_logkappa <- Tps(parameter_vals,
                           SElogkappa_mean,
                           df=10)
