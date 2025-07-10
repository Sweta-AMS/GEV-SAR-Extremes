## Generate plot for the test set:
# rm(list=ls())

## -- Required Package library -- :
source("~/Desktop/Tentative Results-Proj2/R scripts/requiredPackages.R")
source("~/Downloads/LKrigSimSAR.R")
source("~/LKrigSARPareto.R")
source("~/Desktop/Tentative Results-Proj2/R script/LKrigSAREvd.R")
library(reticulate)
library(MASS)
library(survival)
library(fitdistrplus)

## -- Loading the Trained CNN model in R, creating python environment --:
# Ensure that the required Python packages are installed
py_install("tensorflow")

# Import the necessary Python modules using reticulate
tf <- import("tensorflow")
keras <- import("keras")

## Check if the file exists
# model_path <- "~/Downloads/FitSpatialGEVLogNormalCNN30.keras"
model_path <- "~/Downloads/FitSpatialGEVLogNormalCNNLatest (1).keras"
if (file.exists(model_path)) {
  cat("File exists. Proceeding to load the model...\n")
} else {
  stop("File not found. Please check the file path.")
}

# Use an absolute path if needed
model_path_abs <- normalizePath(model_path)

# Load the TensorFlow model
model <- tryCatch({
  keras$models$load_model(model_path_abs)
}, error = function(e) {
  cat("Error loading the model: ", e$message, "\n")
  reticulate::py_last_error()
})

# Print model summary to verify it loaded correctly
if (!is.null(model)) {
  summary <- model$summary()
  cat(summary)
} else {
  cat("Failed to load the model. Check the error message above for details.\n")
}
################################################################################
## Function to use a min-max scaling
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

## Location
M <- 16 # changed from 15
sGrid <- list(x=1:M, y=1:M)
s <- make.surface.grid(sGrid)
plot(s, pch=19)

n <- nrow(s)
print(paste0('Field size: ', n))

## CI based on fitting simple linear regression or just spline: 
true_parameter <- load('~/Downloads/ParameterConfiguration-LatestNov8-Test.RData')
true_parameter <- storeParameterRep
dim(true_parameter)

pred_parameterRep30 <- read.csv('~/Downloads/test_estimates_NGP-5000TestSampleRep30.csv', header=FALSE)
dim(pred_parameterRep30)

pred_parameterRep10 <- read.csv('~/Downloads/test_estimates_NGP-5000TestSampleRep10.csv', header=FALSE)
dim(pred_parameterRep10)


pred_parameterRep1 <- read.csv('~/Downloads/test_estimates_NGP-5000TestSampleRep1.csv', header=FALSE)
dim(pred_parameterRep1)

# RMSE across the true parameter:


setwd("~/Desktop")
png("BoxplotTestSet30.png",
    units="in", 
    width=13,
    height=4,
    res=200)
par(mfrow=c(1,3), 
    mar=c(4,4,3,3) + 0.3,
    oma=c(0.4,2.5,3,8))
bplot.xy(true_parameter[,1],
         pred_parameterRep30[,1]-true_parameter[,1])
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression(xi),
      side=3,
      # las=1,
      line=1)
bplot.xy(0.5*log(true_parameter[,2]-4), 0.5*log(pred_parameterRep30[,2])-0.5*log(true_parameter[,2]-4))
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression('log '~kappa),
      side=3,
      # las=1,
      line=1)
bplot.xy(0.5*log(true_parameter[,3]), 0.5*pred_parameterRep30[,3]-0.5*true_parameter[,3])
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression('log '~tau),
      side=3,
      # las=1,
      line=1)
dev.off()


##Rep 10
setwd("~/Desktop")
png("BoxplotTestSet10.png",
    units="in", 
    width=13,
    height=4,
    res=200)
par(mfrow=c(1,3), 
    mar=c(4,4,3,3) + 0.3,
    oma=c(0.4,2.5,3,8))
bplot.xy(true_parameter[,1],
         pred_parameterRep10[,1]-true_parameter[,1])
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression(xi),
      side=3,
      # las=1,
      line=1)
bplot.xy(0.5*log(true_parameter[,2]-4), 0.5*log(pred_parameterRep10[,2])-0.5*log(true_parameter[,2]-4))
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression('log '~kappa),
      side=3,
      # las=1,
      line=1)
bplot.xy(0.5*log(true_parameter[,3]), 0.5*pred_parameterRep10[,3]-0.5*true_parameter[,3])
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression('log '~tau),
      side=3,
      # las=1,
      line=1)
dev.off()


##Rep 1
setwd("~/Desktop")
png("BoxplotTestSet1.png",
    units="in", 
    width=13,
    height=4,
    res=200)
par(mfrow=c(1,3), 
    mar=c(4,4,3,3) + 0.3,
    oma=c(0.4,2.5,3,8))
bplot.xy(true_parameter[,1],
         pred_parameterRep1[,1]-true_parameter[,1])
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression(xi),
      side=3,
      # las=1,
      line=1)
bplot.xy(0.5*log(true_parameter[,2]-4), 0.5*log(pred_parameterRep1[,2])-0.5*log(true_parameter[,2]-4))
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression('log '~kappa),
      side=3,
      # las=1,
      line=1)
bplot.xy(0.5*log(true_parameter[,3]), 0.5*pred_parameterRep1[,3]-0.5*true_parameter[,3])
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('Bias in estimation',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression('log '~tau),
      side=3,
      # las=1,
      line=1)
dev.off()

## -- Computing conference interval --:
CIshapeFit <- lm(true_parameter[,1]~pred_parameter[,1] + pred_parameter[,2] + pred_parameter[,3])
CIkappa2Fit <- lm(true_parameter[,2]~pred_parameter[,1] + pred_parameter[,2] + pred_parameter[,3])
CItau2Fit <- lm(true_parameter[,3]~pred_parameter[,1] + pred_parameter[,2] + pred_parameter[,3])

# Make predictions with confidence intervals
CIshape <- predict(CIshapeFit, 
                   interval="confidence",
                   level=0.95)

true_shape <- true_parameter[,1]

## CI
# Plotting the band
lines(true_shape, CIshape[,2], col = "blue", lwd = 2)  # Lower curve
plot(true_shape[1:100], pred_parameter[1:100,1])
lines(true_shape,  CIshape[,3], col = "blue", lwd = 2)  # Upper curve





plot(pred_parameter[,1] ,
     pred_shape[,1])

# # # Make predictions and get confidence intervals
# predictions <- predict(pred_function,
#                        newdata=c(pred_parameter),
#                        interval="confidence")

# # Convert predictions to a data frame for easier handling
# predictions_df <- data.frame(pred_parameter,
#                              predictions)


plot(true_parameter[,1])
points(pred_parameter[,1], predictions_df$lwr, col='magenta')
points(pred_parameter[,1], predictions_df$upr, col='green')



# Make predictions with confidence intervals
pred_shape <- predict(pred_function, 
                      interval="confidence",
                      level=0.95)
plot(pred_parameter[,1] ,
     pred_shape[,1])
#pred_kappa2 <- predict(bias_correction_kappa2, interval = "confidence", level = 0.95)

# Combine predictions and confidence intervals into a data frame for plotting
predictions_df_shape <- data.frame(true_parameter[,1], pred_parameter[,1], pred_shape)
# Load necessary library
library(ggplot2)

# Plot for shape
ggplot(predictions_df_shape, aes(x = parameter_shape, y = cnnEst_shape)) +
  geom_point() +
  geom_line(aes(y = fit), color = "red") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(title = "Shape Predictions with 95% CI", x = "Parameter Shape", y = "CNN Estimate Shape") +
  theme_minimal()

















ns <- 1000
AWghtGrid <- 4 + exp(runif(ns, log(0.001), log(1.6)))
lambdaGrid <- exp(runif(ns, log(0.0001), log(0.001)))
shapeGrid <- runif(ns, 0.1, 0.8)

testPars <- cbind(shapeGrid, AWghtGrid, lambdaGrid)
dim(testPars)

# expand.grid(shapeGrid, AWghtGrid)
# cbind(shapeGrid, AWghtGrid, lambdaGrid)
dim(testPars)
summary(testPars)










m <- 30

nc <- nrow(testPars )
print(nc)

parameter_array <- array(NA, dim=c(nc, 3))
dim(parameter_array)

extremalFieldsScaled <- array(NA, dim=c(nc, n, m))
dim(extremalFieldsScaled)

extremalFields <- array(NA, dim=c(nc, n, m))
dim(extremalFields)



tic()  
# Initialize result_list before starting the loop
set.seed(123)
for (i in 1:nc)
{
  cat('GEV parameter loop no', i, '\n')
  
  # Defining GEV parameters:
  shape_val <- testPars[i,1]
  a_wght <- testPars[i,2]
  lambda_val <- testPars[i,3]
  
  print(paste("Loop", i))
  
  # parameter_array[i, ] <- as.numeric(c(shape_val, (a_wght-4), lambda_val))
  #print(paste("Inner Loop", j))
  
  # Setup LKrig
  LKinfo <- LKrigSetup(s,
                       a.wght=a_wght,
                       nlevel=1,
                       nu=1,
                       NC=16,
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
  
  ## Set another seed for the log-normal matrix generation, unique for each iteration
  ## set.seed(111 + i + j)
  LogNormalMat <- exp(matrix(rnorm(n*m,
                                   mean=locParameter,
                                   sd=sqrt(scaleParameter)),
                             n, m))  # n x m matrix
  
  # Scale each column 
  Z <- (ySim * LogNormalMat)
  Z_scaled <- scale_matrix(Z)
  
  # Store: extremalFields array
  extremalFields[i, , ] <- Z 
  extremalFieldsScaled[i, , ] <- Z_scaled # removed Z_scaled
}
toc()  # End timing for the entire process
dim(extremalFieldsScaled)
# time: 342.946 sec elapsed

library(reticulate)
# Load the TensorFlow model
model_path <- "~/Downloads/FitSpatialGEVLogNormalCNNLatest (2).keras"
if (file.exists(model_path)) {
  cat("File exists. Proceeding to load the model...\n")
} else {
  stop("File not found. Please check the file path.")
}
model_path_abs <- normalizePath(model_path)

# Load the TensorFlow model
model <- tryCatch({
  keras$models$load_model(model_path_abs)
}, error = function(e) {
  cat("Error loading the model: ", e$message, "\n")
  reticulate::py_last_error()
})

# Print model summary to verify it loaded correctly
if (!is.null(model)) {
  summary <- model$summary()
  cat(summary)
} else {
  cat("Failed to load the model. Check the error message above for details.\n")
}

# # Create a function to reshape the data correctly
# reshape_data <- function(data) {
#   reshaped_data <- array(data, dim=c(dim(data)[1], 16, 16, 30))
#   return(reshaped_data)
# }

# Check dimensions of the extremalFieldsScaled array
dim(extremalFieldsScaled)

# Initialize cnnEst with correct dimensions
reshaped_data <- array(extremalFieldsScaled, dim=c(dim(extremalFieldsScaled)[1], 16, 16, 30))
dim(reshaped_data)

cnnEst <- model$predict(reshaped_data)
summary(cnnEst)

# # Iterate through the ns and R dimensions to reshape data and predict using the model
# for(i in 1:nc) {
#   data <- extremalFieldsScaled[i, , ]
#   data_array <- reshape_data(data)  # Ensure data is reshaped to (R, 16, 16, 30)
#   dim(data_array)
#   
#   cnnEst[i, ] <- model$predict(data_array)
# }
# Verify the dimensions of cnnEst
dim(cnnEst) 
dim(parameter_array)

# Bias correction for shape
beta0 <- 0.0017
beta1 <- 0.9978

shapeEst_bc <- 0.0017 + 0.9978*(cnnEst[,1])
plot(shapeGrid, cnnEst[,1])

# Bias correction for kappa2
alpha0 <- 0.0124
alpha1 <- 0.9462
kappa2Est_bc <- 0.0124 + 0.9462*(cnnEst[,2])

# ## Bias correction:
# bias_correction_shape <- lm(c(cnnEst[,,1])~c(parameter_array[,,1]))
# beta0 <- 0.09160584    
# beta1 <- 0.74653919 
# 
# bias_correction_kappa2 <- lm(c(cnnEst[,,2])~c(parameter_array[,,2]))
# alpha0 <- 0.09160584  
# alpha1 <-  0.3782559 


bplot.xy(parameter_array[ ,,1], cnnEst[, ,1]-parameter_array[,,1])
abline(h=0, col='magenta')
bplot.xy(parameter_array[ ,,2], cnnEst[, ,2]-parameter_array[,,2])
abline(h=0, col='magenta')
bplot.xy(parameter_array[ ,,3], cnnEst[, ,3]-parameter_array[,,3])
abline(h=0, col='magenta')


## -- RMSE of the estimates -- :
# -- Shape parameter --:
square_error_shape <- (shapeEst_bc -parameter_array[, ,1])^2
dim(square_error_shape) # 125x50

rmse_shape <- sqrt(apply(square_error_shape, 1, mean))
length((rmse_shape))

bplot.xy(testPars[,1], rmse_shape)

square_error_kappa2 <- (kappa2Est_bc-parameter_array[, ,2])^2
dim(square_error_kappa2) # 125x50

rmse_kappa2 <- sqrt(apply(square_error_kappa2, 1, mean))
length((rmse_kappa2))

bplot.xy(testPars[,2], rmse_kappa2)

fit_tps_rmse_shape <- Tps(testPars,
                          rmse_shape,
                          df=10)
fit_tps_rmse_kappa2 <- Tps(testPars,
                           rmse_kappa2,
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
mtext(expression(kappa^2),
      side=2,
      cex=2.5,
      las=1,
      line=4)
mtext(expression('RMSE('~xi~')'),
      side=3,
      cex=2,
      line=3)
surface(fit_tps_rmse_kappa2,
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
mtext(expression(kappa^2),
      side=2,
      cex=2.5,
      las=1,
      line=4)
mtext(expression('RMSE('~kappa^2~')'),
      side=3,
      cex=2,
      line=3)
dev.off()

# plot(cnnEst[,,1], c(parameter_array[, ,1]))
# abline(a=0, b=1, col='magenta')

# plot(cnnEst[,,2], c(parameter_array[, ,2]))
# abline(a=0, b=1, col='magenta')

# plot(cnnEst[,,3], c(parameter_array[, ,3]))
# abline(a=0, b=1, col='magenta')

# ## -- Bias correction --:
# # Bias Correction: IsotonicRegression
# library(isotone)
# dim(cnnEst)
# 
# fit_shape <- gpava(c(cnnEst[, ,1]), c(parameter_array[, ,1]), weights=NULL)
# fit_kappa2 <- gpava(c(cnnEst[, ,2]), c(parameter_array[, ,2]), weights=NULL)
# fit_tau2 <- gpava(c(cnnEst[, ,3]), c(parameter_array[, ,3]), weights=NULL)
# 
# # Extract the fitted valuess
# shape_pred <- fit_shape$x
# kappa2_pred <- fit_kappa2$x
# tau2_pred <- fit_tau2$x
# 
# plot(shape_pred, c(parameter_array[, ,1]))
# abline(a=0, b=1, col='magenta')
# 
# plot(kappa2_pred, c(parameter_array[, ,2]))
# abline(a=0, b=1, col='magenta')
# 
# plot(tau2_pred, c(parameter_array[, ,3]))
# abline(a=0, b=1, col='magenta')

################################################################################
## Bias correction: Quantile mapping in R:
y_obs <- parameter_array[ , ,1]
y_pred <- cnnEst[, ,1]

simple_regression <- lm(c(y_obs)~c(y_pred))
y_pred <- simple_regression$fitted.values

plot(y_pred,y_obs )

# Make predictions and get confidence intervals
predictions <- predict(simple_regression,
                       newdata=c(y_pred), 
                       interval="confidence")

# Convert predictions to a data frame for easier handling
predictions_df <- data.frame(y_pred,
                             predictions)
################################################################################

# Apply Quantile Mapping
cor_predict_shape <- quantile_mapping(TestParameterComb[,1] , pred_shape)
pred_shapeBC <- matrix(cor_predict_shape, nrow=ncomb, ncol=R)


## Bias correction
QuanitleMap <- 
  
  
  
  
  
  ## Latest used
  ncomb <- 4*4*4
cat('Test sample', ncomb, '\n')
B <- 1000
R <- 100 # no of repititon across the parameter configuration for coverage calculation:

## Test Parameter Configuration generated across grid:
AWghtGrid <- 4 + exp(seq(log(0.001), log(2), length.out=(ncomb)^(1/3)))
lambdaGrid <- exp(seq(log(0.0001), log(0.01), length.out=(ncomb)^(1/3)))
shapeGrid <- seq(0.01, 0.9, length.out=(ncomb)^(1/3))

# AWghtGrid <- 4 + exp(runif(n=ncomb, log(0.001), log(1.8)))
# lambdaGrid <- exp(runif(n=ncomb, log(0.0001), log(0.001)))
# shapeGrid <- runif(n=ncomb, 0.01, 0.8)

TestParameterComb <-  expand.grid(shapeGrid, AWghtGrid, lambdaGrid)
dim(TestParameterComb)

# TestParameterComb  <- TestParameterComb[1:ncomb, ]


# For replication size 10:
source("~/Desktop/Tentative Results-Proj2/R script/LKrigSAREvd.R")
m <- 30
R <- 100
extremalFields <- array(NA, dim=c(ncomb, R, n, m))
dim(extremalFields)

extremalFieldsScaled <- array(NA, dim=c(ncomb, R, n, m))
dim(extremalFieldsScaled)

tic()  
# Initialize result_list before starting the loop
for (i in 1:ncomb)
{
  cat('GEV parameter loop no', i, '\n')
  
  # Defining GEV parameters:
  shape_val <- TestParameterComb[i,1]
  a_wght <- TestParameterComb[i,2]
  lambda_val <- TestParameterComb[i,3]
  
  print(paste("Loop", i))
  
  for(j in 1:R)
  {
    ##  Set a seed for reproducibility, unique for each iteration
    ## set.seed(123 + i+ j)
    print(paste("Inner Loop", j))
    
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
    
    ## Set another seed for the log-normal matrix generation, unique for each iteration
    ## set.seed(111 + i + j)
    LogNormalMat <- exp(matrix(rnorm(n*m,
                                     mean=locParameter,
                                     sd=sqrt(scaleParameter)),
                               n, m))  # n x m matrix
    
    # Scale each column 
    Z <- (ySim * LogNormalMat)
    Z_scaled <- scale_matrix(Z)
    
    # Store: extremalFields array
    extremalFields[i,j, , ] <- Z 
    extremalFieldsScaled[i,j, , ] <- Z_scaled # removed Z_scaled
  }
}
toc()  # End timing for the entire process
dim(extremalFieldsScaled)
# time: 19.15948 minutes

save(extremalFieldsScaled, file="extremalFieldsScaled-TestData.RData")
save(TestParameterComb, file="TestParameterComb.RData")

t <- 1
image.plot(matrix(extremalFieldsScaled[t,10,,1], 16, 16))
TestParameterComb[t,]





library(extRemes)
input_data <- array(revd(256*30), dim = c(1, 16, 16, 30))  # Change shape according to your model

# Make predictions using the model
predictions <- model$predict(input_data)

# Print the predictions
print(predictions)
dim( extremalFieldsScaled)

pred_parameter <- array(NA, c(ncomb, R, 3))
dim(pred_parameter)

for(i in 1:ncomb)
{
  data <- array(extremalFieldsScaled[i, , ,], dim = c(R, 16, 16, 30))
  pred_parameter[i, , ] <- model$predict(data)
  
}

pred_shape <- pred_parameter[, ,1]
dim(pred_shape)

pred_kappa2 <- pred_parameter[, ,2]
dim(pred_kappa2)

pred_tau2 <- pred_parameter[, ,3]
dim(pred_tau2)

## Bias correction: Quantile mapping in R:
library(MASS)
library(survival)
library(fitdistrplus)

# Function to perform Quantile Mapping
quantile_mapping <- function(observed, predicted) {
  obs_ecdf <- ecdf(observed)
  pred_ecdf <- ecdf(predicted)
  quantiles <- pred_ecdf(predicted)
  corrected_predicted <- quantile(obs_ecdf, quantiles)
  return(corrected_predicted)
}

# Apply Quantile Mapping
cor_predict_shape <- quantile_mapping(TestParameterComb[,1] , pred_shape)
pred_shapeBC <- matrix(cor_predict_shape, nrow=ncomb, ncol=R)

hist(cor_predict_shape)
hist(TestParameterComb[,1])

cor_predict_kappa2 <- quantile_mapping((TestParameterComb[,2]-4) , pred_kappa2)
pred_kappa2BC <- matrix(cor_predict_kappa2, nrow=ncomb, ncol=R)

hist(cor_predict_kappa2)
hist(TestParameterComb[,2]-4)

cor_predict_tau2 <- quantile_mapping((TestParameterComb[,3]) , pred_tau2)
pred_tau2BC <- matrix(cor_predict_tau2, nrow=ncomb, ncol=R)

hist(cor_predict_tau2)
hist(TestParameterComb[,3])


# # Step 4: Plot the original and corrected predictions against observed values
# plot(observed, predicted, col = alpha("red", 0.5), pch = 19,
#      xlab = "Observed Values", ylab = "Predicted Values",
#      main = "Bias Correction using Quantile Mapping")
# points(observed, corrected_predicted, col = alpha("blue", 0.5), pch = 19)
# legend("topleft", legend = c("Original Predicted", "Corrected Predicted"),
#        col = c("red", "blue"), pch = 19)


## Histogram:
hist(pred_shapeBC[1,])
abline(v=TestParameterComb[1 ,1], col='magenta')

hist(pred_shapeBC[2, ])
abline(v=TestParameterComb[2 ,1], col='magenta')

hist(pred_shapeBC[3, ])
abline(v=TestParameterComb[3 ,1], col='magenta')

hist(pred_shapeBC[4, ])
abline(v=TestParameterComb[4 ,1], col='magenta')

hist(exp(pred_kappa2[1, ]))
abline(v=TestParameterComb[1 ,2]-4, col='magenta')

hist(exp(pred_kappa2[2, ]))
abline(v=TestParameterComb[2 ,2]-4, col='magenta')

hist(exp(pred_kappa2[3, ]))
abline(v=TestParameterComb[3, 2]-4, col='magenta')

hist(exp(pred_kappa2[4, ]))
abline(v=TestParameterComb[4 ,2]-4, col='magenta')


hist(pred_kappa2[5, ])
abline(v=TestParameterComb[5,2]-4, col='magenta')


hist(pred_tau2[6, ])
abline(v=TestParameterComb[6,3], col='magenta')

hist(pred_shape[7, ])
abline(v=TestParameterComb[7,1], col='magenta')

hist(pred_shape[8, ])
abline(v=TestParameterComb[8,1], col='magenta')

dim(pred_parameter) # 27 100   3

## Doing bootstrap: 
# Initialize arrays for bootstrap predictions
source("~/LKsimFunc.R") # for generating bootstrap samples with replication
boot_pred_shape <- matrix(NA, nrow=ncomb, ncol=R)
boot_pred_Kappa2 <- matrix(NA, nrow=ncomb, ncol=R)
boot_pred_tau2 <- matrix(NA, nrow=ncomb, ncol=R)

CI_shape <- array(NA, dim= c(ncomb, R, 2))
CI_Kappa2 <- array(NA, dim= c(ncomb, R, 2))
CI_tau2 <- array(NA, dim= c(ncomb, R, 2))

alpha <-  0.05

tic()
# Perform bootstrap resampling
for (i in 1:ncomb) {
  cat('Loop:', i , "\n")
  tic()
  for (j in 1:R) {
    cat('Inner loop:', j, "\n")
    theta <- pred_parameter[i, j , ] 
    theta[2] <- theta[2]+ 4
    
    # Bootstrap sample: 
    generateBsample <- simulate_spatial_field(s=s, theta=theta, R=B, m=30)
    Bsample <- array(generateBsample$extremalFieldsScaled, c(B, 16, 16, m))
    
    boot_est <- model$predict(Bsample)
    
    # Calculate the confidence interval and coverage probability
    lower_percentile <- alpha/2
    upper_percentile <- 1 - (alpha/2)
    
    lower_bound_shape <-  quantile(boot_est[ ,1], probs=lower_percentile)
    upper_bound_shape <- quantile(boot_est[ ,1], probs=upper_percentile)
    
    lower_bound_Kappa2 <-  quantile(boot_est[ ,2], probs=lower_percentile)
    upper_bound_Kappa2 <-  quantile(boot_est[ ,2], probs=upper_percentile)
    
    lower_bound_tau2 <- quantile(boot_est[ ,3],  probs=lower_percentile)
    upper_bound_tau2 <- quantile(boot_est[ ,3],  probs=upper_percentile)
    
    
    boot_pred_shape[i, j] <- quantile(boot_est[ ,1], 0.5)
    boot_pred_Kappa2[i, j] <- quantile(boot_est[ ,2], 0.5)
    boot_pred_tau2[i, j] <- quantile(boot_est[ ,3], 0.5)
    
    CI_shape[i, j, ] <- c(lower_bound_shape, upper_bound_shape)
    CI_Kappa2[i, j, ] <- c(lower_bound_Kappa2, upper_bound_Kappa2)
    CI_tau2[i, j, ] <- c(lower_bound_tau2, upper_bound_tau2)
  }
  toc()
}
toc()

# Calculate coverage probabilities
coverage_bootShape <- numeric(ncomb)
coverage_bootKappa2 <- numeric(ncomb)
coverage_bootTau2 <- numeric(ncomb)

for (i in 1:ncomb) {
  coverage_bootShape[i] <- mean(CI_shape[i, ,1] < TestParameterComb[i,1] & CI_shape[i, ,2]> TestParameterComb[i, 1])
  coverage_bootKappa2[i] <- mean(CI_Kappa2[i, ,1] < TestParameterComb[i,2] & CI_Kappa2[i, ,2] > TestParameterComb[i, 2])
  coverage_bootTau2[i] <- mean(CI_tau2[i, ,1] < TestParameterComb[i, 3] & CI_tau2[i, ,2] > TestParameterComb[i, 3])
}

# Print the coverage probabilities
print(coverage_boot_loc)
print(coverage_boot_scale)
print(coverage_boot_shape)






