## Simulate the LK fields and compare with matern gaussian field: 
## Working on adding the addition channel for the standardized Gaussian marginal for 
rm(list=ls())

setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")



## -- Evaluate the CNN model on different sample size --:

## -- Sample size 30 --:
# -- Plot for the CNN Estimates Sample Size 30 no replication across the parameter configuration --
TestEstn30 <- read.csv('~/Downloads/testParameterEstSampleSize30-NoRep.csv', 
                       header=FALSE)
dim(TestEstn30)

TrueEstn30 <- read.csv('~/Downloads/trueParameterEstSampleSize30-NoRep.csv',
                       header=FALSE)
dim(TrueEstn30)

setwd("~/Desktop")
png("Boxplotn30-NoRep.png",
    units="in", 
    width=24,
    height=8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(5,5,5,5) + 0.3,
    oma=c(2.5,3,5,10))
bplot.xy(TrueEstn30[,1],
         TestEstn30[,1]-TrueEstn30[,1],
         ylim=c(-.14,.14),
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression(xi),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
bplot.xy(0.5*log(TrueEstn30[,2]-4),
         0.5*(TestEstn30[,2]- log(TrueEstn30[,2]-4)),
         ylim=c(-0.6, 0.6),
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
abline(h=0, col='magenta')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression('log '~kappa),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
mtext('n=30',
      side=3,
      cex=2,
      # las=1,
      line=5)
bplot.xy(0.5*log(TrueEstn30[,3]),
         0.5*(TestEstn30[,3] -log(TrueEstn30[,3])),
         ylim=c(-0.6, 0.6), 
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
abline(h=0, col='magenta')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression('log '~tau),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
dev.off()

## Sample size 10
# -- Plot for the CNN Estimates Sample Size 30 no replication across the parameter configuration --
TestEstn10 <- read.csv('~/Downloads/testParameterEstSampleSize10-NoRep (1).csv', 
                       header=FALSE)
dim(TestEstn10)

TrueEstn10 <- read.csv('~/Downloads/trueParameterEstSampleSize10-NoRep (1).csv',
                       header=FALSE)
dim(TrueEstn10)

setwd("~/Desktop")
png("Boxplotn10-NoRep.png",
    units="in", 
    width=24,
    height=8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(5,5,5,5) + 0.3,
    oma=c(2.5,3,5,10))
bplot.xy(TrueEstn10[,1],
         TestEstn10[,1]-TrueEstn10[,1],
         ylim=c(-.14,.14),
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression(xi),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
bplot.xy(0.5*log(TrueEstn10[,2]-4),
         0.5*(TestEstn10[,2]- log(TrueEstn10[,2]-4)),
         ylim=c(-0.6, 0.6),
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
abline(h=0, col='magenta')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression('log '~kappa),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
mtext('n=10',
      side=3,
      cex=2,
      # las=1,
      line=5)
bplot.xy(0.5*log(TrueEstn10[,3]),
         0.5*(TestEstn10[,3] -log(TrueEstn10[,3])),
         ylim=c(-0.6, 0.6), 
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
abline(h=0, col='magenta')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression('log '~tau),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
dev.off()



## Sample size 1
# -- Plot for the CNN Estimates Sample Size 30 no replication across the parameter configuration --
TestEstn1 <- read.csv('~/Downloads/testParameterEstSampleSize1-NoRep.csv', 
                       header=FALSE)
dim(TestEstn1)

TrueEstn1 <- read.csv('~/Downloads/trueParameterEstSampleSize1-NoRep.csv',
                       header=FALSE)
dim(TrueEstn1)

setwd("~/Desktop")
png("Boxplotn1-NoRep.png",
    units="in", 
    width=24,
    height=8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(5,5,5,5) + 0.3,
    oma=c(2.5,3,5,10))
bplot.xy(TrueEstn1[,1],
         TestEstn1[,1]-TrueEstn1[,1],
         ylim=c(-.14,.14),
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
abline(h=0, col='magenta')
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression(xi),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
bplot.xy(0.5*log(TrueEstn1[,2]-4),
         0.5*(TestEstn1[,2]- log(TrueEstn1[,2]-4)),
         ylim=c(-0.6, 0.6),
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
abline(h=0, col='magenta')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression('log '~kappa),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
mtext('n=1',
      side=3,
      cex=2,
      # las=1,
      line=5)
bplot.xy(0.5*log(TrueEstn1[,3]),
         0.5*(TestEstn1[,3] -log(TrueEstn1[,3])),
         ylim=c(-0.6, 0.6), 
         pch=19,
         cex= 2,
         boxwex=0.4,
         xlab='',
         ylab='',
         axes = FALSE,
         yaxt='n')
abline(h=0, col='magenta')
box()
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext(expression('log '~tau),
      side=3,
      cex=2.4,
      # las=1,
      line=1)
dev.off()



## Not plotting RMSE for no rep but instead 

## -- RMSE of the estimates -- :
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





# ## Location
# M <- 16 # changed from 15
# sGrid <- list(x=1:M, y=1:M)
# s <- make.surface.grid(sGrid)
# plot(s)
# 
# ## Generating the training set    
# ## training grid size 50 and validation is 20, let's keep testing 10
# n <- nrow(s)
# print(paste0('Field size: ', n))
# 
# ## Latest used
# n_full <- 10000
# print(paste0('full size: ', n_full))
# 
# 
# ## -- Parameter configuration --:
# # Parameter Configuration for Training and Validation set:
# # set.seed(222)
# # AWghtGrid <- 4 + exp(runif(n=n_full, log(0.001), log(2))) 
# # summary(AWghtGrid)
# # hist(AWghtGrid)
# # 
# # set.seed(222)
# # lambdaGrid <- exp(runif(n=n_full, log(0.0001), log(0.1)))
# # summary(lambdaGrid)
# # hist(lambdaGrid)
# # 
# # set.seed(222)
# # shapeGrid <- 1- exp(runif(n=n_full, log(0.1), log(0.9)))
# # summary(shapeGrid)
# # hist(shapeGrid)
# 
# # Test parameter configuration
# set.seed(111)
# AWghtGrid <- 4 + exp(runif(n=n_full, log(0.001), log(2)))
# summary(AWghtGrid)
# hist(AWghtGrid)
# 
# set.seed(111)
# lambdaGrid <- exp(runif(n=n_full, log(0.0001), log(0.1)))
# summary(lambdaGrid)
# hist(lambdaGrid)
# 
# set.seed(222)
# shapeGrid <- runif(n=n_full, (0.01), (0.9))
# # 1- exp(runif(n=n_full, log(0.01), log(0.9)))
# summary(shapeGrid)
# hist(shapeGrid)
# 
# ## Parameter Configuration:
# parameterComb <- cbind(shapeGrid, AWghtGrid, lambdaGrid)
# dim(parameterComb)
# 
# 
# ### -- GENERATE TRAINING SET --:
# # Function to use a min-max scaling
# scale_matrix <- function(matrix_temp)
# {
#   # Step 1: Compute column-wise minimum and maximum
#   column_min <- apply(matrix_temp, 2, min)
#   column_max <- apply(matrix_temp, 2, max)
#   
#   # Step 2: Compute the difference between max and min for each column
#   column_max_min <- column_max - column_min
#   
#   # Step 3: Subtract the column-wise minimum from each element in the matrix
#   matrix_temp_scale <- sweep(matrix_temp, 2, column_min, "-")
#   
#   # Step 4: Divide each column by its respective max-min range
#   matrix_temp_scaled <- matrix_temp_scale / matrix(rep(column_max_min, each = nrow(matrix_temp)), nrow = nrow(matrix_temp))
#   
#   # Return the final scaled matrix
#   return(matrix_temp_scaled)
# }
# 
# # For replication size 10:
# source("~/Desktop/Tentative Results-Proj2/R script/LKrigSAREvd.R")
# m <- 30
# extremalFields <- array(NA,
#                         dim=c(n_full, n, m))
# extremalFieldsScl <- array(NA,
#                            dim=c(n_full, n, m))
# dim(extremalFields)
# dim(extremalFieldsScl)
# 
# tic()  # Start timing the entire process
# # Initialize result_list before starting the loop
# set.seed(123)
# for (i in 1:length(shapeGrid))
# {
#   cat('GEV parameter loop no', i, '\n')
#   
#   # Defining GEV parameters:
#   shape_val <- shapeGrid[i]
#   a_wght <- AWghtGrid[i] 
#   lambda_val <- lambdaGrid[i]
#   
#   # Setup LKrig
#   LKinfo <- LKrigSetup(s,
#                        a.wght=a_wght,
#                        nlevel=1,
#                        nu=1,
#                        NC=16,
#                        NC.buffer=4,
#                        fixedFunction=NULL)
#   
#   # Simulate SAR coefficients
#   coefSARList <- LKrigSAREvd(LKinfo,
#                              loc=1,
#                              scale=shape_val,
#                              shape=shape_val,
#                              M=m,
#                              asList=FALSE)
#     
#   coefSAR <- coefSARList$coefSAR
#     
#   # Simulate field
#   PHI <- LKrig.basis(s, LKinfo)
#   ySim <- (PHI %*% coefSAR)
#   cat('Dim of ySim: ', dim(ySim), '\n')
#     
#   # Add measurement error
#   scaleParameter <- log(1 + (lambda_val))  # Variance of the log-normal is tau2
#   locParameter <- (-scaleParameter / 2)
#     
#   # Set another seed for the log-normal matrix generation, unique for each iteration
#   LogNormalMat <- exp(matrix(rnorm(n*m,
#                                    mean=locParameter,
#                                    sd=sqrt(scaleParameter)),
#                               n, m))  # n x m matrix
#   dim(LogNormalMat)
#   
#   # Scale each column (if required by your context)
#   Z <- ySim * LogNormalMat
#   Z_scaled <- scale_matrix(Z)
#   
#   # Store the scaled matrix in the extremalFields array
#   extremalFields[i, , ] <- Z # removed Z_scaled
#   extremalFieldsScl[i, , ] <- Z_scaled
# }
# toc()  # End timing for the entire process
# dim(extremalFields)
# 
# ## Shuffle index in R
# # n_full <- 500
# 
# set.seed(123)
# index_shuffle <- sample(n_full)
# index_shuffle[1:10]
# 
# storeZRep <- extremalFieldsScl[1:n_full,  ,]
# dim(storeZRep)
# 
# 
# ## -- Parameter Configuration --:
# storeParameterRep <- parameterComb[1:n_full, ]
# dim(storeParameterRep)
# 
# summary(storeParameterRep)
# 
# ## -- Shuffling the dataset -- :
# storeZRep <- storeZRep[index_shuffle, , ]
# storeParameterRep <- storeParameterRep[index_shuffle, ]
# 
# ## -- Defining test set --: 
# dim(storeParameterRep)
# save(storeParameterRep, 
#      file="ParameterConfigurationSampleSize30-Dec7-Test-NoRep.RData")
# 
# 
# dim(storeZRep)
# save(storeZRep, file="storeZSampleSize30-Dec7-Test-NoRep.RData")


