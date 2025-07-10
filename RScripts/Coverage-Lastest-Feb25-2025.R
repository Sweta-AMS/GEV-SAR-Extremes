## Load the npy files for RMSE and coverga probability plot:
rm(list=ls())

## -- Load the true parameter configuration along with the dataset --:
setwd("~/Desktop/GEV-SAR/Data/Coverage")
load('PC-cvg-Test-Rep30-Feb22-2025.RData')
dim(parameterCombCvg) #  343   3


## -- CNN estimates --:
file_path_trainI <- path.expand("cnnEst-Cvg-I.npy") # Estimates over the train set
file_path_trainII <- path.expand("cnnEst-Cvg-II.npy") # 
file_path_trainIII <- path.expand("cnnEst-Cvg-III.npy") # 
file_path_trainIV <- path.expand("cnnEst-Cvg-IV.npy") #

## -- Import numpy and load the file -- ##
np <- import("numpy")
matTrainI <- np$load(file_path_trainI)
matTrainII <- np$load(file_path_trainII)
matTrainIII <- np$load(file_path_trainIII)
matTrainIV <- np$load(file_path_trainIV)

## -- Load the NumPy file -- ##
cnnEstTrainI <- as.array(matTrainI)
dim(cnnEstTrainI)  # 100x100x3

cnnEstTrainII <- as.array(matTrainII)
dim(cnnEstTrainII)  # 100x100x3

cnnEstTrainIII <- as.array(matTrainIII)
dim(cnnEstTrainIII)  # 100x100x3

cnnEstTrainIV <- as.array(matTrainIV)
dim(cnnEstTrainIV)  # 100x100x3

## Combined estimates: 
cnnEst <- array(NA, dim=c(343, 500, 3))
cnnEst[1:dim(cnnEstTrainI)[1], , ] <- cnnEstTrainI
cnnEst[(dim(cnnEstTrainI)[1]+1):(dim(cnnEstTrainI)[1]+dim(cnnEstTrainII)[1]), , ] <- cnnEstTrainII
cnnEst[(dim(cnnEstTrainI)[1] + dim(cnnEstTrainII)[1]+1):(dim(cnnEstTrainI)[1]+dim(cnnEstTrainII)[1]+dim(cnnEstTrainIII)[1]), , ] <- cnnEstTrainIII
cnnEst[(dim(cnnEstTrainI)[1] + dim(cnnEstTrainII)[1]+ dim(cnnEstTrainIII)[1] + 1):343, , ] <- cnnEstTrainIV
head(cnnEst)

## -- Spatial Fields -- :
# across each parameter configuration, we have 500 replications
square_error_Xi <- (cnnEst[, ,1] - parameterCombCvg[,1])^2
square_error_lnKappa <- (1/4)*(cnnEst[, ,2] - log(parameterCombCvg[,2]-4))^2
square_error_lnTau <- (1/4)*(cnnEst[, ,3] - log(parameterCombCvg[,3]))^2

rmse_Xi <- sqrt(apply(square_error_Xi, 1, mean))
rmse_lnKappa <- sqrt(apply(square_error_lnKappa, 1, mean))
rmse_lnTau <- sqrt(apply(square_error_lnTau, 1, mean))

## Xi
x <- cbind(log(parameterCombCvg[,2]-4), 
           parameterCombCvg[,1])
fit_tps_rmse_Xi <- Tps(x,
                       rmse_Xi,
                       df=10)

surface(fit_tps_rmse_Xi,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),
        main= 'RMSE OF SHAPE EST - SS 1000') # image plot for the shape parameter

## ln Kappa
x <- cbind(log(parameterCombCvg[,2]-4), 
           parameterCombCvg[,1])
fit_tps_rmse_lnKappa <- Tps(x,
                            rmse_lnKappa,
                            df=10)

surface(fit_tps_rmse_lnKappa,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),
        main= 'RMSE OF Xi - Rep1') # image plot for the shape parameter

## ln Tau
x <- cbind(log(parameterCombCvg[,2]-4), 
           log(parameterCombCvg[,3]))
fit_tps_rmse_lnTau <- Tps(x,
                          rmse_lnTau,
                          df=10)

surface(fit_tps_rmse_lnTau ,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),
        main= 'RMSE OF lnTau - Rep1') # image plot for the shape parameter


## -- Compute the coverage probability --: 
source("~/Desktop/GEV-SAR/RScripts/Proj2-Data2Evaluate.R")
shape_true_train <- trainParameters[ , 1]
kappa2_true_train <- (trainParameters[ , 2]-4)
tau2_true_train <- (trainParameters[ , 3])
length(tau2_true_train)


X_train <- data_frame(shape=shape_pred_train,
                      kappa=0.5*log(kappa2_pred_train),
                      tau=0.5*log(tau2_pred_train))
dim(X_train)

# Quantile levels
quantiles <- c(0.025, 0.975)

# Fit models for each quantile
modelsTau2 <- lapply(quantiles, function(tau){
  rq(0.5*log(tau2_true_train)~., data=X_train, tau = tau)
})
## -------------------------------------------------------------

## log Kappa:
# Fit models for each quantile
modelsKappa2 <- lapply(quantiles, function(tau){
  rq(0.5*log(kappa2_true_train)~., data=X_train, tau = tau)
})
## -------------------------------------------------------------


modelsXi <- lapply(quantiles, function(tau){
  rq((shape_true_train)~., data=X_train, tau = tau)
})
################################################################

## -- Calculate the coverage --:
tauCoverage <- rep(NA)
kappaCoverage <- rep(NA)
xiCoverage <- rep(NA)
for(i in 1:dim(cnnEst)[1])
{
  X_test <- data_frame(shape=cnnEst[i, , 1],
                       kappa=0.5*cnnEst[i, , 2],
                       tau=0.5*cnnEst[i, ,3])
  
  
  ## Kappa
  fittedKappa2Test <- lapply(modelsKappa2, function(model){
    predict(model, newdata=X_test)})
  
  QRKappa2EstTest <- matrix(NA,
                            nrow=length(quantiles),
                            ncol=nrow(X_test))
  for(t in 1:length(quantiles))
  {
    QRKappa2EstTest[t, ] <- fittedKappa2Test[[t]]
  }
  
  
  ## Tau
  fittedTau2Test <- lapply(modelsTau2, function(model){
    predict(model, newdata=X_test)})
  
  QRTau2EstTest <- matrix(NA,
                          nrow=length(quantiles),
                          ncol=nrow(X_test))
  for(t in 1:length(quantiles))
  {
    QRTau2EstTest[t, ] <- fittedTau2Test[[t]]
  }
  
  ## Xi
  fittedXiTest <- lapply(modelsXi, function(model){
    predict(model, newdata=X_test)})
  
  QRXiEstTest <- matrix(NA,
                        nrow=length(quantiles),
                        ncol=nrow(X_test))
  for(t in 1:length(quantiles))
  {
    QRXiEstTest[t, ] <- fittedXiTest[[t]]
  }
  
  tauCoverage[i] <- mean(QRTau2EstTest[1, ] < 0.5*log(parameterCombCvg[i,3]) &
                           QRTau2EstTest[2, ] > 0.5*log(parameterCombCvg[i,3]))
  
  kappaCoverage[i] <- mean(QRKappa2EstTest[1, ] < 0.5*log(parameterCombCvg[i,2]-4) & 
                             QRKappa2EstTest[2, ] > 0.5*log(parameterCombCvg[i,2]-4))
  #kappaCoverageBin[i, ] <- as.integer(QRKappa2EstTest[1, ] < 0.5*log(parameterCombCvg[i,2]-4) & 
  #                                      QRKappa2EstTest[2, ] > 0.5*log(parameterCombCvg[i,2]-4))
  
  xiCoverage[i] <- mean(QRXiEstTest[1, ] < parameterCombCvg[i,1] &
                          QRXiEstTest[2, ] > parameterCombCvg[i,1])
}


##### ------ SURFACE PROBABILITY PLOT ---- #####
True_parameter <- cbind( 0.5*log(parameterCombCvg[,3]), 0.5*log(parameterCombCvg[,2]-4) )
dim(True_parameter) # 100 x 2 

summary( xiCoverage)
summary(kappaCoverage)
summary(tauCoverage)

kappaCoverage[kappaCoverage < quantile(kappaCoverage, 0.25)] <- quantile(kappaCoverage, 0.25) - IQR(kappaCoverage)
tauCoverage[tauCoverage < quantile(tauCoverage, 0.25)] <- quantile(tauCoverage, 0.25) - IQR(tauCoverage)


fit_tps_cvg_shape <- Tps(True_parameter,
                         xiCoverage,
                         df=10)

fit_tps_cvg_logKappa <- Tps(True_parameter,
                            kappaCoverage,
                            df=10)




## Not using right now
fit_tps_cvg_logTau <- Tps(True_parameter,
                          tauCoverage,
                          df=10)


# Color
library(RColorBrewer)
jet.colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

# Generate the desired number of colors from this palette
nbcol <- 50
color <- jet.colors(nbcol)
ramp <- colorRamp(c("blue", "white", "red"))

# colorTable <- rgb( ramp(seq(0, 1, length = 50)), max = 255)
# ramp<- colorRamp(c("white","blue"))
# 
# colorTableII<- rgb( ramp(seq(0, 1,length=50)), max=255)

## Plot Coverage: 
png("CvgSurf.png",
    units="in", 
    width=33, # changed from 33
    height=12,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,8,10.5,8) + 0.3,
    oma=c(2,3,3,3))
surface(fit_tps_cvg_shape,
        zlim=c(0.78,1),
        levels=c(0.8, 0.9, 0.95),
        col= color,
        labcex = 2.5,
        lwd=2,
        #lablwd=2,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.mar = 3,
        legend.width=1.2,
        horizontal=TRUE,
        cex.main=7,
        axis.args=list(cex.axis=5, cex= 5, padj=0.9, lwd=2))
axis(2,
     las=1,
     cex.axis=3.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=3.5,
     lwd=2,
     padj=0.9)
mtext(expression('log('~kappa~')'),
      side=2,
      cex=4,
      line=6,
      las=3)
mtext(expression('log('~tau~')'), #('log('~tau~')')
      side=1,
      cex=3.8,
      line=8)
mtext(expression('Coverage probability of'~xi),
      side=3,
      cex=3.5,
      line=3)

surface(fit_tps_cvg_logKappa,
        zlim=c(0.78,1),
        levels=c(0.8, 0.9, 0.95),
        cex=3,
        labcex = 2.5,
        lwd=2,
        col=color,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        horizontal=TRUE,
        legend.width=1.2,
        cex.main=7,
        axis.args=list(cex.axis=5,padj=0.9, cex=5, lwd=2))
axis(2,
     las=1,
     cex.axis=3.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=3.5,
     lwd=2,
     padj=0.9)
mtext(expression('Coverage probability of ln('~kappa~')'),
      side=3,
      cex=3.5,
      line=3)
mtext(expression('log('~kappa~')'),
      side=2,
      cex=4,
      line=6,
      las=3)
mtext(expression('log('~tau~')'), #('log('~tau~')')
      side=1,
      cex=3.8,
      line=8)

surface(fit_tps_cvg_logTau,
        zlim=c(0.78,1),
        col= color,
        cex=3,
        yaxt='n',
        labcex = 2.5,
        lwd=2,
        xaxt='n',
        xlab='',
        ylab= '',
        levels=c(0.8, 0.9, 0.95),
        legend.width=1.2,
        horizontal=TRUE,
        cex.main=7,
        axis.args=list(cex.axis=5, padj=0.9, cex=5, lwd=2))
axis(2,
     las=1,
     cex.axis=3.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=3.5,
     lwd=2,
     padj=0.9)
mtext(expression('Coverage probability of ln('~tau~')'),
      side=3,
      cex=3.5,
      line=3)
mtext(expression('log('~kappa~')'),
      side=2,
      cex=3.8,
      line=6,
      las=3)
mtext(expression('log('~tau~')'), #('log('~tau~')')
      side=1,
      cex=3.8,
      line=8)
dev.off()
