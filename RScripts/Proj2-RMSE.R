rm(list=ls())

## -- Load the true parameter configuration along with the dataset --:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data")
load('storeZRep-Jan27-Test-Rep30-PC-Rep100-Ltst.RData')
dim(extremalFields) # 1000   60  256   30

load('ParameterConfiguration-Jan27-Test-Rep30-PC-Rep100-Ltst.RData')
dim(parameterComb) # 1000x3

## -- CNN estimates --:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/CNN-Est-forPCwithREP")
file_path_train1 <- path.expand("est_parameter_testRep1-RepAcrossPC60.npy") # Estimates over the train set
file_path_train10 <- path.expand("est_parameter_testRep10-RepAcrossPC60.npy") # 
file_path_train30 <- path.expand("est_parameter_testRep30-RepAcrossPC60.npy") # 

## -- Import numpy and load the file -- ##
np <- import("numpy")
matTrain1 <- np$load(file_path_train1)
matTrain10 <- np$load(file_path_train10)
matTrain30 <- np$load(file_path_train30)

## -- Load the NumPy file -- ##
cnnEstTrain1 <- as.array(matTrain1)
dim(cnnEstTrain1)  # 100x100x3

cnnEstTrain10 <- as.array(matTrain10)
dim(cnnEstTrain10)  # 100x100x3

cnnEstTrain30 <- as.array(matTrain30)
dim(cnnEstTrain30)  # 100x100x3

## -- Spatial Fields wiht Rep 1-- :
# across each parameter configuration, we have 60 replications
square_error_Xi_rep1 <- (cnnEstTrain1[, ,1] - parameterComb[,1])^2
square_error_lnKappa_rep1 <- (1/4)*(cnnEstTrain1[, ,2] - log(parameterComb[,2]-4))^2
square_error_lnTau_rep1 <- (1/4)*(cnnEstTrain1[, ,3] - log(parameterComb[,3]))^2

rmse_Xi_rep1 <- sqrt(apply(square_error_Xi_rep1, 1, mean))
rmse_lnKappa_rep1 <- sqrt(apply(square_error_lnKappa_rep1, 1, mean))
rmse_lnTau_rep1 <- sqrt(apply(square_error_lnTau_rep1, 1, mean))

## Xi
x <- cbind(log(parameterComb[,2]-4), 
           parameterComb[,1])
fit_tps_rmse_Xi_rep1 <- Tps(x,
                            rmse_Xi_rep1,
                            df=10)

surface(fit_tps_rmse_Xi_rep1,
        main= 'RMSE OF SHAPE EST - SS 1000') # image plot for the shape parameter

## ln Kappa
x <- cbind(log(parameterComb[,2]-4), 
           parameterComb[,1])
fit_tps_rmse_lnKappa_rep1 <- Tps(x,
                                 rmse_lnKappa_rep1,
                                 df=10)

surface(fit_tps_rmse_lnKappa_rep1,
        main= 'RMSE OF Xi - Rep1') # image plot for the shape parameter

## ln Tau
x <- cbind(log(parameterComb[,2]-4), 
           log(parameterComb[,3]))
fit_tps_rmse_lnTau_rep1 <- Tps(x,
                               rmse_lnTau_rep1,
                               df=10)

surface(fit_tps_rmse_lnTau_rep1 ,
        main= 'RMSE OF lnTau - Rep1') # image plot for the shape parameter

library(RColorBrewer)
# Comparitive Plot
setwd("~/Desktop")
png("RMSE-Rep1.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5, 8,18,16),
    mar=c(8,13,12.5,10) + 0.3,
    oma=c(0.4,8,6,9))

# Xi
surface(fit_tps_rmse_Xi_rep1 ,
        zlim=c(0.03,0.15),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=4,
      line=6,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( '~xi~' )'),
      side=3,
      cex=3.5,
      line=3.5)



# ln kappa
surface(fit_tps_rmse_lnKappa_rep1,
        zlim=c(0.05,0.6),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=4,
      line=6,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( log('~kappa~') )'),
      side=3,
      cex=3.5,
      line=3.5)

mtext('r=1',
      side=3,
      cex=4,
      line=10)

## ln Tau
surface(fit_tps_rmse_lnTau_rep1,
        zlim=c(0.05,0.6),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression('log('~tau~')'),
      side=2,
      cex=3,
      line=5.5,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( log('~tau~') )'),
      side=3,
      cex=3.5,
      line=3.5)
dev.off()


#### --- Rep 10 ----:
# across each parameter configuration, we have 60 replications
square_error_Xi_rep10 <- (cnnEstTrain10[, ,1] - parameterComb[,1])^2
square_error_lnKappa_rep10 <- (1/4)*(cnnEstTrain10[, ,2] - log(parameterComb[,2]-4))^2
square_error_lnTau_rep10 <- (1/4)*(cnnEstTrain10[, ,3] - log(parameterComb[,3]))^2

rmse_Xi_rep10 <- sqrt(apply(square_error_Xi_rep10, 1, mean))
rmse_lnKappa_rep10 <- sqrt(apply(square_error_lnKappa_rep10, 1, mean))
rmse_lnTau_rep10 <- sqrt(apply(square_error_lnTau_rep10, 1, mean))

## Xi
x <- cbind(log(parameterComb[,2]-4), 
           parameterComb[,1])
fit_tps_rmse_Xi_rep10 <- Tps(x,
                            rmse_Xi_rep10,
                            df=10)

surface(fit_tps_rmse_Xi_rep10,
        main= 'RMSE OF SHAPE EST - SS 1000') # image plot for the shape parameter

## ln Kappa
x <- cbind(log(parameterComb[,2]-4), 
           parameterComb[,1])
fit_tps_rmse_lnKappa_rep10 <- Tps(x,
                                 rmse_lnKappa_rep10,
                                 df=10)

surface(fit_tps_rmse_lnKappa_rep10,
        main= 'RMSE OF Xi - Rep1') # image plot for the shape parameter

## ln Tau
x <- cbind(log(parameterComb[,2]-4), 
           log(parameterComb[,3]))
fit_tps_rmse_lnTau_rep10 <- Tps(x,
                               rmse_lnTau_rep10,
                               df=10)

surface(fit_tps_rmse_lnTau_rep10,
        main= 'RMSE OF lnTau - Rep1') # image plot for the shape parameter

library(RColorBrewer)
# Comparitive Plot
setwd("~/Desktop")
png("RMSE-Rep10.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5, 8,18,16),
    mar=c(8,13,12.5,10) + 0.3,
    oma=c(0.4,8,6,9))

# Xi
surface(fit_tps_rmse_Xi_rep10,
        zlim=c(0.03,0.15),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=4,
      line=6,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( '~xi~' )'),
      side=3,
      cex=3.5,
      line=3.5)



# ln kappa
surface(fit_tps_rmse_lnKappa_rep10,
        zlim=c(0.05,0.6),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=4,
      line=6,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( log('~kappa~') )'),
      side=3,
      cex=3.5,
      line=3.5)

mtext('r=10',
      side=3,
      cex=4,
      line=10)

## ln Tau
surface(fit_tps_rmse_lnTau_rep10,
        zlim=c(0.05,0.6),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression('log('~tau~')'),
      side=2,
      cex=3,
      line=5.5,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( log('~tau~') )'),
      side=3,
      cex=3.5,
      line=3.5)
dev.off()


#### --- Rep 30 ----:
# across each parameter configuration, we have 60 replications
square_error_Xi_rep30 <- (cnnEstTrain30[, ,1] - parameterComb[,1])^2
square_error_lnKappa_rep30 <- (1/4)*(cnnEstTrain30[, ,2] - log(parameterComb[,2]-4))^2
square_error_lnTau_rep30 <- (1/4)*(cnnEstTrain30[, ,3] - log(parameterComb[,3]))^2

rmse_Xi_rep30 <- sqrt(apply(square_error_Xi_rep30, 1, mean))
rmse_lnKappa_rep30 <- sqrt(apply(square_error_lnKappa_rep30, 1, mean))
rmse_lnTau_rep30 <- sqrt(apply(square_error_lnTau_rep30, 1, mean))

## Xi
x <- cbind(log(parameterComb[,2]-4), 
           parameterComb[,1])
fit_tps_rmse_Xi_rep30 <- Tps(x,
                             rmse_Xi_rep30,
                             df=10)

surface(fit_tps_rmse_Xi_rep30,
        main= 'RMSE OF SHAPE EST - SS 1000') # image plot for the shape parameter

## ln Kappa
x <- cbind(log(parameterComb[,2]-4), 
           parameterComb[,1])
fit_tps_rmse_lnKappa_rep30 <- Tps(x,
                                  rmse_lnKappa_rep30,
                                  df=10)

surface(fit_tps_rmse_lnKappa_rep30,
        main= 'RMSE OF Xi - Rep1') # image plot for the shape parameter

## ln Tau
x <- cbind(log(parameterComb[,2]-4), 
           log(parameterComb[,3]))
fit_tps_rmse_lnTau_rep30 <- Tps(x,
                                rmse_lnTau_rep30,
                                df=10)

surface(fit_tps_rmse_lnTau_rep30,
        main= 'RMSE OF lnTau - Rep1') # image plot for the shape parameter

library(RColorBrewer)
# Comparitive Plot
setwd("~/Desktop")
png("RMSE-Rep30.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5, 8,18,16),
    mar=c(8,13,12.5,10) + 0.3,
    oma=c(0.4,8,6,9))

# Xi
surface(fit_tps_rmse_Xi_rep30,
        zlim=c(0.03,0.15),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=4,
      line=6,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( '~xi~' )'),
      side=3,
      cex=3.5,
      line=3.5)



# ln kappa
surface(fit_tps_rmse_lnKappa_rep30,
        zlim=c(0.05,0.6),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=4,
      line=6,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( log('~kappa~') )'),
      side=3,
      cex=3.5,
      line=3.5)

mtext('r=30',
      side=3,
      cex=4,
      line=10)

## ln Tau
surface(fit_tps_rmse_lnTau_rep30,
        zlim=c(0.05,0.6),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        col=colorRampPalette(brewer.pal(11, "RdBu"))(100),  # Corrected
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
mtext(expression('log('~tau~')'),
      side=2,
      cex=3,
      line=5.5,
      las=3)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=3,
      line=7.2)
mtext(expression('RMSE( log('~tau~') )'),
      side=3,
      cex=3.5,
      line=3.5)
dev.off()

