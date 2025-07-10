## Generate plot for the test set:
rm(list=ls())

## Load true parameter configurations:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
load("TruePar-NoNugget-LikComp-01-28-25-Test.RData")
dim(parameterComb) # 32x32

true_Xi <- parameterComb[,1]
true_lnKappa <- 0.5*log(parameterComb[,2]-4)
true_lnKappa <- 0.5*log(parameterComb[,3])

## Load Likelihood Estimates:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")

load("storeMLE-for-Rep30-Sample-32x32grid.RData")
load("storeMLE-for-Rep10-Sample-32x32grid.RData")
load("storeMLE-for-Rep1-Sample-32x32grid.RData")

# load("storeMLE-for-Rep1-Sample-32x32grid-Latest.RData")
# load("storeMLE-for-Rep10-Sample-32x32grid-Latest.RData")
# load("storeMLE-for-Rep30-Sample-32x32grid-Latest.RData")

dim(storeMLEs)
dim(storeMLEs1)
dim(storeMLEsRep10)

source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/LatestLikelihoodFunc-Oct7-2024.R")

## Load CNN Estimates:
# Expand the file path
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE/CNN-Estimates")
file_path_30 <- path.expand("testEst-ML-n-30-Latest-01-28-25.npy")
file_path_10 <- path.expand("testEst-ML-n-10-Latest-01-28-25.npy")
file_path_1 <- path.expand("testEst-ML-n-1-Latest-01-28-25.npy")

# Import numpy and load the file
library(reticulate)
np <- import("numpy")

mat30 <- np$load(file_path_30)
mat10 <- np$load(file_path_10)
mat1 <- np$load(file_path_1)


# Load the NumPy file
cnnRep30 <- as.matrix(mat30)
dim(cnnRep30)


cnnRep10 <- as.matrix(mat10)
dim(cnnRep10)

cnnRep1 <- as.matrix(mat1)
dim(cnnRep1)

setwd("~/Desktop")
png("Bias-ML-Ltst.png",
    units="in",
    width=18, # c
    height=8.8,
    res=200)
par(mfrow=c(1,2),
    mar=c(5,4,4,3) + 0.3,
    oma=c(1.5,3.5,3.5,4))
set.panel(1,2)
bplot.xy((parameterComb[,1]),
         exp(storeMLEs[,1]) -parameterComb[,1],
         ylim=c(-0.2,0.8),
         # ylim=c(0.2, 1),
         cex=1.6,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
abline(h=0,
       col='magenta',
       lwd=3,
       lty=2)
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
      cex=2.4,
      line=5.4,
      las=1)

mtext('Bias',
      side=2,
      cex=2.4,
      las=3,
      line=4.2)
mtext(expression(xi),
      side=3,
      cex=2.6,
      line=1)
mtext(expression('ML Estimates'),
      cex=3,
      side=3,
      adj=1.7,
      line=4.5)

bplot.xy(sqrt((parameterComb[,2]-4)),
         exp(storeMLEs[,2])-sqrt((parameterComb[,2]-4)),
         cex=1.4,
         ylim=c(-0.3,0.62),
        # ylim=c(0, 2),
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
box()
abline(h=0,
       col='magenta',
       lwd=3,
       lty=2)
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
      cex=2.4,
      line=5.4,
      las=1)

mtext('Bias',
      side=2,
      cex=2.4,
      las=3,
      line=4.2)
mtext(expression(kappa),
      side=3,
      cex=2.6,
      line=1)
dev.off()



setwd("~/Desktop")
png("Bias-CNN-Ltst.png",
    units="in",
    width=18, # c
    height=8.8,
    res=200)
par(mfrow=c(1,2),
    mar=c(5,4,4,3) + 0.3,
    oma=c(1.5,3.5,3.5,4))
set.panel(1,2)
bplot.xy(parameterComb[,1],
         cnnRep30[,1]-parameterComb[,1],
         ylim=c(-0.2,0.8),
         cex=1.4,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
abline(h=0,
       col='magenta',
       lwd=3,
       lty=2)
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
      cex=2.4,
      line=5.4,
      las=1)

mtext('Bias',
      side=2,
      cex=2.4,
      las=3,
      line=4.2)
mtext(expression(xi),
      side=3,
      cex=2.6,
      line=1)
mtext(expression('CNN Estimates'),
      cex=2.6,
      side=3,
      adj=1.6,
      line=4.5)
bplot.xy(sqrt((parameterComb[,2]-4)),
         exp(0.5*cnnRep30[,2])-sqrt((parameterComb[,2]-4)),
         cex=1.4,
         ylim=c(-0.3,0.62),
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
box()
abline(h=0,
       col='magenta',
       lwd=3,
       lty=2)
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
      cex=2.4,
      line=5.4,
      las=1)

mtext('Bias',
      side=2,
      cex=2.4,
      las=3,
      line=4.2)
mtext(expression(kappa),
      side=3,
      cex=2.6,
      line=1)
dev.off()



setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Results-Figures")

##  --- Ratio of the ML/CNN Estimates ---: 
setwd("~/Desktop")
png("MLvsCNN-Xi-Ltst.png",
    units="in",
    width=27, 
    height=8.8,
    res=200)
par(mfrow=c(1,3),
    mar=c(5,5,5,5) + 0.3,
    oma=c(2.5, 4.5,4.5,8))
set.panel(1,3)
bplot.xy(parameterComb[,1],
         exp(storeMLEs1[,1])/cnnRep1[,1],
         cex=1.6,
         N=9,
         lwd=2,
         ylim=c(0.4,3),
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
abline(h=1, col='magenta', lwd=3, lty=2)
box()
axis(2,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.7,
      line=6,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2.7,
      #las=2,
      line=5)
mtext('r=1',
      side=3,
      cex=2.8,
      line=1)


bplot.xy(parameterComb[,1],
         exp(storeMLEsRep10[,1])/cnnRep10[,1],
         ylim=c(0.4,3),
         cex=1.6,
         N=9,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
box()
abline(h=1, col='magenta', lwd=3, lty=2)
axis(2,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.7,
      line=6,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2.7,
      #las=2,
      line=5)
mtext('r=10',
      side=3,
      cex=2.8,
      line=1)
mtext(expression('Comparison of MLE/CNN Estimate Ratios '~xi),
      cex=3,
      side=3,
      line=4.5)



bplot.xy(parameterComb[,1],
         exp(storeMLEs[,1])/(cnnRep30[,1]),
         ylim=c(0.4,3),
         cex=1.6,
         N=9,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
abline(h=1, col='magenta', lwd=3, lty=2)
box()
axis(2,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.7,
      line=6,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2.7,
      #las=2,
      line=5)
mtext('r=30',
      side=3,
      cex=2.8,
      line=1)
dev.off()




setwd("~/Desktop")
png("MLvsCNN-Kappa2-Ltst.png",
    units="in",
    width=27, 
    height=8.8,
    res=200)
par(mfrow=c(1,3),
    mar=c(5,5,5,5) + 0.3,
    oma=c(2.5, 4.5,4.5,8))
set.panel(1,3)


bplot.xy(parameterComb[,2]-4,
         exp(2*storeMLEs1[,2])/exp(cnnRep1[,2]),
         cex=1.6,
         N=12,
         lwd=2,
         ylim=c(0.4,3),
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')

abline(h=1, col='magenta', lwd=3, lty=2)
box()
axis(2,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.7,
      line=6,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2.7,
      #las=2,
      line=5)
mtext('r=1',
      side=3,
      cex=2.8,
      line=1)

bplot.xy(parameterComb[,2]-4,
         exp(2*storeMLEsRep10[,2])/exp(cnnRep10[,2]),
         ylim=c(0.4,3),
         cex=1.6,
         N=12,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
box()
abline(h=1, col='magenta', lwd=3, lty=2)
axis(2,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.7,
      line=6,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2.7,
      #las=2,
      line=5)
mtext('r=10',
      side=3,
      cex=2.8,
      line=1)
mtext(expression('Comparison of MLE/CNN Estimate Ratios '~kappa^2),
      cex=3,
      side=3,
      line=4.5)


bplot.xy(parameterComb[,2]-4,
         exp(2*storeMLEs[,2])/exp(cnnRep30[,2]),
         ylim=c(0.4,3),
         cex=1.6,
         lwd=2,
         N=12,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
abline(h=1, col='magenta', lwd=3, lty=2)
box()
axis(2,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.8,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.7,
      line=6,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2.7,
      #las=2,
      line=5)
mtext('r=30',
      side=3,
      cex=2.8,
      line=1)
dev.off()


# 
# 
# 
# 
# 
# 
# 

setwd("~/Desktop")
png("ML-Bias-Xi-Ltst.png",
    units="in",
    width=27, 
    height=8.8,
    res=200)
par(mfrow=c(1,3),
    mar=c(5,5,5,5) + 0.3,
    oma=c(2.5, 4.5,4.5,8))
set.panel(1,3)

bplot.xy(parameterComb[,1],
         exp(storeMLEs1[,1])-parameterComb[,1],
         cex=1.6,
         lwd=2,
         ylim=c(-0.1,1),
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
abline(h=0, col='magenta', lwd=3, lty=2)
box()
axis(2,
     las=1,
     cex.axis=3.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=3.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.8,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.8,
      las=3,
      line=5)
mtext('Rep 1',
      side=3,
      cex=2.8,
      line=1)

bplot.xy(parameterComb[,1],
         exp(storeMLEsRep10[,1])-parameterComb[,1],
         ylim=c(-0.2,1),
         cex=1.6,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
box()
abline(h=0, col='magenta', lwd=3, lty=2)
axis(2,
     las=1,
     cex.axis=3.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=3.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.8,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.8,
      las=3,
      line=5)
mtext('Rep 10',
      side=3,
      cex=2.8,
      line=1)
mtext(expression('MLE for '~xi),
      cex=3.2,
      side=3,
      line=4.5)

bplot.xy(parameterComb[,1],
         exp(storeMLEs[,1])-parameterComb[,1],
         ylim=c(-0.1,1),
         cex=1.6,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
abline(h=0, col='magenta', lwd=3, lty=2)
box()
axis(2,
     las=1,
     cex.axis=3.2,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=3.2,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2.8,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.8,
      las=3,
      line=5)
mtext('Rep 30',
      side=3,
      cex=2.8,
      line=1)
dev.off()


setwd("~/Desktop")
png("MLvsCNN-Kappa2-Ltst.png",
    units="in",
    width=27,
    height=8.8,
    res=200)
par(mfrow=c(1,3),
    mar=c(5,4,4,4) + 0.3,
    oma=c(1.5,3.5,4,8))
set.panel(1,3)

bplot.xy(parameterComb[,2]-4,
         exp(2*storeMLEs1[,2])/exp(cnnRep1[,2]),
         cex=1.6,
         lwd=2,
         ylim=c(0.4,3),
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')

abline(h=1, col='magenta', lwd=3, lty=2)
box()
axis(2,
     las=1,
     cex.axis=2.7,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.7,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext('Rep 1',
      side=3,
      cex=2.2,
      line=1)


bplot.xy(parameterComb[,2]-4,
         exp(2*storeMLEsRep10[,2])/exp(cnnRep10[,2]),
         ylim=c(0.4,3),
         cex=1.6,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
box()
abline(h=1, col='magenta', lwd=3, lty=2)
axis(2,
     las=1,
     cex.axis=2.7,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.7,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext('Rep 10',
      side=3,
      cex=2.2,
      line=1)
mtext(expression('Comparison of MLE/CNN Estimate Ratios '~kappa^2),
      cex=2.4,
      side=3,
      line=4.5)

bplot.xy(parameterComb[,2]-4,
         exp(2*storeMLEs[,2])/exp(cnnRep30[,2]),
         ylim=c(0.4,3),
         cex=1.6,
         lwd=2,
         yaxt='n',
         boxwex=0.6,
         axes = FALSE,
         xlab='',
         ylab= '')
abline(h=1, col='magenta', lwd=3, lty=2)
box()
axis(2,
     las=1,
     cex.axis=2.7,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=2.7,
     lwd=2,
     padj=0.9)
mtext('True',
      side=1,
      cex=2,
      line=5,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=2,
      las=3,
      line=5)
mtext('Rep 30',
      side=3,
      cex=2.2,
      line=1)

dev.off()


