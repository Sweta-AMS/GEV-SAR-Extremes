## Generate plot for the test set:
rm(list=ls())

## Load true parameter configurations:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
load("TruePar-NoNugget-LikComp-01-28-25-Test.RData")
dim(parameterComb)

## Load Likelihood Estimates:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
load("storeMLE-for-Rep30-Sample-32x32grid.RData")
load("storeMLE-for-Rep10-Sample-32x32grid.RData")
load("storeMLE-for-Rep1-Sample-32x32grid.RData")


## Load CNN Estimates:
# Expand the file path
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE/CNN-Estimates")
file_path_30 <- path.expand("testEst-ML-n-30-Latest-01-28-25.npy")
file_path_10 <- path.expand("testEst-ML-n-10-Latest-01-28-25.npy")
file_path_1 <- path.expand("testEst-ML-n-1-Latest-01-28-25.npy")

# Import numpy and load the file
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



## Generating plot:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Results-Figures")

setwd("~/Desktop")
png("MLvsCNN-Xi-Ltst-Rep30.png",
    units="in", 
    width=13,
    height=4,
    res=200)
par(mfrow=c(1,3), 
    mar=c(4,4,3,3) + 0.3,
    oma=c(0.4,2.5,3,8))
set.panel(1,3)
bplot.xy(parameterComb[,1],
         exp(storeMLEs[,1])/(cnnRep30[,1]),
         ylim=c(0.4,3))
abline(h=1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 30',
      side=3,
      line=1)

bplot.xy(parameterComb[,1],
         exp(storeMLEsRep10[,1])/cnnRep10[,1],
         ylim=c(0.4,3))
abline(h=1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 10',
      side=3,
      # las=1,
      line=1)
mtext(expression(xi),
      cex=2,
      side=3,
      line=3)
bplot.xy(parameterComb[,1],
         exp(storeMLEs1[,1])/cnnRep1[,1],
         #N = 16,
         ylim=c(0.4,3))
abline(h=1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 1',
      side=3,
      line=1)
dev.off()

setwd("~/Desktop")
png("MLvsCNN-Kappa2-Ltst-Rep30.png",
    units="in", 
    width=13,
    height=4,
    res=200)
par(mfrow=c(1,3), 
    mar=c(4,4,3,3) + 0.3,
    oma=c(0.4,2.5,3,8))
set.panel(1,3)
bplot.xy(parameterComb[,2]-4,
         exp(storeMLEs[,2])/exp(cnnRep30[,2]),
         ylim=c(0.4,3))
abline(h=1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 30',
      side=3,
      line=1)

bplot.xy(parameterComb[,2]-4,
         exp(storeMLEsRep10[,2])/exp(cnnRep10[,2]),
         ylim=c(0.4,3))
abline(h=1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 10',
      side=3,
      line=1)
mtext(expression(kappa^2),
      cex=2,
      side=3,
      line=3)
bplot.xy(parameterComb[,2]-4,
         exp(storeMLEs1[,2])/exp(cnnRep1[,2]),
         ylim=c(0.4,3))
abline(h=1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE/CNN',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 1',
      side=3,
      line=1)
dev.off()


### Just scatter plot:
setwd("~/Desktop")
png("Boxplot-ML-Xi-Ltst-Rep30.png",
    units="in", 
plot(storeMLEsSmall[1, ,2],
     exp(storeMLEsSmall[1, ,1]),
     xlim=c(min(storeMLEsSmall[1, ,2], cnnRep30[1, ,2]),
            max(storeMLEsSmall[1, ,2], cnnRep30[1, ,2])),
     ylim=c(min(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1]),
            max(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1])),
     pch=19,
     col=alpha("#63B", alpha=0.6))
points(cnnRep30[1, ,2],
       cnnRep30[1, ,1],
       xlim=c(min(storeMLEsSmall[1, ,2], cnnRep30[1, ,2]),
              max(storeMLEsSmall[1, ,2], cnnRep30[1, ,2])),
       ylim=c(min(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1]),
         max(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1])),
       pch=19,
       col=alpha("#8DD641", alpha=0.6))
points(log(parameterComb[1,2]-4),
       parameterComb[1,1],
       xlim=c(min(storeMLEsSmall[1, ,2], cnnRep30[1, ,2]),
              max(storeMLEsSmall[1, ,2], cnnRep30[1, ,2])),
       ylim=c(min(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1]),
              max(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1])),
       col= "#E8336C" ,
       lwd=3,
       pch=3)
bplot.xy(parameterComb[,1],
         exp(storeMLEs[,1]),
         ylim=c(0.2,1.2))
abline(0,1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 30',
      side=3,
      line=1)

bplot.xy(parameterComb[,1],
         exp(storeMLEsRep10[,1]),
         ylim=c(0.2,1.2))
abline(0,1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 10',
      side=3,
      # las=1,
      line=1)
mtext(expression(xi),
      cex=2,
      side=3,
      line=3)
bplot.xy(parameterComb[,1],
         exp(storeMLEs1[,1]),
         ylim=c(0.2,1.2))
abline(0,1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 1',
      side=3,
      line=1)
dev.off()



png("Boxplot-ML-lnKappa-Ltst-Rep30.png",
    units="in", 
    width=13,
    height=4,
    res=200)
par(mfrow=c(1,3), 
    mar=c(4,4,3,3) + 0.3,
    oma=c(0.4,2.5,3,8))
set.panel(1,3)
bplot.xy(log(parameterComb[,2]-4),
         (storeMLEs[,2]))
abline(0,1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 30',
      side=3,
      line=1)

bplot.xy(log(parameterComb[,2]-4),
         (storeMLEsRep10[,2]))
abline(0,1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 10',
      side=3,
      line=1)
mtext(expression('log('~kappa~')'),
      cex=2,
      side=3,
      line=3)
bplot.xy(log(parameterComb[,2]-4),
         (storeMLEs1[,2]))
abline(0,1, col='magenta')
mtext('True',
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext('MLE',
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext('Rep 1',
      side=3,
      line=1)
dev.off()


