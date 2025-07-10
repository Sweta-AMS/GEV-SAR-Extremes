rm(list=ls())

## Required script
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSARGaussian.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd (1).R")

# load('~/Downloads/redseatemperature_extra.rdata')
# load("~/Downloads/redseatemperature.rdata")
# print(dim(loc)) # dim: 16703 2


source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/LKrigSAREvd (1).R")

library( extRemes)
M <- 100 # changed from 15
sGrid <- list(x=1:M, y=1:M)
s <- make.surface.grid(sGrid)

plot(s, pch=19)
L<- 2
a.wght<- 4.1
LKinfo0<- LKrigSetup(s, NC=100, NC.buffer=5,
                     a.wght=4.5, nlevel=1,nu=1)
sGridList<- LKinfo0$latticeInfo$grid[[1]]
sGrid<- make.surface.grid( sGridList)

# the real model
LKinfo<- LKrigSetup(s,
                    NC=100/L,
                    NC.buffer=5,
                    a.wght=a.wght, 
                    nlevel=1,
                    nu=1)


setwd("~/Desktop")
png("SAREVD.png", 
    units = "in",
    width = 8,
    height = 10.875,  #10.875,
    res = 200)
par( oma=c( 0,0,0,4))
set.panel(3,2)
par(mar=c(1,1,2,3))
PHI<- LKrig.basis(sGrid, LKinfo)
set.seed( 222)
look1<- LKrigSARGaussian(LKinfo, M = 1)
look1<- look1/sd( look1)
ySim<- PHI%*%look1
ySim<- (ySim - min(ySim)) / ( max(ySim)- min(ySim))
image( as.surface( sGridList,c(ySim) ),
       axes=FALSE, xlab="",ylab="",
       col=turbo(256)
)
title("Gaussian", cex.main=2)

set.seed( 222)

look1<- LKrigSAREvd(LKinfo, M = 1, loc=0, scale=1, shape=0)
look1 <- look1$coefSAR
look1<- look1/sd(look1)
ySim<- PHI%*%look1
ySim<- (ySim - min(ySim)) / ( max(ySim)- min(ySim))
image( as.surface( sGridList,c(ySim), ),
       axes=FALSE, xlab="", ylab="",
       col=turbo(256))
title("Gumbel", cex.main=2)


set.seed( 222)

look1<- LKrigSAREvd(LKinfo, M = 1, loc=1, scale=0.15,  shape=.15)
look1 <- look1$coefSAR
look1<- look1/sd(look1)
ySim<- PHI%*%look1
ySim<- (ySim - min(ySim)) / ( max(ySim)- min(ySim))
image( as.surface( sGridList,c(ySim) ),
       axes=FALSE, xlab="", ylab="",
       col=turbo(256))
title(xi~'= 0.15' , cex.main=2.5)

set.seed( 222)
look1<- LKrigSAREvd(LKinfo, M = 1, loc=1, scale=0.3, shape=.3)
look1 <- look1$coefSAR
look1<- look1/sd(look1)
ySim<- PHI%*%look1
ySim<- (ySim - min(ySim)) / ( max(ySim)- min(ySim))
image( as.surface( sGridList,c(ySim) ),
       axes=FALSE, xlab="", ylab="",
       col=turbo(256))
title(xi~' = 0.3', cex.main=2.5)


set.seed( 222)
look1<- LKrigSAREvd(LKinfo, M = 1, loc=1, scale=0.5, shape=.5)
look1 <- look1$coefSAR
look1<- look1/sd(look1)
ySim<- PHI%*%look1
ySim<- (ySim - min(ySim)) / ( max(ySim)- min(ySim))
image( as.surface( sGridList,c(ySim) ),
       axes=FALSE, xlab="", ylab="",
       col=turbo(256))
title(xi~" = 0.5", cex.main=2.5)

set.seed( 222)
look1<- LKrigSAREvd(LKinfo, M = 1,loc=1, scale=0.8,  shape= .8)
look1 <- look1$coefSAR
look1<- look1/sd(look1)
ySim<- PHI%*%look1
ySim<- (ySim - min(ySim)) / ( max(ySim)- min(ySim))
image( as.surface( sGridList,c(ySim) ),
       axes=FALSE, xlab="", ylab="",
       col=turbo(256))
title(xi~" = 0.8", cex.main=2.5)

par( oma= c(0,0,0,1))
imagePlot(legend.only=TRUE, zlim=c(0,1), legend.width=3, axis.args=list(cex.axis=2,lwd=2),
          col=turbo(256))

dev.off()


### Now with nugget effect: 10 % 
### i.e., tau^2 = 0.001
tau2_val1 = 0.001
scaleParameter1 <- log(1 + (tau2_val1))  # Variance of the log-normal is tau2 and mean =1
locParameter1 <- (-scaleParameter1/ 2)

# Set another seed for the log-normal matrix generation, unique for each iteration
set.seed(123)
n <- nrow(sGrid)
LogNormalMat1 <- exp(rnorm(n,
                          mean=locParameter1,
                          sd=sqrt(scaleParameter1)))  # n x m matrix

length(LogNormalMat1)
hist(LogNormalMat1)

nuggetTerm0.001 <- LogNormalMat1/sd(LogNormalMat1)
image.plot(as.matrix(nuggetTerm0.001, 100,100))

### i.e., tau^2 = 0.01
tau2_val2 = 0.01
scaleParameter2 <- log(1 + (tau2_val2))  # Variance of the log-normal is tau2 and mean =1
locParameter2 <- (-scaleParameter2/ 2)

# Set another seed for the log-normal matrix generation, unique for each iteration
set.seed(123)
LogNormalMat2 <- exp(rnorm(n,
                           mean=locParameter2,
                           sd=sqrt(scaleParameter2)))  # n x m matrix

length(LogNormalMat2)
hist(LogNormalMat2)

nuggetTerm0.01 <- LogNormalMat2/sd(LogNormalMat2)
image.plot(as.matrix(nuggetTerm0.01, 100,100))

### i.e., tau^2 = 0.1
tau2_val3 = 0.1
scaleParameter3 <- log(1 + (tau2_val3))  # Variance of the log-normal is tau2 and mean =1
locParameter3 <- (-scaleParameter3/ 2)

# Set another seed for the log-normal matrix generation, unique for each iteration
set.seed(123)
LogNormalMat3 <- exp(rnorm(n,
                           mean=locParameter3,
                           sd=sqrt(scaleParameter3)))  # n x m matrix

length(LogNormalMat3)
hist(LogNormalMat3)

nuggetTerm0.1 <- LogNormalMat3/sd(LogNormalMat3)
image.plot(as.matrix(nuggetTerm0.1, 100,100))

### With nugget efffect:
setwd("~/Desktop")
png("SAREvdII+1percentNugget.png", 
    units = "in",
    width =15.875 ,
    height = 5,  #10.875 # changed from 8
    res = 200)
par(mfrow=c(1,3), 
    mar=c(1,1,5,3) + 0.3,
    oma=c( 0.5, 0.5, 2.5, 5))

par()
PHI<- LKrig.basis(sGrid, LKinfo)

set.seed( 222)
look1<- LKrigSAREvd(LKinfo, M = 1, loc=1, scale=0.5, shape=.5)
look1 <- look1$coefSAR
look1<- look1/sd(look1)
ySim<- PHI%*%look1
ySimN <- ySim*nuggetTerm0.001
ySimN <- (ySimN - min(ySimN)) / ( max(ySimN)- min(ySimN))
image( as.surface( sGridList,c(ySimN) ),
       axes=FALSE, xlab="", ylab="",
       col=turbo(256))
title(tau^2~"= 0.001", cex.main=3, line=2.4)


set.seed( 222)
look1<- LKrigSAREvd(LKinfo, M = 1, loc=1, scale=0.5,  shape=.5)
look1 <- look1$coefSAR
look1<- look1/sd(look1)
ySim<- PHI%*%look1
ySimN <- ySim*nuggetTerm0.01
ySimN <- (ySimN - min(ySimN)) / ( max(ySimN)- min(ySimN))
image( as.surface( sGridList,c(ySimN) ),
       axes=FALSE, xlab="", ylab="",
       col=turbo(256))
title(tau^2~"= 0.01", cex.main=3, line=2.4)

mtext(expression(xi~'= 0.5, '~kappa^2~'= 0.1'),
      side=3.5,
      cex=2,
      line=4)

set.seed( 222)
look1<- LKrigSAREvd(LKinfo, M = 1, loc=1, scale=0.5,  shape= .5)
look1 <- look1$coefSAR

look1<- look1/sd(look1)
ySim<- PHI%*%look1
ySimN <- ySim*nuggetTerm0.1
ySimN <- (ySimN - min(ySimN)) / ( max(ySimN)- min(ySimN))
image( as.surface( sGridList,c(ySimN) ),
       axes=FALSE, xlab="", ylab="",
       col=turbo(256))
title(tau^2~"= 0.1", cex.main=3, line=2.4)

par(oma= c(0,2,0,3))
imagePlot(legend.only=TRUE, zlim=c(0,1), legend.width=3, axis.args=list(cex.axis=3,lwd=2),
          col=turbo(256))

dev.off()


