rm(list=ls())

# VSS
# Small Sample Simulation Study
# for comparison of the quantiles
# for different quantile set evaluate at the fixed sample model trained


# -- Package Required --
library(fields)
library(maps)
library(ismev)
library(viridis)
library(ggplot2)
library(tictoc)
library(viridis)
library(colorspace)
library(grDevices)

# Load the test parameter values
# Dimension:  20*20 = 400 grid size
y_true <- read.csv('~/Desktop/Revise results/Store test parameter and test sample/test_parameter_config_version.csv',
                   header=FALSE)

scale_true <- y_true[,2]  # true scale 
shape_true <- (-y_true[,3])  # true scale # previously we where using shape parameter c,
# here we work with Xi


# Just the true shape and scale parameters:
True_parameter <- cbind(scale_true, shape_true) # stacking the parameter vectors together to fit TPS for surface plot
dim(True_parameter) # 400 x 2 

### --- SET1 --- ###
# -- Sample size 30 -- :
shape_est_set1_30 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_shape_30_NN_set1.csv',
                              header=FALSE) # size 400x100

# here 400 is the parameter configuration, 100 is the replication
scale_est_set1_30 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_scale_30_NN_set1.csv',
                              header=FALSE) # size 400x100

square_error_shape_set1_30 <- (-shape_est_set1_30 - shape_true)^2
square_error_scale_set1_30 <- (scale_est_set1_30 - scale_true)^2

rmse_shape_set1_30 <- sqrt(apply(square_error_shape_set1_30, 1, mean))
rmse_scale_set1_30 <- sqrt(apply(square_error_scale_set1_30, 1, mean))

fit_tps_rmse_shape_set1_30 <- Tps(True_parameter,
                                  rmse_shape_set1_30,
                                  df=100)
surface(fit_tps_rmse_shape_set1_30,
        #zlim=c(0.05,0.48),
        main= 'RMSE OF SHAPE EST - SS 30') # image plot for the shape parameter

fit_tps_rmse_scale_set1_30 <- Tps(True_parameter,
                                  rmse_scale_set1_30,
                                  df=100)
surface(fit_tps_rmse_scale_set1_30,
        zlim=c(0,15),
        main= 'RMSE OF SCALE EST - SS 30') # image plot for the scale parameter


# -- Sample size 72 -- :
shape_est_set1_72 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_shape_72_NN_set1.csv',
                              header=FALSE)  # size 400x100

scale_est_set1_72 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_scale_72_NN_set1.csv',
                              header=FALSE) 

square_error_shape_set1_72 <- (-shape_est_set1_72 - shape_true)^2
square_error_scale_set1_72 <- (scale_est_set1_72 - scale_true)^2

rmse_shape_set1_72 <- sqrt(apply(square_error_shape_set1_72, 1, mean))
# min: 0.05386739
# max: 0.2475483
rmse_scale_set1_72 <- sqrt(apply(square_error_scale_set1_72, 1, mean))
# min: 0.01231527
# max: 11.07766
fit_tps_rmse_shape_set1_72 <- Tps(True_parameter,
                                  rmse_shape_set1_72,
                                  df=100)
surface(fit_tps_rmse_shape_set1_72,
        zlim=c(0,0.48),
        main= 'RMSE OF SHAPE EST - SS 72') # image plot for the shape parameter



fit_tps_rmse_scale_set1_72 <- Tps(True_parameter,
                                  rmse_scale_set1_72,
                                  df=100)
surface(fit_tps_rmse_scale_set1_72,
        zlim=c(0,15),
        main= 'RMSE OF SCALE EST - SS 72') # image plot for the scale parameter


# -- Sample size 173 -- :
shape_est_set1_173 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_shape_173_NN_set1.csv',
                               header=FALSE) # size 400x100
scale_est_set1_173 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_scale_173_NN_set1.csv',
                               header=FALSE) 

square_error_shape_set1_173 <- (-shape_est_set1_173 - shape_true)^2
square_error_scale_set1_173 <- (scale_est_set1_173 - scale_true)^2

rmse_shape_set1_173 <- sqrt(apply(square_error_shape_set1_173, 1, mean))
rmse_scale_set1_173 <- sqrt(apply(square_error_scale_set1_173, 1, mean))

fit_tps_rmse_shape_set1_173 <- Tps(True_parameter,
                                   rmse_shape_set1_173,
                                   df=100)
surface(fit_tps_rmse_shape_set1_173,
        #zlim=c(0,0.48),
        main= 'RMSE OF SHAPE EST - SS 173') # image plot for the shape parameter

fit_tps_rmse_scale_set1_173 <- Tps(True_parameter,
                                   rmse_scale_set1_173,
                                   df=100)
surface(fit_tps_rmse_scale_set1_173,
        #zlim=c(0,15),
        main= 'RMSE OF SCALE EST - SS 173')

# -- Sample size 416 -- :
shape_est_set1_416 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_shape_416_NN_set1.csv',
                               header=FALSE) # size 400x100
scale_est_set1_416 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_scale_416_NN_set1.csv',
                               header=FALSE) 

square_error_shape_set1_416 <- (-shape_est_set1_416 - shape_true)^2
square_error_scale_set1_416 <- (scale_est_set1_416 - scale_true)^2

rmse_shape_set1_416 <- sqrt(apply(square_error_shape_set1_416, 1, mean))
rmse_scale_set1_416 <- sqrt(apply(square_error_scale_set1_416, 1, mean))

fit_tps_rmse_shape_set1_416 <- Tps(True_parameter,
                                   rmse_shape_set1_416,
                                   df=100)
surface(fit_tps_rmse_shape_set1_416,
        main= 'RMSE OF SHAPE EST - SS 416') # image plot for the shape parameter

fit_tps_rmse_scale_set1_416 <- Tps(True_parameter,
                                   rmse_scale_set1_416,
                                   df=100)
surface(fit_tps_rmse_scale_set1_416,
        main='RMSE OF SCALE EST - SS 416') # image plot for the scale parameter


# -- Sample size 1000 -- :
shape_est_set1_1000 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_shape_1000_NN_set1.csv',
                                header=FALSE) # size 400x100
scale_est_set1_1000 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-1/VSS_pred_scale_1000_NN_set1.csv',
                                header=FALSE) # size 400x100

square_error_shape_set1_1000 <- (-shape_est_set1_1000 - shape_true)^2
square_error_scale_set1_1000 <- (scale_est_set1_1000 - scale_true)^2

rmse_shape_set1_1000 <- sqrt(apply(square_error_shape_set1_1000, 1, mean))
rmse_scale_set1_1000 <- sqrt(apply(square_error_scale_set1_1000, 1, mean))


fit_tps_rmse_shape_set1_1000 <- Tps(True_parameter,
                                    rmse_shape_set1_1000,
                                    df=100)
surface(fit_tps_rmse_shape_set1_1000,
        main= 'RMSE OF SHAPE EST - SS 1000') # image plot for the shape parameter

fit_tps_rmse_scale_set1_1000 <- Tps(True_parameter,
                                    rmse_scale_set1_1000,
                                    df=100)
surface(fit_tps_rmse_scale_set1_1000,
        main= 'RMSE OF SCALE EST - SS 1000') # image plot for the scale parameter


### --- SET2 --- ###- 
# -- Sample size 30 -- :
shape_est_set2_30 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_shape_30_NN_set2.csv',
                              header=FALSE) # size 400x100

# here 400 is the parameter configuration, 100 is the replication
scale_est_set2_30 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_scale_30_NN_set2.csv',
                              header=FALSE) # size 400x100


square_error_shape_set2_30 <- (-shape_est_set2_30 - shape_true)^2
square_error_scale_set2_30 <- (scale_est_set2_30 - scale_true)^2

rmse_shape_set2_30 <- sqrt(apply(square_error_shape_set2_30, 1, mean))
rmse_scale_set2_30 <- sqrt(apply(square_error_scale_set2_30, 1, mean))


fit_tps_rmse_shape_set2_30 <- Tps(True_parameter,
                                  rmse_shape_set2_30,
                                  df=100)
surface(fit_tps_rmse_shape_set2_30,
        main= 'RMSE OF SHAPE EST - SS 30') # image plot for the shape parameter

fit_tps_rmse_scale_set2_30 <- Tps(True_parameter,
                                  rmse_scale_set2_30,
                                  df=100)
surface(fit_tps_rmse_scale_set2_30,
        main= 'RMSE OF SCALE EST - SS 30') # image plot for the scale parameter


# -- Sample size 72 -- :
shape_est_set2_72 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_shape_72_NN_set2.csv',
                              header=FALSE)  # size 400x100

scale_est_set2_72 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_scale_72_NN_set2.csv',
                              header=FALSE) # size 400x100

square_error_shape_set2_72 <- (-shape_est_set2_72 - shape_true)^2
square_error_scale_set2_72 <- (scale_est_set2_72 - scale_true)^2

rmse_shape_set2_72 <- sqrt(apply(square_error_shape_set2_72, 1, mean))
rmse_scale_set2_72 <- sqrt(apply(square_error_scale_set2_72, 1, mean))

fit_tps_rmse_shape_set2_72 <- Tps(True_parameter,
                                  rmse_shape_set2_72,
                                  df=100)
surface(fit_tps_rmse_shape_set2_72,
        main= 'RMSE OF SHAPE EST - SS 72') # image plot for the shape parameter

fit_tps_rmse_scale_set2_72 <- Tps(True_parameter,
                                  rmse_scale_set2_72,
                                  df=100)
surface(fit_tps_rmse_scale_set2_72,
        main= 'RMSE OF SCALE EST - SS 72') # image plot for the scale parameter


# -- Sample size 173 -- :
shape_est_set2_173 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_shape_173_NN_set2.csv',
                               header=FALSE) # size 400x100
scale_est_set2_173 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_scale_173_NN_set2.csv',
                               header=FALSE) 

square_error_shape_set2_173 <- (-shape_est_set2_173 - shape_true)^2
square_error_scale_set2_173 <- (scale_est_set2_173 - scale_true)^2

rmse_shape_set2_173 <- sqrt(apply(square_error_shape_set2_173, 1, mean))
rmse_scale_set2_173 <- sqrt(apply(square_error_scale_set2_173, 1, mean))

fit_tps_rmse_shape_set2_173 <- Tps(True_parameter,
                                   rmse_shape_set2_173,
                                   df=100)

surface(fit_tps_rmse_shape_set2_173,
        main= 'RMSE OF SHAPE EST - SS 173') # image plot for the shape parameter

fit_tps_rmse_scale_set2_173 <- Tps(True_parameter,
                                   rmse_scale_set2_173,
                                   df=100)
surface(fit_tps_rmse_scale_set2_173,
        main= 'RMSE OF SCALE EST - SS 173') # image plot for the scale parameter

# -- Sample size 416 -- :
shape_est_set2_416 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_shape_416_NN_set2.csv',
                               header=FALSE) # size 400x100
scale_est_set2_416 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_scale_416_NN_set2.csv',
                               header=FALSE) 

square_error_shape_set2_416 <- (-shape_est_set2_416 - shape_true)^2
square_error_scale_set2_416 <- (scale_est_set2_416 - scale_true)^2

rmse_shape_set2_416 <- sqrt(apply(square_error_shape_set2_416, 1, mean))
rmse_scale_set2_416 <- sqrt(apply(square_error_scale_set2_416, 1, mean))


fit_tps_rmse_shape_set2_416 <- Tps(True_parameter,
                                   rmse_shape_set2_416,
                                   df=100)
surface(fit_tps_rmse_shape_set2_416,
        main= 'RMSE OF SHAPE EST - SS 416') # image plot for the shape parameter

fit_tps_rmse_scale_set2_416 <- Tps(True_parameter,
                                   rmse_scale_set2_416,
                                   df=100)
surface(fit_tps_rmse_scale_set2_416,
        main= 'RMSE OF SCALE EST - SS 416') # image plot for the scale parameter


# -- Sample size 1000 -- :
shape_est_set2_1000 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_shape_1000_NN_set2.csv',
                                header=FALSE) # size 400x100
scale_est_set2_1000 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-2/VSS_pred_scale_1000_NN_set2.csv',
                                header=FALSE) 

square_error_shape_set2_1000 <- (-shape_est_set2_1000 - shape_true)^2
square_error_scale_set2_1000 <- (scale_est_set2_1000 - scale_true)^2

rmse_shape_set2_1000 <- sqrt(apply(square_error_shape_set2_1000, 1, mean))
rmse_scale_set2_1000 <- sqrt(apply(square_error_scale_set2_1000, 1, mean))

fit_tps_rmse_shape_set2_1000 <- Tps(True_parameter,
                                    rmse_shape_set2_1000,
                                    df=100)
surface(fit_tps_rmse_shape_set2_1000,
        main= 'RMSE OF SHAPE EST - SS 1000') # image plot for the shape parameter

fit_tps_rmse_scale_set2_1000 <- Tps(True_parameter,
                                    rmse_scale_set2_1000,
                                    df=100)
surface(fit_tps_rmse_scale_set2_1000,
        main= 'RMSE OF SCALE EST - SS 1000') # image plot for the scale parameter


######  -- VSS - SET 3 --  ########
# ramp <- colorRamp(c("blue", "white", "red"))
# colorTable <- rgb( ramp(seq(0, 1, length = 50)), max = 255)

# -- Sample size 30 -- :
shape_est_set3_30 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_shape_30_NN_set3.csv',
                              header=FALSE) # size 400x100

# here 400 is the parameter configuration, 100 is the replication
scale_est_set3_30 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_scale_30_NN_set3.csv',
                              header=FALSE) 

square_error_shape_set3_30 <- (-shape_est_set3_30 - shape_true)^2
square_error_scale_set3_30 <- (scale_est_set3_30 - scale_true)^2

rmse_shape_set3_30 <- sqrt(apply(square_error_shape_set3_30, 1, mean))
rmse_scale_set3_30 <- sqrt(apply(square_error_scale_set3_30, 1, mean))


fit_tps_rmse_shape_set3_30 <- Tps(True_parameter,
                                  rmse_shape_set3_30,
                                  df=100)
surface(fit_tps_rmse_shape_set3_30,
        main= 'RMSE OF SHAPE EST - SS 30')
#zlim=c(0.08,0.36)) # image plot for the shape parameter

fit_tps_rmse_scale_set3_30 <- Tps(True_parameter,
                                  rmse_scale_set3_30,
                                  df=100)
surface(fit_tps_rmse_scale_set3_30,
        main= 'RMSE OF SCALE EST - SS 30')# image plot for the scale parameter


# -- Sample size 72 -- :
shape_est_set3_72 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_shape_72_NN_set3.csv',
                              header=FALSE)  # size 400x100

scale_est_set3_72 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_scale_72_NN_set3.csv',
                              header=FALSE) 

square_error_shape_set3_72 <- (-shape_est_set3_72 - shape_true)^2
square_error_scale_set3_72 <- (scale_est_set3_72 - scale_true)^2

rmse_shape_set3_72 <- sqrt(apply(square_error_shape_set3_72, 1, mean))
rmse_scale_set3_72 <- sqrt(apply(square_error_scale_set3_72, 1, mean))

fit_tps_rmse_shape_set3_72 <- Tps(True_parameter,
                                  rmse_shape_set3_72,
                                  df=100)
surface(fit_tps_rmse_shape_set3_72,
        main='RMSE OF SHAPE EST - SS 72') # image plot for the shape parameter

fit_tps_rmse_scale_set3_72 <- Tps(True_parameter,
                                  rmse_scale_set3_72,
                                  df=100)

surface(fit_tps_rmse_scale_set3_72,
        main= 'RMSE OF SCALE EST - SS 72') # image plot for the scale parameter


# -- Sample size 173 -- :
shape_est_set3_173 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_shape_173_NN_set3.csv',
                               header=FALSE) # size 400x100
scale_est_set3_173 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_scale_173_NN_set3.csv',
                               header=FALSE) # size 400x100

square_error_shape_set3_173 <- (-shape_est_set3_173 - shape_true)^2
square_error_scale_set3_173 <- (scale_est_set3_173 - scale_true)^2

rmse_shape_set3_173 <- sqrt(apply(square_error_shape_set3_173, 1, mean))
rmse_scale_set3_173 <- sqrt(apply(square_error_scale_set3_173, 1, mean))


fit_tps_rmse_shape_set3_173 <- Tps(True_parameter,
                                   rmse_shape_set3_173,
                                   df=100)
surface(fit_tps_rmse_shape_set3_173,
        main= 'RMSE OF SHAPE EST - SS 173') # image plot for the shape parameter

fit_tps_rmse_scale_set3_173 <- Tps(True_parameter,
                                   rmse_scale_set3_173,
                                   df=100)
surface(fit_tps_rmse_scale_set3_173,
        main= 'RMSE OF SCALE EST - SS 173') # image plot for the scale parameter

# -- Sample size 416 -- :
shape_est_set3_416 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_shape_416_NN_set3.csv',
                               header=FALSE) # size 400x100
scale_est_set3_416 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_scale_416_NN_set3.csv',
                               header=FALSE) # size 400x100

square_error_shape_set3_416 <- (-shape_est_set3_416 - shape_true)^2
square_error_scale_set3_416 <- (scale_est_set3_416 - scale_true)^2

rmse_shape_set3_416 <- sqrt(apply(square_error_shape_set3_416, 1, mean))
rmse_scale_set3_416 <- sqrt(apply(square_error_scale_set3_416, 1, mean))

fit_tps_rmse_shape_set3_416 <- Tps(True_parameter,
                                   rmse_shape_set3_416,
                                   df=100)
surface(fit_tps_rmse_shape_set3_416,
        main='RMSE OF SHAPE EST - SS 416') # image plot for the shape parameter

fit_tps_rmse_scale_set3_416 <- Tps(True_parameter,
                                   rmse_scale_set3_416,
                                   df=100)
surface(fit_tps_rmse_scale_set3_416,
        main= 'RMSE OF SCALE EST - SS 416') # image plot for the scale parameter


# -- Sample size 1000 -- :
shape_est_set3_1000 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_shape_1000_NN_set3.csv',
                                header=FALSE) # size 400x100
scale_est_set3_1000 <- read.csv('~/Desktop/Revise results/VSS-Estimated Parameter/Set-3/VSS_pred_scale_1000_NN_set3.csv',
                                header=FALSE) 

square_error_shape_set3_1000 <- (-shape_est_set3_1000 - shape_true)^2
square_error_scale_set3_1000 <- (scale_est_set3_1000 - scale_true)^2

rmse_shape_set3_1000 <- sqrt(apply(square_error_shape_set3_1000, 1, mean))
rmse_scale_set3_1000 <- sqrt(apply(square_error_scale_set3_1000, 1, mean))


fit_tps_rmse_shape_set3_1000 <- Tps(True_parameter,
                                    rmse_shape_set3_1000,
                                    df=100)

surface(fit_tps_rmse_shape_set3_1000,
        main= 'RMSE OF SHAPE EST - SS 1000') # image plot for the shape parameter

fit_tps_rmse_scale_set3_1000 <- Tps(True_parameter,
                                    rmse_scale_set3_1000,
                                    df=100)

surface(fit_tps_rmse_scale_set3_1000,
        main='RMSE OF SCALE EST - SS 1000') # image plot for the scale parameter


# Comparitive Plot
setwd("~/Desktop")
png("VSS-Scale-Set1-SmallSampleCase.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,12.5,10) + 0.3,
    oma=c(0.4,2.5,3,8))

surface(fit_tps_rmse_scale_set1_173,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=5,
      line=8,
      las=2)
mtext(expression(sigma),
      side=1,
      cex=4.5,
      line=7.5)
mtext('n=173',
      side=3,
      cex=4.5,
      line=3)
mtext('VSS-Set1',
      side=3,
      cex=5,
      line=8,
      adj=-0.1)
# n = 416
surface(fit_tps_rmse_scale_set1_416,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext('n=416',
      side=3,
      cex=4.5,
      line=3)
mtext(expression('RMSE('~sigma~')'),
      side=3,
      cex=5,
      line=8)
# n = 1000
surface(fit_tps_rmse_scale_set1_1000,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext('n=1000',
      side=3,
      cex=4.5,
      line=3)
dev.off()

setwd("~/Desktop")
png("VSS-Shape-Set1-SmallSampleCase.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,12.5,10) + 0.3,
    oma=c(0.4,2.5,3,8))

surface(fit_tps_rmse_shape_set1_30,
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=5,
      line=8,
      las=2)
# mtext(expression(sigma),
#       side=1,
#       cex=4.5,
#       line=7.5)
mtext('n=30',
      side=3,
      cex=4.5,
      line=3)
mtext('VSS-Set1',
      side=3,
      cex=5,
      line=8,
      adj=-0.1)
# n = 416
surface(fit_tps_rmse_shape_set1_72,
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))

axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext('n=72',
      side=3,
      cex=4.5,
      line=3)
mtext(expression('RMSE('~xi~')'),
      side=3,
      cex=5,
      line=9)
# n = 1000
surface(fit_tps_rmse_shape_set1_173,
        #zlim=c(0.005,0.18),
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext('n=173',
      side=3,
      cex=4.5,
      line=3)
dev.off()


setwd("~/Desktop/Revise results")
png("VSS-Scale-Set2-Vers1.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,12.5,10) + 0.3,
    oma=c(0.4,2.5,3,8))

surface(fit_tps_rmse_scale_set2_173,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=5,
      line=8,
      las=2)
mtext(expression(sigma),
      side=1,
      cex=4.5,
      line=7.5)
mtext('n=173',
      side=3,
      cex=4.5,
      line=3)
mtext('VSS-Set2',
      side=3,
      cex=5,
      line=8,
      adj=-0.1)

# n = 416
surface(fit_tps_rmse_scale_set2_416,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext('n=416',
      side=3,
      cex=4.5,
      line=3)
# mtext(expression('RMSE('~sigma~')'),
#       side=3,
#       cex=5,
#       line=8)

# n = 1000
surface(fit_tps_rmse_scale_set2_1000,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext('n=1000',
      side=3,
      cex=4.5,
      line=3)
dev.off()

setwd("~/Desktop")
png("VSS-Shape-Set2-SmallSampleCase.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,12.5,10) + 0.3,
    oma=c(0.4,2.5,3,8))

surface(fit_tps_rmse_shape_set2_30,
        #zlim=c(0.005,0.18),
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=5,
      line=8,
      las=2)
# mtext(expression(sigma),
#       side=1,
#       cex=4.5,
#       line=7.5)
mtext('n=30',
      side=3,
      cex=4.5,
      line=3)
mtext('VSS-Set2',
      side=3,
      cex=5,
      line=8,
      adj=-0.1)

# n = 416
surface(fit_tps_rmse_shape_set2_72,
        # zlim=c(0.005,0.18),
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext('n=72',
      side=3,
      cex=4.5,
      line=3)
# mtext(expression('RMSE('~xi~')'),
#       side=3,
#       cex=5,
#       line=8)

# n = 1000
surface(fit_tps_rmse_shape_set2_173,
        # zlim=c(0.005,0.18),
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext('n=173',
      side=3,
      cex=4.5,
      line=3)
dev.off()


setwd("~/Desktop/Revise results")
png("VSS-Scale-Set3-Vers1.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,12.5,10) + 0.3,
    oma=c(0.4,2.5,3,8))

surface(fit_tps_rmse_scale_set3_173,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=5,
      line=8,
      las=2)
mtext(expression(sigma),
      side=1,
      cex=4.5,
      line=7.5)
mtext('n=173',
      side=3,
      cex=4.5,
      line=3)
mtext('VSS-Set3',
      side=3,
      cex=5,
      line=8,
      adj=-0.1)
# n = 416
surface(fit_tps_rmse_scale_set3_416,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(sigma),
      side=1,
      cex=4.5,
      line=7.5)
mtext('n=416',
      side=3,
      cex=4.5,
      line=3)
# mtext(expression('RMSE('~sigma~')'),
#       side=3,
#       cex=5,
#       line=8)
# n = 1000
surface(fit_tps_rmse_scale_set3_1000,
        zlim=c(0,9),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(sigma),
      side=1,
      cex=4.5,
      line=7.5)

mtext('n=1000',
      side=3,
      cex=4.5,
      line=3)
dev.off()

setwd("~/Desktop")
png("VSS-Shape-Set3-SmallSampleCase.png",
    units="in", 
    width=32,
    height=10,
    res=200)
par(mfrow=c(1,3), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,12.5,10) + 0.3,
    oma=c(0.4,2.5,3,8))

surface(fit_tps_rmse_shape_set3_30,
        #zlim=c(0.005,0.18),
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=5,
      line=8,
      las=2)
mtext(expression(sigma),
      side=1,
      cex=4.5,
      line=7.5)
mtext('n=30',
      side=3,
      cex=4.5,
      line=3)
mtext('VSS-Set3',
      side=3,
      cex=5,
      line=8,
      adj=-0.1)

# n = 416
surface(fit_tps_rmse_shape_set3_72,
        # zlim=c(0.005,0.18),
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(sigma),
      side=1,
      cex=4.5,
      line=7.5)
mtext('n=72',
      side=3,
      cex=4.5,
      line=3)
# mtext(expression('RMSE('~xi~')'),
#       side=3,
#       cex=5,
#       line=8)

# n = 1000
surface(fit_tps_rmse_shape_set3_173,
        # zlim=c(0.005,0.18),
        zlim=c(0,0.25),
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=5,
        cex.main=7,
        axis.args=list(cex.axis=3.5,lwd=2))
axis(2,
     las=1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=5,
     lwd=2,
     padj=0.9)
mtext(expression(sigma),
      side=1,
      cex=4.5,
      line=7.5)
mtext('n=173',
      side=3,
      cex=4.5,
      line=3)
dev.off()






