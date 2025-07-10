rm(list=ls())
library(fields)
library(maps)
# Packages required for MLE 
library(ismev)
library(viridis)
library(ggplot2)
library(tictoc)

# do maximum likelihood approach over the Test Set
# load the parameters over the test set
# to be noted: here shape means 'c'
setwd("~/Desktop/Latest Scripts and Datasets - Paper I - March 6/FSS-Data+Results")
parameter_test <- read.csv("parameter_test.csv", 
                          header=FALSE)
colnames(parameter_test) <- c('loc', 'scale', 'shape')

# -- test sample based on  N_test = 10,000 -- #
#sample_test<- read.csv("sample_test.csv", 
#                       header = FALSE)
#colnames(sample_test) <- 1:1000
# ---------------------------------------------

# NN-output over the test set
pred_vals <- read.csv('fixed_sample_predict_NN.csv',
                      header= FALSE)
colnames(pred_vals)<- c('loc', 'scale', 'shape')
# ---------------------------------------------

# true parameter values over the test set
#test_vals<- read.csv('fixed_sample_test_values.csv',
#                     header= FALSE)
#colnames(test_vals) <- c('loc', 'scale', 'shape')
# ---------------------------------------------

# MLEs over the test set
# mle_matrix_over_test <- matrix(NA, nrow=nrow(sample_test), ncol=3)
# se_mle_test <- matrix(NA, nrow=nrow(sample_test), ncol=3)
# 
# tic()
# for(i in 1:nrow(sample_test))
# {
#   #fit<- fevd(as.numeric(sample_test[i,])) # using extReme PACKAGE!
#   #mle_matrix_over_test[i,] <- fit$results$par
#   fit<- gev.fit(as.numeric(sample_test[i,])) # using ismev PACKAGE!
#   mle_matrix_over_test[i,]<- fit$mle
#   se_mle_test[i,] <- fit$se
# }
# toc()
# # Recorder timing: 
# # 637.866 sec close to  10.631 mins
# 
# mle_matrix_over_test <- as.data.frame(mle_matrix_over_test)
# colnames(mle_matrix_over_test) <- c('loc', 'scale', 'shape')
# 
# se_matrix_over_test <- as.data.frame(se_mle_test)
# colnames(se_matrix_over_test) <- c('loc', 'scale', 'shape')
# 
# save(mle_matrix_over_test,
#      file="mle_test_sample.rda")
# save(se_matrix_over_test,
#      file="mle_se_test_sample.rda")

setwd("~/Desktop/Latest Scripts and Datasets - Paper I - March 6/FSS-Data+Results") 
load('mle_test_sample.rda')
load('mle_se_test_sample.rda')

# ---------------------------------------------
# avoid_index <- which(mle_matrix_over_test[,1]<1 |mle_matrix_over_test[,1]>50 )
# mle_vals_loc <- mle_matrix_over_test$loc[-avoid_index]
# test_vals_loc <- test_vals$loc[-avoid_index]
# pred_vals_loc <- pred_vals$loc[-avoid_index]
# test_vals_scale<- test_vals$scale[-which(mle_matrix_over_test[,2]>40)] 
# pred_vals_scale<- pred_vals$scale[-which(mle_matrix_over_test[,2]>40)] 
# mle_vals_scale<- mle_matrix_over_test$scale[-which(mle_matrix_over_test[,2]>40)] 
# index_to_avoid<- which(mle_matrix_over_test$shape> max(pred_vals$shape)|mle_matrix_over_test$shape< (min(pred_vals$shape)))
# index_to_avoid<- which(mle_matrix_over_test$shape> 0.4|mle_matrix_over_test$shape< (-0.4)) # close to 1000 cases
# mle_shape <- mle_matrix_over_test$shape[-(which(mle_matrix_over_test$shape==min(mle_matrix_over_test$shape)))] #xi
# test_shape<- -(test_vals$shape[-(which(mle_matrix_over_test$shape==min(mle_matrix_over_test$shape)))])#xi
# pred_shape<- -(pred_vals$shape[-(which(mle_matrix_over_test$shape==min(mle_matrix_over_test$shape)))]) #xi
# mle_vals_shape <- mle_matrix_over_test$shape
# test_vals_shape<- -(test_vals$shape)
# pred_vals_shape<- -(pred_vals$shape) #xi # getting outlier values
#index_to_avoid<- which(mle_matrix_over_test$shape> max(pred_vals$shape)|mle_matrix_over_test$shape< (min(pred_vals$shape)))
#index_to_avoid<- which(mle_matrix_over_test$shape> 0.4|mle_matrix_over_test$shape< (-0.4)) # close to 1000 cases
#Q_1_err <- quantile(error_loc_NN,
#                    probs=0.25)
#Q_3_err <- quantile(error_loc_NN, 
#                    probs=0.75)
#IQR_err <- Q_3_err - Q_1_err
#lbd <- Q_1_err-2*IQR_err
#ubd <- Q_3_err+2*IQR_err
# ---------------------------------------------

# -- Location: working range (1, 50) --
mle_loc <- mle_matrix_over_test[,1]
test_loc <- parameter_test$loc
pred_loc <- pred_vals$loc

# -- Scale: working range (0.1, 40) --
test_scale <- parameter_test$scale 
pred_scale <- pred_vals$scale 
mle_scale <- mle_matrix_over_test[,2]

# -- Shape: working range (-0.4, 0.4) --
mle_shape  <- mle_matrix_over_test[,3] #xi
test_shape <- -(parameter_test$shape) #xi
pred_shape <- -(pred_vals$shape) #xi

# -- Plots: --
error_shape_ML <- (mle_shape-test_shape)
error_shape_NN <- (pred_shape-test_shape)

error_scale_ML<- log(mle_scale)-log(test_scale)
error_scale_NN <- log(pred_scale)-log(test_scale)

error_loc_ML <- (mle_loc-test_loc)
error_loc_NN <- (pred_loc-test_loc)

# -- Error Plot --
# setwd("~/Desktop/Latest Scripts and Datasets - Paper I - March 6/FSS-Data+Results")
setwd("~/Downloads")
png("Error-ShapeFSS.png",
    units="in",
    width=28,
    height=11,
    res=400)
#par(mfrow=c(1,2), mai=c(2,4,1.2,0.6))
# par(mfrow= c(1,2),
#     mar=c(6, 9, 5, 1),
#     oma=c(4, 6.2, 0.4, 1))
par(mfrow= c(1,2),
    mar=c(6, 9, 5, 1),
    oma=c(4, 9, 0.4, 1))
# NN
bplot.xy(test_shape,
         error_shape_NN,
         xlim=c(-0.4, 0.4),
         # ylim=c(min(error_shape_NN,error_shape_ML),
         #        max(error_shape_NN,error_shape_ML)),
         ylim=c(-0.125, 0.125),
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
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
# mtext('NN',
#       side=3,
#       line=2,
#       cex=4.5)
abline(h=0, col= 'red', lwd=5, lty=2)
mtext(expression(hat(xi)-xi),
      side=2,
      cex=4.5,
      line=9,
      las=2)
mtext(expression(xi),
      side=1,
      cex=4.5,
      line=9)
# ML
bplot.xy(test_shape,
         error_shape_ML,
         xlim=c(-0.4, 0.4),
         # ylim=c(min(error_shape_NN,error_shape_ML),
         #        max(error_shape_NN,error_shape_ML)),
         ylim=c(-0.125, 0.125),
         pch=19,
         cex= 2,
         axes = FALSE,
         boxwex=0.4,
         xlab='',
         ylab='',
         yaxt='n')
box()
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
abline(h=0, col= 'red', lwd=5, lty=2)
mtext(expression(xi),
      side=1,
      cex=4.5,
      line=9)
# mtext('ML',
#       side=3,
#       cex=4.5,
#       line=2)
# Change the plot region color
dev.off()

# -- SCALE --
# Error Plot 
setwd("~/Desktop/Latest Scripts and Datasets - Paper I - March 6/FSS-Data+Results")
setwd("~/Downloads")
png("Error-ScaleFSS.png",
    units="in",
    width=28,
    height=11,
    res=400)
par(mfrow= c(1,2),
    mar=c(6, 9, 5, 1),
    oma=c(4, 9, 0.4, 1))
# par(mfrow= c(1,2),
#     mar=c(6, 9, 5, 1),
#     oma=c(4, 6.2, 0.4, 1))
bplot.xy(log(test_scale),
         error_scale_NN,
         xlim=c(log(0.1), log(40)),
         # ylim=c(min(error_scale_NN, error_scale_ML),
         #      max(error_scale_NN, error_scale_ML)),
         ylim=c(-0.125,0.125),
         pch=19,
         cex= 2,
         #outline=FALSE,
         boxwex=0.4,
         #cex.main= 5,
         xlab='',
         ylab='',
         axes = FALSE,
         #outbg = alpha('black', 0.3),   # Outliers color
         yaxt='n')
box()
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
# mtext('NN',
#       side=3,
#       line=2,
#       cex=4.5)
abline(h=0, col= 'red', lwd=5, lty=2)
mtext(expression('ln'~ hat(sigma)/sigma),
      side=2,
      cex=4,
      line=8.8,las=2)
mtext(expression('ln'~sigma),
      side=1,
      cex=4,
      line=9)
bplot.xy(log(test_scale),
         error_scale_ML,
         xlim=c(log(0.1), log(40)),
         #ylim=c(min(error_scale_NN, error_scale_ML),
         #       max(error_scale_NN, error_scale_ML)),
         ylim=c(-0.125,0.125),
         pch=19,
         cex= 2,
         #outline=FALSE,
         #cex.main= 5,
         axes = FALSE,
         boxwex=0.4,
         #outbg = alpha('black', 0.3),
         xlab='',
         ylab='',
         yaxt='n')

box()
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
# mtext('ML',
#       side=3,
#       line=2,
#       cex=4.5)
abline(h=0, col='red', lwd=5, lty=2)
mtext(expression('ln'~sigma),
      side=1,
      cex=4,
      line=9)
# Change the plot region color
dev.off()
#-----------------------------------

## -- Location -- ##
index_to_avoidI <- which(error_loc_NN> 0.4|error_loc_NN< (-0.4))
index_to_avoidII <- which(error_loc_ML> 0.4|error_loc_ML< (-0.4))
index_to_avoid <- union(index_to_avoidI, index_to_avoidII)
length(index_to_avoid)

# Error Plot 
setwd("~/Desktop/Latest Scripts and Datasets - Paper I - March 6/FSS-Data+Results")
png("Error-LocFSS.png",
    units="in",
    width=28,
    height=11,
    res=400)
#par(mfrow=c(1,2), mai=c(2,4,1.2,0.6))
par(mfrow= c(1,2),
    mar=c(6, 8, 5, 1),
    oma=c(4, 6, 0.4, 1))

bplot.xy(test_loc[-index_to_avoid],
         error_loc_NN[-index_to_avoid],
         xlim=c(1, 50),
         #ylim=c(-2,4),
         ylim=c(-0.41, 0.41),
         pch=19,
         cex= 2,
         #outline=FALSE,
         boxwex=0.4,
         cex.main=4,
         xlab='',
         ylab='',
         axes = FALSE,
         #outbg = alpha('black', 0.3),   # Outliers color
         yaxt='n')
box()
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
abline(h=0, col= 'red', lwd=5, lty=2)
mtext(expression(hat(mu)-mu),
      side=2,
      cex=4.5,
      line=9)
mtext(expression(mu),
      side=1,
      cex=4.5,
      line=9)
mtext('NN',
      side=3,
      line=2,
      cex=4.5)
bplot.xy(test_loc[-index_to_avoid],
         error_loc_ML[-index_to_avoid],
         xlim=c(1, 50),
         # ylim=c(min(error_loc_NN, error_loc_ML),
         #        max(error_loc_NN,error_loc_ML)),
         ylim=c(-0.41, 0.41),
         pch=19,
         #outline=FALSE,
         cex= 2,
         cex.main= 4,
         axes = FALSE,
         boxwex=0.4,
         #outbg = alpha('black', 0.3),
         xlab='',
         ylab='',
         yaxt='n')
box()
axis(1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
axis(2,
     las=1,
     cex.axis=4,
     lwd=2,
     padj=0.9)
abline(h=0, col= 'red', lwd=5, lty=2)
mtext(expression(mu),
      side=1,
      cex=4.5,
      line=9)
mtext('ML',
      side=3,
      cex=4.5,
      line=2)
# Change the plot region color
dev.off()

## hist()
# ##-- Bootstrap --##
# # Bootstrap estimate of first 1000 samples on the test set
# boot_estimate_test<- read.csv("boot_est_test.csv", 
#                               header=FALSE)
# colnames(boot_estimate_test)<- c('loc', 'scale', 'shape')
# 
# boot_estimate_SD_test<- read.csv("boot_sd_test.csv", 
#                                  header=FALSE)
# colnames(boot_estimate_SD_test)<- c('loc', 'scale', 'shape')
# 
# # NN-output over the test set
# pred_vals <- read.csv('fixed_sample_predict_NN.csv',
#                       header= FALSE)
# colnames(pred_vals)<- c('loc', 'scale', 'shape')
# # ---------------------------------------------
# # true parameter values over the test set
# test_vals<- read.csv('fixed_sample_test_values.csv',
#                      header= FALSE)
# colnames(test_vals) <- c('loc', 'scale', 'shape')
# setwd("~/Downloads")
# png("envelopePlot.png",
#     units="in",
#     width=44,
#     heigh=14,
#     res=600) # preferred 800
# par(mfrow=c(1,3), mai=c(1.8, 2, 0.8, 0.6))
# plot(test_vals$loc[1:1000], 
#      pred_vals$loc[1:1000],
#      cex=3,
#      pch=19,
#      xlab='',
#      ylab='',
#      xaxt='n',
#      yaxt='n')
# axis(2,
#      las=2,
#      cex.axis=4.5,
#      lwd=2,
#      padj=0.9)
# axis(1,
#      cex.axis=4.5,
#      lwd=2,
#      padj=0.9)
# lines(boot_estimate_test$loc , (boot_estimate_test$loc - 1.96*boot_estimate_SD_test$loc), 
#       col=alpha('magenta',0.1), 
#       lwd=2,
#       lty=2)
# lines(boot_estimate_test$loc , boot_estimate_test$loc + 1.96*boot_estimate_SD_test$loc,
#       col=alpha('magenta',0.1),
#       lwd=2,
#       lty=2)
# mtext(expression(hat(mu)),
#       side=2,
#       las=2,
#       cex=4.5,
#       line=9)
# mtext(expression(mu),
#       side=1,
#       cex=4.5,
#       line=9)
