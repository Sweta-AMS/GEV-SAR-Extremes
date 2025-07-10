## Generate plot for the test set:
rm(list=ls())
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Proj2-Source-Boxplot-1-10-30.R")


##Rep 30
setwd("~/Desktop")
png("BoxplotTestSet30-L.png",
    units="in", 
    width=27,
    height=8.8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(5,4,5.5,4) + 0.3,
    oma=c(1.5,3.5,4,8))
set.panel(1,3)
bplot.xy(shape_true_test,
         pred_parameterRep30[,1]-shape_true_test- biasXI30,
         ylim=c(-0.6, 0.5),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)  
mtext(expression(xi),
      side=3, 
      cex=3.2,
      line=1)

bplot.xy(0.5*log(kappa2_true_test),
         0.5*(pred_parameterRep30[,2])- 0.5*log(kappa2_true_test) - biaslnKappa30,
         ylim=c(-2.5,2.5),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)  
mtext(expression('log '~kappa),
      side=3,
      cex=3.2,
      line=1)
mtext(expression('r=30'),
      side=3,
      cex=3.2,
      line=5.5)

bplot.xy(0.5*log(tau2_true_test),
         0.5*pred_parameterRep30[,3]- 0.5*log(tau2_true_test)- biaslnTau30,
         ylim=c(-1.5,2),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)   
mtext(expression('log '~tau),
      side=3,
      cex=3.2,
      line=1)
dev.off()



##Rep 10
setwd("~/Desktop")
png("BoxplotTestSet10-L.png",
    units="in", 
    width=27,
    height=8.8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(5,4,5.5,4) + 0.3,
    oma=c(1.5,3.5,4,8))
set.panel(1,3)
bplot.xy(shape_true_test,
         pred_parameterRep10[,1]-shape_true_test- biasXI10,
         ylim=c(-0.6, 0.5),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)  
mtext(expression(xi),
      side=3, 
      cex=3.2,
      line=1)

bplot.xy(0.5*log(kappa2_true_test),
         0.5*(pred_parameterRep10[,2])- 0.5*log(kappa2_true_test) - biaslnKappa10,
         ylim=c(-2.5,2.5),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)  
mtext(expression('log '~kappa),
      side=3,
      cex=3.2,
      line=1)
mtext(expression('r=10'),
      side=3,
      cex=3.2,
      line=5.5)

bplot.xy(0.5*log(tau2_true_test),
         0.5*pred_parameterRep10[,3]- 0.5*log(tau2_true_test)- biaslnTau10,
         ylim=c(-1.5,2),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)   
mtext(expression('log '~tau),
      side=3,
      cex=3.2,
      line=1)
dev.off()



##Rep 1
setwd("~/Desktop")
png("BoxplotTestSet1-L.png",
    units="in", 
    width=27,
    height=8.8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(5,4,5.5,4) + 0.3,
    oma=c(1.5,3.5,4,8))
set.panel(1,3)
bplot.xy(shape_true_test,
         pred_parameterRep1[,1]-shape_true_test- biasXI1,
         ylim=c(-0.6, 0.5),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)   
mtext(expression(xi),
      side=3, 
      cex=3.2,
      line=1)

bplot.xy(0.5*log(kappa2_true_test),
         0.5*(pred_parameterRep1[,2])- 0.5*log(kappa2_true_test) - biaslnKappa1,
         ylim=c(-2.5,2.5),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)   
mtext(expression('log '~kappa),
      side=3,
      cex=3.2,
      line=1)
mtext(expression('r=1'),
      side=3,
      cex=3.2,
      line=5.5)

bplot.xy(0.5*log(tau2_true_test),
         0.5*pred_parameterRep1[,3]- 0.5*log(tau2_true_test)- biaslnTau1,
         ylim=c(-1.5,2),
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
      cex=2.5,
      line=5,
      las=1)
mtext('Bias',
      side=2,
      cex=2.5,
      las=3,
      line=4.8)  
mtext(expression('log '~tau),
      side=3,
      cex=3.2,
      line=1)
dev.off()

# ## -- Computing conference interval --:
# t <- 1
# image.plot(matrix(extremalFieldsScaled[t,10,,1], 16, 16))
# TestParameterComb[t,]
# 
# 
# 
# 
# 
# library(extRemes)
# input_data <- array(revd(256*30), dim = c(1, 16, 16, 30))  # Change shape according to your model
# 
# # Make predictions using the model
# predictions <- model$predict(input_data)
# 
# # Print the predictions
# print(predictions)
# dim( extremalFieldsScaled)
# 
# pred_parameter <- array(NA, c(ncomb, R, 3))
# dim(pred_parameter)
# 
# for(i in 1:ncomb)
# {
#   data <- array(extremalFieldsScaled[i, , ,], dim = c(R, 16, 16, 30))
#   pred_parameter[i, , ] <- model$predict(data)
#   
# }
# 
# pred_shape <- pred_parameter[, ,1]
# dim(pred_shape)
# 
# pred_kappa2 <- pred_parameter[, ,2]
# dim(pred_kappa2)
# 
# pred_tau2 <- pred_parameter[, ,3]
# dim(pred_tau2)
# 
# ## Bias correction: Quantile mapping in R:
# library(MASS)
# library(survival)
# library(fitdistrplus)
# 
# # Function to perform Quantile Mapping
# quantile_mapping <- function(observed, predicted) {
#   obs_ecdf <- ecdf(observed)
#   pred_ecdf <- ecdf(predicted)
#   quantiles <- pred_ecdf(predicted)
#   corrected_predicted <- quantile(obs_ecdf, quantiles)
#   return(corrected_predicted)
# }
# 
# # Apply Quantile Mapping
# cor_predict_shape <- quantile_mapping(TestParameterComb[,1] , pred_shape)
# pred_shapeBC <- matrix(cor_predict_shape, nrow=ncomb, ncol=R)
# 
# hist(cor_predict_shape)
# hist(TestParameterComb[,1])
# 
# cor_predict_kappa2 <- quantile_mapping((TestParameterComb[,2]-4) , pred_kappa2)
# pred_kappa2BC <- matrix(cor_predict_kappa2, nrow=ncomb, ncol=R)
# 
# hist(cor_predict_kappa2)
# hist(TestParameterComb[,2]-4)
# 
# cor_predict_tau2 <- quantile_mapping((TestParameterComb[,3]) , pred_tau2)
# pred_tau2BC <- matrix(cor_predict_tau2, nrow=ncomb, ncol=R)
# 
# hist(cor_predict_tau2)
# hist(TestParameterComb[,3])
# 
# 
# # # Step 4: Plot the original and corrected predictions against observed values
# # plot(observed, predicted, col = alpha("red", 0.5), pch = 19,
# #      xlab = "Observed Values", ylab = "Predicted Values",
# #      main = "Bias Correction using Quantile Mapping")
# # points(observed, corrected_predicted, col = alpha("blue", 0.5), pch = 19)
# # legend("topleft", legend = c("Original Predicted", "Corrected Predicted"),
# #        col = c("red", "blue"), pch = 19)
# 
# 
# 
# dim(pred_parameter) # 27 100   3
# 
# ## Doing bootstrap: 
# # Initialize arrays for bootstrap predictions
# source("~/LKsimFunc.R") # for generating bootstrap samples with replication
# boot_pred_shape <- matrix(NA, nrow=ncomb, ncol=R)
# boot_pred_Kappa2 <- matrix(NA, nrow=ncomb, ncol=R)
# boot_pred_tau2 <- matrix(NA, nrow=ncomb, ncol=R)
# 
# CI_shape <- array(NA, dim= c(ncomb, R, 2))
# CI_Kappa2 <- array(NA, dim= c(ncomb, R, 2))
# CI_tau2 <- array(NA, dim= c(ncomb, R, 2))
# 
# alpha <-  0.05
# 
# tic()
# # Perform bootstrap resampling
# for (i in 1:ncomb) {
#   cat('Loop:', i , "\n")
#   tic()
#   for (j in 1:R) {
#     cat('Inner loop:', j, "\n")
#     theta <- pred_parameter[i, j , ] 
#     theta[2] <- theta[2]+ 4
#     
#     # Bootstrap sample: 
#     generateBsample <- simulate_spatial_field(s=s, theta=theta, R=B, m=30)
#     Bsample <- array(generateBsample$extremalFieldsScaled, c(B, 16, 16, m))
#     
#     boot_est <- model$predict(Bsample)
#     
#     # Calculate the confidence interval and coverage probability
#     lower_percentile <- alpha/2
#     upper_percentile <- 1 - (alpha/2)
#     
#     lower_bound_shape <-  quantile(boot_est[ ,1], probs=lower_percentile)
#     upper_bound_shape <- quantile(boot_est[ ,1], probs=upper_percentile)
#     
#     lower_bound_Kappa2 <-  quantile(boot_est[ ,2], probs=lower_percentile)
#     upper_bound_Kappa2 <-  quantile(boot_est[ ,2], probs=upper_percentile)
#     
#     lower_bound_tau2 <- quantile(boot_est[ ,3],  probs=lower_percentile)
#     upper_bound_tau2 <- quantile(boot_est[ ,3],  probs=upper_percentile)
#     
#     
#     boot_pred_shape[i, j] <- quantile(boot_est[ ,1], 0.5)
#     boot_pred_Kappa2[i, j] <- quantile(boot_est[ ,2], 0.5)
#     boot_pred_tau2[i, j] <- quantile(boot_est[ ,3], 0.5)
#     
#     CI_shape[i, j, ] <- c(lower_bound_shape, upper_bound_shape)
#     CI_Kappa2[i, j, ] <- c(lower_bound_Kappa2, upper_bound_Kappa2)
#     CI_tau2[i, j, ] <- c(lower_bound_tau2, upper_bound_tau2)
#   }
#   toc()
# }
# toc()
# 
# # Calculate coverage probabilities
# coverage_bootShape <- numeric(ncomb)
# coverage_bootKappa2 <- numeric(ncomb)
# coverage_bootTau2 <- numeric(ncomb)
# 
# for (i in 1:ncomb) {
#   coverage_bootShape[i] <- mean(CI_shape[i, ,1] < TestParameterComb[i,1] & CI_shape[i, ,2]> TestParameterComb[i, 1])
#   coverage_bootKappa2[i] <- mean(CI_Kappa2[i, ,1] < TestParameterComb[i,2] & CI_Kappa2[i, ,2] > TestParameterComb[i, 2])
#   coverage_bootTau2[i] <- mean(CI_tau2[i, ,1] < TestParameterComb[i, 3] & CI_tau2[i, ,2] > TestParameterComb[i, 3])
# }






