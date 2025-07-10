## Simulate the LK fields and compare with matern gaussian field: 
## Working on adding the addition channel for the standardized Gaussian marginal for 
rm(list=ls())

## Required script
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")


## CNN est: Calling the training estimates:
file_path_trainC <- path.expand("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/Coverage/cnnEst-Cvg.npy")
##cnnEst-Cvg-fixedTau2-0.01.npy") 

## -- Import numpy and load the file -- ##
np <- import("numpy")
matTrainC <- np$load(file_path_trainC)

## -- Load the NumPy file -- ##
cnnEstC <- matTrainC
dim(cnnEstC)  # 50x1000x3


## -- Defining test set --: 

## True parameter
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/Coverage")
# save(parameterCombCvg,
#      file="PC-cvg-Test-Rep30-PC-Rep1000-fixedTau2-0.001.RData")
# save(extremalFieldsCvg, file="storeZRep-cvg-Test-Rep30-PC-Rep1000-fixedTau2-0.001.RData")

# load("PC-cvg-Test-Rep30-PC-Rep1000-fixedTau2-0.01.RData")

setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/Coverage")
load("PC-cvg-Test-Rep30-PC-Rep1000.RData")
dim(parameterCombCvg) # 50x3

# load("storeZRep-cvg-Test-Rep30-PC-Rep1000.RData")
# dim(extremalFieldsCvg)
# index <- (parameterCombCvg[,1] > 0.4 & 0.6 > parameterCombCvg[,1]) 
# slc_par <- parameterCombCvg[index,]
# slc_fields <- extremalFieldsCvg[index,  , ,]
# image.plot(matrix(slc_fields[1, 1, ,1], 16, 16))

plot(cnnEstC[3, ,1],
     cnnEstC[3, ,2],
     pch=19,
     col=alpha('grey', alpha=0.3))
points(parameterCombCvg[3,1],
       log(parameterCombCvg[3,2]-4),
       pch=19,
       col='magenta')

####################################################################
## Load the fitted model for CI computation:
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Proj2-Data2Evaluate.R")
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
kappaCoverageBin <- matrix(NA, nrow=50, ncol=1000)
xiCoverage <- rep(NA)
for(i in 1:dim(cnnEstC)[1])
{
  X_test <- data_frame(shape=cnnEstC[i, , 1],
                       kappa=0.5*cnnEstC[i, , 2],
                       tau=0.5*cnnEstC[i, ,3])
  
  
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
  kappaCoverageBin[i, ] <- as.integer(QRKappa2EstTest[1, ] < 0.5*log(parameterCombCvg[i,2]-4) & 
                            QRKappa2EstTest[2, ] > 0.5*log(parameterCombCvg[i,2]-4))
  
  xiCoverage[i] <- mean(QRXiEstTest[1, ] < parameterCombCvg[i,1] &
                        QRXiEstTest[2, ] > parameterCombCvg[i,1])
}

t_test_result <- t.test(c(kappaCoverageBin[10,]), mu = 0.95)
t_test_result


setwd("~/Desktop")
png("CoveragePlot-Tau2-0.01.png",
    units="in", 
    width=24.5,
    height=8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(6,5,4,5) + 0.3,
    oma=c(0.4,2.3,2.8,6))
set.panel(1,3)
plot(parameterCombCvg[,1],
     xiCoverage,
     pch=19,
     cex=4.5,
     ylim=c(0.6, 1),
     xlab='',
     ylab='',
     xaxt='n',
     yaxt='n',
     main='',
     col=alpha("#63B", 0.5))
text(x = min(parameterCombCvg[,1]) + 0.5,  # Adjust X position dynamically
     y = 0.94,  # Slightly above the line for better visibility
     labels = "95% Coverage",
     col = 'grey3',
     cex = 2,   # Adjust font size
     font = 2)  # Bold text
abline(h=0.95, 
       col='grey3',
       lty=2,
       lwd=5)
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
mtext('True parameter',
      side=1,
      cex=2.2,
      line=4.8,
      las=1)
mtext('Coverage probability',
      side=2,
      cex=2.2,
      las=3,
      line=4.8)
mtext(expression(xi),
      side=3,
      cex=3.3,
      line=1.2)

plot(0.5*log(parameterCombCvg[,2]-4),
     kappaCoverage,
     pch=19,
     cex=4.5,
     xlab='',
     ylab='',
     xaxt='n',
     yaxt='n',
     ylim=c(0.6, 1),
     main='',
     col=alpha("#63B", 0.5))
text(x = min(0.5*log(parameterCombCvg[,2]-4)) + 2.5,  # Adjust X position dynamically
     y = 0.94,  # Slightly above the line for better visibility
     labels = "95% Coverage",
     col = 'grey3',
     cex = 2,   # Adjust font size
     font = 2)  # Bold text
abline(h=0.95, 
       col='grey3',
       lty=2,
       lwd=5)
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
mtext('True parameter',
      side=1,
      cex=2.2,
      line=4.8,
      las=1)
mtext('Coverage probability',
      side=2,
      cex=2.2,
      las=3,
      line=4.8)
mtext(expression('log('~kappa~')'),
      side=3,
      cex=3,
      line=1.2)

plot(0.5*log(parameterCombCvg[,3]),
     tauCoverage,
     ylim=c(0.6, 1),
     pch=19,
     cex=4.5,
     xlab='',
     ylab='',
     xaxt='n',
     yaxt='n',
     main='',
     col=alpha("#63B", 0.5))
abline(h=0.95, 
       col='grey3',
       lty=2,
       lwd=5)
# Label for the reference line
text(x = min(0.5 * log(parameterCombCvg[,3])) + 3.6,  # Adjust X position dynamically
     y = 0.94,  # Slightly above the line for better visibility
     labels = "95% Coverage",
     col = 'grey3',
     cex = 2,   # Adjust font size
     font = 2)  # Bold text
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
mtext('True parameter',
      side=1,
      cex=2.2,
      line=4.8,
      las=1)
mtext('Coverage probability',
      side=2,
      cex=2.2,
      las=3,
      line=4.8)
mtext(expression('log('~tau~')'),
      side=3,
      cex=3,
      line=1.2)
dev.off()

#########################################################################################################

##### ------ SURFACE PROBABILITY PLOT ---- #####
### -- Surface plot of coverage probability--:
# True_parameter <- cbind(0.5*log(parameterCombCvg[,3]), 0.5*log(parameterCombCvg[,2]-4))
True_parameter <- cbind(0.5*log(parameterCombCvg[,2]-4), (parameterCombCvg[,1]))
dim(True_parameter) # 100 x 2 

fit_tps_cvg_shape <- Tps(True_parameter,
                         xiCoverage,
                         df=10)

# True_parameterKappa <- cbind(0.5*log(parameterCombCvg[,3]),
#                              0.5*log(parameterCombCvg[,2]-4))
# dim(True_parameterKappa) # 100 x 2 

fit_tps_cvg_logKappa <- Tps(True_parameter,
                            kappaCoverage,
                            df=10)


# True_parameterTau <- cbind(0.5*log(parameterCombCvg[,2]-4),
#                            0.5*log(parameterCombCvg[,3]))
# dim(True_parameterTau) # 100 x 2 

## Not using right now
# fit_tps_cvg_logTau <- Tps(True_parameter,
#                           tauCoverage,
#                           df=10)


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
setwd("~/Desktop")
png("CoverageSurfs-Tau0.01.png",
    units="in", 
    width=22, # changed from 33
    height=12,
    res=200)
par(mfrow=c(1,2), 
    mai=c(18.5,5,15,16),
    mar=c(8,10,10.5,8) + 0.3,
    oma=c(2,3,3,5))
surface(fit_tps_cvg_shape,
        zlim=c(0.78,1),
        col= color,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        legend.width=1.2,
        horizontal=TRUE,
        cex.main=7,
        axis.args=list(cex.axis=3.5,padj=0.9, lwd=2))
axis(2,
     las=1,
     cex.axis=3.5,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=3.5,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2,
      cex=4,
      line=8,
      las=2)
mtext(expression('log('~kappa~')'), #('log('~tau~')')
      side=1,
      cex=3.5,
      line=6)
mtext(expression('Coverage probability of'~xi),
      side=3,
      cex=3.5,
      line=3)

surface(fit_tps_cvg_logKappa,
        zlim=c(0.78,1),
        col=color,
        yaxt='n',
        xaxt='n',
        xlab='',
        ylab= '',
        horizontal=TRUE,
        legend.width=1.2,
        cex.main=7,
        axis.args=list(cex.axis=3.5,padj=0.9, lwd=2))
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
mtext(expression(xi),
      side=2,
      cex=4,
      line=8,
      las=2)
mtext(expression('log('~kappa~')'), #('log('~tau~')')
      side=1,
      cex=3.5,
      line=6)
dev.off()



# Ensure XB is properly initialized
# XB <- matrix(NA, 
#              nrow = nrow(parameterCombCvg), 
#              ncol = ncol(XB))
# dim(XB)
# for(i in 1:nrow(XB))
# {
#   # true_p <- 0.5 * log(parameterCombCvg[i, 2] - 4)
#   for(j in 1:ncol(XB))
#   {
#     XB[i, j] <- as.integer(QRKappaEst[i,j,1] < 0.5*log(parameterCombCvg[i,2]-4) & 
#                    QRKappaEst[i,j,2] > 0.5*log(parameterCombCvg[i,2]-4))
#       # as.integer(QRKappaEst[i, j, 1] < true_p & true_p < QRKappaEst[i, j, 2])
#   }
# }








###########################################################################################################
### -- Fitting logistic regression to see the effect true parameter values: xi, log Kappa, log Tau -- ###
# library(MASS)
# 
# # Define response variable: Coverage probability (Binary Indicator)
# # Assuming coverage probabilities are stored in tauCoverage, kappaCoverage, and xiCoverage
# coverage_binary_tau <- ifelse(tauCoverage >= 0.90, 1, 0)  # 1 if coverage ≥ 95%, else 0
# coverage_binary_kappa <- ifelse(kappaCoverage >= 0.90, 1, 0)
# coverage_binary_xi <- ifelse(xiCoverage >= 0.90, 1, 0)
# 
# # Define predictor variables (True parameters)
# xi <- parameterCombCvg[, 1]
# log_kappa <- 0.5 * log(parameterCombCvg[, 2] - 4)
# log_tau <- 0.5 * log(parameterCombCvg[, 3])
# 
# # Create a data frame
# data_tau <- data.frame(coverage_logTau=coverage_binary_tau,
#                        xi,
#                        log_kappa,
#                        log_tau)
# data_kappa <- data.frame(coverage_logKappa=coverage_binary_kappa,
#                          xi,
#                          log_kappa, 
#                          log_tau)
# data_xi <- data.frame(coverage_Xi=coverage_binary_xi, 
#                       xi,
#                       log_kappa,
#                       log_tau)
# 
# # Fit Logistic Regression for each coverage
# model_tau <- glm(coverage_logTau ~log_tau,
#                  data=data_tau,
#                  family=binomial)
# model_kappa <- glm(coverage_logKappa ~ log_kappa,
#                    data=data_kappa,
#                    family=binomial)
# model_xi <- glm(coverage_Xi ~ xi + log_kappa + log_tau, 
#                 data=data_xi, 
#                 family=binomial)
# 
# # Summarize results
# summary(model_tau)
# summary(model_kappa)
# summary(model_xi)
# 
# # Model diagnostics: Checking goodness-of-fit
# pseudo_R2_tau <- 1 - (model_tau$deviance / model_tau$null.deviance)
# pseudo_R2_kappa <- 1 - (model_kappa$deviance / model_kappa$null.deviance)
# pseudo_R2_xi <- 1 - (model_xi$deviance / model_xi$null.deviance)
# 
# cat("Pseudo R^2 for tau model:", pseudo_R2_tau, "\n")
# cat("Pseudo R^2 for kappa model:", pseudo_R2_kappa, "\n")
# cat("Pseudo R^2 for xi model:", pseudo_R2_xi, "\n")
# 
# # Visualizing the effect of true parameters
# par(mfrow = c(1, 3))
# plot(xiCoverage,
#      fitted(model_xi),
#      ylim=c(0.65,1),
#      main = expression(xi),
#      xlab = "True xi", 
#      ylab = "Predicted Coverage", 
#      col = "blue", pch = 19)
# plot(log_kappa,
#      fitted(model_kappa),
#      main = expression(log(kappa)), 
#      xlab = "True log(kappa)",
#      ylab = "Predicted Coverage",
#      col = "red",
#      pch = 19)
# plot(log_tau, 
#      fitted(model_tau), 
#      main = expression(log(tau)),
#      xlab = "True log(tau)", 
#      ylab = "Predicted Coverage", 
#      col = "green",
#      pch = 19)
# 
# 



# setwd("~/Desktop")
# png("CoverageSurfs-Tau0.01.png",
#     units="in", 
#     width=22, # changed from 33
#     height=12,
#     res=200)
# par(mfrow=c(1,2), 
#     mai=c(18.5,5,15,16),
#     mar=c(8,10,10.5,8) + 0.3,
#     oma=c(0.6,3,3,5))
# surface(fit_tps_cvg_shape,
#         #zlim=c(0.6,1),
#         col= color,
#         yaxt='n',
#         xaxt='n',
#         xlab='',
#         ylab= '',
#         legend.width=1.2,
#         horizontal=TRUE,
#         cex.main=7,
#         axis.args=list(cex.axis=3.5,padj=0.9, lwd=2))
# axis(2,
#      las=1,
#      cex.axis=5,
#      lwd=2,
#      padj=0.9)
# axis(1,
#      cex.axis=5,
#      lwd=2,
#      padj=0.9)
# mtext(expression('log('~kappa~')'),
#       side=2,
#       cex=4,
#       line=8,
#       las=3)
# mtext(expression(xi), #('log('~tau~')')
#       side=1,
#       cex=4,
#       line=8.2)
# mtext(expression('Coverage probability of'~xi),
#       side=3,
#       cex=3.5,
#       line=3)
# 
# surface(fit_tps_cvg_logKappa,
#         #zlim=c(0.6,1),
#         col=color,
#         yaxt='n',
#         xaxt='n',
#         xlab='',
#         ylab= '',
#         horizontal=TRUE,
#         legend.width=1.2,
#         cex.main=7,
#         axis.args=list(cex.axis=3.5,padj=0.9, lwd=2))
# axis(2,
#      las=1,
#      cex.axis=5,
#      lwd=2,
#      padj=0.9)
# axis(1,
#      cex.axis=5,
#      lwd=2,
#      padj=0.9)
# mtext(expression('Coverage probability of ln('~kappa~')'),
#       side=3,
#       cex=3.5,
#       line=3)
# mtext(expression(xi), #('log('~tau~')')
#       side=1,
#       cex=4,
#       line=8.2)
# mtext(expression('log('~kappa~')'),
#       side=2,
#       cex=4,
#       line=7.5)
# 
# surface(fit_tps_cvg_logTau,
#         col=color,
#         yaxt='n',
#         xaxt='n',
#         xlab='',
#         ylab= '',
#         legend.width=1.2,
#         cex.main=7,
#         horizontal=TRUE,
#         axis.args=list(cex.axis=3.5,padj=0.9, lwd=2))
# axis(2,
#      las=1,
#      cex.axis=5,
#      lwd=2,
#      padj=0.9)
# axis(1,
#      cex.axis=5,
#      lwd=2,
#      padj=0.9)
# # mtext(expression('log('~tau~')'),
# #       side=1,
# #       cex=4,
# #       line=8.2)
# # mtext(expression('log('~kappa~')'),
# #       side=2,
# #       cex=4,
# #       line=7.5)
# mtext(expression('log('~tau~')'),
#       side=1,
#       cex=4,
#       line=8.2)
# mtext(expression('log('~kappa~')'),
#       side=2,
#       cex=4,
#       line=7.5)
# mtext(expression('Coverage probability of ln('~tau~')'),
#       side=3,
#       cex=3.5,
#       line=3)
# dev.off()




